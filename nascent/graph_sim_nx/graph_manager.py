# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This file is part of Nascent.
##
##    Nascent is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Affero General Public License as
##    published by the Free Software Foundation, either version 3 of the
##    License, or (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Affero General Public License for more details.
##
##    You should have received a copy of the GNU Affero General Public License
##    along with this program. If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=C0103,C0301,W0511,W0611,W0142

## TODO: Enable all pylint checks

"""

Module with all logic related to the structural graph representing the system.
This is used for calculation of activity, etc.


"""

from __future__ import absolute_import, print_function, division

import itertools
import math
import os
import pdb
from collections import deque
from functools import reduce
from itertools import groupby
from math import log as ln
from operator import itemgetter, mul

import networkx as nx
from networkx.algorithms.shortest_paths import shortest_path

from .system_graphs import InterfaceGraph
from .nx_utils import draw_graph_and_save
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        AVOGADRO_VOLUME_NM3,
                        HELIX_XOVER_DIST, HELIX_WIDTH, HELIX_STACKING_DIST,
                        chain_model_gamma_exponent as gamma)
from .domain import DomainEnd


## Module constants:
SEGMENT_STIFFNESS, SEGMENT_LENGTHS, SEGMENT_CONTOUR_LENGTH, SEGMENT_LENGTH_SQUARED = 0, 1, 0, 1


def prod(sequence):
    """ Return the product of all elements in sequence. """
    return reduce(mul, sequence, 1)



def determine_interaction(source, target, edge):
    """
    :source:, :target: DomainEnd instances connected by edge.
    :edge: NetworkX edge connecting :source: and :target:.
    """
    # Fallback method: Manually determine what type of interaction we have based on source, target:
    print("WARNING: ends5p3p_graph edge (%s - %s) has no interaction attribute: %s" % (
        source, target, edge))
    if target == source.stack_partner:
        interaction = STACKING_INTERACTION
        length = 1
        length_sq = 1
    elif target == source.hyb_partner:
        interaction = HYBRIDIZATION_INTERACTION
        length = 1 # one nm from one helix to the next. We really don't know for sure because it turns.
        length_sq = 1
    # Check phosphate backbone up- and downstream (5' and 3').
    elif target in (source.pb_upstream, source.pb_downstream):
        if source.hyb_partner and \
            ((source.end == "5p" and target == source.pb_downstream) or
             (source.end == "3p" and target == source.bp_upstream)):
            # Above could probably be done with an XOR:
            # (source.end == "5p") != (target == source.pb_downstream) # boolean xor
            # (source.end == "5p") ^ (target == source.pb_downstream)  # bitwise xor
            # We have a stacked duplex:
            interaction = STACKING_INTERACTION
            length = source.domain.ds_dist_ee_nm
            length_sq = source.domain.ds_dist_ee_sq
        else:
            # We have a single-stranded domain:
            interaction = PHOSPHATEBACKBONE_INTERACTION
            length = source.domain.ss_dist_ee_nm     # mean end-to-end length; not contour length
            length_sq = source.domain.ss_dist_ee_sq
    else:
        raise ValueError("Could not determine interaction between %s and %s" % (source, target))
    return interaction


def determine_end_end_distance(source, target, edge_attrs, interaction):
    """
    Returns a two-tuple of (length_sq, length_nm).
    WAIT: We *need* this info to already be in the graph in order to efficiently calculate the shortest path!
    """
    stiffness = 0    # 0 for ss, 1 for ds, optionally at some point 2 for helix-bundle.
    edge_required_attrs = ("dist_ee_sq", "dist_ee_nm", "stiffness")
    if 'dist_ee_sq' in edge_attrs:
        length_sq = edge_attrs['dist_ee_sq']
        length_nm = edge_attrs.get('dist_ee_nm') or math.sqrt(length_sq) # short-circut; don't calculate default value.
        stiffness = edge_attrs.get('stiffness')
        if all(at in edge_attrs for at in edge_required_attrs):
            return length_nm, length_sq, stiffness
    print("\n\nWARNING: Edge (%s, %s) with attrs %s is missing attributes %s! (It should have to calculate shortest path!)\n"
          % (source, target, edge_attrs, [at for at in edge_required_attrs if at not in edge_attrs]))

    # We need to determine length, length_sq and stiffness manually... (Mostly for reference...)
    if interaction is None:
        interaction = determine_interaction(source, target, edge_attrs)

    ## TODO: Move this logic the the graph edge creation location -- ComponentMgr (hybridize, stack, etc).
    ## Note that we use multigraphs. We can easily have
    if interaction == PHOSPHATEBACKBONE_INTERACTION:
        if source.domain == target.domain:
            ## "Long 5p-3p edge", connecting the 5p end to the 3p end of the same domain.
            # Length depends on whether domain is hybridized:
            if source.domain.partner is not None:
                # hybridized - consider as rigid stack:
                length = source.domain.ds_dist_ee_nm
                length_sq = source.domain.ds_dist_ee_sq
                stiffness = 1
            else:
                length = source.domain.ss_dist_ee_nm     # mean end-to-end length; not contour length
                length_sq = source.domain.ss_dist_ee_sq
        else:
            # Cross-over (or similar):
            length = HELIX_XOVER_DIST
            length_sq = HELIX_XOVER_DIST**2
    elif interaction == HYBRIDIZATION_INTERACTION:
        length = HELIX_WIDTH
        length_sq = length**2
        stiffness = 1
        # Regardless of what type of connection (stacking/hybridization/phosphate), if
        # source.domain != target.domain, then length is always 1 (or maybe it should be 2...)
        # Edit: No.
        # if length is None:
        #     length = 1
        # if length_sq is None:
        #     length_sq = 1
    elif interaction == STACKING_INTERACTION:
        length = HELIX_STACKING_DIST
        length_sq = length**2
        stiffness = 1
    else:
        raise ValueError("UNKNOWN INTERACTION: %r" % (interaction, ))
    return length, length_sq, stiffness





class GraphManager(object):
    """
    Class for approximating the system structure using a graph.
    The GraphManager uses four graphs to represent connections in the system:
        DomainEnds graph
        InterfaceNodes graph
        Domains graph
        Strands graph
    """

    def __init__(self, strands=None, ends5p3p_graph=None, domain_graph=None, strand_graph=None):
        ### System graphs ###
        # Not sure whether we should have system-level graphs or a graph for each complex.
        # Domain-level graphs can be created for each graph, in cases where we need cheap reduced graphs.
        # Otherwise, system-level graphs are probably cheaper and easier.
        # Why have graphs at multiple levels?
        # - Some operations, e.g. connected_strands() are much cheaper if we have strand-level graph.
        # MultiGraph vs Graph?
        # - Should be a multigraph to support stacking
        # - domains and ends can connect to the same partner via both backbone and stacking interactions
        # - hairpins without loops has hybridization to same domain as phosphate backbone
        #    (but we can disallow this and say we must always have a loop in hairpins).
        # - Not sure what the performance penalty of multi graphs are vs regular graphs?
        # - We could have a separate graphs without stacking, but not sure when those would be useful?
        self.ends5p3p_graph = ends5p3p_graph or nx.MultiGraph()
        self.domain_graph = domain_graph or nx.MultiGraph()
        self.strand_graph = strand_graph or nx.MultiGraph()
        self.interface_graph = InterfaceGraph()
        if strands:
            # TODO: Add strand-node attrs to strand-graph:
            # TODO: Add "Label" attr to all nodes (for Gephi)
            self.strand_graph.add_nodes_from(strands)
            for strand in strands:
                self.domain_graph.add_nodes_from(strand.nodes(data=True))
                self.domain_graph.add_edges_from(strand.edges(keys=True, data=True))
                self.ends5p3p_graph.add_nodes_from(strand.ends5p3p_graph.nodes(data=True))
                self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph.edges(keys=True, data=True))

                # self.interface_graph.add_nodes_from(
                #     end for end5p in strand.ends5p3p_graph.nodes() if isinstance(end5p, DomainEnd5p)
                #     for end in (end5p, end5p.bp_downstream) if end is not None)
                # self.interface_graph.add_nodes_from(
                #     InterfaceNode((end5p, end5p.bp_downstream) if end5p.bp_downstream else (end5p, ))
                #     for end5p in strand.ends5p3p_graph.nodes() if isinstance(end5p, DomainEnd5p))
                ## Edit: InterfaceNodes are now a persistent part of DomainEnds
        self.interface_graph.add_nodes_from(end.ifnode for end in self.ends5p3p_graph.nodes())
        domain_end_edges = self.ends5p3p_graph.edges(data=True, keys=True)
        # We must start from a graph without hybridization/stacking.
        # If we have hybridization/stacking, we have to perform node-delegation on those nodes
        # and we don't want to implement that right now...
        assert all(key == PHOSPHATEBACKBONE_INTERACTION for source, target, key, data in domain_end_edges)
        self.interface_graph.add_edges_from((source.ifnode, target.ifnode, data)
                                            for source, target, key, data in domain_end_edges)

        ### Other attributes: ###
        self.fnprefix = ""

        ### Loops: ###
        # We currently just use domain objects to track loops. Could probably be domain species.
        self.enable_loop_tracking = True


    def interfaces_shortest_path(self, ifnode1, ifnode2):
        """
        TODO: This should certainly be cached.
        """
        if isinstance(ifnode1, DomainEnd):
            ifnode1, ifnode2 = ifnode1.ifnode.top_delegate(), ifnode2.ifnode.top_delegate()
        return shortest_path(self.interface_graph, ifnode1, ifnode2)


    def group_interfaces_path_by_stiffness(self, path):
        """
        Returns an iterator of structural elements based on a interface-level path (list of InterfaceGraph nodes),
        where neighboring edges with the same stiffness has been grouped:
        path_edges = [(stiffness, [(length, length_sq, stiffness, source, target), ...]),
                      (stiffness, [(length, length_sq, stiffness, source, target), ...]),
                      (stiffness, [(length, length_sq, stiffness, source, target), ...]),
                      ...]
        That is, if the path is:
            4 nm ss-backbone + 3 nm ss-backbone + 5 nm ds+backbone + 1 nm ss-backbone + 2 nm ds-backbone + 8 nm ds-backbone,
        The result is:
            [(0, [(4 nm, 16 nm2, source, target), (3 nm, 9 nm2, source, target)]),
             (1, [(5 nm, 25 nm2, source, target)]),
             (0, [(1 nm, 1 nm2, source, target)]),
             (1, [(2 nm, 4 nm2, source, target), (8 nm, 64 nm2, source, target)])]

        Where stiffness=0 indicates a list of single-stranded edges,
        stiffness=1 indicates a hybridized duplex edge, and stiffness=3 is used for helix-bundles.
        """
        ## First create a generator of (length, length_sq, stiffness, source, target)
        ## Then use itertools.groupby to group elements by stiffness
        ## https://docs.python.org/3/library/itertools.html#itertools.groupby
        path_source_target_eattr = ((source, target, self.interface_graph[source][target])
                                    for source, target in zip(path, path[1:]))
        # We really should have all three lengths: len_contour, dist_ee_nm, and dist_ee_sq.
        path_tuples = ((edge_attrs['len_contour'], edge_attrs['dist_ee_sq'], edge_attrs['stiffness'], source, target)
                       for source, target, edge_attrs in path_source_target_eattr)
        stiffness_key = itemgetter(2) # Because stiffness = path_tup[2] for path_tup i path_tuples
        return itertools.groupby(path_tuples, key=stiffness_key)


    def interfaces_path_segments(self, path):
        """
        Arguments:
            :path:  A list of ifnodes describing the path.
        Return a list of (stiffness, [length, sum_length_squared]) tuples
        where each tuple represents a path segment with the same stiffness.
        This basically reduces the path to its constituent parts, which are usually all
        that is needed to calculate end-to-end probability distribution functions.
        """
        grouped_path_edges = self.group_interfaces_path_by_stiffness(path)
        # Trim away stiffness, source, target from (length, length_sq, stiffness, source, target) tuples:
        # zip(*[(1,1), (2, 4), (3, 9), (4, 16)]) = [(1, 2, 3, 4), (1, 4, 9, 16)]
        return [(stiffness, [sum(lengths) for lengths in zip(*[tup[:2] for tup in elem_tuples])])
                for stiffness, elem_tuples in grouped_path_edges]

    def grouped_path_edges_to_segments(self, grouped_path_edges):
        """
        Reference function, should be inlined.
        """
        return [(stiffness, [sum(lengths) for lengths in zip(*[tup[:2] for tup in elem_tuples])])
                for stiffness, elem_tuples in grouped_path_edges]


    def join_two_edge_groups(self, part_a, part_b, simulate_reaction, reaction_elems=(0, None)):
        """
        Join two edges groups.
            Assumes that the reaction is between the start of part_a and the end of part_b, i.e.:
                [part_a_last_segment, ..., part_a_first_segment] <::> [part_b_last_segment, ..., part_b_first_segment]

        Will join the first segment/group of part_a a with the last segment/group of part b:
        The join is done in-place by updating part_b. Make sure to make a copy if you want to retain the original!
        For stacking reactions, this is quite simple.
        For hybridization reactions, this requires a bit more thought..

        :part_a: and :part_b: are lists of segment edges groups with the form of:
            [(segment_stiffness, [segment_edge_attr_tuples]), ...],
        e.g.
            [(1, [(e1_length, e1_len_sq, e1_stiffness, e1_source, e1_target), (e2 attr tuple), ...])
             (0, [edge tuples for next segment with stiffness=0]),
             ...]
        """
        ## Adjust path segments as they would be *after* reaction.
        part_a_first_group = part_a[0]
        part_b_last_group = part_b[-1]
        if simulate_reaction is STACKING_INTERACTION:
            # Here we can assume that stacking *goes through* the duplexes:
            # TODO: Not sure this assertion is correct!
            # One case where it is not the case is the backbone edge connecting a stacked ifnode (self-loop)
            assert part_a[0][SEGMENT_STIFFNESS] > 0
            assert part_b[-1][SEGMENT_STIFFNESS] > 0
            # First join the first segment of a to the last segment of b: (using extend in-place)
            part_b[-1][1].extend(part_b[0][1])
            # Then append the rest of the segments from part a to part b:
            part_b.extend(part_a[1:])
        ## TODO: Determine this for hybridization reactions.
        elif simulate_reaction is HYBRIDIZATION_INTERACTION:
            assert part_a[0][SEGMENT_STIFFNESS] == 0
            assert part_b[-1][SEGMENT_STIFFNESS] == 0
            # Remove the domain edge from part_a and part_b:
            domain_edge_a = part_a[0][1].pop(0)
            domain_edge_b = part_b[-1][1].pop()
            ## part_a    d1end2    d1end1
            ##  <------o-o---------o
            ##            \         \       part_b
            ##             o---------o-------<-<-<--
            ##             d2end2    d2end1
            # length, len_sq, stiffness, source, target
            d1end1_ifnode = domain_edge_a[3]
            d1end2_ifnode = domain_edge_a[4]
            d2end1_ifnode = domain_edge_b[3]
            d2end2_ifnode = domain_edge_b[4]
            d1 = d1end1_ifnode.domain_end.domain
            d2 = d2end1_ifnode.domain_end.domain
            assert d1end1_ifnode != d1end2_ifnode != d2end1_ifnode != d2end2_ifnode
            assert d1end2_ifnode.domain_end.domain == d1
            assert d2end2_ifnode.domain_end.domain == d2
            assert d1end1_ifnode.domain_end.end == d2end2_ifnode.domain_end.end
            assert d2end1_ifnode.domain_end.end == d1end2_ifnode.domain_end.end
            # Use the ifnode guaranteed to be connected in interface_graph:
            duplex_tup = (d1.ds_len_contour, d1.ds_dist_ee_sq, 1, d2end1_ifnode, d1end2_ifnode)
            # Add new 1-element segment group to part_b:
            part_b.append((1, [duplex_tup]))
            # Append modified part_a to part_b:
            part_b.extend(part_a[1:])
            # NB: the newly hybridized duplex is still connected by flexible phosphate backbone links,
            # since it has not had any oppertunity to stack, so the surrounding segments are still flexible,
            # and I don't have to worry about joining the stiff duplex edge with neighboring segments.
        else:
            raise NotImplementedError("Only STACKING_INTERACTION and HYBRIDIZATION_INTERACTION currently supported.")
        return part_b



    def join_edges_groups(self, edge_groups, simulate_reaction=None, reaction_elems=(0, None)):
        """
        Returns
            closed_loop_groups:  A list of segments making up the loop after the loop has been closed by the given reaction.
                                 Each segment is a list of edges with the same stiffness.
        Note that edges_groups are slightly different from path segments:
            segments are tuples of    (stiffness, [segment_contour_length, segment_dist_ee_sq])
            path groups are tuples of (stiffness, [(length, length_sq, stiffness, source, target), ...])
            i.e. each path group has a full list of path edges for each stiffness group,
            while each segment only has information about the segment's contour length and sum-of-squared-links.
        Raises:
            IndexError if the path is already only a single, stiff duplex.
        Arguments:
        :loop_path: List of ifnodes describing the full loop path.
            (InterfaceGraph is not a multigraph so adj[source][target] is unambiguous.)
        :simulate_reaction: Assume that the first and last path elements have undergone this type of reaction.
            Only valid argument is None or STACKING_INTERACTION.
        :reaction_elems: If the elements undergoing the reaction described by simulate_reaction,
            add them here. (NOT CURRENTLY SUPPORTED)
        Alternative method names:
            join_path_groups, join_path_edge(s)_groups
            join_edge_segment_groups, join_segment_edge_groups
            close_loop_edges_groups
        """
        # edge_groups = self.group_interfaces_path_by_stiffness(path)
        # Grouped path edges is a list of:
        # path_edges_groups = [(stiffness, [(length, length_sq, stiffness, source, target), ...]), ...]

        ## Adjust path segments as they would be *after* reaction.
        if simulate_reaction is STACKING_INTERACTION:
            # A duplex end can either be connected through its own duplex (the domain backbone),
            # or throught the backbone connection to the next domain.
            # We have both cases for both duplex ends.
            first_group = edge_groups[0]
            last_group = edge_groups[-1]
            if last_group[SEGMENT_STIFFNESS] > 0 and first_group[SEGMENT_STIFFNESS] > 0:
                if len(edge_groups) == 1:
                    # The path may or may not go through the duplexed double-helix for the ends being stacked:
                    if first_group[SEGMENT_STIFFNESS] > 0:
                        # Trying to stack opposite ends of duplexed double-helix:
                        raise IndexError("Path consists of a single, stiff element, cannot circularize")
                else: # if len(edge_groups) > 1:
                    ## Both duplex ends are connected via their own (stiff) duplex;
                    ## we should put their segments "together" as they would be when stacked:
                    # Remove the last segment and apply it onto the first.
                    # Even though the group is a tuple, the edge attrs in group_tup[1] is still a mutable list, i.e:
                    #       first_group = (first_group_stiffness, [list of edge tuples]
                    first_group[1] += last_group[1]   # Extend first group edges with edges from last group
                    edge_groups.pop()  # Remove the last segment and append it to the first
                    # Edit: You should not sum the square, but re-calculate the square based on actual length:
                    # first_group[1][1] = first_group[1][0]**2
                    # We re-calculate sum-of-squares for stiff groups/segments later, not here.
            # else:
            # Ends are connected by one or more flexible segments, either a single phosphate-backbone link,
            # or through one or more unhybridized domains.
        ## TODO: Consider hybridization reactions, c.f. join_two_edge_groups(...) method.
        return edge_groups


    def closed_loop_segments(self, loop_path, simulate_reaction=None, reaction_elems=(0, None)):
        """
        Like closed_loop_path_groups, but returns segments rather than path groups.
        It does this a little different, though:
            First get segments (not edge groups), then:
            * if STACKING reaction and first and last segments are both stiff, merge first and last segment.
        This is the order used by the original calculate_loop_activity.
        """
        # path_segments is a list of (stiffness, [length, sum_length_squared]) tuples
        # where stiffness is (0 for ssDNA, 1 for dsDNA and 2 for helix-bundles).
        path_segments = self.interfaces_path_segments(loop_path)

        ## Adjust path segments as they would be *after* reaction.
        if simulate_reaction is STACKING_INTERACTION:
            # A duplex end can either be connected through its own duplex (the domain backbone),
            # or throught the backbone connection to the next domain.
            # We have both cases for both duplex ends.
            first_segment = path_segments[0]
            if len(path_segments) > 1:
                # Remove the last segment and apply it onto the first.
                last_segment = path_segments[-1]
                # SEGMENT_STIFFNESS, SEGMENT_LENGTHS, SEGMENT_CONTOUR_LENGTH, SEGMENT_LENGTH_SQUARED = 0, 1, 0, 1
                if last_segment[SEGMENT_STIFFNESS] > 0 and first_segment[SEGMENT_STIFFNESS] > 0:
                    ## Both duplex ends are connected via their own (stiff) duplex;
                    ## we should put their segments "together" as they would be when stacked:
                    # stiffness = first_segment[0] if first_segment[0] >= last_segment[0] else last_segment[0]
                    # Even though segment is a tuple, the segment lengths in segment[1] is still a mutable list:
                    #path_segments[0] = (first_segment[0], []
                    path_segments.pop()  # Remove the last segment and append it to the first
                    first_segment[1][0] += last_segment[1][0]   # Add segment contour lengths
                    # first_segment[1][1] += last_segment[1][1]   # Add segment squared lengths
                    # Edit: You should not sum the square, but re-calculate the square based on actual length:
                    first_segment[1][1] = first_segment[1][0]**2
                # else:
                    # If the duplexes are not connected via their own stiff/duplexed domains, then downstream
                    # calculation using LRE/SRE should produce correct result...
            else:
                # Only a single path segment:
                if first_segment[SEGMENT_STIFFNESS] > 0:
                    # We have a single, fully-stacked segment; this cannot ever stack back upon it self.
                    # return 0
                    raise IndexError("Path consists of a single, stiff element, cannot circularize")
                else:
                    # A single, flexible connection...
                    # This must be just a single phosphate connection between upstream 3' and downstream 5' ends, right?
                    # Nah, it's the same for duplex hybridization of duplex domains separated by T-loop.
                    # Thus, we also need to check for loops, etc... Let's just let the general routine take its course.
                    pass
                    # mean_sq_ee_dist = first_segment[1][1]
                    # effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2)
                    # activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
                    # return 2.0
        return path_segments



    def calculate_loop_activity(self, loop_path, simulate_reaction=None, reaction_elems=(0, None)):
        """
        Returns
            activity - the "probability" that the loop ends are within a critical radius (relative to a 1 M solution).
        Arguments:
        :loop_path: List of ifnodes describing the full loop path.
            (InterfaceGraph is not a multigraph so adj[source][target] is unambiguous.)
        :simulate_reaction: Assume that the first and last path elements have undergone this type of reaction.
            Only valid argument is None or STACKING_INTERACTION.
        :reaction_elems: If the elements undergoing the reaction described by simulate_reaction,
            add them here. (NOT CURRENTLY SUPPORTED)
        """
        # path_segments is a list of (stiffness, [length, sum_length_squared]) tuples
        # where stiffness is (0 for ssDNA, 1 for dsDNA and 2 for helix-bundles).
        try:
            path_segments = self.closed_loop_segments(loop_path, simulate_reaction, reaction_elems)
        except IndexError:
            return 0

        ## Find LRE and SRE parts
        ## TODO: Replace "element" by "segment", which better describes that a segment is composed of one or more
        ## edges with same stiffness.
        try:
            # LRE: "Longest Rigid Element"; "SRE": "Sum of Remaining Elements".
            _, LRE_len, LRE_len_sq, LRE_idx = max((stiffness, elem_length, elem_len_sq, i)
                                                  for i, (stiffness, (elem_length, elem_len_sq))
                                                  in enumerate(path_segments)
                                                  if stiffness > 0)
            # No need to check if LRE_len <= HELIX_WIDTH - InterfaceGraph only includes backbone connections.
        except ValueError:
            # No stiff elements: (summarizing segment_len_sq is OK for non-stiff, single-stranded segments)
            LRE_len, LRE_len_sq = 0, 0
            SRE_len = sum(elem_length for (stiffness, (elem_length, elem_len_sq)) in path_segments if stiffness == 0)
            SRE_len_sq = sum(elem_len_sq for (stiffness, (elem_length, elem_len_sq)) in path_segments if stiffness == 0)
        else:
            # Exclude LRE when calculating SRE length:
            # Don't sum segment_len_sq for stiff segments; You must, instead, square the total segment_length.
            if len(path_segments) > 1:
                SRE_lengths_and_sq = [(elem_length, elem_length**2 if stiffness > 0 else elem_len_sq)
                                      for sub_path in (path_segments[:LRE_idx], path_segments[LRE_idx+1:])
                                      for stiffness, (elem_length, elem_len_sq) in sub_path]
                SRE_len, SRE_len_sq = [sum(vals) for vals in zip(*SRE_lengths_and_sq)]
            else:
                # SRE_lengths, SRE_sq_lengths = [], []
                SRE_len, SRE_len_sq = 0, 0

        ## Minor corrections to SRE_len and SRE_len_sq:
        # SRE_len -= 0.25
        # SRE_len_sq -= 0.35

        ## Check if the contour-length is long enough for the elements to "touch":
        if LRE_len > SRE_len:
            # The domains cannot reach each other. (Unless we allow helical bending which is not yet implemented)
            return 0

        # If LRE is significant, then we cannot use the simple end-end-distance < rc approximation. Instead, we must
        # consider the end-to-end PDF of the SRE part at radius r=LRE_len and within critical volume v_c.
        # Can this be approximated by a simple LRE factor?
        # LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq)) if LRE_len_sq > 0 else 1  - Uh, yes.
        # But only relevant for LRE_len_sq > SRE_len_sq/4:

        if LRE_len_sq > 0 and LRE_len_sq*4 > SRE_len_sq:
            # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
            # Use end-to-end distance PDF for the of the SRE segments using LRE as radius
            activity = (3/(2*math.pi*SRE_len_sq))**gamma * math.exp(-3*LRE_len_sq/(2*SRE_len_sq)) * AVOGADRO_VOLUME_NM3
        else:
            mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm
            activity = (3/(2*math.pi*mean_sq_ee_dist))**gamma * AVOGADRO_VOLUME_NM3

        ## Calculate mean end-to-end squared distance between the two domains, aka E_r_sq:
        # We already have the squared length, Nᵢbᵢ², so we just need to sum LRE and SRE:
        ## EDIT, no, wait - if we have stacked, stiff double-helices/helix-bundles, then we CANNOT
        ## just add the squared lengths together. We have to square the entire element length.
        ## We can do it for ss-helices, since the "squared length" is already just N_nt*ss_rise_per_nt*ss_kuhn_length,
        ## but that isn't valid for fully-stacked continuous segments of ds-helices.
        ## Done: Re-calculate SRE_len_sq based on the summarized segment lengths, not the path edges.

        # There shouldn't be any need to test for if mean_sq_ee_dist <= HELIX_XOVER_DIST**2 in the InterfaceGraph

        ## TODO: (Optimization) Factor out and pre-calculate (3/(2*math.pi))**gamma * AVOGADRO_VOLUME_NM3
        ##       such that activity = PREFACTOR * mean_sq_ee_dist^(-gamma)

        ## Note: "effective_volume" is just an informal term for P_loop(rc) / P_v0(rc) x AVOGADRO_VOLUME_NM3
        # effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**gamma
        # activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
        # if gamma_corr > 1:
        #     activity = activity**gamma_corr

        return activity




    def loop_formation_activity(self, elem1, elem2, reaction_type):
        r"""
        Returns
            :intracomplex_activity:
        between domain1 and domain2, so that
            c_j = k_j * intracomplex_activity
        The intracomplex activity is basically just:
            activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
        where NA is Avogadro's constant, 6.022e23/mol.

        """
        activity, loop_info = self.loop_formation_effects(elem1, elem2, reaction_type)
        return activity


    def loop_breakage_effects(self, elem1, elem2, reaction_type):
        r"""
        Returns tuple with:
            loop_info
        Where
          * loop_info is a dict describing the loop(s) that will be changed (and, perhaps, split)
            by deleting the connection between elem1 and elem2.
        Arguments:
            :elem1:, :elem2:, :reaction_type: are the same as for intercomplex_activity().

        """
        if reaction_type is HYBRIDIZATION_INTERACTION:
            domain1, domain2 = elem1, elem2
            cmplx = domain1.strand.complex
            d1end5p, d2end5p = domain1.end5p, domain2.end5p
            d1end3p, d2end3p = domain1.end3p, domain2.end3p
            # path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
            # Domains are currently NOT hybridized and thus also NOT stacked. There should be NO delegation!
            assert d1end5p.ifnode.top_delegate() == d1end5p.ifnode
            assert d1end3p.ifnode.top_delegate() == d1end3p.ifnode
            assert d2end5p.ifnode.top_delegate() == d2end5p.ifnode
            assert d2end3p.ifnode.top_delegate() == d2end3p.ifnode
            reactant_nodes = (d1end5p.ifnode, d1end3p.ifnode, d2end5p.ifnode, d2end3p.ifnode)
            ## TODO: (OPTIMIZATION) Instead of finding the shortest path, try to
            ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP.
            ## If they are, then the existing loop should give the shortest path.
            ## (For now, just assert that IF the two loops are part of an existing, non-obsolete loop,
            ## then the shortest path is in that loop.)
            elem1_ifnode, elem2_ifnode = d1end5p.ifnode, d2end5p.ifnode
            # new_loop_path = shortest_path(self.interface_graph, d1end5p.ifnode, d2end5p.ifnode)
            # slice_start = 1 if new_loop_path[1] == d1end3p.ifnode else 0
            # slice_end = -1 if new_loop_path[-2] == d2end3p.ifnode else None
            # pdb.set_trace()
            # if slice_start == 1 or slice_end is not None:
            #     new_loop_path = new_loop_path[slice_start:slice_end]
        elif reaction_type is STACKING_INTERACTION:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
            cmplx = h1end5p.domain.strand.complex
            #  DUPLEX/DOUBLE-HELIX 1 (d: domain, h: Single helix, dh: double-helix / duplex-ends)
            #             dh1, d1         dh2, d4
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            #             dh1, d2         dh2, d3
            elem1_ifnode, elem2_ifnode = h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate()
            # new_loop_path = shortest_path(self.interface_graph, dh1_delegate, dh2_delegate)
        else:
            assert isinstance(elem1, DomainEnd) and isinstance(elem2, DomainEnd)
            d1end5p, d2end5p = elem1, elem2
            cmplx = d1end5p.domain.strand.complex
            elem1_ifnode, elem2_ifnode = elem1.ifnode.top_delegate(), elem2.ifnode.top_delegate()

        # We assume elem1 and elem2 has been split apart, so they should no longer be represented by a single ifnode.
        assert elem1_ifnode != elem2_ifnode

        # We assume connection is already broken between elem1 and elem2.
        if len(cmplx.loops) == 0:
            return None

        # TODO, FIX: While stacking/unstacking always creates/deletes a loop, that is not necessarily the case
        # TODO, FIX: for HYBRIDIZATION reactions.
        # TODO, FIX: Actually that is not the problem. The problem is that loops can be affected in different ways.
        # TODO, FIX: We assumed that the "shared path" was broken. But that may not be true for duplex dehybridization.
        # TODO, FIX: This should be fixed when looping over affected loops.

        # 1. Find all affected loops:
        # cmplx.ifnode_loopids_index[elem1_ifnode] => {set of loopids for loops touching ifnode}
        # Wait: If elem1 and elem2 are currently hybridized or stacked, wouldn't they have the
        # same top_delegate and thus share the exact same loops?
        # No: At this point, the breaking reaction has occured, so the two elem no longer has the same ifnode.
        # However, since we haven't updated cmplx.ifnode_loopids_index index yet, one of them should be empty
        # and the other should be full.
        if reaction_type is HYBRIDIZATION_INTERACTION:
            d1end5p_loops = cmplx.ifnode_loopids_index[d1end5p.ifnode]
            d1end3p_loops = cmplx.ifnode_loopids_index[d1end3p.ifnode]
            d2end5p_loops = cmplx.ifnode_loopids_index[d2end5p.ifnode]
            d2end3p_loops = cmplx.ifnode_loopids_index[d2end3p.ifnode]
            if len(d1end5p_loops) > 0:
                assert len(d2end3p_loops) == 0
                end1_loops = d1end5p_loops
                end1_former_delegate = d1end5p.ifnode
                end1_new_delegate = d2end3p.ifnode
            else:
                # assert len(d2end3p_loops) > 0  # There is no guarantee we have any loops
                assert len(d1end5p_loops) == 0
                end1_loops = d2end3p_loops
                end1_former_delegate = d2end3p.ifnode
                end1_new_delegate = d1end5p.ifnode
            if len(d2end5p_loops) > 0:
                assert len(d1end3p_loops) == 0
                end2_loops = d2end5p_loops
                end2_former_delegate = d2end5p.ifnode
                end2_new_delegate = d1end3p.ifnode
            else:
                # assert len(d1end3p_loops) > 0   # There is no guarantee we have any loops
                assert len(d2end5p_loops) == 0
                end2_loops = d1end3p_loops
                end2_former_delegate = d1end3p.ifnode
                end2_new_delegate = d2end5p.ifnode
            previous_top_delegates = ends_former_top_delegates = (end1_former_delegate, end2_former_delegate)
            newly_undelegated_ifnodes = ends_new_delegates = (end1_new_delegate, end2_new_delegate)
            affected_loopids = set.union(*[cmplx.ifnode_loopids_index[ifnode] for ifnode in reactant_nodes])
        else:
            # For STACKING (and probably backbone ligation/nicking) we can do it simpler:
            elem1_loops = cmplx.ifnode_loopids_index[elem1_ifnode]
            elem2_loops = cmplx.ifnode_loopids_index[elem2_ifnode]
            if len(elem1_loops) > 0:
                previous_top_delegate = elem1_ifnode    # Or maybe there just wasn't any loops?
                previous_top_delegates = [elem1_ifnode]    # Or maybe there just wasn't any loops?
                newly_undelegated_ifnode = elem2_ifnode
                newly_undelegated_ifnodes = [elem2_ifnode]
                assert len(elem2_loops) == 0
                affected_loopids = elem1_loops
            else:
                previous_top_delegate = elem2_ifnode
                newly_undelegated_ifnode = elem1_ifnode
                previous_top_delegates = [elem2_ifnode]
                newly_undelegated_ifnodes = [elem1_ifnode]
                assert len(elem1_loops) == 0
                affected_loopids = elem2_loops
        if len(affected_loopids) == 0:
            return None
        # Make sure not to modify affected_loopids - we use that later! (This is just for debugging)
        # unprocesses_affected_loopids = affected_loopids.copy()
        intact_path_loops = set()
        broken_path_loops = set()

        # 2. Update all loops to include the ifnode split.
        # Edit: TODO: Should this be described genericly? Otherwise, maybe move this function to Complex
        #       to indicate that this is actually modifying the loops, not just "calculating loop effects" !
        #       OTOH, I'm making use of self.interface_graph which is not available within Complex.
        # We need to have the updated ifnodes in the paths because we query the edge length in
        # find_alternative_shortest_path using self.interface_graph[source][target]['len_contour']
        # changed_loopids_lst = sorted(affected_loopids, key=lambda loopid: cmplx.loops[loopid]['activity'], reverse=True)
        # changed_loopids_deque = deque(changed_loopids_lst)
        affected_loop_tuples = sorted([(cmplx.loops[loopid]['activity'],
                                        cmplx.loops[loopid]['loop_hash'],
                                        cmplx.loops[loopid]['path_spec'],
                                        cmplx.loops[loopid]['path'],
                                        loopid)
                                       for loopid in affected_loopids], reverse=True)
        # affected_loops = [cmplx.loops[loopid] for loopid in changed_loopids_lst]
        # for loopid, affected_loop in zip(changed_loopids_lst, affected_loops):
        print("Pre-processing affected loop paths:")
        print("\n".join(str(tup) for tup in affected_loop_tuples))
        for activity, loop_old_hash, path_spec, path, loopid in affected_loop_tuples:
            # path = affected_loop['path']
            cmplx.loops[loopid]['debug_info']['path_before_breakage_preprocessing'] = tuple(path)  # TODO: Remove debug
            path_node_index = {ifnode: idx for idx, ifnode in enumerate(path)}
            # TODO, FIX: For duplex dehybridization, we should perhaps pre-process the path such that if
            # the previous_top_delegate ifnode is no longer connected to the next ifnode in the path,
            # then check if the newly_delegated_ifnode is and if so, swap them out.
            # Then, if both duplex ends are still connected to the rest of the path, then
            # no further path changes are necessary. You can add effects caused by changes in contour length,
            # so perhaps re-calculate loop activity, but you don't need to change the path ifnodes.
            if reaction_type is HYBRIDIZATION_INTERACTION:
                ends_path_idxs = [path_node_index.get(ifnode) for ifnode in ends_former_top_delegates]
                end_on_path = [idx is not None for idx in ends_path_idxs]
                assert any(end_on_path) # path should go through duplex somehow.
                if all(end_on_path):
                    # Both duplex ends are part of the path. Figure out which is up/down stream in the path's direction
                    assert abs(ends_path_idxs[0] - ends_path_idxs[1]) == 1
                    # path_is_parallel = ends_path_idxs[0] < ends_path_idxs[1]
                    end_downstream_idx, end_upstream_idx = ((ends_path_idxs[0], ends_path_idxs[1])
                                                            if ends_path_idxs[0] < ends_path_idxs[1] else
                                                            (ends_path_idxs[1], ends_path_idxs[0]))
                    neighbor_idxs = (end_downstream_idx-1, # OK, even if idx=0
                                     (end_upstream_idx+1) % len(path))  # modulu list len to wrap around
                    neighbors = [path[idx] for idx in neighbor_idxs]
                    for (neighbor, top_delegate, new_delegate, idx) in zip(
                            neighbors, ends_former_top_delegates, ends_new_delegates, ends_path_idxs):
                        if top_delegate not in self.interface_graph.adj[neighbor]:
                            # I can't think of any situations where neither delegates are connected in the path:
                            assert new_delegate in self.interface_graph.adj[neighbor]
                            # Swap ifnodes:
                            print("Swapping old top delegate %s with new delegate %s on loop %s path %s" %
                                  (top_delegate, new_delegate, loopid, path))
                            path[idx] = new_delegate
                        else:
                            pdb.set_trace()
                    # Check if the two path-connected delegates are connected to each other.
                    # TODO: This check belongs in the "find-new-loop-path-maybe" loop below!
                    # Path is updated to use the best-connected ifnode.
                    path_connected_ifnodes = [path[idx] for idx in ends_path_idxs]
                    if path_connected_ifnodes[0] in self.interface_graph.adj[path_connected_ifnodes[1]]:
                        print("Loop %s: Path-connected Ifnodes %s are still connected to each other" % (
                            loopid, path_connected_ifnodes))
                        intact_path_loops.add(loopid)
                    else:
                        broken_path_loops.add(loopid)
                else:
                    # If not all(ends_on_path), then I don't think we need to do any pre-processing.
                    # The loop path must surely change, it is just a question of how.
                    # (It could only happen if the end is stacked, but then it couldn't dehybridize..
                    # Also, we don't allow duplexes to be directly
                    ifnode_idx, former_top_delegate = next((idx, ifnode) for idx, ifnode in zip(
                        ends_path_idxs, ends_former_top_delegates) if idx is not None)
                    neighbor_before, neighbor_after = path[ifnode_idx-1], path[(ifnode_idx+1) % len(path)]
                    if former_top_delegate in self.interface_graph.adj[neighbor_before]:
                        assert neighbor_before in self.interface_graph.adj[former_top_delegate]
                        assert neighbor_after not in self.interface_graph.adj[former_top_delegate]
                        assert former_top_delegate not in self.interface_graph.adj[neighbor_after]
                    else:
                        assert neighbor_after in self.interface_graph.adj[former_top_delegate]
                        assert former_top_delegate in self.interface_graph.adj[neighbor_after]
                        assert neighbor_before not in self.interface_graph.adj[former_top_delegate]
                        assert former_top_delegate not in self.interface_graph.adj[neighbor_before]
                    broken_path_loops.add(loopid)


            else:
                # STACKING_INTERACTION
                # TODO: If unstacking backbone-connected domain-ends, then after ifnode-split path updating below,
                # the path is essentially intact. True, we will also find the path with find_alternative_shortest_path,
                # but that search isn't needed. Instead, just check if after processing the ifnode split, the
                # two ifnodes are (a) still connected to each other, and (b) still connected to the path.
                # if so, add loop to intact_path_loops.
                if len(path) == 1:
                    # Special case, backbone-connected stacked duplex ends;
                    path.append(newly_undelegated_ifnode)
                    print(("Inserted newly_undelegated_ifnode %s AFTER previous_top_delegate %s in loop %s path") % (
                        newly_undelegated_ifnode, previous_top_delegate, loopid))
                    broken_path_loops.add(loopid)
                    print(" - loop %s = %s marked as broken (length 1)." % (loopid, path))
                    continue
                ifnode_idx = path.index(previous_top_delegate)
                ifnode_before = path[ifnode_idx-1] # = (ifnode_idx % len(path))-1
                ifnode_after = path[(ifnode_idx+1) % len(path)]
                if len(path) == 2:
                    # For path of length 2, we should still be connected..  A--X --> A-B--X <=> A--X--B
                    # We can still short-circuit the "where to insert newly_undelegated_ifnode (B)" check, but
                    # we must still do the following "are A and B still connected or is loop broken" test below.
                    path.append(newly_undelegated_ifnode)
                    print(("Appending newly_undelegated_ifnode %s AFTER previous_top_delegate %s at end of loop %s."
                           " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
                elif previous_top_delegate in self.interface_graph.adj[ifnode_before]:
                    # former_top_delegate is still connected to ifnode_before
                    # newly_undelegated_ifnode should be inserted *after* former_top_delegate in the path.
                    # Not sure about these assertions, verify if you encounter AssertionErrors:
                    # One reason why it may be correct is that ifnode_idx is not at the start or end of the path,
                    # so the path must be 3 or more elements long.
                    assert ifnode_before in self.interface_graph.adj[previous_top_delegate]
                    assert previous_top_delegate not in self.interface_graph.adj[ifnode_after]
                    assert ifnode_after not in self.interface_graph.adj[previous_top_delegate]
                    assert newly_undelegated_ifnode in self.interface_graph.adj[ifnode_after]
                    assert newly_undelegated_ifnode not in self.interface_graph.adj[ifnode_before]
                    # TODO: Remove the excessive assertions above.
                    # Insert newly_undelegated_ifnode AFTER previous_top_delegate in the path at position idx+1:
                    if ifnode_idx+1 == len(path): # former_top_delegate at end of list; we can just append.
                        path.append(newly_undelegated_ifnode)
                        print(("Appending newly_undelegated_ifnode %s AFTER previous_top_delegate %s at end of loop %s."
                               " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
                    else:
                        path.insert(ifnode_idx+1, newly_undelegated_ifnode)
                        print(("Inserted newly_undelegated_ifnode %s AFTER previous_top_delegate %s in loop %s path at "\
                               "position %s+1") % (newly_undelegated_ifnode, previous_top_delegate, loopid, ifnode_idx))
                else:
                    # former_top_delegate is NO longer connected to ifnode_before (but should be to ifnode_after)
                    # newly_undelegated_ifnode should be connected to ifnode_before
                    # newly_undelegated_ifnode should be inserted *before* former_top_delegate in the path.
                    assert ifnode_before not in self.interface_graph.adj[previous_top_delegate]
                    assert previous_top_delegate in self.interface_graph.adj[ifnode_after]
                    assert ifnode_after in self.interface_graph.adj[previous_top_delegate]
                    assert newly_undelegated_ifnode not in self.interface_graph.adj[ifnode_after]
                    assert newly_undelegated_ifnode in self.interface_graph.adj[ifnode_before]
                    # TODO: Remove the excessive assertions above.
                    # Insert newly_undelegated_ifnode *before* previous_top_delegate in the path:
                    if ifnode_idx == 0:
                        # previous_top_delegate is at start of path; newly_undelegated_ifnode before = at end.
                        path.append(newly_undelegated_ifnode)
                        print(("Appending newly_undelegated_ifnode %s BEFORE previous_top_delegate %s at end of loop %s."
                               " path: %s ") % (newly_undelegated_ifnode, previous_top_delegate, loopid, path))
                    else:
                        path.insert(ifnode_idx, newly_undelegated_ifnode)
                        print(("Inserted newly_undelegated_ifnode %s BEFORE previous_top_delegate %s in loop %s path at "\
                               "position %s yielding the path: %s") % (
                            newly_undelegated_ifnode, previous_top_delegate, loopid, ifnode_idx, path))
                # Check if loop path is broken:
                if newly_undelegated_ifnode in self.interface_graph.adj[previous_top_delegate]:
                    # Loop path still valid, re-use:
                    assert previous_top_delegate in self.interface_graph.adj[newly_undelegated_ifnode]
                    intact_path_loops.add(loopid)
                    print(" - loop %s = %s marked as INTACT." % (loopid, path))
                else:
                    broken_path_loops.add(loopid)
                    print(" - loop %s = %s marked as BROKEN." % (loopid, path))




        # 3. Find the loop with the highest activity. This is the one we want to remove.
        # Concern: Do we need the path ifnode-split processing of this loop? Well, yeah, because we use it
        # as alternative path when piecing together new paths for broken loops.
        # Concern: Do we really need a deque? Or can we just use an iterable? Do we need to append new loops to the
        # list of affected loops? I don't think so.
        # del_loop_id = changed_loopids_deque.popleft()
        # del_loop_hash = cmplx.loops[del_loop_id]['loop_hash']
        affected_loop_tups_iter = iter(affected_loop_tuples)
        activity, del_loop_hash, del_loop_path_spec, del_loop_path, del_loop_id = next(affected_loop_tups_iter)

        # changed_loops = {del_loop_id: None}
        changed_loops_by_hash = {}
        changed_loops_hash_by_id = {}
        alternative_loop_paths = [del_loop_path]
        # reaction loop_effects directive = the collection of loop changes;
        # loop_update info = how to update each affected loop
        loop_effects = {
            'is_forming': False,
            'del_loop_id': del_loop_id,
            'del_loop_hash': del_loop_hash,
            'del_loop_path_spec': del_loop_path_spec, # for debugging
            'del_loop_activity': activity,
            # 'changed_loops': changed_loops,
            'changed_loops_by_hash': changed_loops_by_hash,
            'newly_undelegated_ifnodes_specs': {ifnode.state_fingerprint() for ifnode in newly_undelegated_ifnodes},
            'previous_top_delegates_specs': {ifnode.state_fingerprint() for ifnode in previous_top_delegates},
            # It is nice to permanently save complex state fingerprints to ensure that the complex
            # that tries to apply this loop_effects change directive has the same state
            'complex_loop_ensemble_fingerprint': cmplx.loop_ensemble_fingerprint,
            'complex_state_fingerprint': cmplx.state_fingerprint(),
            'directive_debug_info': {
                # Instance-specific debug info and other unused values
                'changed_loops_hash_by_id': changed_loops_hash_by_id,
                # 'previous_top_delegate': previous_top_delegate,         # Do not use, for debugging only
                # 'newly_undelegated_ifnode': newly_undelegated_ifnode,   # Do not use, for debugging only
                # 'previous_top_delegate_spec': previous_top_delegate.state_fingerprint(),
            },
            'description': "loop_breakage_effects directive",
        }

        # unprocesses_affected_loopids.remove(loop1_id) # changed_loopids

        # System:   Λ₀           .------.
        #            .-------,--´------. \
        #           /  Λa   /          \  \ c
        #         a \   eᵤᵥ/:/   Λb    /   \  Λc
        #            \      /       b /     \
        #             `-.--´---------´      |
        #                `-----------------´
        # We have three affected loops: Λa (aᵢ+e₁ᵢ), Λb (bᵢ+e₂ᵢ), and Λc (cᵢ+e₃ᵢ).
        # Note that the edge designation is a little weird, since it is basically a matrix:
        # Λa vs Λb:  Λa=a₂+e₂₁, Λb=b₁+e₁₂,  e₁₂==e₂₁
        # Λa vs Λc:  Λa=a₃+e₃₁, Λc=c₁+e₁₃,  e₁₃==e₃₁
        # Λb vs Λc:  Λb=b₃+e₃₂, Λc=c₂+e₂₃,  e₂₃==e₃₂

        # That is, the "shared subpath" depends on what loops we are comparing.
        # If we assume that two loops will only ever have one shared sub-path,
        # then any shared path should be eᵤᵥ, the path that is no longer usable.

        # All loops goes through eᵤᵥ which is broken apart, so all loops are affected.
        # We must: (1) Delete one loop, and (2) find new paths for the remaining loops.
        # When finding new loops, care must be taken to make the new loops (a) unique and
        # (b) use the shortest-path available.
        # In the case above:
        # affected_loopids = {Λa, Λb, Λc}
        # del_loop = Λa
        # changed_loopids_deque = (Λb, Λc)
        # alternative_loop_paths = [Λa]
        # Expected (if both Λb and Λc opts to use sub-path a for the new shortest path)
        #   Λb:  b₁+e₁₂ --> b₁+a₂
        #   Λc:  c₁+e₁₃ --> a₃+c₁
        # while loop 1:
        #   # loop1 is the loop before
        #   loop1 = Λb
        #   path2_


        # Propagate the change iteratively same way as you determine loop formation efffects:
        # while changed_loopids_deque:
        #     # For the first round, the loop considered is Λ₁ with loop_id None.
        #     # It is not appropriate to use "loop0" to describe this; you can use "loop1" or "new_loop"
        #     loop1_id = changed_loopids_deque.popleft()
        for old_activity, loop_old_hash, loop_old_path_spec, path, loopid in affected_loop_tups_iter:
            # looping over iter because we don't want to include the first to-be-deleted loop.
            if loopid in intact_path_loops:
                # Just need to update loop's activity using the updated loop path:
                path_is_broken = False
                new_loop_path = path # We can use old path
                a2 = self.calculate_loop_activity(path, simulate_reaction=None)
            else:
                # loop1_info = cmplx.loops[loopid]
                # loop_old_hash = loop1_info['loop_hash']
                # path = loop1_info['path']
                # loop0_path is the alternative_loop_path used to create the new_loop_path.
                assert loopid in broken_path_loops
                path_is_broken = True
                new_path_length, new_loop_path, e1, e2, e3, loop0_path = self.find_alternative_shortest_path(
                    path, alternative_loop_paths)
                a2 = self.calculate_loop_activity(new_loop_path, simulate_reaction=None)
            # [ifnode.state_fingerprint() for ifnode in new_loop_path]
            new_loop_path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
            new_loop_hash = cmplx.calculate_loop_hash(new_loop_path_spec)
            # new_loop_path_spec = tuple(new_loop_path_spec) # TODO: Remove cast
            loop_change_factor = a2/old_activity
            loop_update_info = {
                'loop_hash': new_loop_hash,
                'path_spec': tuple(new_loop_path_spec), # We need the path-spec to re-create re-built loop
                'activity': a2,
                'dS': ln(a2), # loop_entropy in units of R
                # Values for the original loop (mostly for debugging)
                # 'a0': old_activity,
                # debug info for the last loop change:
                'loop_update_debug_info': {
                    'path_is_broken': path_is_broken,
                    'description': 'loop_breakage_effects() %s loop' % ("rebuilt" if path_is_broken else "modified"),
                    'new_path_tuple': tuple(new_loop_path), # Only for debugging; the result is cached and uses between complexes
                    'new_path_spec': tuple(new_loop_path_spec),
                    'loop_change_factor': loop_change_factor, # = a2/a0
                    'old_activity': old_activity,
                    'old_loop_hash': loop_old_hash, # so we can relate this to the original/parent loop.
                    'old_path_spec': loop_old_path_spec,
                    'source_loop_id': loopid,
                },
            }
            # changed_loops[loopid] = loop_update_info
            changed_loops_by_hash[loop_old_hash] = loop_update_info
            changed_loops_hash_by_id[loopid] = loop_old_hash
            # unprocesses_affected_loopids.remove(loop1_id)
        # assert len(unprocesses_affected_loopids) == 0  # TODO: Remove debug check
        return loop_effects


    def find_alternative_shortest_path(self, loop1_path, alternative_loop_paths):
        r"""
        Again,
            loop1 is the "old" loop that no longer works because it is broken.
            alternative_loops are other loop paths that are also broken, but should be usable...
        The hypothesis is currently that we should be able to piece together a new shortest path
        using the previous (broken-up and processed) paths.
        ## Consider the system:   .------.
        ##            .-------,--´------. \
        ##           /  Λa   /          \  \ c = e1
        ##         a \   eᵤᵥ/:/   Λb    /   \  Λc = loop1
        ##            \      /       b /     \
        ##             `-.--´---------´      |
        ##                `-----------------´
        We are considering loop1_path = Λc = c+eᵤᵥ
        against alternative_loop_paths = [Λa, Λb]     "loop0"
        for each loop0 in alternative_loop_paths
        we find the shared part eᵤᵥ aka e3 and the not-shared-but-possibly-new-shorter sub-path e2.
        Thus, for loop1_path=Λc and alternative_loop_path=Λa, we have the following loop aliases:
            e3 = eᵤᵥ = e₁₃ = e₃₁
            e2 = a
            e1 = c
            loop2_path = e1 + e2 = c + a
        # (Yes, there was lots of changes to the nomenclature, but I think it is fixed now:)
        # EDIT EDIT EDIT: e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.

        # Note: When forming, e3 is the "possibly shorter" path and e1 the shared path part;
        # when breaking, e3 is the "broken part" (shared) and e2 is the "possibly shorter" part.
        Previous considerations:
            # Or should we make it the same: e3 is the "possibly shorter" path, e1 is the shared (but broken) path.
            # e3 is "possibly shorter", e1 is the shared part.
            # edit edit: e1 is the shared part, e3 is the "unique to loop1" part, e2 is the "unique to loop0" part.
            # edit edit edit: e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.

        Discussion: Where do the "intersection nodes" go?
        - in find_maybe_shorter_e1_e3_subpaths() we added the intersection nodes to the two candidate sub-paths
            because we wanted to compare the two candidate sub-paths directly
        We don't have that concern here: we are comparing the full loop path length, since the different
        alternatives have different shared/overlapping sub-paths.


        """
        alternatives = []
        loop1_nodes = set(loop1_path)
        for loop0_path in alternative_loop_paths:
            loop0_nodes = set(loop0_path)
            e2_nodes = loop0_nodes - loop1_nodes # set(e3)
            e3_nodes = loop0_nodes & loop1_nodes # shared nodes
            # This is probably an easier and more reliable way:
            loop1_e3 = [ifnode for ifnode in loop1_path if ifnode in e3_nodes]
            loop0_e3 = [ifnode for ifnode in loop0_path if ifnode in e3_nodes]
            if loop1_e3 == loop0_e3:
                e3_is_parallel = True
            else:
                e3_is_parallel = False
                assert loop1_e3 == loop0_e3[::-1]
            loop1_e1 = [ifnode for ifnode in loop1_path if ifnode not in e3_nodes]
            loop0_e2 = [ifnode for ifnode in loop0_path if ifnode not in e3_nodes]
            # intersection_nodes = (loop1_e3[0], loop1_e3[-1])
            assert len(loop1_e3) >= 2
            assert loop1_e3[0] != loop1_e3[-1]
            loop1_e3_idxs = (loop1_path.index(loop1_e3[0]), loop1_path.index(loop1_e3[-1]))
            loop1_e3_idxs = (loop1_path.index(loop0_e3[0]), loop1_path.index(loop0_e3[-1]))
            if e3_is_parallel:
                # We need to reverse loop0_e2 to make it fit.
                loop2_path = loop1_e1 + [loop1_e3[0]] + loop0_e2[::-1] + [loop1_e3[-1]]
                # It is generally faster to extend+append instead of concatenate, but no difference for small lists.
            else:
                loop2_path = loop1_e1 + [loop1_e3[0]] + loop0_e2 + [loop1_e3[-1]]
            path2_length = sum(self.interface_graph[source][target]['len_contour']
                               for source, target in zip(loop2_path, loop2_path[1:]))
            alternatives.append((path2_length, loop2_path, loop1_e1, loop0_e2, loop1_e3, loop0_path))


            #
            # ## Old groupby approach:
            # # e3 is the shared part, e1 is the "unique to loop1" part, e2 is the "unique to loop0" part.
            # e1_or_e3_sub_paths = [(is_shared, list(group_iter)) for is_shared, group_iter in
            #                       groupby(loop1_path, lambda ifnode: ifnode not in loop0_nodes)] # pylint: disable=W0640
            # assert len(e1_or_e3_sub_paths) <= 3
            # if len(e1_or_e3_sub_paths) == 1:
            #     # We only have one sub-path: make sure it is e1
            #     assert e1_or_e3_sub_paths[0][0] is True  # verify that e₁ is *not* on Λ₀
            #     e1 = e1_or_e3_sub_paths[0][1]
            #     e3a, e3b = [e1[0]], [e1[-1]]
            # elif len(e1_or_e3_sub_paths) == 2:
            #     if e1_or_e3_sub_paths[0][0] is True:
            #         # We have e3b > 0; e3a starts at e1 (has zero length)
            #         e1 = e1_or_e3_sub_paths[0][1]
            #         e3a, e3b = [e1[0]], [e1[-1]]
            #         e3b.extend(e1_or_e3_sub_paths[1][1])
            #     else:
            #         # We have e3a > 0; e3b starts and ends on e1 (has zero length)
            #         e3a, e1 = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1]
            #         e3a, e3b = [e1[0]], [e1[-1]]
            #         e3a.append(e1[0])
            # else:
            #     # We cannot know whether loop1_path starts on e1 or e3, but we can easily check:
            #     if e1_or_e3_sub_paths[0][0] is True: # True means "on e1, not on e3".
            #         # loop1_path starts on e1
            #         assert tuple(tup[0] for tup in e1_or_e3_sub_paths) == (True, False, True)
            #         e1a, e3, e1b = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1], e1_or_e3_sub_paths[2][1]
            #         e1 = e1b + e1a
            #         e3 = [e1[-1]] + e3 + [e1[0]] # e3 starts with the last node in e1 and ends with the first e1 node
            #     else:
            #         # loop1_path starts on e3:
            #         assert tuple(tup[0] for tup in e1_or_e3_sub_paths) == (False, True, False)
            #         e3a, e1, e3b = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1], e1_or_e3_sub_paths[2][1]
            #         e3 = [e1[-1]] + e3b + e3a + [e1[0]] # e3 starts with last node in e1 and end in the first e1 node
            #
            # ## Unlike when forming new loops, here we are not very interested in e1.
            # ## However, we *are* interested in e2 and the total loop path length.
            #
            # # Find e2:
            # # loop1 and loop0 can potentially go in opposite directions. We need to be aware of that:
            # e3_is_parallel = loop0_path.index(e3[0]) < loop0_path.index(e3[-1])
            # e2_or_e3 = [(is_shared, list(group_iter)) for is_shared, group_iter
            #             in groupby(loop0_path, lambda node: node in e2_nodes)]  # pylint: disable=W0640
            # # e2_subpaths = [sub_path_group[1] for sub_path_group in e2_or_e3 if sub_path_group[0]]
            # assert 1 <= len(e2_or_e3) <= 3
            # if len(e2_or_e3) == 3:
            #     if e2_or_e3[0][0] is True: # we start in e2
            #         e2a, e3s, e2b = e2_or_e3[0][1] + e2_or_e3[1][1] + e2_or_e3[2][1]
            #         e2 = e2b + e2a
            #     else:
            #         e3a_, e2, e3b_ = e2_or_e3[0][1] + e2_or_e3[1][1] + e2_or_e3[2][1]
            #         e3_ = e3b_ + e3a_  # Using underscore to distinguish from e3 calculated using loop1_path
            # elif len(e2_or_e3) == 2:
            #     if e2_or_e3[0][0] is True: # we start in e2
            #         e2, e3_ = e2_or_e3[0][1] + e2_or_e3[1][1]
            #     else:
            #         e3_, e2 = e2_or_e3[0][1] + e2_or_e3[1][1]
            # else:
            #     # I don't think this should happen since e3 includes the "intersection nodes" (as does e1).
            #     assert len(e3) == 0
            #     assert e2_or_e3[0][0] is True # we start in e2
            #     e2 = e2_or_e3[0][1]
            # # alternatively, you could have looked at whether loop0_path[0] and loop0_path[-1]
            # # are both in e2_nodes
            #
            # ## Determine how e2 and e1 should be combined to form a new shortest candidate:
            # # Note: e1 and e3 e3a/e3b share the "intersection nodes", but e2 does not.
            # # Again: e3 is the shared/broken sub-path, e1 is "unique to loop1" sub-path, e2 is "unique to loop0".
            # if e3_is_parallel:
            #     # we must reverse e2 when concatenating with e1:
            #     assert e2[-1] in self.interface_graph.adj[e1[-1]] # the last nodes are next to each other
            #     assert e2[0] in self.interface_graph.adj[e1[0]]   # and the first nodes are next to each other
            #     loop2_path = e1 + e2[::-1]                        # so definitely need to reverse
            # else:
            #     assert e2[0] in self.interface_graph.adj[e1[-1]]
            #     assert e2[-1] in self.interface_graph.adj[e1[0]]
            #     loop2_path = e1 + e2
            # # loop2_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
            # #                 in self.group_interfaces_path_by_stiffness(loop2_path)]
            # # We only evaluate shortest paths based on length, we do not calculate activity or anything like that:
            # # path_source_target_eattr = ((source, target, self.interface_graph[source][target])
            # #                             for source, target in zip(loop2_path, loop2_path[1:]))
            # # # We really should have all three lengths: len_contour, dist_ee_nm, and dist_ee_sq.
            # # path_tuples = ((edge_attrs['len_contour'], edge_attrs['dist_ee_sq'], edge_attrs['stiffness'], source, target)
            # #                for source, target, edge_attrs in path_source_target_eattr)
            # path2_length = sum(self.interface_graph[source][target]['len_contour']
            #                    for source, target in zip(loop2_path, loop2_path[1:]))
            # alternatives.append((path2_length, loop2_path, e1, e2, e3, loop0_path))
        # end for loop0_path in alternative_loop_paths
        return min(alternatives)



    def loop_formation_effects(self, elem1, elem2, reaction_type):
        r"""
        Returns tuple with:
            activity, loop_info
        Where
          * Activity is a numeric value describing the probability that the two elements are
            in close enough proximity to react, relative to the probability when the two elements
            are free in solution at a concentration of 1 M, and
          * loop_info is a dict describing the loop(s) that will be created (and, perhaps, split)
            by connecting elem1 and elem2.
        Arguments:
            :elem1:, :elem2:, :reaction_type: are the same as for intercomplex_activity().

        Making a new connection has two effects:
            1: Creating/splitting loops:
                1a: Creating a new loop.
                1b: Splitting an existing loop in two.
                1c: A split that is then merged into a single Ifnode (e.g. stacking of neighboring duplexes).
            2: Altering existing loops:
                e.g. stacking causes two helices to fuse.
        We must consider the actiity by evaluating the effect *after* performing the reaction.
        We should keep track of the effects:
            Create a new/primary loop.
            Split an existing loop in two.
            Ifnode merge.
            Joined helices.
            Other affected loops (side-effects).

        We typically want to get activity for the final state, i.e. for hybridization we want to calculate
        what the activity is with the duplex formed.
        Which strategy is better:
            (1) Perform the hybridization/stacking, calculate activity, then revert when done.
                Actually connecting the graph could make it hard to find shortest path (they are directly connected),
                but you could make other changes, e.g. adjust edge length/rigidity to what it would be when
                stacked/hybridized.
            (2) Perform calculations on unhybridized/unstacked state, but make adjustments as-needed
                to give results that are equivalent to the formed state.
                I'm currently doing this (2).

        """
        print("\nCalculating loop_formation_effects for '%s' reaction between %s and %s" % (reaction_type, elem1, elem2))

        if reaction_type is HYBRIDIZATION_INTERACTION:
            domain1, domain2 = elem1, elem2
            cmplx = domain1.strand.complex
            d1end5p, d2end5p = domain1.end5p, domain2.end5p
            d1end3p, d2end3p = domain1.end3p, domain2.end3p
            # path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
            # Domains are currently NOT hybridized and thus also NOT stacked. There should be NO delegation!
            assert d1end3p.ifnode.top_delegate() == d1end3p.ifnode
            assert d2end3p.ifnode.top_delegate() == d2end3p.ifnode
            assert d1end5p.ifnode.top_delegate() == d1end5p.ifnode
            assert d2end5p.ifnode.top_delegate() == d2end5p.ifnode
            reactant_nodes = {d1end3p.ifnode, d2end3p.ifnode, d1end5p.ifnode, d2end5p.ifnode}
            ## TODO: (OPTIMIZATION) Instead of finding the shortest path, try to
            ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP.
            ## If they are, then the existing loop should give the shortest path.
            ## (For now, just assert that IF the two loops are part of an existing, non-obsolete loop,
            ## then the shortest path is in that loop.)
            new_loop_path = shortest_path(self.interface_graph, d1end5p.ifnode, d2end5p.ifnode)
            slice_start = 1 if new_loop_path[1] == d1end3p.ifnode else 0
            slice_end = -1 if new_loop_path[-2] == d2end3p.ifnode else None
            # pdb.set_trace()
            if slice_start == 1 or slice_end is not None:
                new_loop_path = new_loop_path[slice_start:slice_end]
        elif reaction_type is STACKING_INTERACTION:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
            cmplx = h1end5p.domain.strand.complex
            # if h1end3p.domain.name in ("e", "C"):
            # pdb.set_trace()
            # Nomenclature: d: domain, h: Single helix, dh: double-helix
            # (in this case, duplex ends and double-helix ends are equivalent)
            #  DUPLEX/DOUBLE-HELIX 1
            #             dh1, d1         dh2, d4
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            #             dh1, d2         dh2, d3
            # So, we could use alternative names:
            # dh1end3p = d1end3p = h1end3p
            # dh1end5p = d2end5p = h2end5p
            # dh2end3p = d3end3p = h2end3p
            # dh2end5p = d4end5p = h1end5p
            # and elem1, elem2 would then be:
            # assert (dh1end3p, dh1end5p), (dh2end3p, dh2end5p) == (elem1, elem2) # double-helix ends
            # assert (d1end3p, d2end5p), (d3end3p, d4end5p) == (elem1, elem2)     # domain ends

            dh1_delegate, dh2_delegate = h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate()
            ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
            assert dh1_delegate, dh2_delegate == (h2end5p.ifnode.top_delegate(), h1end5p.ifnode.top_delegate())

            new_loop_path = shortest_path(self.interface_graph, dh1_delegate, dh2_delegate)

            ## Adjusting new_loop_path shouldn't be required when using InterfaceGraph/Nodes:
            # slice_start = 1 if new_loop_path[1] == h2end5p else 0
            # slice_end = -1 if new_loop_path[-2] == h1end5p else None
            # if slice_start == 1 or slice_end is not None:
            #     new_loop_path = new_loop_path[slice_start:slice_end]
            # However, you MUST put the start and end nodes together as they would be when stacked,
            # (unless you choose to resolve this by increasing unstacking rate (together with throttle)..
            # - This is done AFTER getting the new_loop_path segments/elements...
        else:
            assert isinstance(elem1, DomainEnd) and isinstance(elem2, DomainEnd)
            d1end5p, d2end5p = elem1, elem2
            cmplx = d1end5p.domain.strand.complex
            # domain1, domain2 = d1end5p.domain, d2end5p.domain # only used for debugging
            # new_loop_path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
            new_loop_path = shortest_path(self.interface_graph,
                                 d1end5p.ifnode.top_delegate(),
                                 d2end5p.ifnode.top_delegate())
            slice_start, slice_end = 0, None


        # activity for shortest_path aka "a1" aka "new_loop_activity":
        activity = self.calculate_loop_activity(new_loop_path, simulate_reaction=reaction_type)
        if activity == 0:
            return 0, None

        new_loop_nodes = set(new_loop_path)
        processed_secondary_loopids = set()
        all_affected_loops = set()
        # changed_loops = {}   # old_id => [new_loop1_dict, new_loop2_dict, ...]
        changed_loop_paths = {} # # lookup of loopid: new_loop_path for all changed loops
        changed_loops_by_hash = {} # old_loop_hash: {dict with new loop info.}
        # It would be nice to also have the activities for the new loops created (not just the shortest-loop)
        # Maybe have a 'changed_loops' old_loopid => [list of dicts describing new loops]
        loop_split_factors = []
        side_effects_factor = 1
        side_effects_activities = []
        side_effects_factors = []
        path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
        ## TODO: Remove the checks below >>
        ifnode_by_hash_index_before = set(cmplx.ifnode_by_hash.items())
        cmplx.rebuild_ifnode_by_hash_index()
        cmplx.rebuild_loopid_by_hash_index()
        ifnode_by_hash_index_after = set(cmplx.ifnode_by_hash.items())
        assert path_spec == cmplx.calculate_loop_path_spec(new_loop_path)
        assert ifnode_by_hash_index_before == ifnode_by_hash_index_after
        print("new loop (shortest_path):", new_loop_path)
        print(" - new_loop_path_spec:", path_spec)
        ## TODO: Remove until here <<
        new_loop = {
            'loop_hash': cmplx.calculate_loop_hash(path_spec),
            'path_spec': tuple(path_spec),
            'activity': activity,
            'dS': ln(activity),
            'new_loop_debug_info': {
                'path': tuple(new_loop_path), # Use immutable tuples to avoid mutation bugs
                'description': "loop_formation_effects new loop"
            },
        }
        print(" - new_loop info: %s" % (new_loop,))
        changed_loop_paths[None] = new_loop_path

        changed_loops_hash_by_id = {}
        new_stacking_loops = {}
        loop_effects = {
            # total loop formation activity; Product of shortest-path activity and any/all loop_change_factors:
            'activity': 0,
            'is_forming': True,
            # 'shortest_path': tuple(new_loop_path), # use ['new_loop']['path'] instead.
            # 'shortest_path_spec': path_spec, # use ['new_loop']['path_spec'] instead.
            # 'shortest_path_activity': activity, # use ['new_loop']['activity'] instead.
            'new_loop': new_loop,
            # 'changed_loops': changed_loops, # old_loopid => [list of new loop dicts]
            'changed_loops_by_hash': changed_loops_by_hash, # old_loop_hash:  {dict with new loop info.}
            'directive_debug_info': {
                # Instance-specific debug info and other unused values
                'changed_loops_hash_by_id': changed_loops_hash_by_id,
            },
            # 'changed_loops_hash_by_id': changed_loops_hash_by_id,
            # changed_loops is a list of dicts, each representing the new loop formed, with keys 'path' and activity.
            # 'loops_considered': processed_secondary_loopids, # for debugging
            # 'loops_affected': all_affected_loops,  # Obsoleted by changed_loops
            'loop_split_factors': loop_split_factors,
            # loops not reached by "neighbor-loop-propagation" algorithm but which should still be considered.
            # There may be additional loops that are relevant due to stacking (network-connected ends)
            # These are collected in new_stacking_loops:
            'new_stacking_loops': new_stacking_loops,
            'stacking_side_effects_activities': side_effects_activities,
            'stacking_side_effects_factors': side_effects_factors,
            'complex_loop_ensemble_fingerprint': cmplx.loop_ensemble_fingerprint,
            'complex_state_fingerprint': cmplx.state_fingerprint(),
            'description': "loop_formation_effects directive",
        }

        if not cmplx.loops:
            # Complex has no existing loops, just return already:
            return activity, loop_effects


        ## Useful NetworkX functions:
        ## networkx.algorithms.cycles.cycle_basis

        ## TODO: I'm still contemplating if it's easier to check for side-effects if we
        ## temporarily modify the interface_graph to look as it would if the reactants were reacted
        ## and then undo the change afterwards. Right now we are simulating this in the methods:
        ## * join_edges_groups, join_two_edge_groups, and
        ## * closed_loop_segments
        ## These may be called repeatedly for each secondary loop that is checked for side-effects.
        ## ...but we are usually modifying the already-calculated shortest_path, so no.



        ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
        ## TODO: CHECK FOR SECONDARY LOOPS!
        ## We can easily do this simply by whether the two nodes are already in an existing loop.
        ## NOTE: If both reactant nodes are part of an existing loop, we don't strictly have to
        ## spend resources calculating the shortest path: The shortest path should be a sub-path
        ## of the existing loop. And, we also need the "other" part of the loop:
        ## We calculate reaction activity for "internal" loop node reactions as:
        ##      ab̌c = ab̌ + b̌c - ac
        ## where ab̌c is the internal loop closure, ac is closing the outer loop,
        ## and ab and bc is closing each of the inner loops.
        ## There might also be a term describing the change that e.g. hybridization
        ## would have on the existing loop (secondary side-effects).
        ## Note: The above expression is for shape energies. If we are talking activities, then
        ## we simply multiply (instead of add) and divide (instead of subtract).

        ## Q: When joining two ifnodes, are we sure that at most ONE existing loop is dividied into two?
        ## I would expect this to be the case, I can't currently imagine how joining two nodes
        ## should be able to split more than one existing loop...

        ## Q: Should cmplx.loops_by_interface contain *all* loops or only the current effective loops?
        ## A: Only the effective; determining other "outer" loops could be ambigious, I think.

        ## TODO: Short-circuit shortest_path calculation by checking if reactants are part of existing loop first
        ## (right now we are using it for assertions)

        ## TODO: There are cases when we are splitting an existing loop, but not by direct "pinching":
        ## In fact, IN MOST CASES, path[0] and path[-1] will NOT BOTH be part of an existing loop, even though
        ## they are effectively splitting a loop in two. See notes on Wang-Uhlenbeck entropy.
        # if False and (path[0] in cmplx.ifnode_loopids_index and len(cmplx.ifnode_loopids_index[path[0]]) > 0
        #     and path[-1] in cmplx.ifnode_loopids_index and len(cmplx.ifnode_loopids_index[path[-1]]) > 0):
        #     ## Uh, OK so both elems are part of one or more existing loops. But we don't know if they are
        #     ##  part of the SAME loop. Not the best way to do it.. In any case, it is probably better
        #     ##  to handle this case like any other, as it is done below.
        #     ## For the record, this is what you should have done:
        #     loop_with_both_nodes = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[0]]
        #     ## TODO: Need to test this - e.g. case 1(a) or 1(b)
        #     ## TODO: I feel like we should still cosider *all* shared loops and not just break out by this one case ^^
        #     ## TODO: THIS CASE NEEDS TO BE COMPLETELY RE-WORKED!
        #     ## Perfect "pinching" of existing loop (rare)
        #     # We have a secondary loop.
        #     # Find loops containing both the first and last path interface nodes:
        #     shared_loopids = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[-1]]
        #     assert len(shared_loopids) == 1           ## Not sure about this...
        #     loop0id = next(iter(shared_loopids))
        #     loop0 = cmplx.loops[loop0id]
        #     loop0_path = loop0['path'] # list of ifnodes from when the loop was made.
        #     # IF the reaction/node-merge splits more than ONE loop, we want to know!
        #     # Can we always use the smallest loop, or is there a particular order that must be obeyed?
        #     # Also, is the number of nodes "len(loop['nodes'])" the way to go?
        #     # It is not the same as the "shortest" path method..
        #     ## TODO: Check/fix this code:  - EDIT: Replaced by len(shared_loops) == 1 above.
        #     # loop0 = min(shared_loops, key=lambda loop: len(loop['nodes']))  # NOT GUARANTEED
        #     shortest_path_nodes = set(path)
        #     assert all(node in loop0['ifnodes'] for node in path)             # NOT GUARANTEED
        #     assert shortest_path_nodes <= loop0['ifnodes']
        #     # Using set operations is faster (path is subset of existing loop0)
        #     # Set operators: &: intersection, -: difference), ^: symmetric_difference, <=: is subset, >=: superset.
        #     path2_nodes = loop0['ifnodes'] - shortest_path_nodes
        #     loop0_path_idxs = (loop0_path.index(path[0]), loop0_path.index(path[-1]))
        #     if loop0_path_idxs[0] > loop0_path_idxs[1]:
        #         loop0_path_idxs = loop0_path_idxs[::-1]
        #     path1 = e1 = loop0_path[loop0_path_idxs[0]:loop0_path_idxs[1]+1]
        #     path2 = e2 = loop0_path[loop0_path_idxs[1]:] + loop0_path[:loop0_path_idxs[0]+1]
        #     # path is the shortest path and should equal either e1 or e2:
        #     assert all(node in path1 for node in path) or all(node in path2 for node in path)
        #     # a1 = self.calculate_loop_activity(path1, simulate_reaction=reaction_type)
        #     a2 = self.calculate_loop_activity(path2, simulate_reaction=reaction_type)
        #     a0 = loop0["activity"]
        #     ## THIS SHOULD JUST BE APPENDED TO loop_split_factors:
        #     loop_split_factors.append(a2/a0)
        #     processed_secondary_loopids.add(loop0id)
        #     pdb.set_trace()
        #     # If hybridization, special case where we account for the length of the newly-formed duplex.
        #     # TODO: This should not preclude processing of other shared loops, but rather be just an optimization.

        # elif cmplx.ifnode_loopids_index:

        ## TODO, WIP: Add bredth-first algorithm to search for all changed loops.
        # Algorithm:
        # Strategy:
        # * Initialize simply by the new loop changed_loops and new_loop_id=None to changed_loops_deque,
        # * Then loop over changed_loops_deque and process just like you are already doing.

        # Initialize bredth-first search for changed loops:
        changed_loops_deque = deque()
        changed_loops_deque.append(None) # None is used to temporarily identify the newly formed loop

        ## Note:
        ## After forming the new loop, Λ₁, if e3 < e2, then the old Λ₀ is updated to be the Λ₂ loop (aka Λ₀').
        ## If e3 is longer than e2, then we don't do anything. This is arguably an in-complete description
        ## of the shape entropy. But then, shape entropy should be described by the edges, not loops,
        ## since loops are not actually additive.

        ## Check for non-pinching loop-splitting cases: (partial loop sharing)
        ##           Λ₀              .------.
        ##            .----------,--´------. \
        ##           /   Λ₁     /          \  \ e₄ (not part of the original discussion)
        ##        e₁ \      e₃ /:/   Λ₂    /   \
        ##            \         /         / e₂  \
        ##             `----.--´---------´      |
        ##                   `-----------------´
        ## Where symbols Λ₀ Λ₁ Λ₂ e₁ e₂ e₃:
        ## Λ₁ = e₁+e₃ is the shortest path loop, for the reaction in question,
        ## Λ₀ = e₁+e₂ is an existing loop, overlapping the reaction shortest loop path,
        ## Λ₂ = e₂+e₃ is the path for the other loop, formed by splitting Λ₀.
        ## e₁ is the part shared by both Λ₁ and Λ₀ but not in Λ₂
        ## e₃ is the path in both Λ₁ and Λ₂ but not in Λ₀,
        ## e₂ is the part shared by both Λ₂ and Λ₀ but not in Λ1
        ## e2 > e1 because by shortest-path definition e₁+e₃ < e₂+e₃.
        ## The Wang-Uhlenbeck matrix can be written as:
        ##      ┌ e1+e3   e3  ┐   ┌ e1+e3   e1  ┐   ┌ e1+e2   e2  ┐
        ##  W = │             │ = │             │ = │             │
        ##      └  e3   e2+e3 ┘   └  e1   e1-e2 ┘   └  e2   e3+e2 ┘
        ## In general we have Wang-Uhlenbeck determinant (where dS = R α ln(det(W)))
        ##      Before forming e3: det(W₀) = Λ₀ = e₁+e₂
        ##      After  forming e3: det(W₁) = e₁e₂ + e₂e₃ + e₃e₁
        ##                         = (e1+e3) (e2+e3) - e3^2  = Λ₁ Λ₂ - e3^2
        ##                         = (e3+e1) (e2+e1) - e1^2  = Λ₁ Λ₀ - e1^2
        ## For e₃ == 0, we have the special case:
        ##  det(W₀) = e₁e₂       and det(W₀) = Λ₀ = e₁+e₂
        ##      a = a₁*a₂/a₀     where aᵢ is the activity calculated for loop Λᵢ alone.
        ## What about the case where e3 is a single, rigid element?
        ##  * This would be the stacking-reaction side-effects considered above..
        ##  * Means that Λ₁ and Λ₂ are independent. Same as for e3 == 0.
        ## Maybe all cases where Λ₁ and Λ₂ are independent can be calculated by multiplying a₁ with a₂/a₀.
        ## Can we consider the a₂/a₀ a "side-effect"?
        ## Thus, just go over all shared loops and ask:
        ## * Is e3 splitting the loop into two independent loops? - Calculate factor a₂/a₀.
        ## The criteria for e3 yielding independent loop splitting is:
        ## * If e3 is a single, rigid segment,
        ## * If e1 >> 0 and e3 < e1 (e1 < e2 by definition).
        ## For e1 ≈ 0:
        ##  det(W₁) = Λ₁ Λ₀ - e1^2 ≈ Λ₁ Λ₀
        ##  a ≈ a₁*a₀/a₀ = a₁
        ## This latter condition should also be applicable if e3 >> e2 > e1
        ## However, these case approximations assume e1, e2, e3 are made up of a big number of flexible kuhn links.
        ##
        ## Is there a way to look at these purely from a
        ##
        ## Can we view this as "the shortest-path activity modulated by the
        ## This is really non-obvious, but let's give it a try...
        ## First see if there is any overlap between the shortest path and existing loops:
        loop_split_factor = 1
        while changed_loops_deque:
            # For the first round, the loop considered is Λ₁ with loop_id None.
            # It is not appropriate to use "loop0" to describe this; you can use "loop1" or "new_loop"
            # pdb.set_trace()
            loop1_id = changed_loops_deque.popleft()
            # loop1_info = changed_loops[loop1_id]
            # loop1_info = cmplx.loops[loop1_id]
            # loop1_shortest_path = tuple(loop1_info['path'])
            loop1_shortest_path = tuple(changed_loop_paths[loop1_id])


            # loops_with_both_nodes = cmplx.ifnode_loopids_index[path[0]] & cmplx.ifnode_loopids_index[path[-1]]
            loop1_path_nodes = set(loop1_shortest_path)
            # Quick way to get all nodes on the shortest path that are already on *ANY* existing loop.
            shared_nodes = loop1_path_nodes.intersection(cmplx.ifnode_loopids_index.keys())
            if shared_nodes:
                shared_loopids = set.union(*[
                    cmplx.ifnode_loopids_index[ifnode] for ifnode in shared_nodes
                ]).difference(processed_secondary_loopids)
            else:
                shared_loopids = {}
            # assert loops_with_both_nodes <= shared_loopids
            # assert len(loops_with_both_nodes) <= 1
            ## We have a few special cases:
            ## If the new edge e₃ is completely rigid, then e₂ has no influence on Λ₁, nor does e₁ influence on Λ₂.
            ## - The "e₃ == 0" aka "loop pinching" above can be seen as a special case of this.
            ## Perhaps just simplify as:
            ## * If e₃ has fewer segments than e₁, e₂ then calculate activity as:
            ##      a = a₁*a₂/a₀, where a₁ = f(e₁+e₃)e₂, a₂ = f(e₂+e₃), and a₀ = f(e₁+e₂) is the outer loop activity.
            ##      and f is calculate_loop_activity function.
            ##   How is this different from what is done in the side-effects?
            ##   * Currently we only search for side-effects from stacking.
            ## * If e₃ has more segments than e₁, e₂ then calculate activity simply using the shortest path:
            ##      a = a₁ if e₁ < e₂ else a₂
            ## Symbols: Λ₀ Λ₁ Λ₂ e₁ e₂ e₃

            # if shared_loopids:
            # shared_loops = [ for loop0id in shared_loopids]
            ## TODO: Check that this does not overlap with the "side-effects" factors calculation..
            ## It actually seems like this could supplant the side-effects calculation.
            # Find nodes that are not part of any loops:
            # fully_unshared_nodes = shortest_path_nodes.difference(cmplx.ifnode_loopids_index.keys())
            for loop0_id in shared_loopids:
                loop0 = cmplx.loop_by_loopid[loop0_id]
                print("Processing possibly-affected loop0 %s: %s" % (loop0_id, loop0))
                # path is shortest-path loop Λ₁, shared loop is Λ₀ aka "loop0".
                # Makes lookup faster for long loop0 but with constant overhead.
                ## TODO (optimization): If path[0] and path[-1] are both in the loop, then this can be optimized,
                ##  because you don't have to group, but can simply calculate the e1, e2 paths directly.
                loop0_path = tuple(loop0['path'])
                # When stacking backbone-connected duplex-ends, we form (self-)loops with a single node and no edge.
                # Do not consider these for optimization.
                if len(loop0_path) < 2:
                    continue
                # TODO: Remove loop_hash assertion:
                loop0_path_spec = cmplx.calculate_loop_path_spec(loop0_path)
                loop0_hash = cmplx.calculate_loop_hash(loop0_path_spec)
                assert loop0['loop_hash'] == loop0_hash  # TODO: Remove recalculation check of loop0_hash
                loop0_nodes = set(loop0_path)
                assert len(loop0_nodes) > 1  # If path has more than 1 elems, the set should as well.
                # Casting group_iter to list: (TODO: Check if that is really needed...)

                # Note: loop1_shortest_path (e1+e3) is *NOT* a candidate path for loop0.
                # loop1_shortest_path is used to find the overlap (e1) between loop0 and loop1 and determine if the
                # "non-overlapping subpath" (e3) is a shorter route for loop0.
                e1, e3a, e3b, e3_groups, use_e3 = self.find_maybe_shorter_e1_e3_subpaths(
                    loop0_nodes, loop1_shortest_path, reaction_type)

                # Calculate new activity for new loop0 path
                if use_e3:
                    # New loops Λ₁ and Λ₂ are independent, multiply activity by factor a₁/a₀:
                    # TODO: Account for cases where e3a == [] or e3b == [] or both (e3==0) !
                    a0 = loop0['activity']

                    # Find e2:
                    e2_nodes = loop0_nodes - set(e1)
                    e2_or_e1 = [(is_shared, list(group_iter)) for is_shared, group_iter
                                in groupby(loop0_path, lambda node: node in e2_nodes)]
                    e2_subpaths = [sub_path_group[1] for sub_path_group in e2_or_e1 if sub_path_group[0]]
                    # If Λ₀ loop path starts inside e2, we will have two sub-paths, otherwise only 1:
                    if len(e2_subpaths) > 1:
                        e2 = e2_subpaths[1] + e2_subpaths[0]
                    else:
                        e2 = e2_subpaths[0]
                    # alternatively, you could have looked at whether loop0_path[0] and loop0_path[-1]
                    # are both in e2_nodes
                    # Figure out how e2 and e3 should be combined: e2+e3 or e3+e2.
                    # Note: e1, e2 and e3a/e3b all share the node at the intersection, right? Edit: No, only e1, e3
                    # Note: We are only concerned with the intersection nodes that are e3a[-1] and e3b[0]
                    # It does not make sense to use e3a[0] or e3b[-1] for this
                    # (which are the reaction pair and shortest_path start/end)
                    if e2[0] in self.interface_graph[e3a[-1]]:
                        # Connection between first part of e3 and start of e2, and last part of e3 and end of e2:
                        assert e3a[-1] in self.interface_graph[e2[0]]
                        assert e2[-1] in self.interface_graph[e3b[0]]
                        assert e3b[0] in self.interface_graph[e2[-1]]
                        new_loop_path = e3a + e2 + e3b
                    else:
                        # Connection between last node of first part of e3 and the last  node of e2,
                        # and between        first node of last part of e3 and the first node of e2:
                        assert e3a[-1] in self.interface_graph[e2[-1]] # last node of e3a + last of e2
                        assert e2[-1] in self.interface_graph[e3a[-1]]
                        assert e3b[0] in self.interface_graph[e2[0]]   # first node of e3a + first of e2
                        assert e2[0] in self.interface_graph[e3b[0]]
                        # Need to reverse direction of e3 or e2 before joining them:
                        # e3 is e3b extended by e3a (list performs best by extending at the end)
                        new_loop_path = e3a + e2[::-1] + e3b
                    print(" - The new loop produced a shorter path for existing loop0 %s, new_loop_path = %s" % (
                        loop0_id, new_loop_path))
                    a2 = self.calculate_loop_activity(new_loop_path, simulate_reaction=reaction_type)
                    if a2 == 0:
                        return 0, None
                    loop_split_factor = a2/a0
                    ## Question: Do we need to add loop1 = e1+e3 as replacement loop?
                    ## Current answer: No, loop1 will always be created as the shortest_path loop.
                    ## The 'e1' and 'e2' edges may be different for different loop0s,
                    ## but loop1 should always be the same (PROVIDED THAT THE SHORTEST PATH IS ON loop0)
                    ## Note that this is NOT the always the case below where I'm considering loops on the full ds helix.
                    # [ifnode.state_fingerprint() for ifnode in new_loop_path]
                    new_loop_path_spec = cmplx.calculate_loop_path_spec(new_loop_path)
                    # DEBUG: Test that the path_spec doesn't change if updating ifnode_by_hash index:
                    cmplx.rebuild_ifnode_by_hash_index()
                    cmplx.rebuild_loopid_by_hash_index()
                    assert new_loop_path_spec == cmplx.calculate_loop_path_spec(new_loop_path) # TODO: Remove assertion
                    assert cmplx.calculate_loop_hash(new_loop_path_spec) == cmplx.calculate_loop_hash(tuple(new_loop_path_spec))
                    # TODO: Remove index re-calculation and assertions above
                    loop_update_info = {
                        'loop_hash': cmplx.calculate_loop_hash(new_loop_path_spec), # Used to identify the loop instance
                        'path_spec': tuple(new_loop_path_spec), # Is used to re-create the loop path ifnodes
                        'activity': a2,
                        'dS': ln(a2), # in units of R
                        # Values for the original loop (mostly for debugging)
                        'loop_update_debug_info': {
                            'description': "loop_formation_effects updated loop found by breath-first search starting at new loop.",
                            'new_path_tuple': tuple(new_loop_path), # Only for debugging; the result is cached and uses between complexes
                            'loop_change_factor': loop_split_factor, # = a2/a0
                            'old_activity': a0,
                            'old_loop_hash': loop0_hash, # so we can relate this to the original/parent loop.
                            'old_path_spec': loop0_path_spec,
                            'source_loop_id': loop0_id,
                        },
                    }
                    print(" - New info (path, hash, activity) for loop0 %s: %s\n" % (loop0_id, loop_update_info))
                    assert loop0_id not in changed_loop_paths
                    # changed_loops[loop0_id] = loop_update_info
                    changed_loop_paths[loop0_id] = new_loop_path
                    assert loop0_hash not in changed_loops_by_hash
                    changed_loops_by_hash[loop0_hash] = loop_update_info
                    changed_loops_hash_by_id[loop0_id] = loop0_hash
                    loop_split_factors.append(loop_split_factor)
                    # Add loop0_id to deque so we can check if loop0 has any affected neighboring loops:
                    changed_loops_deque.append(loop0_id)
                # end if use_e3 (loop0 should be updated)
            # end for loop0_id in shared_loopids:
            ## TODO: THIS IS NOT DONE YET!
            # Union-assignment operator should modify in-place and be equilvalent to myset.update(other_set)
            processed_secondary_loopids |= shared_loopids
            processed_secondary_loopids.add(loop1_id)
        # end while changed_loops_deque


        #### Checking for further secondary side-effects: ####
        ## TODO: Make sure side-effects does not overlap with primary loop calculation
        ## Even if the reaction form a new loop, the reaction can still have structural effects
        ## that alters the energy of existing loops. The presense of these loops would equally
        ## influence the activity for forming this new loop.
        ## E.g. if stacking of two helices would require excessive stretching.
        ## Evaluate secondary effects of the reaction for existing loops:
        if reaction_type is STACKING_INTERACTION:
            # It should be enough to consider existing loops; if there are connections not part of
            # an existing loop, that connection should be part of the new loop created by this reaction.
            # Any loop with nodes on dh1/dh2:
            # What are double-helices? DomainEnds? Or IfNodes? Or an object instance with both?
            # dh1 = cmplx.helix_by_domainend[h1end3p]  # double-helix, duplexed and fully-stacked
            # dh2 = cmplx.helix_by_domainend[h2end3p]
            # dh1_ifnodes, dh2_ifnodes = ({end.ifnode for end in dh1}, {end.ifnode for end in dh2})
            # Actually, it's probably easier to calculate it as-needed using:
            # dh1_upstream_ifnodes = list(h1end3p.upstream_stacked_top_ifnodes_generator())
            # dh2_upstream_ifnodes = list(h2end3p.upstream_stacked_top_ifnodes_generator())
            # Treat the two stacking double-helices as arms, pivoting at the stacking junctiong:
            # Both are in direction inside-out away from the junction and includes the first reactant.
            dh1_upstream_ifnodes = h1end3p.upstream_stacked_top_ifnodes_list([dh1_delegate])
            dh2_upstream_ifnodes = h2end3p.upstream_stacked_top_ifnodes_list([dh2_delegate])
            # dh_ifnodes_list = (dh1_upstream_ifnodes[::-1]
            #                    + [top_ifnode for top_ifnode in
            #                       (end.ifnode.top_delegate() for end in (h1end3p, h2end3p))
            #                       if (top_ifnode != dh1_upstream_ifnodes[0] and
            #                           top_ifnode != dh2_upstream_ifnodes[0])]
            #                    + dh2_upstream_ifnodes)
            # dh_arm_idx_by_ifnode = {sign*i: ifnode
            #                         for sign, nodes in ((1, dh1_upstream_ifnodes), (-1, dh2_upstream_ifnodes))
            #                         for i, ifnode in enumerate(nodes)}
            dh1_ifnodes = set(dh1_upstream_ifnodes)
            dh2_ifnodes = set(dh2_upstream_ifnodes)
            dh_ifnodes = dh1_ifnodes | dh2_ifnodes
            dh1_nodes_in_loops = dh1_ifnodes.intersection(cmplx.ifnode_loopids_index.keys())
            dh2_nodes_in_loops = dh2_ifnodes.intersection(cmplx.ifnode_loopids_index.keys())
            # Edit: Loops are currently stored as mutable dicts, you cannot make a set of loops.
            # You must either use loop IDs or create an object for each loop.
            dh1_loopids = {loopid for ifnode in dh1_nodes_in_loops if ifnode in cmplx.ifnode_loopids_index
                           for loopid in cmplx.ifnode_loopids_index[ifnode]}
            dh2_loopids = {loopid for ifnode in dh2_nodes_in_loops if ifnode in cmplx.ifnode_loopids_index
                           for loopid in cmplx.ifnode_loopids_index[ifnode]}
            # Alternatively, go through all loops and check if any of the loop's nodes are on the doublehelix:
            # dh1_loops = {loopid for loopid in cmplx.loops.keys() if any(...) }
            ## NOTE: We should not consider loops that are only on one of the dh arms and not the other,
            ## hence the two set comprehensions which are then followed by the intersection.
            ## Affected loops are loops which goes through both double-helix 1 and 2:
            affected_loops = dh1_loopids & dh2_loopids
            all_affected_loops |= affected_loops  # equivalent to in-place update(other_set)
            ## By definition, the shortest path for loops with nodes on both dh1 and dh2 should now
            ## go through the created stacked junction.
            ## In general, a stacking interaction will cause affected loops to split into two loops.
            ## However, there is a special but frequent case where one of the new loops will simply be the
            ## phosphate backbone between two adjacent duplexes. In that case we don't have to split the loop in two,
            ## but simply re-calculate it.
            ## Loop over all affected loops and re-calculate the activity for them.
            ## We then take the product: side_effects_factor = ∏(a_new/a_old for loop in affected loops)
            ## However, we must either modify the path to ensure it goes through the stacked junction,
            ## or re-calculate the shortest path for the loop, expecting it to go through the junction.
            # Optimization:
            affected_loops -= processed_secondary_loopids
            for loop0id in affected_loops:
                print("After updating the initial search for changed loops, this should never happen anymore!")
                pdb.set_trace()
                if loop0id in processed_secondary_loopids:
                    # pass
                    continue
                ## TODO: Add another (optional; expensive) check to look for "network-connected helix stacking branches"

                loop0 = cmplx.loop_by_loopid[loop0id]
                # What is stored in a loop? Let's just start by letting loop be a dict with whatever info we need:
                loop0_path = loop0["path"]   # List of ifnodes
                assert loop0_path == [ifnode.top_delegate() for ifnode in loop0_path] # make sure we have top delegates
                loop0_path_spec = cmplx.calculate_loop_path_spec(loop0_path)
                loop0_hash = cmplx.calculate_loop_hash(loop0_path_spec)
                assert loop0['loop_hash'] == loop0_hash  # TODO: Remove recalculation check of loop0_hash
                loop0_nodes = set(loop0_path)
                assert len(loop0_path) == len(loop0_nodes)
                # Make sure the shortest path is not a sub-path of loop0, this should've been processed above.
                assert not new_loop_nodes <= loop0_nodes
                # new_paths = []
                a0 = loop0['activity']
                a1a2 = []
                # changed_loops[loop0id] = {
                #     'new_paths': new_paths,
                #     'a1a2': a1a2,   # Primary loop closing activity for each new path
                #     'a0': a0,       # Store this, for good measure...
                #     #'loop_pinching_activity': loop_pinching_activity  # Combined effect of splitting old loop in two
                # }
                replacement_loop1 = {}
                # list of dicts. Two dicts, unless the two ends are currently adjacent.
                # Note: Any loop with nodes present on the shortest-path should be considered in the check above.
                # If the two ends are adjacent in the interface graph, then they are certainly on the shortest path!
                if dh1_delegate in loop0["path"] and dh2_delegate in loop0["path"]: #and \
                    ## TODO: MAKE SURE WE ARE NOT DOUBLE-COUNTING; C.F. CALCULATION USING SECONDARY LOOP SEARCH BELOW.
                    print("WEIRD: Reaction ifnodes %s and %s are both on an existing loop, but somehow this loop "
                          "was not included when processing shared_loops??")
                    # abs(loop0["path"].index(dh1_delegate) - loop0["path"].index(dh2_delegate)) == 1:
                    # We can just modify the path, joining it together at the new junction:
                    idx1, idx2 = loop0["path"].index(dh1_delegate), loop0["path"].index(dh2_delegate)
                    if idx1 > idx2:
                        idx1, idx2 = idx2, idx1
                    new_path1 = loop0["path"][idx2:] + loop0["path"][:idx1+1]  # +1 to include the node.
                    # new_paths.append(new_path1)  # new_paths = [new_path1, new_path2] aka [e1, e2]
                    a1 = self.calculate_loop_activity(new_path1, simulate_reaction=reaction_type)
                    if a1 == 0:
                        return a1, None
                    a1a2.append(a1)
                    replacement_loop1['path'] = new_path1
                    replacement_loop1['path_spec'] = cmplx.calculate_loop_path_spec(new_path1)
                    replacement_loop1['loop_hash'] = cmplx.calculate_loop_hash(replacement_loop1['path_spec'])
                    replacement_loop1['a0'] = a0
                    replacement_loop1['loop0_hash'] = loop0_hash
                    replacement_loop1['old_loop0_id'] = loop0id
                    replacement_loop1['activity'] = a1
                    replacement_loop1['loop_change_factor'] = a1/a0
                    replacement_loop1['description'] = "loop_formation_effects updated loop (stacked helices).",

                    if idx2-idx1 == 1:
                        # Case 1a: Special case: The two joined interface nodes are adjacent; just re-calculate current loop0.
                        # No extra paths
                        new_stacking_loops[loop0_hash] = [replacement_loop1]
                        # Case 0 is shortest_path is a complete subset of loop0, i.e. e3=0.
                        # Case 1 is loop0 is partially overlapping with shortest path,
                        # Case 3 is no overlap with shortest_path; loop0 detected through stacked double-helix.
                        replacement_loop1['description'] = ("Case 0: shared loop affected by duplex stacking "
                                                            "of adjacent ifnodes.")
                        print("WEIRD: ifnode %s and %s are adjacent for STACKING reaction, "
                              "yet they were not processed when considering loops on the shortest path!?")
                        pdb.set_trace()
                    else:
                        # Case 1b: Split the existing path into two new loops:
                        # We have already split out one of the loops, just need the other:
                        new_path2 = loop0["path"][idx1:] + loop0["path"][:idx2+1]
                        # new_paths.append(new_path2)
                        a2 = self.calculate_loop_activity(new_path2, simulate_reaction=reaction_type)
                        if a2 == 0:
                            return 0, None
                        a1a2.append(a2)
                        replacement_loop2 = {}
                        replacement_loop2['path'] = new_path2
                        replacement_loop2['path_spec'] = cmplx.calculate_loop_path_spec(new_path2)
                        replacement_loop2['loop_hash'] = cmplx.calculate_loop_hash(replacement_loop2['path_spec'])
                        replacement_loop2['a0'] = a0
                        replacement_loop2['loop0_hash'] = loop0_hash
                        replacement_loop2['old_loop0_id'] = loop0id
                        replacement_loop2['activity'] = a2
                        replacement_loop2['loop_change_factor'] = a2/a0
                        replacement_loop2['description'] = "Case 0: shared loop splitted by duplex stacking."
                        new_stacking_loops[loop0_hash] = [replacement_loop1, replacement_loop2]
                else:
                    # Case 2: Stacking forms a genuine new loop paths: loop0 -> loop1, loop2
                    # (I think this will also take care of case 1b and even 1a)
                    # We have to generate two new loops, partially composed of the old path
                    # and partially of the newly-stacked double-helix:
                    # OBS OBS: One of the two loops you are creating will already have been
                    # First, split the existing loop0 in two on the nodes that touch the newly-stacked double-helix:
                    grouped_path = [(is_in_dh, list(group_iter)) for is_in_dh, group_iter in
                                    itertools.groupby(loop0_path, lambda ifnode: ifnode in dh_ifnodes)]
                    # if len(grouped_path) in (2, 3): print("Case 1a")
                    # should be either 4 or 5 elements.
                    # 4 if the loop0 path happens to start/end right between the double-helix and the other parts.
                    assert 4 <= len(grouped_path) <= 5
                    if len(grouped_path) == 5:
                        # Remove the first element and append it to the last:
                        grouped_path[-1].extend(grouped_path.pop(0))
                        assert grouped_path[0][0][0] == grouped_path[0][-1][0]
                    assert len(grouped_path) == 4
                    ## TODO: Having done the grouping above, we could probably treat all cases (1a/1b/2a/2b) the same
                    if grouped_path[0][0]: # The first group is on the double-helix
                        e1, t1, e2, t2 = [path_group for is_in_dh, path_group in grouped_path]
                    else:
                        t1, e2, t2, e1 = [path_group for is_in_dh, path_group in grouped_path]

                    # Note: The variable nomenclature here is different from the "loops on shortest path" search above;
                    # e1, e2 are parts of the newly-formed arm, i.e. e3a and e3b above.
                    new_stacking_loops[loop0id] = []
                    for dh_before, t_group, dh_after in ((e1, t1, e2), (e2, t2, e1)):
                        # double-helix before the t-loop, then t-loop, then the double-helix after.
                        # ARGH: One of the loops just formed will already have been formed when processing/"splitting"
                        # the previous loop.
                        # These are signed, + for dh1 and - for dh2:
                        # e1_start, e1_end, e2_start, e2_end = (
                        #     dh_arm_idx_by_ifnode[dh_before[0]], dh_arm_idx_by_ifnode[dh_before[-1]],
                        #     dh_arm_idx_by_ifnode[dh_after[0]], dh_arm_idx_by_ifnode[dh_after[-1]])
                        # dh_start_idx = dh_ifnodes_list.index(dh_before[-1]) # dh_arm_idx_by_ifnode[dh_before[-1]]
                        # dh_end_idx = dh_ifnodes_list.index(dh_after[0]) # dh_arm_idx_by_ifnode[dh_after[0]]
                        if dh_before[-1] in dh1_ifnodes:
                            assert dh_after[0] in dh2_ifnodes
                            arm1, arm2 = dh1_upstream_ifnodes, dh2_upstream_ifnodes
                        else:
                            assert dh_before[-1] in dh2_ifnodes
                            assert dh_after[0] in dh1_ifnodes
                            arm1, arm2 = dh2_upstream_ifnodes, dh1_upstream_ifnodes
                        arm1_idx = arm1.index(dh_before[-1])
                        arm2_idx = arm2.index(dh_after[0])
                        # if arm1_idx == 0 and arm2_idx == 0: print("Case 1(b)")
                        # Revert arm2 so order is outside-in:
                        new_path = arm1[:arm1_idx+1] + t_group + arm2[arm2_idx::-1]  # start-at:stop-before:step-by
                        assert not set(new_path) == new_loop_nodes
                        # This "new_path" is already considered as the shortest path.
                        # Note that this MUSTN'T happen: The present loop0 should have been included in
                        # processed_secondary_loopids when checking shared_loops above.
                        # new_paths.append(new_path)
                        print("Adding new stacking-induced path:\n", new_path)
                        #pdb.set_trace()
                        # ... so, is the idea to just take the product of all of them ?
                        # There is a chance that one of the loops is going to be the primary Λ₁ loop with activity a₁.
                        # So if you do it this way, then make sure you are not including this twice. ₀₁₂
                        # Also, usually I would say a = a₁ * prod(a₂/a₀ for a₀, a₂ in independent loops)
                        # However, here I am not including a₀ in my side-effect factor,
                        # but only including it at the end. Not sure that is right.
                        a2 = self.calculate_loop_activity(new_path, simulate_reaction=reaction_type)
                        if a2 == 0:
                            return 0, None
                        # side_effects_new_loop_activities.append(a2)
                        side_effects_factors.append(a2/a0)
                        a1a2.append(a2)
                        replacement_loop = {}
                        replacement_loop['path'] = new_path
                        replacement_loop['a0'] = a0
                        replacement_loop['old_loop_id'] = loop0id
                        replacement_loop['activity'] = a2
                        replacement_loop['loop_change_factor'] = a2/a0
                        replacement_loop['description'] = "Case 1+stack: shared loop on shortest_path from stacking."
                        new_stacking_loops[loop0id].append(replacement_loop)
                    side_effects_activities.append(a1a2[0]*a1a2[1]/a0)

                # pdb.set_trace()
                # if len(new_paths) == 1: # Optimization for this special case (excessive optimization)
                #     path_activities = [self.calculate_loop_activity(new_paths[0], simulate_reaction=reaction_type)]
                #     loop_pinching_activity = path_activities[0]/loop0["activity"]
                # else:
                # path_activities = [self.calculate_loop_activity(new_path, simulate_reaction=reaction_type)
                #                    for new_path in new_paths]
                #loop_pinching_activity = prod(path_activities)/loop0["activity"]
                #side_effects_factors.append(loop_pinching_activity)
                # For the moment, just: old_loopid => dict with new loop paths
                # changed_loops[loop0id] = {
                #     'new_paths': new_paths,
                #     'a1a2': a1a2,   # Primary loop closing activity for each new path
                #     'a0': a0,       # Store this, for good measure...
                #     #'loop_pinching_activity': loop_pinching_activity  # Combined effect of splitting old loop in two
                # }
                processed_secondary_loopids.add(loop0id)
            # end for loop in affected_loops

            # side_effects_factors = [(self.calculate_loop_activity(loop["path"], simulate_reaction=reaction_type)/
            #                          loop["activity"]) for loop in affected_loops]

            # if any(factor == 0 for factor in side_effects_factors):  # Edit: Is checked as early as possible
            #     return 0
            side_effects_factor = prod(side_effects_factors)


        ## TODO: Compare stacking side_effect_factors and loop_split_factors
        # if any(factor == 0 for factor in loop_split_factors):
        #     return 0
        if loop_split_factors:
            loop_split_factor = prod(loop_split_factors)
            print("Shortest-loop activity: %0.05g" % activity)
            print("Loop split factors:",
                  " * ".join("%0.02g" % v for v in loop_split_factors),
                  " = %0.03g" % loop_split_factor)
            print("Activity incl. loop split factors: %0.02g * %0.02g = %0.03g" %
                  (activity, loop_split_factor, activity*loop_split_factor))
        if side_effects_factors:
            # Stacking side-effect factors include a0 for all loops; a = a1/a0 * a2/a0
            print("Stacking side-effects activities: ",
                  " * ".join("%0.02g" % v for v in side_effects_activities),
                  " = %0.03g" % prod(side_effects_activities))
            print("Stacking side-effects factors:    ",
                  " * ".join("%0.02g" % v for v in side_effects_factors),
                  " = %0.03g" % side_effects_factor)
            print("Note: Stacking side-effect factors include a0 for all loops; a = a1/a0 * a2/a0. "
                  "I.e. they divide by a0 once too many.")
            #pdb.set_trace()


        """
        ### ALTERNATIVE IDEA: ###
        Assume you are always gonna form exactly one new loop.
        Everything else is side-effects.
        * You form the primary loop (shortest path).
        * Then you look at secondary loops = loops touching the primary loop / shortest path.
        * For secondary loops which were modified:
        * Then you look at tertiary loops = loops touching the modified secondary loops
        * For the tertiary loops that were modified:
        * Repeat, etc.
        At every level you define:
            e1 = the edge which this loop have in common with the previous loop.
            e3 = the edge which split the previous loop in two (e.g. a stacking duplex).
            e2 = the edge in this loop that is not e1

        This should probably be done with a breath-first algorithm using a deque or similar.

        For loops between two stacking double helices, since all existing loops will be split, this should work.
        What if you have more complex networks, where the loops are not just


        ### TODO - a better way to determine the effect of stacking double-helices: ###
        For stacking interactions of domains, there is no guarantee that looking at (primary/shortest loops)
        will detect if the helices being stacked are connected by other means upstream of the stacking ends.
        A better approach would be to detect branch-points:
            1. Find nodes on each arm that are branching off, and order them by distance to the stacking ends,
                with the nodes furthest away first.
            2. For each pair of furthest away from the stacking ends, find the shortest path between the two nodes.
                If the shortest-path goes through a node downstream (closer to the stacking ends),
                make that node the current node and discart the more distant nodes.
            3. Calculate loop activity for path = (arm1 + current_shortest_path + arm2)
            4. Proceed to the next pair in the list until the list is empty.

        """



        # loop_effects = {
        #     'activity': 0,
        #     'shortest_path': path,
        #     'shortest_path_spec': path-of-ifnode-fingerprints
        #     'shortest_path_activity': activity,
        #     'changed_loops': changed_loops,  # old_id => [new_loop1_dict, (new_loop2_dict]
        #     'changed_loops_specs': changed_loops_specs, # loop_state_hash => [path of ifnodes_state_fingerprints]
        #     'loops_considered': processed_secondary_loopids,
        #     'loops_affected': all_affected_loops,
        #     'loop_split_factors': loop_split_factors,
        #     'stacking_side_effects_activities': side_effects_activities,
        #     'stacking_side_effects_factors': side_effects_factors,
        # }

        # "activity" variable is loop activity for the newly-created loop.
        # How to get the actual path and side-effects when performing the reaction?
        # Return activity *and* loop info with path and side_effects (dict?) - for caching?
        if loop_split_factors:
            activity *= loop_split_factor
        loop_effects['activity'] = activity
        return activity, loop_effects


    def find_maybe_shorter_e1_e3_subpaths(self, loop0_nodes, loop1_path, reaction_type):
        r"""
        Arguments:
            loop0_nodes: Nodes in a current loop for which we are considering there might be shorter path.
            loop1_path: A recently created or updated loop forming a new shortest path.
        Note: loop1_path is *NOT* the "possibly shorter candidate path for loop0.
        It is a path that overlaps with loop0 and which is the shortest path for another
        new or recently updated loop, but it is not a path that we are considering for loop0.
        Instead, what we are considering is a new loop path for loop0 formed partially by
        loop0 and partially by loop1_path, specifically the part of loop1 which does not overlap with loop0.
        If you want to compare a current loop with a candidate loop, there is another method for doing that.

        Consider the system below. It currently has a single loop, Λ₀. We have are forming a new
        edge connection, e₃, producing a new loop (Λ₁) with shortest path e₁+e₃.
        We want to know if Λ₀ (e₂+e₁) can be made shorter by using e₃ instead of e₁, forming the updated
        shortest-path e₂+e₃ (Λ₂).
        To do that, we must first determine what e₁ and e₃ sub-paths actually are (and by subtraction e₂).
            Λ₀  .-------,--------.
               /   Λ₁  /  e₃     \  e₂
            e₁ \      /:/    Λ₂  /
                `------´--------´
        This method will find e1 and e3 sub-paths for loop1, where
        e1 is the sub-path that is on Λ₀ and Λ₁ but not Λ₂,
        e3 is the sub-path that is on Λ₁ and Λ₂ but not Λ₀,
        e2 is the sub-path that is on Λ₀ and Λ₂ but not Λ₁.

        Note: The nomenclature of e1, e3 is relative to the loop being considered defined as Λ₀ and the
        recently updated shortest-path for another loop Λ₁.

        That is, after creating the initial "new loop" (= Λ₁) and determining that we should update the old loop,
        Λ₀ --> Λ₂, then we may want to look at a another existing loop, Λ₃ = (e₂+e₄), and determine if this loop
        should use the possibly-shorter path e₃+e₄ instead.
        We would then call find_maybe_shorter_loop2_subpaths(loop0_path=Λ₃=e₂+e₄, loop1_path=Λ₂=e₂+e₃, ...)
        Within this method, Λ₃=e₂+e₄ --> Λ₀ (e₂+e₁)  and  Λ₂=e₂+e₃ --> Λ₁=e₁+e₃.

        Question: Should the two "intersection nodes" where e1, e2 and e3 meet be part of e1, e2, e3?
        It makes sense that they are part of e2 because obviously they are shared.
        However, when comparing the original e1 subpath length with the new e3 subpath, it is nice if both sub-paths
        contain the intersection nodes because we need to include the edge-length to/from the intersection nodes.
        So, although it is not 100% intuitive, we define the intersection nodes
        to be part of e1 and e3 and NOT PART of e2.

        Also note that loop1_path must start and end on e3, i.e. at the new connection being formed.

        """
        ### 1. Find e1 and e3 subpaths: ###
        ### TODO: DOES THIS WORK if loop1_path does not start and end on e3??
        ### This may be the case when working through all the "possibly changed" loops..
        ### Although I think for changed_loops, we always set the path to start at e3.

        e1_or_e3_sub_paths = [(is_shared, list(group_iter)) for is_shared, group_iter in
                              groupby(loop1_path, lambda ifnode: ifnode in loop0_nodes)]
        # grouped_path consists of the first part of e3, then e1, then the last part of e3.
        # Do we ever encounter e2? Well, not when we are grouping shortest path (e3+e1 by definition)..
        # Except if e3a or e3b have zero length (full or half pinching).
        # If e3 == 0 (e3_nodes = Ø), then we'll only have a single group!
        # This can happen both for stacking and hybridization.
        # There are also cases where len(e1_or_e3_sub_paths) == 2!
        assert len(e1_or_e3_sub_paths) <= 3
        if len(e1_or_e3_sub_paths) == 1:
            # If we only have one sub-path, it is e1; e3 is just the two nodes at the start/end of e1.
            assert e1_or_e3_sub_paths[0][0] is True  # verify that e₁ is on Λ₀
            e1 = e1_or_e3_sub_paths[0][1]
            e3a, e3b = [e1[0]], [e1[-1]]
        elif len(e1_or_e3_sub_paths) == 2:
            if e1_or_e3_sub_paths[0][0] is True:
                # We have e3b > 0; e3a starts at e1 (has zero length)
                e1 = e1_or_e3_sub_paths[0][1]
                e3a, e3b = [e1[0]], [e1[-1]]
                e3b.extend(e1_or_e3_sub_paths[1][1])
            else:
                # We have e3a > 0; e3b starts and ends on e1 (has zero length)
                e3a, e1 = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1]
                e3a, e3b = [e1[0]], [e1[-1]]
                e3a.append(e1[0])
        else:
            assert tuple(tup[0] for tup in e1_or_e3_sub_paths) == (False, True, False)
            # Note: grouped_path is grouped by whether node is on Λ₀ or not, not stiffness:
            # The "intersection" nodes are usually considered part of both e1, e2 and e3:
            # intersect_nodes = (e1[0], e1[-1])
            e3a, e1 = e1_or_e3_sub_paths[0][1], e1_or_e3_sub_paths[1][1]
            e3a.append(e1[0])
            # e3b starts with the last node in e1 and then extends the rest of the sub-path:
            e3b = [e1[-1]]
            e3b.extend(e1_or_e3_sub_paths[2][1])

        ### 2. Compare the lengths of e1 and e3 subpaths: ###
        # TODO: Use a Blist to store path-nodes (you will be doing lots of arbitrary inserts/deletions)
        # grouped by stiffness, i.e. similar to segments
        e1_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                     in self.group_interfaces_path_by_stiffness(e1)]
        e3a_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                      in self.group_interfaces_path_by_stiffness(e3a)] # First part of e3
        e3b_groups = [(stiffness, list(group_iter)) for stiffness, group_iter
                      in self.group_interfaces_path_by_stiffness(e3b)] # Last part of e3
        ## If e3a and e3b are just single nodes, we won't have any edges!
        ## Maybe it is better to join the groups first before joining the edges?
        ## Or maybe we need to include the intersection node in both e1 and e3a/e3b?

        if e3a_groups and e3b_groups:
            e3_groups = self.join_two_edge_groups(e3a_groups, e3b_groups, simulate_reaction=reaction_type)
        # elif e3a_groups:
        #     # start of e3a stacks directly with end of e1
        #     e3_groups = self.join_two_edge_groups(e3a_groups, e1, simulate_reaction=reaction_type)
        else:
            # This is new (but not unexpected; will happen for empty e3a/e3b) - seems to work.
            # TODO: CHECK THIS MORE.
            e3_groups = [group for segment in (e3a_groups, e3b_groups) if segment for group in segment]
        # Uh, it would, perhaps be nice to be able to join/react un-grouped paths,
        # and then group afterwards...
        # edge groups is a list of (stiffness, [(length, len_sq, stiffness, source, target), ...]) tuples.
        # Use edge groups if you need to process something.
        # (it is also slightly cheaper because it does not calculate segment length/length_sq sums).
        if len(e3_groups) <= 1:
            use_e3 = True
        else:
            # Just summing length and see if e3 is shorter than e1.
            e1_length = sum(edge_tup[0] for stiffness, edge_tuples in e1_groups for edge_tup in edge_tuples)
            e3_length = sum(edge_tup[0] for stiffness, edge_tuples in e3_groups for edge_tup in edge_tuples)
            # summing segment-length squared:
            e1_len_sq = sum(sum(etup[0] for etup in edge_tuples)**2 if stiffness > 0 else
                            sum(etup[1] for etup in edge_tuples)
                            for stiffness, edge_tuples in e1_groups)
            e3_len_sq = sum(sum(etup[0] for etup in edge_tuples)**2 if stiffness > 0 else
                            sum(etup[1] for etup in edge_tuples)
                            for stiffness, edge_tuples in e3_groups)
            if e3_length < e1_length and e3_len_sq < e1_len_sq:
                use_e3 = True
            else:
                use_e3 = False
        return e1, e3a, e3b, e3_groups, use_e3




    ### OLD/OBSOLETE METHODS: ###

    def domains_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        """
        return shortest_path(self.domain_graph, domain1, domain2)

    def ends5p3p_shortest_path(self, end5p3p1, end5p3p2):
        """ :end5p3p1:, :end5p3p2: DomainEnd nodes (either End5p or End3p),
        calculates path from node1 to node2 using ends5p3p_graph. The returned path is a list
        starting with node1 and ending with node2, and has the form: [node1, <all nodes in between, ...>, node2]
        TODO: This should certainly be cached.
        TODO: Verify shortest path for end3p as well?
        """
        return shortest_path(self.ends5p3p_graph, end5p3p1, end5p3p2, weight='dist_ee_sq')

    def ends5p3p_path_partial_elements(self, path, length_only=False, summarize=False):
        """
        Returns a list of structural elements based on a 5p3p-level path (list of DomainEnd nodes).

        For this, I will experiment with primitive list-based representations rather than making full
        element objects.
        path_edges = [(1, [(length, length_sq, source, target), ...]),
                      (3, [(length, length_sq, source, target), ...]),
                      (2, [(length, length_sq, source, target)]
        Where 1 indicates a list of single-stranded edges,
        2 indicates a hybridization edge, and 3 indicates a list of stacked edges.
        Since only the stiffness type and lenghts are important, maybe just
        path_edges = [(1, [length, length, ...]), (3, [length, length, ...], ...)]
        Edit: Instead of 1/2/3, use the standard INTERACTION constants values.
        """
        path_edges = []
        stiffness = None
        last_stiffness = None
        stiffness_group = None
        ## CONGRATS: YOU JUST RE-INVENTED itertools.groupby()
        for i in range(len(path)-1):
            # source, target are DomainEnd instances, either end5p, end3p or reversed.
            source, target = path[i], path[i+1]
            # Method 1: Use ends5p3p_graph edge attributes (which must include stacking edges!):
            edge_attrs = self.ends5p3p_graph[source][target]
            if self.ends5p3p_graph.is_multigraph():
                # edge_attrs is actually a dict with one or more edge_key => edge_attrs pairs.
                # Select edge with lowest R squared: ('dist_ee_sq' should always be present for 5p3p graphs)
                key, edge_attrs = min(edge_attrs.items(), key=lambda kv: kv[1]['dist_ee_sq'])
            # for multigraphs, graph[u][v] returns a dict with {key: {edge_attr}} containing all edges between u and v.
            # interaction = edge.get('interaction')
            # length = edge_attrs.get('length')     # TODO: Reach consensus on "length" vs "len" vs "weight" vs ?
            # length_sq = edge_attrs.get('weight')     # TODO: Reach consensus on "length" vs "len" vs "weight" vs ?
            interaction = edge_attrs.get('interaction')
            length, length_sq, stiffness = determine_end_end_distance(source, target, edge_attrs, interaction)

            # If stiffness of the current path element is different from the last, then
            # add the old stiffness group to path_edges, and
            # create a new stiffness group to put this element in...
            if stiffness != last_stiffness:
                if last_stiffness is not None:
                    path_edges.append((last_stiffness, stiffness_group))
                stiffness_group = []
                last_stiffness = stiffness
            # stiffness group is a group of elements with the same stiffness type,
            # path_edges is a list of (stiffness, stiffness groups) tuples, e.g.
            #   [('b', [elements with ss backbone interaction]), ('s', [stacked elms]), ('b', [backbone elms]), ...]
            if length_only:
                if length_only == 'sq':
                    stiffness_group.append(length_sq)
                elif length_only == 'both':
                    stiffness_group.append((length, length_sq))
                else:
                    stiffness_group.append(length)
            else:
                stiffness_group.append((length, length_sq, source, target))
        # end iteration over all path elements

        # Append the last stiffness group to path_edges:
        path_edges.append((stiffness, stiffness_group))

        # printd("Path edges: %s" % path_edges)
        #printd("stiffness groups: %s" % stiffness_group)

        if summarize and length_only:
            if length_only == 'both':
                # Return a list of (stiffness, (length, length_squared)) tuples:
                return [(stiffness, [sum(lengths) for lengths in zip(*lengths_tup)]) # pylint: disable=W0142
                        for stiffness, lengths_tup in path_edges]
            else:
                return [(stiffness, sum(lengths)) for stiffness, lengths in path_edges]
        return path_edges


    # def domain_path_elements(self, path):
    #     """
    #     (CURRENTLY NOT USED)
    #     Returns a list of structural elements based on a domain-level path (list of domains).
    #     """
    #     #path_set = set(path)
    #     remaining_domains = path[:] # deque(path)
    #     elements = [] # list of structural elements on path
    #     while remaining_domains:
    #         domain = remaining_domains.pop(0) # .popleft() # use popleft if using a deque
    #         if not domain.partner:
    #             elem = SingleStrand(domain)
    #         else:
    #             elem = DsHelix(domain)
    #             # Determine if DsHelix is actually a helix bundle:
    #             # Uhm...
    #             # Would it be better to have a "complete" helix structure description of the complex
    #             # at all time?
    #             # Whatever... for now
    #         elements.append(elem)
    #         if not remaining_domains:
    #             break
    #         i = 0
    #         # while i < len(remaining_domains) and remaining_domains[i] in elem.domains:
    #         for i, domain in remaining_domains:
    #             if domain not in elem.domains:
    #                 break
    #         else:
    #             remaining_domains = []
    #         if i > 0:
    #             remaining_domains = remaining_domains[i:]



    ## KEPT HERE BECAUSE OF ALL THE NOTES ETC:
    # def intracomplex_activity_old(self, elem1, elem2, reaction_type=HYBRIDIZATION_INTERACTION):
    #     r"""
    #     Returns
    #         :intracomplex_activity:
    #     between domain1 and domain2, so that
    #         c_j = k_j * intracomplex_activity
    #     The intracomplex activity is basically just:
    #         activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
    #     where NA is Avogadro's constant, 6.022e23/mol.
    #
    #     The activity has the same value as the unitless (P_loop/P_v0) just multiplied with "× M" to get unit of M.
    #     Thus, the activity returned by this function can be interpreted as a relative probability ratio
    #     denoting the probability that two domains/reactants will be in sufficient proximity to react,
    #     relative to the probability of two reactants confined within an Avogadro volume, v0 will react, with
    #         v0 = 1/(NA M) = 1/(6.022e23 mol⁻¹ mol L⁻¹) = 1.6e-27 m³ = 1.6 nm³
    #     Two molecules confined in v0 is equivalent to a solution of reactants with a concentration of 1 M
    #     (i.e. standard conditions and a standard activity of 1). The reason we want to specify a relative collision
    #     probability / activity (instead of an absolute) is that determining reaction probabilities from
    #     first principles is very hard. Calculating a relative probability/activity allows us to use empirical data
    #     such as k_on, and energies (ΔH, ΔS) at standard conditions.
    #
    #     For a walkthrough of this argument, see Dannenberger et al, 2015.
    #
    #     In addition to determining the stochastic rate constant, activity (or volume or inverse volume)
    #     can also be used to calculate the loop energy:
    #     dG  = R T ln(effective_volume/avogadro_volume)      # avogadro_volume = 1/(NA × M) = 1.6 nm³
    #         = R T ln((NA × M)/molecular_activity)           # molecular_activity = 1/effective_volume
    #         = - R T ln(loop_activity × M⁻¹)                 # activity = 1/(effective_volume × NA)
    #
    #     Initially, I considered returning inverse_volume (L-1) or mean-root-square-radius (m2), but I think
    #     activity is the most straight-forward to use and it can always be converted to either of the other values.
    #
    #     Regarding names:
    #     - "activity" -  should be have unit of M=mol/L...
    #                     although at the single molecule level it could be argued to be just 1/L.
    #     - volume    -   Unit of L, obviously. (Or maybe m³ or nm³.)
    #     - probability   Should be unit-less. Although we can just say that we return a relative probability factor
    #                     so that c_j = k_j × rel_prob_factor × M.
    #
    #     Alternative names:
    #     - loop_activity
    #     - intracomplex_activity             [could be return value as 1/L or mol/L=M]
    #     - intracomplex_stochastic_activity
    #     - molar_collision_probability  (except it is relative)
    #     - spatial_overlap_factor
    #     - localization_cost  (except we usually associate "cost" with energy)
    #
    #     ## IMPLEMENTATION: ##
    #
    #     a: Flow-chart style:
    #     1. Are the two domains connected by a single ds helix?
    #     2. (...)
    #
    #     b: More direct approach:
    #     1. Determine one (or multiple?) path(s) connecting domain 1 and 2.
    #     2. Determine the structural elements that make up this path:
    #         Flexible, single-stranded connections.
    #         Semi-rigid double-stranded helices.
    #         Rigid, multi-helix bundles.
    #          * Question: How about stacked/rigid interface between a ds-helix and multi-helix bundle?
    #             - Consider as a single, hard-to-calculate element?
    #     4. Determine the length of longest rigid element (LRE),
    #         and the length of all other elements plus half of the length of domain1/2.
    #         (sum remaining elements, SRE)
    #     5. If the SRE is as long or longer than the LRE, then the domains can immediately hybridize,
    #         under loop energy EX
    #     6. If the SRE is shorter than the LRE, but the LRE is a single ds helix, then the
    #         domains can hybridize under helix bending energy EB.
    #     7. Otherwise, if the LRE is longer than the SRE and the LRE is a rigid multi-bundle element,
    #         then for now we assume that the domains cannot bind.
    #
    #     ## Regarding "effective volume" vs P_loop/P_v0: ##
    #     Using "effective volume" is more intuitive than P^rc_loop/P^rc_v0.
    #     However, stricly speaking we ARE considering reactant proximity/localization/collision probability
    #     relative to standard conditions (1 M) using probability distribution function of the loop-connected reactants.
    #     This happens to reduce to something that looks like P_loop/P_v0 = v_0/v_eff,
    #     where we can use the mean squared end-to-end distance (mean_sq_ee_dist aka E[r²]) to calculate v_eff.
    #     However, this is just an approximation valid only for E[r²] >> rc (rc = critical interaction distance),
    #     which happens to reduce to a simple expression for P_loop and an expression for P_loop/P_v0 that does not
    #     depend on rc. In general, especially if rc is significant, P_loop/P_v0 could very well depend on rc!
    #     There is also no general guarantee that the expression for P_rel = P_loop/P_v0 will reduce to something with
    #     an easily-identifiable v_eff subexpression. In that case we can simply define v_eff as
    #       v_eff = v_0 * P_v0(rc) / P_loop(rc) = v_0 / P_rel(rc)       # P_rel < 1 (unless they are really close)
    #     although in practice we would not need to calculate v_eff, we would just use
    #       activity = P_v0(rc) / P_loop(rc) × M
    #
    #     Regarding mixed-level optimization (APPROXIMATION!):
    #     1. Get the path at the strand level
    #     2. Get domain-level subgraph only for domains with strand in the strand-level path
    #     3. Get domain-level shortest path using the subgraph from (2).
    #     4. Get 5p3p-level subgraph with 5p3p ends whose domain is in the domain-level shortest path.
    #     5. Get 5p3p-level shortest path using the subgraph from (4).
    #     Critizism:
    #     * May be far from the "proper" 5p3p shortest path.
    #     * The strand-level graph has a hard time evaluating edge distances (weights).
    #         For instance, all staples are connected to the scaffold. What is the distance between two staples?
    #
    #     Edit: I really feel the only way to properly represent the path is to use 5p3p representation.
    #     Domain-level representation is simply not sufficient.
    #     For instance, in the "bulge" structure below, domains C and c are certainly within reach:
    #                                         _ _ _C_ _ _  3'
    #     5'------------A---------------B----/
    #     3'------------a----------
    #                              \ _ _ _c_ _ _ 5'
    #     However, the domain level representation:
    #         A -- B -- C
    #         |
    #         a
    #           \- c
    #     Is equivalent to the would-be circular structure:
    #             3'_ _ _C_ _ _
    #                          \----B----------------A----------5'
    #                                    ------------a----------
    #                                                           \ _ _ _c_ _ _ 3'
    #     Where, depending in helix Aa, domains C and c may not be within reach.
    #     The domain-level graph looses some detail.
    #     This *can* be recovered via edge attributes, e.g. edge A-B can be directed or otherwise
    #     inform that B is on the same side of A as c is. But that might be just as much work
    #     as just using the 5p3p-ends graph.
    #
    #     Note: Previously I also intended this function to determine whether domain1 and domain2 can hybridize
    #     and what the energy penalty is, i.e. loop energy, helix-bending energy, etc.
    #
    #
    #     ## GLOBAL vs LOCAL model: Using only the minimum loop (shortest path) vs all paths ##
    #
    #     Consider the following two model cases:
    #      ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
    #      |           ⁞‾‾‾‾‾‾‾‾‾‾‾               -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾
    #      |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃                      |           ⁞___________
    #      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'
    #      ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
    #      |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|  -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|
    #      |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃   ₃__________⌡       |           ⁞___________ ₃__________⌡
    #      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'
    #
    #     In both cases, when connecting (B3p-E5p), we would only consider the shortest path, (B3p-A5p-A3p-D5p-D3p-E5p).
    #     However, the second-shortest path (B3p-B5p-C3p-C5p-F3p-F5p-E3p-E5p) would clearly also have an influence
    #     on the PDF overlap (or, activity/effective_volume) of domain B and E.
    #
    #     To make the energy calculation more precise, you could do a search for secondary and tertiary loops.
    #     If these are present, that would increase the activity.
    #
    #     More refs:
    #     * https://en.wikipedia.org/wiki/Loop_entropy
    #
    #     TODO: Check for secondary loops
    #     TODX: Saving secondary loops/paths might be be useful when determining if a dehybridization will split
    #           up a complex. Although that is pretty much the reverse process and can't really be used. Nevermind.
    #
    #     TODO: How does this perform for small loops, e.g. this:     --------------------------------
    #           which is just 0.67 nm.                                ---------------.________________
    #           effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2) = 1.66 nm³
    #           activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3 = 1.66 nm³ / 1.66 nm³ = 1.
    #
    #     TODO: Add option to calculate shape energy contribution from ssDNA-dsDNA during duplex formation.
    #           This energy is entropic and can be viewed in two ways:
    #             (1) After duplex is formed, the effective activity is slightly different
    #                 than it was before the domain was hybridized, because the duplex length and rigidity is changed.
    #                 (long domains become longer; higher rigidity).
    #             (2) Duplex formation changes the shape energy (loop entropy) of the complex by reducing
    #                 the number of available micro-states (the conformational space).
    #           Should duplex-formation-based loop energy be applied to the off-rate or on-rate?
    #             If the duplex is actually formed before the crossover is made, then an argument could be made
    #             that the duplex-formation-shape-energy should be applied to the the on-rate (k_on).
    #             However, this argument is not very convincing: The duplex cannot be formed without at least
    #             partially closing the loop (if the furthest end is hybridized, and there is a new loop with size equal
    #             to the ss length of the domain..).
    #           Thus, it is probably most realistic to apply duplex-formation-shape-dS to the off-rate.
    #           This "post-hybridization duplex-formation-shape-dG" calculation could go hand-in-hand with a general
    #           "post-reaction-processing" where we determine not only adjustments to k_off but also:
    #             * Does the reaction (hybridization/stacking) produce an impossible structure?
    #                 E.g. if
    #
    #
    #     """
    #     #path = self.domains_shortest_path(domain1, domain2)
    #     #path_segments = self.domain_path_elements(path)
    #     # NOTE: The path does not have to span the full length of the element!
    #     # The element could be really long, e.g. a long DsHelix or the full ss scaffold,
    #     # with the path only transversing a fraction of the element.
    #     # Also, for domain-level helices, it is hard to know if the path traverses
    #     # the domain, or just uses it for the hybridization, traversed at one end only.
    #     #      To here
    #     #         |
    #     # -------´ ---------
    #     # ------------------ ------from here
    #     #LRE_len, LRE = max((elem.length_nm, elem) for elem in path_segments if not isinstance(elem, SingleStrand))
    #
    #     ## 5p3p-level shortest path:
    #     # domain1, domain2
    #     if isinstance(elem1, Domain):
    #         """
    #         For domain-hybridization, what path exactly are we considering?
    #         - depends on zippering vs ??
    #
    #         """
    #         domain1, domain2 = elem1, elem2
    #         cmplx = domain1.strand.complex
    #         d1end5p, d2end5p = domain1.end5p, domain2.end5p
    #         d1end3p, d2end3p = domain1.end3p, domain2.end3p
    #         path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
    #         slice_start = 1 if path[1] == d1end3p else 0
    #         slice_end = -1 if path[-2] == d2end3p else None
    #         if slice_start == 1 or slice_end is not None:
    #             path = path[slice_start:slice_end]
    #     elif isinstance(elem1, DomainEnd):
    #         d1end5p, d2end5p = elem1, elem2
    #         cmplx = d1end5p.domain.strand.complex
    #         # domain1, domain2 = d1end5p.domain, d2end5p.domain # only used for debugging
    #         path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
    #         slice_start, slice_end = 0, None
    #     else:
    #         # Assume stacking interaction:
    #         # TODO: Consolidate this for domain hybridization and end stacking
    #         (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
    #         cmplx = h1end5p.domain.strand.complex
    #         # import pdb; pdb.set_trace()
    #         path = self.ends5p3p_shortest_path(h1end3p, h2end3p)
    #         slice_start = 1 if path[1] == h2end5p else 0
    #         slice_end = -1 if path[-2] == h1end5p else None
    #         if slice_start == 1 or slice_end is not None:
    #             path = path[slice_start:slice_end]
    #     ## TODO: Check if d1end3p or d2end3p is in the shortest_path and adjust accordingly.
    #     # printd("5p3p graph shortest path from %s to %s:" % (d1end5p, d2end5p))
    #     # pprintd(path)
    #     path_elements = self.ends5p3p_path_partial_elements(path, length_only='both', summarize=True)
    #     # printd("5p3p graph path elements:")
    #     # pprintd(path_elements)
    #     # import pdb; pdb.set_trace()
    #
    #     ## TODO: CHECK FOR SECONDARY LOOPS!
    #     ## We can easily do this simply by whether the two nodes are already in an existing loop.
    #     if path[0] in cmplx.loops_by_domainend and len(cmplx.loops_by_domainend[path[0]]) > 0 \
    #         and path[-1] in cmplx.loops_by_domainend and len(cmplx.loops_by_domainend[path[-1]]) > 0:
    #         # We have a secondary loop.
    #         shared_loops = cmplx.loops_by_domainend[path[0]] & cmplx.loops_by_domainend[path[-1]]
    #         assert len(shared_loops) > 0
    #         # Can we always use the smallest loop, or is there a particular order that must be obeyed?
    #         smallest_loop = min(shared_loops, key=lambda loop: len(loop['nodes']))
    #         shortest_path_nodes = set(path)
    #         assert all(node in smallest_loop['nodes'] for node in path)
    #         assert shortest_path_nodes <= smallest_loop['nodes']  # Using set operations is faster
    #         # Set operators: & = intersection, - = difference, ^ = symmetric_difference, <= subset, >= superset.
    #         path2_nodes = smallest_loop['nodes'] - shortest_path_nodes
    #
    #     # path_elements is a list of tuples as: [(stiffness, (total-length, sum-of-squared-lengths)), ...]
    #     # For rigid, double-helical elements, element length, l, is N_bp * 0.34 nm.
    #     # For single-stranded elements, we estimate the end-to-end distance by splitting the strand into
    #     #   N = N_nt*0.6nm/1.8nm segments, each segment being the Kuhn length 1.8 nm of ssDNA,
    #     # where 0.6 nm is the contour length of ssDNA and N_nt is the number of nucleotides in the strand.
    #     #   E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (1.8 nm)²
    #     #         = round(N_nt*0.6nm/1.8nm) (1.8 nm)² = N_nt * 0.6/1.8 * 1.8*1.8 nm² = N_nt 0.6*1.8 nm²
    #     #         = N_nt * lˢˢ * λˢˢ = N_nt * 1.08 nm²
    #     # Why use stiffness as first maximum criteria?? Shouldn't it just be stiffness > 0? Meh.
    #     try:
    #         # LRE: "Longest Rigid Element"; "SRE": "Sum of Remaining Elements".
    #         _, LRE_len, LRE_len_sq, LRE_idx = max((stiffness, elem_length, elem_len_sq, i)
    #                                               for i, (stiffness, (elem_length, elem_len_sq))
    #                                               in enumerate(path_elements)
    #                                               if stiffness > 0)
    #         if LRE_len <= HELIX_WIDTH:
    #             # LRE is a helix; this is within the "interaction radius" and is ignored.
    #             # TODO: This is really just a quick hack; more thought should go into this.
    #             LRE_len, LRE_len_sq = 0, 0
    #     except ValueError:
    #         # No stiff elements:
    #         LRE_len, LRE_len_sq = 0, 0
    #         SRE_len = sum(elem_length for (stiffness, (elem_length, elem_len_sq)) in path_elements if stiffness == 0)
    #         SRE_len_sq = sum(elem_len_sq for (stiffness, (elem_length, elem_len_sq)) in path_elements if stiffness == 0)
    #     else:
    #         #LRE = path_elements[LRE_idx] # .pop(LRE_idx)
    #         # Exclude LRE when calculating SRE length:
    #         if len(path_elements) > 1:
    #             # SRE_lengths, SRE_sq_lengths = ...
    #             SRE_lengths_and_sq = [(elem_length, elem_len_sq)
    #                                   for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
    #                                   for stiffness, (elem_length, elem_len_sq) in sub_path]
    #             # SRE_len_sq = sum(elem_len_sq for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
    #             #               for stiffness, (elem_length, elem_len_sq) in sub_path)
    #             ## HACK: If path is [ss-backbone-connection, duplex-hybridization-width, ss-backbone-connection] then
    #             ## we (probably?) have two domains on the same helix connected by hybridization of neighboring domains.
    #             ## In that case, we need to get an activity of 1, either by increasing the interaction radius, or
    #             ## by decreasing the path. For now, I will just decrease the path.
    #             ## NOTE that this may give a very inacurate activity for hair-pins!
    #             if len(SRE_lengths_and_sq) > 1:
    #                 if SRE_lengths_and_sq[-1][0] == HELIX_XOVER_DIST:
    #                     SRE_lengths_and_sq.pop()
    #                 elif SRE_lengths_and_sq[0] == HELIX_XOVER_DIST:
    #                     SRE_lengths_and_sq.pop(0)
    #             SRE_len = sum(tup[0] for tup in SRE_lengths_and_sq)
    #             SRE_len_sq = sum(tup[1] for tup in SRE_lengths_and_sq)
    #         else:
    #             # SRE_lengths, SRE_sq_lengths = [], []
    #             SRE_len, SRE_len_sq = 0, 0
    #
    #     # Comparing mean-end-to-end-squared values vs mean end-to-end lengths vs full contour length?
    #     # There is a difference that sum of squares does not equal the square of sums, so
    #     # even if LRE_len_sq is > SRE_len_sq, LRE_len could be less than SRE_len.
    #     # Also, while ds duplexes has full contour length equal to mean end-to-end length,
    #     # this is not true for ssDNA. Indeed if calculating whether two domains *can* reach,
    #     # it is probably better to use domain.ds_length.
    #
    #     if LRE_len > SRE_len:
    #         # printd("LRE_len > SRE_len for path between %r and %r" % (domain1, domain2))
    #         # pprintd(path)
    #         # pprintd(path_elements)
    #         # The domains cannot reach each other.
    #         # Hybridization requires helical bending; Not implemented yet; just returning 0 meaning "impossible".
    #         # TODO: Implement hybridization via helical bending.
    #         #   Persistance length 50 nm (physiological salt)
    #         #   - Depends on ionic strength and cationic valency
    #         # TODO: Look at formulas for k_on and k_off rates under stress.
    #         #   For DNA, there is certainly a difference between axial "ripping" and perpendicular "zipping".
    #         #   - Zippering occours at about 10-15 pN (sequence dependent).
    #         #   -
    #         return 0
    #     ## There is probably some profound relation between the elements and the gamma factor.
    #     ## E.g. if the contour length is long enough for the domains to reach, but the
    #     ## SRE mean squared end-to-end distance is less than the LRE, then the SRE will rarely
    #     ## be sufficiently extended for the domains to hybridize. This decrease in spatial pdf
    #     ## can be considered equivalent to an increase in effective volume.
    #     ## Another, more complex case, is when the SRE has only (a) one, or (b) a few links,
    #     ## in which case the mean squared end-to-end distance is not a good measure of spatial pdf.
    #     ## In the case where SRE is a rigid 1-element chain of same length as the LRE, the pdf
    #     ## is essentially a sphere centered at the joint between LRE and SRE. (Similar case when
    #     ## the LRE is flanked by two rigid elements.)
    #
    #     # Example 1: SRE_len = 10 * 4 nm = 40 nm; SRE_len_sq = 10 * (4 nm)**2 = 160 nm2.
    #     #            LRE_len = 20 nm,             LRE_len_sq = (20 nm)**2 = 400 nm2.
    #     #            LRE_len < SRE_len, but LRE_len_sq > SRE_len_sq
    #     #            SRE_len/LRE_len = 2 -- higher => lower gamma_corr.
    #     #            LRE_len_sq/SRE_len_sq = 2.5 -- higher => higher gamma_corr.
    #     # We could, for instance, say: gamma_corr = 1 + ln(LRE_len_sq/SRE_len_sq)
    #     # Hmm... probably need to do some further analysis of different examples and try to figure out
    #     # a proper relationship between link-elements and gamma_corr... And it might not be as simple
    #     # as a simple exponential correction to (P_loop/P_v0).
    #
    #     # If LRE_len_sq > SRE_len_sq, then the approximation assumption "we many links of length l_i"
    #     # is certainly not valid (we only have 1 link of length LRE_len).
    #     # Instead of considering P_loop(r<rc)/P_v0(r<rc), we have to consider P_SRE(r=LRE+/-rc)/V(r=LRE)/P_v0(r<rc).
    #     # that is, the probability of the SRE end-end distance equaling LRE length, normalized by the shell
    #     # volume at r=LRE.
    #     # This gives us a factor that seems to be:
    #     # 1/(4 π LRE_len_sq) exp(-3*LRE_len_sq / (2*SRE_len_sq))
    #     # although the first part would give us a non-unitless factor which is not acceptable. It probably has to be
    #     # normalized in some way, but for now just use the exponential part.
    #     # Edit: Actually, the first part is probably something like LRE_len/rc
    #     #
    #     # For LRE_len_sq == SRE_len_sq, this gives us exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 1/e = 0.22.
    #     # For example 1, this will give us:
    #     # exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 0.02.
    #     LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq)) if LRE_len_sq > 0 else 1
    #
    #     gamma_corr = 1
    #     if LRE_len_sq > SRE_len_sq:
    #         # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
    #         # Should probably be approximated with some continuous function.
    #         # gamma = (3/2)*gamma_corr;
    #         gamma_corr += 0.5
    #     # Mean end-to-end squared distance between the two domains, aka E_r_sq:
    #     # Mean end-to-end squared distance.
    #     # We already have the squared length, Nᵢbᵢ², so we just need to sum:
    #
    #     mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm
    #     if mean_sq_ee_dist <= HELIX_XOVER_DIST**2:
    #         # Ensure that activity is not too high:
    #         return 1
    #
    #     ## Regarding "effective volume" vs P_loop/P_v0: ##
    #     # Using "effective volume" is more intuitive than P^rc_loop/P^rc_v0.
    #     # However, stricly speaking we ARE considering reactant proximity/localization/collision probability
    #     # relative to standard conditions (1 M) using probability distribution function of the loop-connected reactants.
    #     # This happens to reduce to something that looks like P_loop/P_v0 = v_0/v_eff,
    #     # where we can use the mean squared end-to-end distance (mean_sq_ee_dist aka E[r²]) to calculate v_eff.
    #     # However, this is just an approximation valid only for E[r²] >> rc (rc = critical interaction distance),
    #     # which happens to reduce to a simple expression for P_loop and an expression for P_loop/P_v0 that does not
    #     # depend on rc. In general, especially if rc is significant, P_loop/P_v0 could very well depend on rc!
    #     # There is also no general guarantee that the expression for P_rel = P_loop/P_v0 will reduce to something with
    #     # an easily-identifiable v_eff subexpression. In that case we can simply define v_eff as
    #     #   v_eff = v_0 * P_v0(rc) / P_loop(rc) = v_0 / P_rel(rc)       # P_rel < 1 (unless they are really close)
    #     # although in practice we would not need to calculate v_eff, we would just use
    #     #   activity = P_v0(rc) / P_loop(rc) × M
    #
    #     effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2)
    #     #effective_volume = (2/3*math.pi*mean_sq_ee_dist)**gamma * 1e-24 # 1e-24 to convert nm3 to L.
    #     #activity = (1/N_AVOGADRO)*(1/effective_volume)
    #     # Using AVOGADRO_VOLUME_NM3 to avoid redundant conversions:
    #     activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
    #     if gamma_corr > 1:
    #         activity = activity**gamma_corr
    #
    #     ## When/where to apply extra gamma? ##
    #     # Note: The extra gamma really should be applied to the activity, not just the effective volume:
    #     # That is, we would always use exponent of 3/2 for calculating effective volume from mean_sq_ee_dist,
    #     # and apply the extra gamma to the unitless activity (v0/effective_volume) as:
    #     #   activity = (v0/effective_volume) ** gamma_corr
    #     # where gamma_corr = 1+x for gamma_org = 3/2+x
    #     # This is also what is effectively done in (Dannenberger, 2015), where γ is outside the expression
    #     #     ΔG = - R T γ ln(C/E[r2])
    #     # TODO: Check that this is also how Dannenberger actually does it when calculating k₊ in the java code.
    #     # To keep ΔG "reasonable" for large E[r2], (Dannenberger et al, 2015) adjusts C approximately as:
    #     #   C = 2.2e-18 m² γ - 2.7e-18 m² = 3.34 C0 γ - 4.0 C0,  C(γ=1.5) = C0
    #     # Where C0 = 3/(2π) v0**(2/3) = 6.7e-19 m²
    #     # This corresponds to increasing Avogadro volume, v0.
    #
    #     # Note: The "activity" as calculated above appears on paper to be a unitless ratio (P_loop/P_v0).
    #     # However, since it is relative to standard conditions (1 M), we just have to implicitly multiply
    #     # the unitless ratio with "× 1 M" to get a proper molar activity.
    #
    #     # TODO: Currently not accounting for bending or zipping energy.
    #     # TODO: Account for secondary (and tertiary?) loops
    #     if activity < 1:
    #         print("Activity %s less than 1 - printing locals():" % activity)
    #         pprint(locals())
    #
    #     return activity



    def draw_and_save_graphs(self, directory=None, fnfmt="{prefix}{graph}_{n}.png",
                             n=1, prefix=None,
                             layout="graphviz", apply_attrs=True):
        # sysgraphs = {"strand_graph": self.strand_graph,
        #              "domain_graph": self.domain_graph,
        #              "ends5p3p_graph": self.ends5p3p_graph}
        # sysgraphs.update({"cmplx-%s" % c.cuid: c for c in self.complexes})
        sysgraphs = {"domain_graph": self.domain_graph}
        prefix = prefix if prefix is not None else self.fnprefix
        ## As of python 3.5 you can do:
        ## allgraphs = {**sysgraphs, **cmplxgraphs}
        if directory is None:
            directory = self.params.get("working_directory", ".")
        for graph_name, graph in sysgraphs.items():
            if len(graph) < 1:
                print("Graph %s contains %s nodes, skipping..." % (graph_name, len(graph)))
                continue
            print("Plotting graph %s ..." % graph_name)
            fn = os.path.join(directory, fnfmt.format(graph=graph_name, n=n, prefix=prefix))
            if apply_attrs:
                nodes, node_attrs = list(zip(*graph.nodes(data=True)))
                ## For some reason, nx_pylab will not allow mixed tuples and strings as colors:
                node_colors = [at.get('_color', (1.0, 0.0, 0.0)) for at in node_attrs]
                if len(graph.edges()) > 0:
                    edges, edge_attrs = list(zip(*[((u, v), data) for u, v, data in graph.edges(data=True)]))
                    edge_colors = [at.get('_color', (0.0, 0.0, 0.0)) for at in edge_attrs]
                else:
                    edges = graph.edges()
                    edge_attrs = None
                    edge_colors = None
                node_labels = {node: node.name for node in graph.nodes()}
                # print("node_attrs:")
                # pprint(node_attrs)
                # print("edge_attrs:")
                # pprint(edge_attrs)
                # print("node colors:")
                # pprint(node_colors)
                # print("edge colors:")
                # pprint(edge_colors)
                ## NOTE: GraphViz doesn't play nice with unknown colors.
                ## If a node/edge has an unknown color, execution will fail.
                ## (which is why I've switched to '_color' for now...)
                draw_graph_and_save(graph, outputfn=fn, layout=layout,
                                    labels=node_labels,
                                    nodes=nodes,
                                    edges=edges,
                                    node_color=node_colors,
                                    edge_color=edge_colors,
                                    alpha=0.7)
            else:
                draw_graph_and_save(graph, outputfn=fn, layout=layout,
                                    labels=node_labels, alpha=0.7)
