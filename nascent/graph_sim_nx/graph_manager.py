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
import os
import networkx as nx
from networkx.algorithms.shortest_paths import shortest_path
# import numpy as np
import math
import itertools
from itertools import chain, groupby
from pprint import pprint
from operator import itemgetter, mul
from functools import reduce
import pdb

from .system_graphs import InterfaceGraph
from .nx_utils import draw_graph_and_save
from .structural_elements.strand import SingleStrand
from .structural_elements.helix import DoubleHelix
from .structural_elements.bundle import HelixBundle
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        N_AVOGADRO, AVOGADRO_VOLUME_NM3,
                        HELIX_XOVER_DIST, HELIX_WIDTH, HELIX_STACKING_DIST)
from .debug import printd, pprintd
from .domain import Domain, DomainEnd


## Module constants:
SEGMENT_STIFFNESS, SEGMENT_LENGTHS, SEGMENT_CONTOUR_LENGTH, SEGMENT_LENGTH_SQUARED = 0, 1, 0, 1

def prod(sequence):
    """ Return the product of all elements in sequence. """
    return reduce(mul, sequence)



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


    def group_interfaces_path_by_stiffness(self, path):
        """
        Returns a list of structural elements based on a interface-level path (list of InterfaceGraph nodes),
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
        path_tuples = ((edge_attrs['dist_ee_nm'], edge_attrs['dist_ee_sq'], edge_attrs['stiffness'], source, target)
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



    def domains_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        """
        return shortest_path(self.domain_graph, domain1, domain2)

    def interfaces_shortest_path(self, ifnode1, ifnode2):
        """
        TODO: This should certainly be cached.
        """
        if isinstance(ifnode1, DomainEnd):
            ifnode1, ifnode2 = ifnode1.ifnode.top_delegate(), ifnode2.ifnode.top_delegate()
        return shortest_path(self.interface_graph, ifnode1, ifnode2)

    def ends5p3p_shortest_path(self, end5p3p1, end5p3p2):
        """ :end5p3p1:, :end5p3p2: DomainEnd nodes (either End5p or End3p),
        calculates path from node1 to node2 using ends5p3p_graph. The returned path is a list
        starting with node1 and ending with node2, and has the form: [node1, <all nodes in between, ...>, node2]
        TODO: This should certainly be cached.
        TODO: Verify shortest path for end3p as well?
        """
        return shortest_path(self.ends5p3p_graph, end5p3p1, end5p3p2, weight='dist_ee_sq')


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

        # TODO: Rename to "path_summary" or "path_segments":
        path_segments = self.interfaces_path_segments(loop_path)
        # path_summary is a list of (stiffness, [length, sum_length_squared]) tuples
        # where stiffness is (0 for ssDNA, 1 for dsDNA and 2 for helix-bundles).

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
                    first_segment[1][1] += last_segment[1][1]   # Add segment squared lengths
                # else:
                    # If the duplexes are not connected via their own stiff/duplexed domains, then downstream
                    # calculation using LRE/SRE should produce correct result...
            else:
                # Only a single path segment:
                if first_segment[SEGMENT_STIFFNESS] > 0:
                    # We have a single, fully-stacked segment; this cannot ever stack back upon it self.
                    return 0
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
        # LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq)) if LRE_len_sq > 0 else 1

        gamma = 3/2
        gamma_corr = 1
        if LRE_len_sq > SRE_len_sq:
            # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
            # Should probably be approximated with some continuous function.
            # gamma = (3/2)*gamma_corr;
            gamma_corr += 0.5
            gamma = 5/3

        ## Calculate mean end-to-end squared distance between the two domains, aka E_r_sq:
        # We already have the squared length, Nᵢbᵢ², so we just need to sum LRE and SRE:
        ## EDIT, no, wait - if we have stacked, stiff double-helices/helix-bundles, then we CANNOT
        ## just add the squared lengths together. We have to square the entire element length.
        ## We can do it for ss-helices, since the "squared length" is already just N_nt*ss_rise_per_nt*ss_kuhn_length,
        ## but that isn't valid for fully-stacked continuous segments of ds-helices.
        ## Done: Re-calculate SRE_len_sq based on the summarized segment lengths, not the path edges.


        mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm
        # There shouldn't be any need to test for if mean_sq_ee_dist <= HELIX_XOVER_DIST**2 in the InterfaceGraph

        ## Note: "effective_volume" is just an informal term for P_loop(rc) / P_v0(rc) x AVOGADRO_VOLUME_NM3
        effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**gamma
        activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
        if gamma_corr > 1:
            activity = activity**gamma_corr

        return activity




    def intracomplex_activity(self, elem1, elem2, reaction_type=HYBRIDIZATION_INTERACTION):
        r"""
        Returns
            :intracomplex_activity:
        between domain1 and domain2, so that
            c_j = k_j * intracomplex_activity
        The intracomplex activity is basically just:
            activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
        where NA is Avogadro's constant, 6.022e23/mol.

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
            ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
            path = shortest_path(self.interface_graph, d1end5p.ifnode, d2end5p.ifnode)
            slice_start = 1 if path[1] == d1end3p.ifnode else 0
            slice_end = -1 if path[-2] == d2end3p.ifnode else None
            # pdb.set_trace()
            if slice_start == 1 or slice_end is not None:
                path = path[slice_start:slice_end]
        elif reaction_type is STACKING_INTERACTION:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
            cmplx = h1end5p.domain.strand.complex
            if h1end3p.domain.name in ("e", "C"):
                pdb.set_trace()
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
            dh1end3p = d1end3p = h1end3p
            dh1end5p = d2end5p = h2end5p
            dh2end3p = d3end3p = h2end3p
            dh2end5p = d4end5p = h1end5p
            # and elem1, elem2 would then be:
            assert (dh1end3p, dh1end5p), (dh2end3p, dh2end5p) == (elem1, elem2) # double-helix ends
            assert (d1end3p, d2end5p), (d3end3p, d4end5p) == (elem1, elem2)     # domain ends

            dh1_delegate, dh2_delegate = h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate()
            ## DETECT IF THE TWO DOMAINS ARE ALREADY PART OF AN EXISTING LOOP:
            assert dh1_delegate, dh2_delegate == (h2end5p.ifnode.top_delegate(), h1end5p.ifnode.top_delegate())

            path = shortest_path(self.interface_graph, dh1_delegate, dh2_delegate)

            ## Adjusting path shouldn't be required when using InterfaceGraph/Nodes:
            # slice_start = 1 if path[1] == h2end5p else 0
            # slice_end = -1 if path[-2] == h1end5p else None
            # if slice_start == 1 or slice_end is not None:
            #     path = path[slice_start:slice_end]
            # However, you MUST put the start and end nodes together as they would be when stacked,
            # (unless you choose to resolve this by increasing unstacking rate (together with throttle)..
            # - This is done AFTER getting the path segments/elements...
        else:
            assert isinstance(elem1, DomainEnd) and isinstance(elem2, DomainEnd)
            d1end5p, d2end5p = elem1, elem2
            cmplx = d1end5p.domain.strand.complex
            # domain1, domain2 = d1end5p.domain, d2end5p.domain # only used for debugging
            # path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
            path = shortest_path(self.interface_graph,
                                 d1end5p.ifnode.top_delegate(),
                                 d2end5p.ifnode.top_delegate())
            slice_start, slice_end = 0, None


        ## Useful NetworkX functions:
        ## networkx.algorithms.cycles.cycle_basis


        #### Checking for secondary side-effects: ####
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
            dh1_nodes_in_loops = dh1_ifnodes.intersection(cmplx.loopids_by_interface.keys())
            dh2_nodes_in_loops = dh2_ifnodes.intersection(cmplx.loopids_by_interface.keys())
            # Edit: Loops are currently stored as mutable dicts, you cannot make a set of loops.
            # You must either use loop IDs or create an object for each loop.
            dh1_loopids = {loopid for ifnode in dh1_nodes_in_loops if ifnode in cmplx.loopids_by_interface
                           for loopid in cmplx.loopids_by_interface[ifnode]}
            dh2_loopids = {loopid for ifnode in dh2_nodes_in_loops if ifnode in cmplx.loopids_by_interface
                           for loopid in cmplx.loopids_by_interface[ifnode]}
            # Alternatively, go through all loops and check if any of the loop's nodes are on the doublehelix:
            # dh1_loops = {loopid for loopid in cmplx.loops.keys() if any(...) }
            ## NOTE: We should not consider loops that are only on one of the dh arms and not the other,
            ## hence the two set comprehensions which are then followed by the intersection.
            ## Affected loops are loops which goes through both double-helix 1 and 2:
            affected_loops = list(dh1_loopids & dh2_loopids)
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
            new_loops = {}   # old_id => [new_loop1_path, new_loop2_path, ...]
            side_effects_factors = []
            for loopid in affected_loops:
                loop = cmplx.loops[loopid]
                old_path = loop["path"]
                new_paths = []
                if dh1_delegate in loop["path"] and dh2_delegate in loop["path"]: #and \
                    # abs(loop["path"].index(dh1_delegate) - loop["path"].index(dh2_delegate)) == 1:
                    # We can just modify the path, joining it together at the new junction:
                    idx1, idx2 = loop["path"].index(dh1_delegate), loop["path"].index(dh2_delegate)
                    if idx1 > idx2:
                        idx1, idx2 = idx2, idx1
                    new_path1 = loop["path"][idx2:] + loop["path"][:idx1+1]  # +1 to include the node.
                    new_paths.append(new_path1)
                    if idx2-idx1 == 1:
                        # Case 1a: Special case: The two joined interface nodes are adjacent; just re-calculate current loop.
                        # No extra paths
                        pass
                    else:
                        # Case 1b: Split the existing path into two new loops:
                        # We have already split out one of the loops, just need the other:
                        new_path2 = loop["path"][idx1:] + loop["path"][:idx2+1]
                        new_paths.append(new_path2)
                else:
                    # Case 2: Stacking forms a genuine new loop paths:
                    # (I think this will also take care of case 1b and even 1a)
                    # We have to generate two new loops, partially composed of the old path
                    # and partially of the newly-stacked double-helix:
                    # First, split the existing loop in two on the nodes that touch the newly-stacked double-helix:
                    grouped_path = [(is_in_dh, list(group_iter)) for is_in_dh, group_iter in
                                    itertools.groupby(old_path, lambda ifnode: ifnode in dh_ifnodes)]
                    # if len(grouped_path) in (2, 3): print("Case 1a")
                    assert 4 <= len(grouped_path) <= 5 # should be either 4 or 5 elements.
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
                    for dh_before, t_group, dh_after in ((e1, t1, e2), (e2, t2, e1)):
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
                        new_paths.append(arm1[:arm1_idx+1] + t_group + arm2[:arm2_idx+1:-1])
                new_loops[loopid] = new_paths


                    #new_path = shortest_path(dh1_delegate,
                new_activity = self.calculate_loop_activity(loop["path"], simulate_reaction=reaction_type)
                side_effects_factors.append(new_activity/loop["activity"])

            side_effects_factors = [(self.calculate_loop_activity(loop["path"], simulate_reaction=reaction_type)/
                                     loop["activity"]) for loop in affected_loops]
            if any(factor == 0 for factor in side_effects_factors):
                return 0
            side_effects_factor = prod(side_effects_factors)


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
        if (path[0] in cmplx.loops_by_interface and len(cmplx.loops_by_interface[path[0]]) > 0
            and path[-1] in cmplx.loops_by_interface and len(cmplx.loops_by_interface[path[-1]]) > 0):
            ## Perfect "pinching" of existing loop (rare)
            # We have a secondary loop.
            # Find loops containing both the first and last path interface nodes:
            shared_loops = cmplx.loops_by_interface[path[0]] & cmplx.loops_by_interface[path[-1]]
            assert len(shared_loops) == 1           ## Not sure about this...
            outer_loop = next(iter(shared_loops))
            outer_path = outer_loop['path'] # list of ifnodes from when the loop was made.
            # IF the reaction/node-merge splits more than ONE loop, we want to know!
            # Can we always use the smallest loop, or is there a particular order that must be obeyed?
            # Also, is the number of nodes "len(loop['nodes'])" the way to go?
            # It is not the same as the "shortest" path method..
            ## TODO: Check/fix this code:  - EDIT: Replaced by len(shared_loops) == 1 above.
            # outer_loop = min(shared_loops, key=lambda loop: len(loop['nodes']))  # NOT GUARANTEED
            path1_nodes = set(path)
            assert all(node in outer_loop['ifnodes'] for node in path)             # NOT GUARANTEED
            assert path1_nodes <= outer_loop['ifnodes']
            # Using set operations is faster (path is subset of existing outer_loop)
            # Set operators: &: intersection, -: difference), ^: symmetric_difference, <=: is subset, >=: superset.
            path2_nodes = smallest_loop['ifnodes'] - path1_nodes
            outer_path_idxs = (outer_path.index(path[0]), outer_path.index(path[-1]))
            if outer_path_idxs[0] > outer_path_idxs[1]:
                outer_path_idxs = outer_path_idxs[::-1]
            path1 = outer_path[outer_path_idxs[0]:outer_path_idxs[1]]
            path2 = outer_path[outer_path_idxs[1:]] + outer_path[outer_path_idxs[:0]]
            assert all(node in path1 for node in path) or all(node in path2 for node in path)
            a1 = self.calculate_loop_activity(path1, simulate_reaction=reaction_type)
            a2 = self.calculate_loop_activity(path2, simulate_reaction=reaction_type)
            a0 = outer_loop["actiity"]
            return a1*a2/a0
            # If hybridization, special case where we account for the length of the newly-formed duplex.
        else:
            ## Check for non-pinching loop-splitting cases: (partial loop sharing)
            ##                           .------.
            ##            .----------,--´------. \
            ##           /   Λ₁     /          \  \ e₄ (not part of the original discussion)
            ##        e₁ \      e₃ /:/   Λ₂    /   \
            ##            \         |         / e₂  \
            ##             `----.--´---------´      |
            ##                   `-----------------´
            ## This is really non-obvious, but let's give it a try...
            ## First see if there is any overlap between the shortest path and existing loops:
            path1_nodes = set(path)
            shared_nodes = path1_nodes.intersection(cmplx.loops_by_interface.keys())
            shared_loops = {cmplx.loops_by_interface[ifnode] for ifnode in shared_loops}
            ## We have a few special cases:
            ## If the new edge e₃ is completely rigid, then e₂ has no influence on Λ₁, nor does e₁ influence on Λ₂.
            ## - The "e₃ == 0" aka "loop pinching" above can be seen as a special case of this.
            ## Perhaps just simplify as:
            ## * If e₃ has fewer segments than e₁, e₂ then calculate activity as:
            ##      a = a₁*a₂/a₀, where a₁ = f(e₁+e₃)e₂, a₂ = f(e₂+e₃), and a₀ = f(e₁+e₂) is the outer loop activity.
            ## * If e₃ has more segments than e₁, e₂ then calculate activity simply using the shortest path:
            ##      a = a₁ if e₁ < e₂ else a₂
            if shared_loops:
                # Find nodes that are not part of any loops:
                unshared_nodes = path1_nodes.difference(cmplx.loops_by_interface.keys())
                # for loop in shared_loop:
            else:
                ## NO SHARED LOOPS: PATH IS NEITHER PARTIALLY OR FULLY PART OF ANY CURRENT LOOPS:
                activity = self.calculate_loop_activity(path, simulate_reaction=reaction_type)








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
    #         path1_nodes = set(path)
    #         assert all(node in smallest_loop['nodes'] for node in path)
    #         assert path1_nodes <= smallest_loop['nodes']  # Using set operations is faster
    #         # Set operators: & = intersection, - = difference, ^ = symmetric_difference, <= subset, >= superset.
    #         path2_nodes = smallest_loop['nodes'] - path1_nodes
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
