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


import os
import networkx as nx
from networkx.algorithms.shortest_paths import shortest_path
import numpy as np
import math
import itertools
from pprint import pprint


from .nx_utils import draw_graph_and_save
from .structural_elements.strand import SingleStrand
from .structural_elements.helix import DsHelix
from .structural_elements.bundle import HelixBundle
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        N_AVOGADRO, AVOGADRO_VOLUME_NM3,
                        HELIX_XOVER_DIST, HELIX_WIDTH, HELIX_STACKING_DIST)
from .debug import printd, pprintd
from .domain import Domain, DomainEnd


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



class GraphManager():

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
        # - hairpins without loops has hybridization to same domain as phosphate backbone.
        # - Not sure what the performance penalty of multi graphs are vs regular graphs?
        # - We could have a separate graphs without stacking, but not sure when those would be useful?
        self.ends5p3p_graph = ends5p3p_graph or nx.MultiGraph()
        self.domain_graph = domain_graph or nx.MultiGraph()
        self.strand_graph = strand_graph or nx.MultiGraph()
        if strands:
            # TODO: Add strand-node attrs to strand-graph:
            # TODO: Add "Label" attr to all nodes (for Gephi)
            self.strand_graph.add_nodes_from(strands)
            for strand in strands:
                self.domain_graph.add_nodes_from(strand.nodes(data=True))
                self.domain_graph.add_edges_from(strand.edges(keys=True, data=True))
                self.ends5p3p_graph.add_nodes_from(strand.ends5p3p_graph.nodes(data=True))
                self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph.edges(keys=True, data=True))

        ### Other attributes: ###
        self.fnprefix = ""



    def domain_path_elements(self, path):
        """
        (CURRENTLY NOT USED)
        Returns a list of structural elements based on a domain-level path (list of domains).
        """
        #path_set = set(path)
        remaining_domains = path[:] # deque(path)
        elements = [] # list of structural elements on path
        while remaining_domains:
            domain = remaining_domains.pop(0) # .popleft() # use popleft if using a deque
            if not domain.partner:
                elem = SingleStrand(domain)
            else:
                elem = DsHelix(domain)
                # Determine if DsHelix is actually a helix bundle:
                # Uhm...
                # Would it be better to have a "complete" helix structure description of the complex
                # at all time?
                # Whatever... for now
            elements.append(elem)
            if not remaining_domains:
                break
            i = 0
            # while i < len(remaining_domains) and remaining_domains[i] in elem.domains:
            for i, domain in remaining_domains:
                if domain not in elem.domains:
                    break
            else:
                remaining_domains = []
            if i > 0:
                remaining_domains = remaining_domains[i:]

    def ends5p3p_path_partial_elements(self, path, length_only=False, summarize=False):
        """
        Returns a list of structural elements based on a 5p3p-level path (list of 5p3p ends).

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

        printd("Path edges: %s" % path_edges)
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


    def ends5p3p_shortest_path(self, end5p3p1, end5p3p2):
        """
        TODO: This should certainly be cached.
        TODO: Verify shortest path for end3p as well?
        """
        return shortest_path(self.ends5p3p_graph, end5p3p1, end5p3p2, weight='dist_ee_sq')


    def intracomplex_activity(self, elem1, elem2):
        r"""
        Returns
            :intracomplex_activity:
        between domain1 and domain2, so that
            c_j = k_j * intracomplex_activity
        The intracomplex activity is basically just:
            activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
        where NA is Avogadro's constant, 6.022e23/mol.

        The activity has the same value as the unitless (P_loop/P_v0) just multiplied with "× M" to get unit of M.
        Thus, the activity returned by this function can be interpreted as a relative probability ratio
        denoting the probability that two domains/reactants will be in sufficient proximity to react,
        relative to the probability of two reactants confined within an Avogadro volume, v0 will react, with
            v0 = 1/(NA M) = 1/(6.022e23 mol⁻¹ mol L⁻¹) = 1.6e-27 m³ = 1.6 nm³
        Two molecules confined in v0 is equivalent to a solution of reactants with a concentration of 1 M
        (i.e. standard conditions and a standard activity of 1). The reason we want to specify a relative collision
        probability / activity (instead of an absolute) is that determining reaction probabilities from
        first principles is very hard. Calculating a relative probability/activity allows us to use empirical data
        such as k_on, and energies (ΔH, ΔS) at standard conditions.

        For a walkthrough of this argument, see Dannenberger et al, 2015.

        In addition to determining the stochastic rate constant, activity (or volume or inverse volume)
        can also be used to calculate the loop energy:
        dG  = R T ln(effective_volume/avogadro_volume)      # avogadro_volume = 1/(NA × M) = 1.6 nm³
            = R T ln((NA × M)/molecular_activity)           # molecular_activity = 1/effective_volume
            = - R T ln(loop_activity × M⁻¹)                 # activity = 1/(effective_volume × NA)

        Initially, I considered returning inverse_volume (L-1) or mean-root-square-radius (m2), but I think
        activity is the most straight-forward to use and it can always be converted to either of the other values.

        Regarding names:
        - "activity" -  should be have unit of M=mol/L...
                        although at the single molecule level it could be argued to be just 1/L.
        - volume    -   Unit of L, obviously. (Or maybe m³ or nm³.)
        - probability   Should be unit-less. Although we can just say that we return a relative probability factor
                        so that c_j = k_j × rel_prob_factor × M.

        Alternative names:
        - loop_activity
        - intracomplex_activity             [could be return value as 1/L or mol/L=M]
        - intracomplex_stochastic_activity
        - molar_collision_probability  (except it is relative)
        - spatial_overlap_factor
        - localization_cost  (except we usually associate "cost" with energy)

        ## IMPLEMENTATION: ##

        a: Flow-chart style:
        1. Are the two domains connected by a single ds helix?
        2. (...)

        b: More direct approach:
        1. Determine one (or multiple?) path(s) connecting domain 1 and 2.
        2. Determine the structural elements that make up this path:
            Flexible, single-stranded connections.
            Semi-rigid double-stranded helices.
            Rigid, multi-helix bundles.
             * Question: How about stacked/rigid interface between a ds-helix and multi-helix bundle?
                - Consider as a single, hard-to-calculate element?
        4. Determine the length of longest rigid element (LRE),
            and the length of all other elements plus half of the length of domain1/2.
            (sum remaining elements, SRE)
        5. If the SRE is as long or longer than the LRE, then the domains can immediately hybridize,
            under loop energy EX
        6. If the SRE is shorter than the LRE, but the LRE is a single ds helix, then the
            domains can hybridize under helix bending energy EB.
        7. Otherwise, if the LRE is longer than the SRE and the LRE is a rigid multi-bundle element,
            then for now we assume that the domains cannot bind.

        Regarding mixed-level optimization (APPROXIMATION!):
        1. Get the path at the strand level
        2. Get domain-level subgraph only for domains with strand in the strand-level path
        3. Get domain-level shortest path using the subgraph from (2).
        4. Get 5p3p-level subgraph with 5p3p ends whose domain is in the domain-level shortest path.
        5. Get 5p3p-level shortest path using the subgraph from (4).
        Critizism:
        * May be far from the "proper" 5p3p shortest path.
        * The strand-level graph has a hard time evaluating edge distances (weights).
            For instance, all staples are connected to the scaffold. What is the distance between two staples?

        Edit: I really feel the only way to properly represent the path is to use 5p3p representation.
        Domain-level representation is simply not sufficient.
        For instance, in the "bulge" structure below, domains C and c are certainly within reach:
                                            _ _ _C_ _ _  3'
        5'------------A---------------B----/
        3'------------a----------
                                 \ _ _ _c_ _ _ 5'
        However, the domain level representation:
            A -- B -- C
            |
            a
              \- c
        Is equivalent to the would-be circular structure:
                3'_ _ _C_ _ _
                             \----B----------------A----------5'
                                       ------------a----------
                                                              \ _ _ _c_ _ _ 3'
        Where, depending in helix Aa, domains C and c may not be within reach.
        The domain-level graph looses some detail.
        This *can* be recovered via edge attributes, e.g. edge A-B can be directed or otherwise
        inform that B is on the same side of A as c is. But that might be just as much work
        as just using the 5p3p-ends graph.

        Note: Previously I also intended this function to determine whether domain1 and domain2 can hybridize
        and what the energy penalty is, i.e. loop energy, helix-bending energy, etc.


        ## GLOBAL vs LOCAL model: Using only the minimum loop (shortest path) vs all paths ##

        Consider the following two model cases:
         ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
         |           ⁞‾‾‾‾‾‾‾‾‾‾‾               -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾
         |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃                      |           ⁞___________
         ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'
         ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
         |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|  -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|
         |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃   ₃__________⌡       |           ⁞___________ ₃__________⌡
         ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'

        In both cases, when connecting (B3p-E5p), we would only consider the shortest path, (B3p-A5p-A3p-D5p-D3p-E5p).
        However, the second-shortest path (B3p-B5p-C3p-C5p-F3p-F5p-E3p-E5p) would clearly also have an influence
        on the PDF overlap (or, activity/effective_volume) of domain B and E.

        To make the energy calculation more precise, you could do a search for secondary and tertiary loops.
        If these are present, that would increase the activity.

        More refs:
        * https://en.wikipedia.org/wiki/Loop_entropy

        TODO: Check for secondary loops
        TODX: Saving secondary loops/paths might be be useful when determining if a dehybridization will split
              up a complex. Although that is pretty much the reverse process and can't really be used. Nevermind.

        TODO: How does this perform for small loops, e.g. this:     --------------------------------
              which is just 0.67 nm.                                ---------------.________________
              effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2) = 1.66 nm³
              activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3 = 1.66 nm³ / 1.66 nm³ = 1.
        """
        #path = self.domains_shortest_path(domain1, domain2)
        #path_elements = self.domain_path_elements(path)
        # NOTE: The path does not have to span the full length of the element!
        # The element could be really long, e.g. a long DsHelix or the full ss scaffold,
        # with the path only transversing a fraction of the element.
        # Also, for domain-level helices, it is hard to know if the path traverses
        # the domain, or just uses it for the hybridization, traversed at one end only.
        #      To here
        #         |
        # -------´ ---------
        # ------------------ ------from here
        #LRE_len, LRE = max((elem.length_nm, elem) for elem in path_elements if not isinstance(elem, SingleStrand))

        ## 5p3p-level shortest path:
        # domain1, domain2
        if isinstance(elem1, Domain):
            domain1, domain2 = elem1, elem2
            d1end5p, d2end5p = domain1.end5p, domain2.end5p
            d1end3p, d2end3p = domain1.end3p, domain2.end3p
        else:
            d1end5p, d2end5p = elem1, elem2
            domain1, domain2 = d1end5p.domain, d2end5p.domain
        path = self.ends5p3p_shortest_path(d1end5p, d2end5p)
        printd("5p3p graph shortest path from %s to %s:" % (d1end5p, d2end5p))
        pprintd(path)
        path_elements = self.ends5p3p_path_partial_elements(path, length_only='both', summarize=True)
        printd("5p3p graph path elements:")
        pprintd(path_elements)

        ## TODO: Check for secondary loops!

        # list of [(stiffness, total-length), ...]
        # For rigid, double-helical elements, element length, l, is N_bp * 0.34 nm.
        # For single-stranded elements, we estimate the end-to-end distance by splitting the strand into
        #   N = N_nt*0.6nm/1.8nm segments, each segment being the Kuhn length 1.8 nm of ssDNA,
        # where 0.6 nm is the contour length of ssDNA and N_nt is the number of nucleotides in the strand.
        #   E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (1.8 nm)²
        #         = round(N_nt*0.6nm/1.8nm) (1.8 nm)² = N_nt * 0.6/1.8 * 1.8*1.8 nm² = N_nt 0.6*1.8 nm²
        #         = N_nt * lˢˢ * λˢˢ = N_nt * 1.08 nm²
        # Why use stiffness as first maximum criteria??
        try:
            _, LRE_len, LRE_len_sq, LRE_idx = max((stiffness, elem_length, elem_len_sq, i)
                                                  for i, (stiffness, (elem_length, elem_len_sq))
                                                  in enumerate(path_elements)
                                                  if stiffness > 0)
        except ValueError:
            # No stiff elements:
            LRE_len = 0
            LRE_len_sq = 0
            SRE_len = sum(elem_length for (stiffness, (elem_length, elem_len_sq)) in path_elements if stiffness == 0)
            SRE_len_sq = sum(elem_len_sq for (stiffness, (elem_length, elem_len_sq)) in path_elements if stiffness == 0)
        else:
            #LRE = path_elements[LRE_idx] # .pop(LRE_idx)
            # Exclude LRE when calculating SRE length:
            if len(path_elements) > 1:
                SRE_lengths, SRE_sq_lengths = zip(*[(elem_length, elem_len_sq)
                                                    for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
                                                    for stiffness, (elem_length, elem_len_sq) in sub_path])
                # SRE_len_sq = sum(elem_len_sq for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
                #               for stiffness, (elem_length, elem_len_sq) in sub_path)
                SRE_len, SRE_len_sq = sum(SRE_lengths), sum(SRE_sq_lengths)
            else:
                # SRE_lengths, SRE_sq_lengths = [], []
                SRE_len, SRE_len_sq = 0, 0

        # Comparing mean-end-to-end-squared values vs mean end-to-end lengths vs full contour length?
        # There is a difference that sum of squares does not equal the square of sums, so
        # even if LRE_len_sq is > SRE_len_sq, LRE_len could be less than SRE_len.
        # Also, while ds duplexes has full contour length equal to mean end-to-end length,
        # this is not true for ssDNA. Indeed if calculating whether two domains *can* reach,
        # it is probably better to use domain.ds_length.

        if LRE_len > SRE_len:
            print("LRE_len > SRE_len for path between %r and %r" % (domain1, domain2))
            pprint(path)
            pprint(path_elements)
            # The domains cannot reach each other.
            # Hybridization requires helical bending; Not implemented yet; just returning 0 meaning "impossible".
            # TODO: Implement hybridization via helical bending.
            #   Persistance length 50 nm (physiological salt)
            #   - Depends on ionic strength and cationic valency
            # TODO: Look at formulas for k_on and k_off rates under stress.
            #   For DNA, there is certainly a difference between axial "ripping" and perpendicular "zipping".
            #   - Zippering occours at about 10-15 pN (sequence dependent).
            #   -
            return 0
        ## There is probably some profound relation between the elements and the gamma factor.
        ## E.g. if the contour length is long enough for the domains to reach, but the
        ## SRE mean squared end-to-end distance is less than the LRE, then the SRE will rarely
        ## be sufficiently extended for the domains to hybridize. This decrease in spatial pdf
        ## can be considered equivalent to an increase in effective volume.
        ## Another, more complex case, is when the SRE has only (a) one, or (b) a few links,
        ## in which case the mean squared end-to-end distance is not a good measure of spatial pdf.
        ## In the case where SRE is a rigid 1-element chain of same length as the LRE, the pdf
        ## is essentially a sphere centered at the joint between LRE and SRE. (Similar case when
        ## the LRE is flanked by two rigid elements.)

        # Example 1: SRE_len = 10 * 4 nm = 40 nm; SRE_len_sq = 10 * (4 nm)**2 = 160 nm2.
        #            LRE_len = 20 nm,             LRE_len_sq = (20 nm)**2 = 400 nm2.
        #            LRE_len < SRE_len, but LRE_len_sq > SRE_len_sq
        #            SRE_len/LRE_len = 2 -- higher => lower gamma_corr.
        #            LRE_len_sq/SRE_len_sq = 2.5 -- higher => higher gamma_corr.
        # We could, for instance, say: gamma_corr = 1 + ln(LRE_len_sq/SRE_len_sq)
        # Hmm... probably need to do some further analysis of different examples and try to figure out
        # a proper relationship between link-elements and gamma_corr... And it might not be as simple
        # as a simple exponential correction to (P_loop/P_v0).

        # If LRE_len_sq > SRE_len_sq, then the approximation assumption "we many links of length l_i"
        # is certainly not valid (we only have 1 link of length LRE_len).
        # Instead of considering P_loop(r<rc)/P_v0(r<rc), we have to consider P_SRE(r=LRE+/-rc)/V(r=LRE)/P_v0(r<rc).
        # that is, the probability of the SRE end-end distance equaling LRE length, normalized by the shell
        # volume at r=LRE.
        # This gives us a factor that seems to be:
        # 1/(4 π LRE_len_sq) exp(-3*LRE_len_sq / (2*SRE_len_sq))
        # although the first part would give us a non-unitless factor which is not acceptable. It probably has to be
        # normalized in some way, but for now just use the exponential part.
        # Edit: Actually, the first part is probably something like LRE_len/rc
        #
        # For LRE_len_sq == SRE_len_sq, this gives us exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 1/e = 0.22.
        # For example 1, this will give us:
        # exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 0.02.
        LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq)) if LRE_len_sq else 1

        gamma_corr = 1
        if LRE_len_sq > SRE_len_sq:
            # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
            # Should probably be approximated with some continuous function.
            # gamma = (3/2)*gamma_corr;
            gamma_corr += 0.5
        # Mean end-to-end squared distance between the two domains, aka E_r_sq:
        # Mean end-to-end squared distance.
        # We already have the squared length, Nᵢbᵢ², so we just need to sum:

        mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm

        ## Regarding "effective volume" vs P_loop/P_v0: ##
        # Using "effective volume" is more intuitive than P^rc_loop/P^rc_v0.
        # However, stricly speaking we ARE considering reactant proximity/localization/collision probability
        # relative to standard conditions (1 M) using probability distribution function of the loop-connected reactants.
        # This happens to reduce to something that looks like P_loop/P_v0 = v_0/v_eff,
        # where we can use the mean squared end-to-end distance (mean_sq_ee_dist aka E[r²]) to calculate v_eff.
        # However, this is just an approximation valid only for E[r²] >> rc (rc = critical interaction distance),
        # which happens to reduce to a simple expression for P_loop and an expression for P_loop/P_v0 that does not
        # depend on rc. In general, especially if rc is significant, P_loop/P_v0 could very well depend on rc!
        # There is also no general guarantee that the expression for P_rel = P_loop/P_v0 will reduce to something with
        # an easily-identifiable v_eff subexpression. In that case we can simply define v_eff as
        #   v_eff = v_0 * P_v0(rc) / P_loop(rc) = v_0 / P_rel(rc)       # P_rel < 1 (unless they are really close)
        # although in practice we would not need to calculate v_eff, we would just use
        #   activity = P_v0(rc) / P_loop(rc) × M

        effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2)
        #effective_volume = (2/3*math.pi*mean_sq_ee_dist)**gamma * 1e-24 # 1e-24 to convert nm3 to L.
        #activity = (1/N_AVOGADRO)*(1/effective_volume)
        # Using AVOGADRO_VOLUME_NM3 to avoid redundant conversions:
        activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
        if gamma_corr > 1:
            activity = activity**gamma_corr

        ## When/where to apply extra gamma? ##
        # Note: The extra gamma really should be applied to the activity, not just the effective volume:
        # That is, we would always use exponent of 3/2 for calculating effective volume from mean_sq_ee_dist,
        # and apply the extra gamma to the unitless activity (v0/effective_volume) as:
        #   activity = (v0/effective_volume) ** gamma_corr
        # where gamma_corr = 1+x for gamma_org = 3/2+x
        # This is also what is effectively done in (Dannenberger, 2015), where γ is outside the expression
        #     ΔG = - R T γ ln(C/E[r2])
        # TODO: Check that this is also how Dannenberger actually does it when calculating k₊ in the java code.
        # To keep ΔG "reasonable" for large E[r2], (Dannenberger et al, 2015) adjusts C approximately as:
        #   C = 2.2e-18 m² γ - 2.7e-18 m² = 3.34 C0 γ - 4.0 C0,  C(γ=1.5) = C0
        # Where C0 = 3/(2π) v0**(2/3) = 6.7e-19 m²
        # This corresponds to increasing Avogadro volume, v0.

        # Note: The "activity" as calculated above appears on paper to be a unitless ratio (P_loop/P_v0).
        # However, since it is relative to standard conditions (1 M), we just have to implicitly multiply
        # the unitless ratio with "× 1 M" to get a proper molar activity.

        # TODO: Currently not accounting for bending or zipping energy.
        # TODO: Account for secondary (and tertiary?) loops
        if activity < 1:
            print("Low activity, %s - printing locals():" % activity)
            pprint(locals())
        return activity



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
