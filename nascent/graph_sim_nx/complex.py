# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=C0103

"""

Module for

"""


import itertools
from collections import deque
import networkx as nx
from networkx.algorithms.shortest_paths import shortest_path
import numpy as np

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)
from .structural_elements import strand, helix, bundle
from .structural_elements.strand import SingleStrand
from .structural_elements.helix import DsHelix
from .structural_elements.bundle import HelixBundle

# Module-level constants and variables:
make_sequential_id = sequential_number_generator()



class Complex(nx.Graph):
    """
    This class represents a graph of connected domains.

    The edges represent ALL types of connections:
     1. Phosphate backbone connections between connected domains in the same strand.
     2. Hybridization connections.
     3. Stacking connections.


    It is not for representing a graph of connected strands,
    although it can generate such a graph using strands_graph().
    Graph representations that can be dynamically generated includes:
        5p3p graph
        strands graph   Graph with strands. Can be generated using blockmodeling.
        helixI graph.   Graph with helices and interfaces:
                        A vertice represents the interface between two domains,
                        or four domais if they are stacked.
        junctionI graph: Like helixI, but represents N-way junctions, e.g. a Holliday junction.
    """
    def __init__(self, strands=None, origin="o"):
        super().__init__(sid=next(make_sequential_id),
                        )
        if strands is None:
            strands = {}
        for strand in strands:
            strand.graph['complex'] = self
        # Should strands be a set or list? Set is most natural, but a list might play better with a adjacency matrix.
        # Then, you could just have a dict mapping strand -> index of the matrix
        self.strands = set(strands)
        self.domains = self.nodes
        #self._domains = itertools.chain(s.domains for s in strands)
        self.add_nodes_from(self.domains_gen()) # Add all domains as nodes

        # Distances between domains.
        # If we have N domains, then we have sum(1..(N-1)) = N**2 - N possible distances?
        # Alternatively, use a proper matrix. Yes, it will be degenerate and include zero-distance, but
        # it might offer better performance and more natural calculations.
        self.domain_distances = {} # {frozenset(d1, d2): dist}
        # Distances between complementary domains (not neseccarily hybridized)
        # Is a subset of domain_distances: {d1d2: v for d1d2 in domain_distances if d1.Name}
        self.compl_domain_distances = {} # {frozenset(d1, d2): dist}

        # Having a list of strands indexed by strand name is useful:
        self.strands_by_name = {}
        for strand in strands:
            if strand.name not in self.strands_by_name:
                self.strands_by_name[strand.name] = [strand]
            else:
                self.strands_by_name[strand.name].append(strand)

        # State fingerprint is used for cache lookup
        # To reduce recalculation overhead, it is separated into three parts:
        # state = strands + hybridizations + stacking fingerprints.
        self._state_fingerprint = None
        self._strands_fingerprint = None
        self._hybridization_fingerprint = None
        self._stacking_fingerprint = None

        # Graph with only hybridization connections between strands. Might be useful for some optimizations.
        self.strand_hybridization_graph = nx.Graph()
        self.stacking_graph = nx.Graph()      # Graph with only stacking connections.
        # (^^ Note: I'm not yet sure whether stacking interactions are part of the main complex graph...)
        # Not sure if each complex should have a 5p3p graph. After all, we already have a system-level 5p3p graph in the simulator.
        self.end5p3p_graph = nx.Graph()
        self.hybridized_domains = set() # set of frozenset({dom1, dom2}) sets, specifying which domains are hybridized.
        self.stacked_domains = set()    # set of tuples. (5p-domain, 3p-domain)
        # Note that stacking is directional.
        # Also, "5p_domain" is the domain at the 5' end of the interface, which is STACKING USING ITS 3' END.
        self.cuid = next(make_sequential_id)
        self.uuid = next(sequential_uuid_gen)


    def state_fingerprint(self):
        """
        Make a unique fingerprint or hash for cache lookup.

        What does the state of a complex depend on?
         1. Strands, either as a list or set {(strand-specie, count)}
            (vertices)
         2. Hybridization connections.
         3. Stacking connections.

        [NOTE: I am currently accepting hybridization degeneracies since I don't want to test for graph isomorphisms]

        Regarding connections:
        We really have several different types of connections (edges):
         1. phosphate-backbone connections within each strand.
         2. hybridizations between domains.
         3. stacking between the ends of domains.

        There are many different ways to uniquely represent a complex, but the requirements are:
         * Must be at the specie level, not instance. I.e. if during the simulation,
                two complexes are in the same state, they should produce the same fingerprint.
        Additionally, a fingerprint should if possible
         * Prevent degeneracies. I.e. if two identical strands are part of the same complex,
                even a sorted-by-strand-domains adjacency matrix could have permutable rows/columns.

        Note: a fingerprint is different from a "persistent representation" in that:
         1. It only has to be one-way. We don't have to be able to re-generate the full complex
            from the fingerprint, we just have to ensure that the fingerprint is unique.
         2. A persistent representation can have degeneracies that produce the same complex.
            A state fingerprint should avoid these degeneracies if possible.

        How to uniquely specify a single strand within a complex?
         - Sometimes, the reason we need a complex state is for caching, where we are indexing together with
            a strand/domain, e.g.
            (strand-domain-specie, complex-state-fingerprint)
         - However, if we have multiple copies of the same strand, they will both have the
            same complex-state-fingerprint. How to know them apart?
            The best is probably to probe the local environment.

        """
        if not self._state_fingerprint:
            self._state_fingerprint = hash((
                self.strands_fingerprint(),
                self.hybridization_fingerprint(),
                self.stacking_fingerprint()
                ))

    def domains_gen(self):
        return (domain for strand in self.strands for domain in strand.domains)

    def strands_species_count(self):
        species_counts = {}
        for strand in self._strands:
            if strand.name not in species_counts:
                species_counts[strand.name] = 1
            else:
                species_counts[strand.name] += 1
        #return species_counts
        return frozenset(species_counts.items())

    def strands_fingerprint(self):
        if not self._strands_fingerprint:
            self._strands_fingerprint = hash(self.strands_species_count())
        return self._strands_fingerprint

    def hybridization_fingerprint(self):
        """
        Return a hash based on the current hybridizations in the complex.
        (As opposed to other factors determining the complex structure such as domains,
        backbone and stacking interactions).
        Challenge: How to
        (1) Ensure that domain connections are precisely specified. E.g. if we have two copies of strand sA and
            two copies of strand sB, and we have the connection: {sA.da, sB.dA}
            How do we know which if the strands sA are connected to which strand sB?
            (Also, specifying sA.da may not be unique if sA has two da domains...)
            We can specify this by adding a unique or sequential id to the strands. Then we have
            connections of type {sA#1.da#1, sB#3.dA#4}.
        (2) But that produces a new issue: How to avoid degeneracies.
            Since the two sA strands are interchangeable, the following connection sets give the same complex:
                edges1 = {{sA#1.da, sB.dA}, {sA#2.db, sB.dB}}
                edges2 = {{sA#2.da, sB.dA}, {sA#1.db, sB.dB}}
            However, they are only interchangeable if all sA#1 connections in edges1 maps to sA#2 in edges2.
            (a kind of isomorphism)

        (3) What about domains that are their own palindrome? E.g. 5'-AGAGCTCT-3' can bind to it self. dP:dP
                then the edges set becomes {sA#1.dP#1, sA#1.dP#2} - not a problem, each species have unique id.

        Problem description and possible solutions:
            We need to determine if edges1 is isomorphic to edges2. This is quite concretely a graph isomorphism
            problem, and determining if the two graphs are isomorphic can be approached using graph theory,
            where the two sets (edges1 and edges2) are modelled as a biparite graph (mapping edges1 to edges2).
            "The graph isomorphism problem is in NP, but it is not known whether it is NP complete.
            In practice, GI can be solved quickly for hundreds of vertices with e.g., NAUTY."
            "Isomorphic graphs have the same degree sequence. However, two graphs with the same degree sequence
                are not necessarily isomorphic."
            Networkx uses the vf2 algorithm to determine if two graphs are isomorphic, based on
                "An Improved Algorithm for Matching Large Graphs." (Cordella et al, )
            However, we don't want to determine whether all elements are isomorphic, we only want to test
            whether sA#1.da can be fully switched with sA#2.da.
            More isomorphism refs:
            * http://www.naftaliharris.com/blog/groupiso/
            * http://math.stackexchange.com/questions/90923/isomorphism-of-sets
            * http://www.python-course.eu/graphs_python.php
            * https://en.wikipedia.org/wiki/Graph_canonization

        I will defer the issue of isomorphic degeneracies to later.
        """
        if not self._hybridization_fingerprint:
            hyb_edges = self.hybridization_edges()
            edgesfs = frozenset(frozenset(d.domain_strand_specie for d in edge) for edge in hyb_edges)
            self._hybridization_fingerprint = hash(edgesfs)
        return self._hybridization_fingerprint


    def hybridization_edges(self):
        """
        How to get a list of domain hybridization?
        * self.edges() will generate a list of all connections: backbone, hybridization and stacking.

        """
        # you can loop over all edges and filter by type:
        #hyb_edges = [frozenset(d1, d2) for d1, d2, cnxtype in self.edges(data='type')
                     #if cnxtype == DUPLEX_HYBRIDIZATION]
        # or maybe it is easier to keep a dict with hybridization connections:
        hyb_edges = self.hybridized_domains
        return hyb_edges

    def stacking_fingerprint(self):
        """
        Return a stacking fingerprint for use with caching.
        This is just hash(frozenset(tuple(d.duid for d in edge) for edge in self.stacked_domains))
        """
        if not self._stacking_fingerprint:
            # This is using (d1, d2) tuple rather than {d1, d2} frozenset, since directionality matters.
            edgesfs = frozenset(tuple(d.duid for d in edge) for edge in self.stacked_domains)
            self._stacking_fingerprint = hash(edgesfs)
        return self._stacking_fingerprint

    def stacking_edges(self):
        """
        How to get a list of domain hybridization?
        Note: stacking is directional and cannot be saved within an undirected graph.
        Perhaps I should just switch to using a directed graph?

        """
        # you can loop over all edges and filter by type:
        # But need directionality to know which end of the domain is pairing.
        #hyb_edges = [frozenset(d1, d2) for d1, d2, cnxtype in self.edges(data='type')
                     #if cnxtype == STACKING_INTERACTION]
        # For now, I have to keep a dict with hybridization connections:
        return self.stacked_domains


    def stacked_subgraph(self):
        """
        Return a graph only with stacked, double-helical domains.
        This graph has the following edges (interactions):
        * Stacking (up/downstream)
        * Hybridization
        As always, we have to determine whether to use the domain-level or 5p3p-level graph.
        This function returns the domain-level subgraph.

        Do we need to filter the domains? We could include the ss domains as non-connected nodes.
        """
        hybridized_domains = [domain for domain in self.domains if domain.partner is not None]
        edges = [(s, t) for s, t, attrs in self.edges(data=True)
                 if attrs['interaction'] in (HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)]
        subg = self.subgraph(hybridized_domains)
        subg.add_edges_from(edges)
        return subg


    def domain_path_elements(self, path):
        """
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

    def end5p3p_path_partial_elements(self, path, length_only=False, summarize=False):
        """
        Returns a list of structural elements based on a 5p3p-level path (list of 5p3p ends).

        For this, I will experiment with primitive list-based representations rather than making full
        element objects.
        path_edges = [(1, [(length, length_sq, source, target), ...]),
                      (3, [(length, length_sq, source, target), ...]),
                      (2, [(length, length_sq, source, target)]
        Where 1 indicates a list of single-stranded edges,
        2 indicates a hybridization edge, and 3 indicates a list of stacked edges.
        Since only the interaction type and lenghts are important, maybe just
        path_edges = [(1, [length, length, ...]), (3, [length, length, ...], ...)]
        Edit: Instead of 1/2/3, use the standard INTERACTION constants values.
        """
        path_edges = []
        last_interaction = None
        interaction_group = None
        for i in range(len(path)-1):
            source, target = path[i], path[i+1]
            # Method 1: Use the networkx graph API and edge attributes:
            edge = self.end5p3p_graph[source][target]
            length = edge['length']
            interaction = edge['interaction']
            # Method 2: Manually determine what type of interaction we have based on source, target:
            if target == source.stack_partner:
                interaction = STACKING_INTERACTION
                length = 1
                length_sq = 1
            elif target == source.hyb_partner:
                interaction = HYBRIDIZATION_INTERACTION
                length = 1 # one nm from one helix to the next. We really don't know for sure because it turns.
                length_sq = 1
            elif target in (source.bp_upstream, source.pb_downstream):
                if source.hyb_partner and \
                    ((source.end == "5p" and target == source.pb_downstream) or
                     (source.end == "3p" and target == source.bp_upstream)):
                    # Above could probably be done with an XOR:
                    # (source.end == "5p") != (target == source.pb_downstream) # boolean xor
                    # (source.end == "5p") ^ (target == source.pb_downstream)  # bitwise xor
                    # We have a stacked duplex:
                    interaction = STACKING_INTERACTION
                    length = source.domain.ds_length_nm
                    length_sq = source.domain.ds_length_sq
                else:
                    # We have a single-stranded domain:
                    interaction = PHOSPHATEBACKBONE_INTERACTION
                    length = source.domain.ss_length_nm     # mean end-to-end length; not contour length
                    length_sq = source.domain.ss_length_sq
            else:
                raise ValueError("Could not determine interaction between %s and %s" % (source, target))

            if interaction == HYBRIDIZATION_INTERACTION:
                # We only distinguish between ds-helix vs single-strand; use STACKING_INTERACTION to indicate ds-helix:
                interaction = STACKING_INTERACTION

            if interaction != last_interaction:
                path_edges.append((last_interaction, interaction_group))
                interaction_group = []
                last_interaction = interaction
            if length_only:
                if length_only == 'sq':
                    interaction_group.append(length_sq)
                elif length_only == 'both':
                    interaction_group.append((length, length_sq))
                else:
                    interaction_group.append(length)
            else:
                interaction_group.append((length, length_sq, source, target))
        if summarize and length_only:
            if length_only == 'both':
                # Return a list of (interaction, (length, length_squared)) tuples:
                return [(interaction, (sum(lengths) for lengths in zip(*lengths_tup))) # pylint: disable=W0142
                        for interaction, lengths_tup in path_edges]
            else:
                return [(interaction, sum(lengths)) for interaction, lengths in path_edges]
        return path_edges



    def domains_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        """
        return shortest_path(self, domain1, domain2)

    def ends5p3p_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        TODO: Verify shortest path for end3p as well?
        """
        return shortest_path(self.end5p3p_graph, domain1.end5p, domain2.end5p)


    def intra_complex_hybridization(self, domain1, domain2):
        """
        Determine whether domain1 and domain2 can hybridize and what the energy penalty is,
        i.e. loop energy, helix-bending energy, etc.
        Returns:
            :activity:
        so that
            c_j = k_j * activity * N_A⁻¹

        Flow-chart style:
        1. Are the two domains connected by a single ds helix?
        2. (...)

        More direct approach:
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
        path = self.ends5p3p_shortest_path(domain1, domain2)
        path_elements = self.end5p3p_path_partial_elements(path, length_only='both', summarize=True)
        # list of [(interaction, total-length), ...]
        # For rigid, double-helical elements, element length, l, is N_bp * 0.34 nm.
        # For single-stranded elements, we estimate the end-to-end distance by splitting the strand into
        #   N = N_nt*0.6nm/1.8nm segments, each segment being the Kuhn length 1.8 nm of ssDNA,
        # where 0.6 nm is the contour length of ssDNA and N_nt is the number of nucleotides in the strand.
        #   E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (1.8 nm)²
        #         = round(N_nt*0.6nm/1.8nm) (1.8 nm)² = N_nt * 0.6/1.8 * 1.8*1.8 nm² = N_nt 0.6*1.8 nm²
        #         = N_nt * lˢˢ * λˢˢ = N_nt * 1.08 nm²
        _, LRE_len, LRE_len_sq, LRE_idx = max((interaction, elem_length, elem_len_sq, i)
                                              for i, (interaction, (elem_length, elem_len_sq))
                                              in enumerate(path_elements)
                                              if interaction in (STACKING_INTERACTION, HYBRIDIZATION_INTERACTION))
        LRE = path_elements[LRE_idx] # .pop(LRE_idx)
        # Exclude LRE when calculating SRE length:
        SRE_len_sq = sum(elem_len_sq for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
                      for interaction, elem_len_sq in sub_path)

        # Comparing mean-end-to-end-squared values vs mean end-to-end lengths vs full contour length?
        # There is a difference that sum of squares does not equal the square of sums, so
        # even if LRE_len_sq is > SRE_len_sq, LRE_len could be less than SRE_len.
        # Also, while ds duplexes has full contour length equal to mean end-to-end length,
        # this is not true for ssDNA. Indeed if calculating whether two domains *can* reach,
        # it is probably better to use domain.ds_length.


        if LRE_len_sq > SRE_len_sq:
            # The domains cannot reach each other.
            # Hybridization requires helical bending; TODO: Not implemented yet.
            return 0
        # TODO: Implement loop energy or activity calculation...
        # Currently not accounting for bending or zipping energy.
        # Mean end-to-end squared distance.
        # We already have the squared length, Nᵢbᵢ², so we just need to sum:
        E_r_sq = sum()







    def fqdn(self):
        return "C[%s]" % (self.graph['sid'])

    def __repr__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "Complex[%s] at %s" % (self.graph['sid'], hex(id(self)))

    def __str__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "C[%s]" % (self.graph['sid'])
