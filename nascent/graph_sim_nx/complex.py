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


"""

Module for

"""


import itertools
import networkx as nx
import numpy as np

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen,
                    PHOSPHATE_BACKBONE, DUPLEX_HYBRIDIZATION, STACKING_INTERACTION)

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
            strands = []
        for strand in strands:
            strand.graph['complex'] = self
        self._strands = strands
        self._domains = itertools.chain(s.domains for s in strands)
        self.add_nodes_from(self._domains) # Add all domains as nodes

        # Distances between domains.
        # If we have N domains, then we have sum(1..(N-1)) possible distances?
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



    def fqdn(self):
        return "C[%s]" % (self.graph['sid'])

    def __repr__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "Complex[%s] at %s" % (self.graph['sid'], hex(id(self)))

    def __str__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "C[%s]" % (self.graph['sid'])
