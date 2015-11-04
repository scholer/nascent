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

# pylint: disable=C0103

"""

Module for

"""


from collections import defaultdict #, deque
import networkx as nx

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        N_AVOGADRO, AVOGADRO_VOLUME_NM3)

# Module-level constants and variables:
make_sequential_id = sequential_number_generator()



class Complex(nx.MultiGraph):
    """
    This class represents a graph of connected domains.

    The edges represent ALL types of connections:
     1. Phosphate backbone connections between connected domains in the same strand.
     2. Hybridization connections.
     3. Stacking connections.

    I believe we need a MultiGraph to represent this, since
    two domains can be edge-connected with both a backbone interaction AND
    a stacking interaction, or, if we allow zero-loop hairpins, both
    backbone interaction AND hybridization interaction.

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
    def __init__(self, data=None, strands=None, origin="o"):
        self.cuid = next(make_sequential_id)
        self.uuid = next(sequential_uuid_gen)   # Universally unique id; mostly used for debugging.
        super().__init__(data=data, cuid=self.cuid)
        if strands is None:
            strands = []
        for strand in strands:
            strand.complex = self
        self.domains = self.nodes  # alias
        # Should strands be a set or list? Set is most natural, but a list might play better with a adjacency matrix.
        # Then, you could just have a dict mapping strand -> index of the matrix
        self.strands = set(strands)
        #self._domains = itertools.chain(s.domains for s in strands)
        if data is None:
            #self.add_nodes_from(self.domains_gen()) # Add all domains as nodes
            # Use the (linear) strand graphs to build an initial graph with phosphate backbones:
            for strand in strands:
                self.add_edges_from(strand)

        # Distances between domains.
        # If we have N domains, then we have sum(1..(N-1)) = N**2 - N possible distances?
        # Alternatively, use a proper matrix. Yes, it will be degenerate and include zero-distance, but
        # it might offer better performance and more natural calculations.
        self.domain_distances = {} # {frozenset((d1, d2)): dist}
        # Distances between complementary domains (not neseccarily hybridized)
        # Is a subset of domain_distances: {d1d2: v for d1d2 in domain_distances if d1.name}
        self.compl_domain_distances = {} # {frozenset((d1, d2)): dist}

        # Having a list of strands indexed by strand name is useful (e.g. for fingerprinting)
        self.strands_by_name = defaultdict(set)
        for strand in strands:
            self.strands_by_name[strand.name].add(strand)

        # Domains, indexed by name:
        # Should it be a set or a list?
        # We do have random access and modification, so maybe a list is better?
        self.domains_by_name = defaultdict(set)
        for domain in self.domains():
            self.domains_by_name[domain.name].add(domain)

        self.hybridized_domains = set() # set of frozenset({dom1, dom2}) sets, specifying which domains are hybridized.
        self.stacked_domains = set()    # set of tuples. (5p-domain, 3p-domain)

        # State fingerprint is used for cache lookup
        # To reduce recalculation overhead, it is separated into three parts:
        # state = strands + hybridizations + stacking fingerprints.
        self._state_fingerprint = None
        self._strands_fingerprint = None
        self._hybridization_fingerprint = None
        self._stacking_fingerprint = None

        ### Graphs: ###

        # TODO: If you have a strand_graph, then strands set is not needed
        # TODO: If you are going to use multi-level per-complex graphs, they must be compatible with initial data!
        # Edit: Breaking and merging per-complex domain-level graphs is already super tedious.
        #       I don't want to do that for strand-level and ends5p3p as well.
        #       Probably better to make use of system-level graphs as much as possible.
        #self.strand_graph = nx.MultiGraph() # Only add if you really know it is needed.
        # Not sure if each complex should have a 5p3p graph.
        # After all, we already have a system-level 5p3p graph in the simulator.
        #self.ends5p3p_graph = nx.Graph()

        # Graph with only hybridization connections between strands. Might be useful for some optimizations.
        #self.strand_hybridization_graph = nx.Graph()
        #self.stacking_graph = nx.Graph()      # Graph with only stacking connections.
        # (^^ Note: I'm not yet sure whether stacking interactions are part of the main complex graph...)
        # Note that stacking is directional.
        # Also, "5p_domain" is the domain at the 5' end of the interface, which is STACKING USING ITS 3' END.


    def add_domain(self, domain, update_graph=False):
        """
        Update domains_by_name.
        update_graph defaults to False, because graph-related stuff
        is usually handled externally for domains.
        """
        #self.domains.add(strand)  # domains are tracked by the graph, self.nodes()
        self.domains_by_name[domain.name].add(domain)
        if update_graph:
            # strand is also a (linear) graph of domains:
            self.add_node(domain)
            # edge_attrs = {'weight': len(domain), 'len': len(domain), 'type': PHOSPHATEBACKBONE_INTERACTION}
            # self.ends5p3p_graph.add_path(domain.end5p, domain.end3p, edge_attrs)

    def remove_domain(self, domain, update_graph=False):
        """ Remove a single domain, updating self.domains_by_name. """
        self.domains_by_name[domain.name].remove(domain)
        if update_graph:
            # strand is also a (linear) graph of domains:
            self.remove_node(domain)
            # self.ends5p3p_graph.remove_nodes_from([domain.end5p, domain.end3p])


    def add_strand(self, strand, update_graph=False):
        """ We keep track of strands for use with fingerprinting, etc. """
        self.strands.add(strand)
        self.strands_by_name[strand.name].add(strand)
        strand.complex = self
        if update_graph:
            # strand is also a (linear) graph of domains:
            self.add_edges_from(strand)
            # self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph)
            # self.strand_graph.add_node(strand)

    def remove_strand(self, strand, update_graph=False):
        """ We keep track of strands for use with fingerprinting, etc. """
        self.strands.remove(strand)
        self.strands_by_name[strand.name].add(strand)
        if strand.complex == self:
            strand.complex = None
        if update_graph:
            # strand is also a (linear) graph of domains:
            self.remove_nodes_from(strand)
            # self.ends5p3p_graph.add_nodes_from(strand.ends5p3p_graph)
            # self.strand_graph.remove_node(strand)

    def add_strands(self, strands, update_graph=False):
        """ Strands must be a set. """
        for strand in strands:
            self.strands_by_name[strand.name].add(strand)
            strand.complex = self
            if update_graph:
                # strand is also a (linear) graph of domains:
                self.add_edges_from(strand)
                # self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph)
        if not isinstance(strands, set):
            strands = set(strands)
        self.strands |= strands

    def remove_strands(self, strands, update_graph=False):
        """ Strands must be a set. """
        for strand in strands:
            self.strands_by_name[strand.name].remove(strand)
            if strand.complex == self:
                strand.complex = None
            if update_graph:
                # strand is also a (linear) graph of domains:
                self.remove_nodes_from(strand)
                # self.ends5p3p_graph.remove_nodes_from(strand.ends5p3p_graph)
                # self.strand_graph.remove_node(strand)
        if not isinstance(strands, set):
            strands = set(strands)
        self.strands -= strands



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
            self._state_fingerprint = (str(self), hash((
                self.strands_fingerprint(),
                self.hybridization_fingerprint(),
                self.stacking_fingerprint()
                )))
        return self._state_fingerprint

    # def domains_gen(self):
    #     return (domain for strand in self.strands for domain in strand.domains)

    def strands_species_count(self):
        """ Count the number of strand species. Used as part of complex state finger-printing. """
        species_counts = {}
        for strand in self.strands:
            if strand.name not in species_counts:
                species_counts[strand.name] = 1
            else:
                species_counts[strand.name] += 1
        #return species_counts
        # Used for hashing, so return a hashable frozenset(((specie1, count), ...))
        return frozenset(species_counts.items())

    def strands_fingerprint(self):
        """ Create a finger-print of the current strands (species). """
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
        #hyb_edges = [frozenset((d1, d2)) for d1, d2, cnxtype in self.edges(data='type')
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
        #hyb_edges = [frozenset((d1, d2)) for d1, d2, cnxtype in self.edges(data='type')
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



    def fqdn(self):
        """ Return a fully-qualified name. """
        return "C[%s]" % (self.cuid)

    def __repr__(self):
        # return "%s[%s]" % (self.name, self.ruid % 100)
        return "Complex[%s] at %s" % (self.cuid, hex(id(self)))

    def __str__(self):
        # return "%s[%s]" % (self.name, self.ruid % 100)
        # String representation should be invariant through the life-time of an object:
        return "C[%s]" % (self.cuid)
