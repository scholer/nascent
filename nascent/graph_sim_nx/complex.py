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

Module for Complexes and SuperComplexes (complex of complexes).

Reminder: Why do we need complexes when we have system-graphs?
Because we need state fingerprints. And state fingerprints
requires complexes with a particular state.

Q: Is SuperComplexes a good idea?
The idea was that stacking interactions are often very brief,
whereas hybridization interaction will eventually "solidify",
thus making hybridized complexes more stable.
Using super-complexes would prevent having to merge and break
complexes continuously.
However, at least for now, this seems much too advanced a solution,
and I will revert back to just having a single "Complex" class,
which can be held together by either hybridization OR STACKING
interactions.


"""

from __future__ import absolute_import, print_function, division
import sys
from collections import defaultdict, Counter, deque
import networkx as nx
from pprint import pprint
from functools import wraps
import pdb
# import numpy as np

# Relative imports
from .connected_multigraph import ConnectedMultiGraph
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .nx_utils import draw_graph_and_save
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        N_AVOGADRO, AVOGADRO_VOLUME_NM3)
from .debug import printd, pprintd

# Module-level constants and variables:
make_sequential_id = sequential_number_generator()
# A generator, get next sequential id with next(make_helix_sequential_id)
helix_sequential_id_generator = sequential_number_generator()

# supercomplex_sequential_id_gen = sequential_number_generator()


def state_should_change(method):
    """ Decorator for methods that are expected to change the Complex state fingerprint.
    Will save the current-soon-to-be-obsolete fingerprint, invoke the method to change the state,
    and finally ensure that the the method actually changed the Complex' state by
    asserting that the new state fingerprint is not the same as the old fingerprint."""
    @wraps(method)
    def replacement_func(self, *args, **kwargs):
        old_fingerprint = self.state_fingerprint()
        ret = method(self, *args, **kwargs)
        self.reset_state_fingerprint()
        new_fingerprint = self.state_fingerprint()
        assert new_fingerprint != old_fingerprint
        return ret
    return replacement_func




# connected_component_subgraphs is not implemented for directed graphs.
# But we can just use system graphs instead...

class Complex(nx.MultiDiGraph):
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

    Discussion: What graph representation do we need for complexes?
      - I thought we would use the complex graph to calculate intra-complex activity,
        but that was never the case: We use system-level graphs, either ends5p3p_graph or interface_graph
      - So, what graph, if any, do we need the complexes to be?
      - It *is* nice to have per-complex graphs, if only to visualize the different complexes.
      - Do we need the per-complex graphs to be directed graphs or undirected?
    Discussion: How to generate proper state fingerprints?
    """
    def __init__(self, data=None, strands=None, origin="o", reaction_deque_size=5): #, reaction_spec_pair=None):
        self.cuid = next(make_sequential_id)
        self.uuid = next(sequential_uuid_gen)   # Universally unique id; mostly used for debugging.
        super(Complex, self).__init__(data=data, cuid=self.cuid)
        if strands is None:
            strands = []
        for strand in strands:
            strand.complex = self
        self.domains = self.node.keys  # alias.. but perhaps better to have a separate set?
        # Should strands be a set or list? Set is most natural, but a list might play better with a adjacency matrix.
        # Then, you could just have a dict mapping strand -> index of the matrix
        self.strands = set(strands)

        # If we know how many of each domain species we have, we can make the problem of graph isomorphism
        # a little easier. (In particular the case where we only have *one* of each domain species.)
        # For domains with specie count >1, instead of using just domain specie, you could use
        # domain.in_complex_identifier. It will not be 100% ceretain/accurate, but it should be close enough.
        self.domain_species_counter = Counter()
        # We are already calculating a strand-species count, might as well make it permanent running counter:
        self.strand_species_counter = Counter()

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

        # self.hybridized_domains = set() # set of frozenset({dom1, dom2}) sets, specifying which domains are hybridized.
        self.hybridized_pairs = set() # set of domain pairs: frozenset({dom1, dom2})
        self.stacked_pairs = set()    # set of domain two-tuples. (5p-domain, 3p-domain)

        # State fingerprint is used for cache lookup
        # To reduce recalculation overhead, it is separated into three parts:
        # state = strands + hybridizations + stacking fingerprints.
        self._state_fingerprint = None
        self._strands_fingerprint = None
        self._hybridization_fingerprint = None
        self._stacking_fingerprint = None

        ## DEBUG attributes: ##
        self._historic_strands = []
        self._historic_fingerprints = [0]
        # edit: deque are not good for arbitrary access, using a list
        self.reaction_deque = deque([0], maxlen=reaction_deque_size)
        #self.reaction_deque = [0]    # deque(maxlen=reaction_deque_size)
        # self.reaction_deque_size = reaction_deque_size
        self.reaction_invocation_count = Counter()  # maps reaction_key => Number of times that reaction has occured.
        self.reaction_throttle_cache = {}
        self.previous_reaction_attrs = {}
        self.N_strand_changes = 0
        self.icid_radius = 5
        # Use domain instance instead of in_complex_identifier. Caching efficiency will decrease.
        # If Complex.icid_use_instance has been set to True, fall back to using domain instances as
        # in_complex_identifier (icid) for *all* domains. If icid_use_instance is not True but is
        # boolean True, it is assumed to be a set of domains for which to use domain instance as icid.
        self.icid_use_instance = False or True

        self.history = deque(maxlen=100) # []

        ## Load initial strands:
        if data is None:
            #self.add_nodes_from(self.domains_gen()) # Add all domains as nodes
            # Use the (linear) strand graphs to build an initial graph with phosphate backbones:
            for strand in strands:
                if self.is_directed() and not strand.is_directed():
                    strand = nx.MultiDiGraph(strand)
                self.add_nodes_from(strand.nodes(data=True))
                self.add_edges_from(strand.edges(data=True, keys=True))

        for strand in strands:
            strand.complex = self
            self.strands_by_name[strand.name].add(strand)
            self.strand_species_counter[strand.name] += 1

        # Domains, indexed by name:
        # Should it be a set or a list?
        # We do have random access and modification, so maybe a list is better?
        self.domains_by_name = defaultdict(set)
        for domain in self.domains():
            self.domains_by_name[domain.name].add(domain)
            self.domain_species_counter[domain.name] += 1

        ### Helix bundling: ###
        # See also modules   helix_duplex  and   structural_elements.helix   !
        # For now, a "helix" is just a set of fully hubridized and stacked Domain (or DomainEnd) nodes.
        # Actually: What happens when merging two complexes? Keeping track of helix indices would be a pain.
        # Probably best to have helix objects which can remain unique between complex merges.
        # But: Making a new object for every hybridization reaction might be expensive (and un-cacheable).
        # Better, then, to use helix_sequential_id_generator to generate unique numbers to identify helices with.
        # What are double-helices? DomainEnds? Or IfNodes? Or an object instance with both?
        # - If we know a node is part of a helix, getting the full helix is pretty easy, since all
        #   DomainEnds are linked. Just use DomainEnd.stacked_upstream() and DomainEnd.stacked_downstream().
        #
        self.helices = set()
        # For a pair of helices hi, hj: A list of paths connecting bp pos (or DomainEnd nodes) between hi hj:
        # (hi, hj): [(hi_bpidx, hj_bpidx, [path]), ...]   or maybe   {(DomainEnd1, DomainEnd2): [path], ...} dict?
        self.helix_connections = {}
        # For each pair of helices ha, hb: For each pair of crossover paths ti, tj:
        # A value denoting how "connected" the two helices are because of those two crossovers:
        # B_ij = (len(t_i) + len(t_j))/(d1 + d2),   where d1 and d2 are distances between crossovers on ha and hb.
        # This is effectively a matrix with indices same as for helix_connections, e.g.
        # connection t_i in [t0, t1, ..., tj, ...] is index B[j,:] and B[:,j] in matrix B for helix-pair hi, hj.
        # If connections don't have a regular index, consider using a numpy.array?
        # Or just a dict {(ti, tj): b_ij}
        self.helix_connectivities = {}
        # For each pair of helices (ha, hb), the minimum b_ij of all b in matrix B:
        self.helix_minimum_freedom_factor = {}
        # When a new connection is made, instead of completely re-calculating B, just re-calculate for the new
        # connection vs all existing connections and then see if any of these are less than the minimum b_ij.
        # When a connection is broken, we must find the new minimum of the B matrix.
        # Note that if a connection is broken due to a duplex-dehybridization reaction,
        # you probably have to split up the matrix. Annoying, but should be doable.
        # This is always the case when a helix-duplex is broken (unstacked, unhybridized),
        # so the data structures should be small/light and easy to break down and re-build.
        # It is probably nice to know given a DomainEnd node, whether it is part of a helix-duplex and which one:
        self.helix_by_domainend = {}  # (duplexed, fully-stacked double-helices)
        self.helix_by_ifnode = {}
        # if domain_end in self.helix_by_domainend: helix = self.helix_by_domainend[domain_end]; ...
        self.helix_bundles = []
        self.helix_bundle_by_helix_idx = {}  # helix_idx:

        ### Complex energy: ###
        # Includes: Volume energies (for every bi-molecular reaction), Shape/Loop energy and of course,
        # stacking- and hybridization energies. For each, we have both dH and dS.
        # self.energies_dHdS = {'volume': [0, 0], 'shape': [0, 0],
        #                       'stacking': [0, 0], 'hybridization': [0, 0]}
        # Enthalpies in self.energies_dHdS[0], entropies in self.energies_dHdS[1]
        self.energies_dHdS = [{'volume': 0.0, 'shape': 0.0, 'stacking': 0.0, 'hybridization': 0.0},
                              {'volume': 0.0, 'shape': 0.0, 'stacking': 0.0, 'hybridization': 0.0}]
        self.energy_total_dHdS = [0, 0]
        # energies_dHdS = np.array([(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0)],
        #     dtype=[('volume', float), ('shape', float), ('stacking', float), ('hybridization', float)])
        # Edit: numpy's "structured array" is almost useless. Use ndarray if you want speed.


        ### Keeping track of loops: ###
        # Dict with loops. Keyed by LoopID ? Maybe a set instead?
        # No, loop is a dict (mutable). We don't have frozendict in Python, but there is the MappingProxyType obj class.
        # Having loops as dicts indexed by an ID is nearly identical to having loops being objects in Python.
        # In Python, objects can be a little slow because they are not just data containers but also have methods etc.
        # In e.g. Julia, having a separate "Loop" type would be faster than Dict. (Julia is ALL about types).
        self.loops = {}  # LoopID => loop dict
        # For each loop, we have a dict with entries (suggested):
        #     path_list: The original path with the nodes forming the loop in the ends.
        #     ifnodes: list of InterfaceNodes. Should this always be the top delegate?
        #           What if we have a duplex in a loop and the duplex dehybridises but the loop is not broken.
        #           Would we update the ifnodes list?
        #     nodes: frozenset of all nodes.
        #     edges:    frozenset of all edges.
        #     edges_specs: frozenset of all edges as Domain-species (rather than objects) - for caching.
        #     edges_properties: a dict with {edge_spec: (length, len_sq, stiffness) tuples for each edge}
        #     loop_entropy: a measure of the loops entropy in units of kb or R.
        # The 'edges' entry is a frozenset of "edge_pairs", each edge_pair a tuple with
        #   (frozenset(source, target), interaction_type)
        #    -- hashing these gives "state-independent" loop fingerprint.
        #    But if one strand/domain is replaced by another (but otherwise identical) strand, this hash changes.
        # The 'edges_specs' is like edges, but uses domain species (e.g. 'A:5p', 'b:3p', etc).
        # This should give the same hash for identical strands.
        # The loop is keyed by the hash: loop_elms_hash = hash(edges_specs)
        # If we want a state-dependent fingerprint, we could use:
        #   hash(tuple(edge_spec+edges_properties[edge_spec] for edge_spec in edges_specs))

        # We also want a way to quickly determine whether a DomainEnd node is part of a loop.
        # self.loops_by_interface = defaultdict(set)  # InterfaceNode => {loop1, loop3, ...}
        # Since loops are mutable we cannot have a set of loops, so use a LoopID and
        self.loopids_by_interface = defaultdict(set)  # InterfaceNode => {loop1_id, loop3_id, ...}
        # Should we use a set or list? We probably have random deletion of loops, so set is arguably better.


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


    # def add_domain(self, domain, update_graph=False):
    #     """
    #     Update domains_by_name.
    #     update_graph defaults to False, because graph-related stuff
    #     is usually handled externally for domains.
    #     """
    #     #self.domains.add(strand)  # domains are tracked by the graph, self.nodes()
    #     self.domains_by_name[domain.name].add(domain)
    #     if update_graph:
    #         # strand is also a (linear) graph of domains:
    #         self.add_node(domain)
    #         # edge_attrs = {'weight': len(domain), 'len': len(domain), 'type': PHOSPHATEBACKBONE_INTERACTION}
    #         # self.ends5p3p_graph.add_path(domain.end5p, domain.end3p, edge_attrs)
    #
    # def remove_domain(self, domain, update_graph=False):
    #     """ Remove a single domain, updating self.domains_by_name. """
    #     self.domains_by_name[domain.name].remove(domain)
    #     if update_graph:
    #         # strand is also a (linear) graph of domains:
    #         self.remove_node(domain)
    #         # self.ends5p3p_graph.remove_nodes_from([domain.end5p, domain.end3p])

    # @property
    # def domains(self):
    #     return list(self.node.keys())


    def print_history(self, history=None, level=0, indent_str="    ", search_str=None, limit=20, totlimit=100,
                      reverse=True, last_first=True):
        print("History for %r (%s, %s)" % (self, "reversed" if reverse else "", "last-first" if last_first else ""))
        if history is None:
            history = self.history
        if last_first ^ reverse: # last_first != reverse, XOR
            entries = self.gen_history_records(history=history, level=level, search_str=search_str,
                                               limit=limit, totlimit=totlimit, reverse=reverse)
        else:
            entries = reversed(list(self.gen_history_records(
                history=history, level=level, search_str=search_str,
                limit=limit, totlimit=totlimit, reverse=reverse)))
        for (level, entry) in entries:
            print(indent_str*level + entry)

    def history_str(self, history=None, level=0, indent_str="    ", search_str=None, sep="\n",
                    limit=20, totlimit=100, reverse=False):
        return sep.join((indent_str*level + entry) for level, entry
            in self.gen_history_records(history=history, level=level, search_str=search_str,
                                        limit=limit, totlimit=totlimit, reverse=reverse))

    def gen_history_records(self, history=None, level=0, search_str=None, limit=20, totlimit=1000, reverse=False):
        if history is None:
            history = self.history
        if limit and limit < len(history):
            org_length = len(history)
            history = history[-limit:] # Makes a slice copy
            history[0] = "(...history truncated to %s of %s entries...)" % (len(history), org_length)
        if reversed:
            history = reversed(history)
        for entry in history:
            totlimit -= 1
            if totlimit < 0:
                yield "(...totlimit reached, breaking of here...)"
                break
            if isinstance(entry, str):
                if search_str is None or search_str in entry:
                    yield (level, entry)
            else:
                # Returning a final value from a generator and obtaining it with "yield from" is
                # a new feature of python 3.3. Not available in pypy yet.
                nextgen =  self.gen_history_records(history=entry, level=level+1,
                                                    search_str=search_str,
                                                    limit=limit-1, totlimit=totlimit)
                # if sys.version_info > (3, 2):
                # totlimit = yield from nextgen
                # else:
                for val in nextgen:
                    yield val
        # if sys.version_info > (3, 2):
        #     return totlimit

    # @state_should_change
    def add_strand(self, strand, update_graph=False):
        """ We keep track of strands for use with fingerprinting, etc. """
        # printd("%r: Adding strand %r..." % (self, strand))
        self.strands.add(strand)
        self.strands_by_name[strand.name].add(strand)
        self.strand_species_counter[strand.name] += 1
        strand.complex = self
        for domain in strand.domains:
            self.domains_by_name[domain.name].add(domain)
            self.domain_species_counter[domain.name] += 1
        if update_graph:
            # strand is also a (linear) graph of domains:
            if self.is_directed() and not strand.is_directed():
                strand = nx.MultiDiGraph(strand)
            self.add_nodes_from(strand.nodes(data=True))
            self.add_edges_from(strand.edges(data=True, keys=True))
            # self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph)
            # self.strand_graph.add_node(strand)
        self._historic_strands.append(sorted(self.strands, key=lambda s: (s.name, s.suid)))
        self.N_strand_changes += 1
        # # self.history.append("add_strand: Adding strand %r (update_graph=%s)" % (strand, update_graph))
        self.reset_state_fingerprint()


    # @state_should_change
    def remove_strand(self, strand, update_graph=False, update_edge_pairs=True):
        """ We keep track of strands for use with fingerprinting, etc. """
        # printd("%r: Removing strand %r..." % (self, strand))
        all_removed_hybridization_pairs, all_removed_stacking_pairs = set(), set()
        self.strands.remove(strand)
        self.strands_by_name[strand.name].remove(strand)
        self.strand_species_counter[strand.name] -= 1
        if self.strand_species_counter[strand.name] == 0:
            del self.strand_species_counter[strand.name]
        if strand.complex == self:
            strand.complex = None
        for domain in strand.domains:
            self.domains_by_name[domain.name].remove(domain)
            self.domain_species_counter[domain.name] -= 1
            if self.domain_species_counter[domain.name] < 1:
                del self.domain_species_counter[domain.name]
            if update_edge_pairs:
                obsolete_hybridization_pairs = {pair for pair in self.hybridized_pairs
                                                if domain in pair}
                self.hybridized_pairs -= obsolete_hybridization_pairs
                obsolete_stacking_pairs = {pair for pair in self.stacked_pairs if domain in pair}
                self.stacked_pairs -= obsolete_stacking_pairs
                all_removed_hybridization_pairs |= obsolete_hybridization_pairs
                all_removed_stacking_pairs |= obsolete_stacking_pairs
        if update_graph:
            # strand is also a (linear) graph of domains:
            self.remove_nodes_from(strand)
            # self.ends5p3p_graph.add_nodes_from(strand.ends5p3p_graph)
            # self.strand_graph.remove_node(strand)
        self._historic_strands.append(sorted(self.strands, key=lambda s: (s.name, s.suid)))
        self.N_strand_changes += 1
        # # self.history.append("remove_strand: Removing strand %r (update_graph=%s)" % (strand, update_graph))
        return all_removed_hybridization_pairs, all_removed_stacking_pairs


    # @state_should_change
    def add_strands(self, strands, update_graph=False):
        """ Strands must be a set. """
        # printd("%r: Adding strands %s..." % (self, strands))
        for strand in strands:
            self.strands_by_name[strand.name].add(strand)
            self.strand_species_counter[strand.name] += 1
            strand.complex = self
            for domain in strand.domains:
                self.domains_by_name[domain.name].add(domain)
                self.domain_species_counter[domain.name] += 1
            if update_graph:
                if self.is_directed() and not strand.is_directed():
                    strand = nx.MultiDiGraph(strand) # to_directed() returns a DEEP copy of the data. We DO NOT want that.
                # strand is also a (linear) graph of domains:
                self.add_nodes_from(strand.nodes(data=True))
                self.add_edges_from(strand.edges(data=True, keys=True))
                # self.ends5p3p_graph.add_edges_from(strand.ends5p3p_graph)
        if not isinstance(strands, set):
            strands = set(strands)
        self.strands |= strands
        self._historic_strands.append(sorted(self.strands, key=lambda s: (s.name, s.suid)))
        self.N_strand_changes += 1
        # # self.history.append("add_strands: Adding strands %s (update_graph=%s)" % (strands, update_graph))


    # @state_should_change
    def remove_strands(self, strands, update_graph=False, update_edge_pairs=True):
        """ Strands must be a set. """
        # printd("%r: Removing strands %s..." % (self, strands))
        # # self.history.append("remove_strands: Removing strands %s (update_graph=%s)" % (strands, update_graph))
        all_removed_hybridization_pairs, all_removed_stacking_pairs = set(), set()
        for strand in strands:
            self.strands_by_name[strand.name].remove(strand)
            self.strand_species_counter[strand.name] -= 1
            if self.strand_species_counter[strand.name] < 1:
                del self.strand_species_counter[strand.name]
            if strand.complex == self:
                strand.complex = None
            for domain in strand.domains:
                self.domains_by_name[domain.name].remove(domain)
                self.domain_species_counter[domain.name] -= 1
                if self.domain_species_counter[domain.name] < 1:
                    del self.domain_species_counter[domain.name]
                if update_edge_pairs:
                    obsolete_hybridization_pairs = {pair for pair in self.hybridized_pairs
                                                    if domain in pair}
                    obsolete_stacking_pairs = {pair for pair in self.stacked_pairs if domain in pair}
                    self.hybridized_pairs -= obsolete_hybridization_pairs
                    self.stacked_pairs -= obsolete_stacking_pairs
                    all_removed_hybridization_pairs |= obsolete_hybridization_pairs
                    all_removed_stacking_pairs |= obsolete_stacking_pairs
            if update_graph:
                # strand is also a (linear) graph of domains:
                self.remove_nodes_from(strand)
                # self.ends5p3p_graph.remove_nodes_from(strand.ends5p3p_graph)
                # self.strand_graph.remove_node(strand)
        if not isinstance(strands, set):
            strands = set(strands)
        self.strands -= strands
        self._historic_strands.append(sorted(self.strands, key=lambda s: (s.name, s.suid)))
        self.N_strand_changes += 1
        return all_removed_hybridization_pairs, all_removed_stacking_pairs


    # @state_should_change
    def add_hybridization_edge(self, domain_pair):
        # self.history.append("add_hybridization_edge: domain_pair = %s" % (domain_pair,))
        domain1, domain2 = domain_pair
        self.add_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        if self.is_directed():
            self.add_edge(domain2, domain1, key=HYBRIDIZATION_INTERACTION)
        self.hybridized_pairs.add(frozenset(domain_pair))
        self._hybridization_fingerprint = None
        self._state_fingerprint = None

    # @state_should_change
    def remove_hybridization_edge(self, domain_pair):
        # self.history.append("add_hybridization_edge: domain_pair = %s" % (domain_pair,))
        domain1, domain2 = domain_pair
        self.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        if self.is_directed():
            self.remove_edge(domain2, domain1, key=HYBRIDIZATION_INTERACTION)
        self.hybridized_pairs.remove(frozenset(domain_pair))
        self._hybridization_fingerprint = None
        self._state_fingerprint = None

    # @state_should_change
    def add_stacking_edge(self, stacking_pair):
        """
        Stacking pair must be tuple ((h1end3p, h2end5p), (h2end3p, h1end5p))
        or frozenset((h1end3p, h2end5p), (h2end3p, h1end5p)).
        """
        # self.history.append("add_stacking_edge: stacking_pair = %s" % (stacking_pair,))
        (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
        self.add_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION)
        self.add_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION)
        self.stacked_pairs.add((h1end3p.domain, h1end5p.domain))
        self.stacked_pairs.add((h2end3p.domain, h2end5p.domain))
        self._stacking_fingerprint = None
        self._state_fingerprint = None

    # @state_should_change
    def remove_stacking_edge(self, stacking_pair):
        # self.history.append("remove_stacking_edge: stacking_pair = %s" % (stacking_pair,))
        (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
        self.remove_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION)
        self.remove_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION)
        self.stacked_pairs.remove((h1end3p.domain, h1end5p.domain))
        self.stacked_pairs.remove((h2end3p.domain, h2end5p.domain))
        self._stacking_fingerprint = None
        self._state_fingerprint = None



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
        ## TODO: Add check that no two domains in the complex has the same in_complex_identifier:
        #in_complex_id_counts = Counter([d.in_complex_identifier() for d in self.nodes()])
        nonzero_icids = [icid for icid in [d.in_complex_identifier() for d in self.nodes()] if icid != 0]
        in_complex_id_counts = Counter(nonzero_icids)
        # if 0 in in_complex_id_counts:
        #     # icid of 0 means "this domain is the only of its kind, no need to calculate an icid."
        #     del in_complex_id_counts[0]
        if any(count > 1 for icid, count in in_complex_id_counts.items()):
            self.adjust_icid_radius_or_use_instance(in_complex_id_counts)
        ## TODO: Add fall-back to using domain instances (rather than state species) for fingerprint.
        ## This will not be useful for caching between complexes, but might still be used within the same
        ## complex, as long as we don't have strand swapping.
        if self._state_fingerprint is None:
            ## TODO: Re-enable hashing when I'm done debugging
            ## Must not include anything unique to this complex instance such as str(self)
            self._state_fingerprint = hash((
                self.strands_fingerprint(),
                self.hybridization_fingerprint(),
                self.stacking_fingerprint()  ## TODO: Implement stacking
                )) % 100000  # TODO: Remove modulus when done debugging.
            # print("\nRe-calculated state fingerprint for %r: %s" % (self, self._state_fingerprint))
            #pdb.set_trace()
            # historic fingerprints are now added only through assert_state_change:
            # self._historic_fingerprints.append((self._state_fingerprint,))  # TODO: Remove when done debugging.
            # self.history.append("state_fingerprint: Calculated fingerprint: %r" % (self._state_fingerprint,))
        return self._state_fingerprint

    def adjust_icid_radius_or_use_instance(self, in_complex_id_counts=None, n_tries=3):
        if in_complex_id_counts is None:
            # icid of 0 means "this domain is the only of its kind, no need to calculate an icid."
            # We can have as many of these as we'd like.
            nonzero_icids = [icid for icid in [d.in_complex_identifier() for d in self.nodes()] if icid != 0]
            in_complex_id_counts = Counter(nonzero_icids)
        while n_tries > 0:
            # if 0 in in_complex_id_counts:
            #     del in_complex_id_counts[0]
            self.icid_radius *= 2 # Double the range (up to n_tries=4 times: 2, 4, 8, 16 times original)
            for domain in self.nodes():
                domain.state_change_reset(reset_complex=False)
            nonzero_icids = [icid for icid in [d.in_complex_identifier() for d in self.nodes()] if icid != 0]
            in_complex_id_counts = Counter(nonzero_icids)
            n_tries -= 1
            if all(count < 2 for icid, count in in_complex_id_counts.items()):
                print("adjust_icid_radius_or_use_instance: Unique icid found at icid_radius %s for %r" %
                      (self.icid_radius, self))
                pdb.set_trace()
                break
        else:
            # Increasing icid range did not help; fall back to using domain instances...
            # printd("Complex: Increasing icid range did not help; falling back to using domain instances...")
            # If Complex.icid_use_instance has been set to True, fall back to using domain instances as
            # in_complex_identifier (icid) for *all* domains. If icid_use_instance is not True but is
            # boolean True, it is assumed to be a set of domains for which to use domain instance as icid.
            self.icid_use_instance = True
            for domain in self.nodes():
                domain.state_change_reset(reset_complex=False)
            #pdb.set_trace()

    def assert_state_change(self, reaction_spec_pair, reaction_attr, reset=False,
                            expected_state_fingerprint=None):
        """
        Will
        1) Make a new state_fingerprint and add it to historic_fingerprints list,
            asserting that the state really has changed.
        2) Append reaction_pair to reaction_deque and check for reaction cycles.

        Arguments:
        :reaction_spec_pair: The reaction_spec_pair for the reaction that caused the change to this new state.

        Returns:
            new_state_fingerprint, microcycle

        If a micro-cycle containing reaction_pair is found, microcycle is a slice of reaction_deque with reaction_pairs:
            [reaction_pair, pair2, pair3, ..., reaction_pair]
        where reaction_pair is the same reaction pair as the one provided as first arg, :reaction_pair:
        Otherwise, microcycle will be None.

        ### Micro-cycles ###
        Regarding reaction micro-cycles, do we record complex states:
            state1 -> state2 -> state1
        or reaction_pairs:
            pair1 -> pair2 -> pair1
        or both:
                   pair1        pair2
            state1 ----> state2 ----> state1
        I think reaction_pairs are the most useful, since that is what we want to throttle.
        This method is invoked after every reaction. Performance is important.
        Consider using a heap or blist instead of standard list:
            http://kmike.ru/python-data-structures/
        ### Edit: We no longer detect micro-cycles, only reaction invocations ###

        """
        if reset:
            ## We could reset state fingerprints (strands, hybridizations, stackings);
            ## we can use reaction_attr.reaction_type to know whether to reset hyb or stack,
            ## but we would need to probe result['case'] in order to know whether to reset
            ## strands/backbone fingerprint.
            self.reset_state_fingerprint()
        assert reaction_spec_pair != self.reaction_deque[-1]
        state_fingerprint = self.state_fingerprint()
        assert state_fingerprint != self._historic_fingerprints[-1]
        if expected_state_fingerprint is not None:
            assert state_fingerprint == expected_state_fingerprint
        self._historic_fingerprints.append(state_fingerprint)
        # print("assert_state_change invoked for %r with reaction_spec_pair %s" % (self, reaction_spec_pair))
        # pprint(locals())
        # print("historic fingerprints:", self._historic_fingerprints)
        # print("reaction deque: ", self.reaction_deque)

        ## Check reaction_attrs:
        if reaction_spec_pair in self.previous_reaction_attrs:
            assert self.previous_reaction_attrs[reaction_spec_pair] == reaction_attr
        else:
            self.previous_reaction_attrs[reaction_spec_pair] = reaction_attr

        ## Increment reaction count:
        self.reaction_invocation_count[reaction_spec_pair] += 1

        # pdb.set_trace()
        ## Detect reaction cycle:  (obsolete)
        # if reaction_spec_pair in self.reaction_deque[-self.reaction_deque_size:]:
        #     # Using x in list is about 40% faster than list.index(x) - note: deque don't have an index because deques are not for random access.
        #     # This is a little complex since we want to make sure we have the highest index.
        #     # Note: list.index(x) is very slow if x is at the end of a large list - time complexity O(n)
        #     #highest_idx = self.reaction_deque[-self.reaction_deque_size:-1:-1]
        #     # change slice start to -2 if :reaction_pair: has already been added!
        #     highest_idx = self.reaction_deque[-1:-self.reaction_deque_size-1:-1].index(reaction_spec_pair) - 1
        #     # While loop solution is much slower especially if several loops are required.
        #     highest_idx2 = self.reaction_deque_size
        #     while True:
        #         try:
        #             highest_idx2 = self.reaction_deque.index(reaction_spec_pair, highest_idx2)
        #         except ValueError:
        #             break
        #     assert highest_idx2 == len(self.reaction_deque) - highest_idx
        #     assert self.reaction_deque[highest_idx2] == self.reaction_deque[highest_idx] == reaction_spec_pair
        #     self.reaction_deque.append(reaction_spec_pair)
        #     return self.reaction_deque[highest_idx:]
        # else:
        #     self.reaction_deque.append(reaction_spec_pair)
        #     return None
        return state_fingerprint, None


    def get_all_fingerprints(self):
        return (self._state_fingerprint, self._strands_fingerprint, self._stacking_fingerprint)


    def reset_state_fingerprint(self, reset_strands=True, reset_hybridizations=True, reset_stacking=False,
                                reset_domains=True):
        ## TODO: Make sure reset_stacking is properly set when required!
        # self.history.append("reset_state_fingerprint: Unsetting fingerprints: %r" % (locals(),))
        self._state_fingerprint = None
        if reset_strands:
            self._strands_fingerprint = None
        if reset_hybridizations:
            self._hybridization_fingerprint = None
        if reset_stacking:
            self._stacking_fingerprint = None
        if reset_domains:
            for domain in self.domains():
                domain.state_change_reset(reset_complex=False)

    # def domains_gen(self):
    #     return (domain for strand in self.strands for domain in strand.domains)

    def strands_species_count(self):
        """ Count the number of strand species. Used as part of complex state finger-printing. """
        species_counts = {}
        domain_species_counts = {}
        for strand in self.strands:
            if strand.name not in species_counts:
                species_counts[strand.name] = 1
            else:
                species_counts[strand.name] += 1
            for domain in strand.domains:
                if domain.name not in domain_species_counts:
                    domain_species_counts[domain.name] = 1
                else:
                    domain_species_counts[domain.name] += 1
        # Remove entries with zero count:
        depleted = [sname for sname, count in self.strand_species_counter.items() if count < 1]
        for sname in depleted:
            del self.strand_species_counter[sname]

        assert species_counts == self.strand_species_counter

        # Remove entries with zero count:
        depleted = [name for name, count in self.domain_species_counter.items() if count < 1]
        for name in depleted:
            del self.domain_species_counter[name]
        assert domain_species_counts == self.domain_species_counter
        #return species_counts
        # Used for hashing, so return a hashable frozenset(((specie1, count), ...))
        return frozenset(species_counts.items())

    def strands_fingerprint(self):
        """ Create a finger-print of the current strands (species). """
        if not self._strands_fingerprint:
            ## TODO: Re-add hashing when I'm done debugging
            # self._strands_fingerprint = hash(self.strands_species_count())
            self._strands_fingerprint = self.strands_species_count()
            # self.history.append("strands_fingerprint: Calculated strands fingerprint: %r" % (self._strands_fingerprint,))
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
            ## TODO: Hey, uhm, does this really work if we have multiple copies of the same domain specie?
            ## If we have multiple, I think we might need to COUNT. And even then,
            ## If we have hybridizations A#1:a#2 and A#3:a#4,
            ## how do we tell that apart from A#1:a#4 and A#3:a#2, if all we have is
            ## ({A, a}, 2) = "There is two A:a connections".
            ## Use (domain.specie, domain.in_complex_identifier) instead?
            edge_pairs = [frozenset((d.domain_strand_specie, d.in_complex_identifier()) for d in edge)
                          for edge in self.hybridization_edges()]
            edgesfs = frozenset(edge_pairs)
            # If complex is a directed graph, we add edges both back and forth.
            # Thus, we get two edges between each pair of domains.
            # We expect to get half the size if we make a set:
            if len(edgesfs) != len(edge_pairs)/2:
                print("Domain.hybridization_fingerprint: WARNING! len(edgesfs) != len(edge_pairs) (%s vs %s)" %
                      (len(edgesfs), len(edge_pairs)))
                if len(edge_pairs) < 10:
                    print("edge_pairs:")
                    pprint(edge_pairs)
                    print("edgesfs:")
                    pprint(edgesfs)
            ## TODO: Re-add hashing when I'm done debugging
            self._hybridization_fingerprint = edgesfs # hash(edgesfs)
            # self.history.append("hybridization_fingerprint: Calculated hybridization fingerprint: %r" % (self._hybridization_fingerprint,))
        return self._hybridization_fingerprint


    def hybridization_edges(self):
        """
        How to get a list of domain hybridization?
        * self.edges() will generate a list of all connections: backbone, hybridization and stacking.

        """
        # you can loop over all edges and filter by type:
        edges = [(d1, d2) for d1, d2, key, interaction in self.edges(keys=True, data='interaction')
                 if key == HYBRIDIZATION_INTERACTION or interaction == HYBRIDIZATION_INTERACTION]
        # or maybe it is easier to keep a dict with hybridization connections:
        # hyb_edges = self.hybridized_domains
        # if self.is_directed():
        assert set(frozenset(tup) for tup in edges) == self.hybridized_pairs
        # else:
        #     assert set(frozenset(tup) for tup in edges) == self.hybridized_pairs
        return edges


    def stacking_fingerprint(self):
        """
        Return a stacking fingerprint for use with caching.
        A stacked_pair is: {(h1end3p, h2end5p), (h2end3p, h1end5p)}
        """
        if not self._stacking_fingerprint:
            # This is using (d1, d2) tuple rather than {d1, d2} frozenset, since directionality matters.
            # edgesfs = frozenset(frozenset((e.domain.domain_strand_specie, e.end) for e in edge)
            #                     for edge in self.stacking_edges())
            # NOTE: STACKING EDGES HAVE DIRECTIONALITY! h1: d1end3p -> d2.end5p - DO NOT USE FROZENSET.
            edge_pairs = [tuple((d.domain_strand_specie, d.in_complex_identifier()) for d in edge)
                          for edge in self.stacking_edges()]
            edgesfs = frozenset(edge_pairs)
            if len(edgesfs) != len(edge_pairs):
                # Primitive check - not all edges are unique. Should not happen. But is not enough to avoid problems!
                print("Domain.hybridization_fingerprint: WARNING! len(edgesfs) != len(edge_pairs) (%s vs %s)" %
                      (len(edgesfs), len(edge_pairs)))
                assert len(edgesfs) == len(edge_pairs)
            ## TODO: Re-enable hashing when I'm done debugging:
            self._stacking_fingerprint = edgesfs # hash(edgesfs)
            # self.history.append("stacking_fingerprint: Calculated stacking fingerprint: %r" % (self._stacking_fingerprint,))
        return self._stacking_fingerprint


    def stacking_edges(self):
        """
        stacking_edges vs stacking_ends:
         - stacking_edge is a tuple: (h1end3p, h1end5p)
         - stacking pair is a pair of stacking edges:
            {(h1end3p, h2end5p), (h2end3p, h1end5p)}

        How to get a list of domain hybridization?
        Note: stacking is directional and cannot be saved within an undirected graph.
        Perhaps I should just switch to using a directed graph?
        """
        # you can loop over all edges and filter by type:
        # But need directionality to know which end of the domain is pairing.
        # New: Complex is now a directed graph. It won't work for analysis,
        # but we have system graphs for that, so let's try.
        stack_edges = [(d1, d2) for d1, d2, key, interaction in self.edges(keys=True, data='interaction')
                       if key == STACKING_INTERACTION or interaction == STACKING_INTERACTION]
        # Only if Complex is a DiGraph:
        assert set(stack_edges) == self.stacked_pairs and len(set(stack_edges)) == len(self.stacked_pairs)
        # For now, I have to keep a dict with hybridization connections:
        return stack_edges
        #return self.stacked_pairs


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
        edges = [(s, t, key, attrs) for s, t, key, attrs in self.edges(data=True, keys=True)
                 if attrs['interaction'] in (HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)]
        subg = self.subgraph(hybridized_domains)
        subg.add_edges_from(edges)
        return subg

    def draw_graph_and_save(self, outputfn, node_labels=True, **kwargs):
        if node_labels is True:
            node_labels = {node: node.instance_name for node in self.nodes()}
        draw_graph_and_save(self, outputfn, labels=node_labels, **kwargs)


    def fqdn(self):
        """ Return a fully-qualified name. """
        return "C[%s]" % (self.cuid)
        # return "[C#%s]" % (self.cuid)

    def __repr__(self):
        # return "%s[%s]" % (self.name, self.ruid % 100)
        return "Complex[%s] at %s" % (self.cuid, hex(id(self)))

    def __str__(self):
        # return "%s[%s]" % (self.name, self.ruid % 100)
        # String representation should be invariant through the life-time of an object:
        return "C[%s]" % (self.cuid)
        #return "[C#%s]" % (self.cuid)



# class SuperComplex(ConnectedMultiGraph):
#     """
#     Each node is a complex.
#     """
#     def __init__(self, **kwargs):
#         super().__init__(**kwargs)
#         self.scuid = next(supercomplex_sequential_id_gen)
#         # Use this or self.nodes() ??
#         self.complexes = set()
#         self.strands = set() # TODO: Add support for non-complexed strands in supercomplex
#         # Use this or self.edges() ??
#         # No, each node is a complex, stacking_pairs are pairs of:
#         # {(h1end3p, h2end5p), (h2end3p, h1end5p)}
#         self.stacking_pairs = set()
#         # Edge could be:
#         # c1, c2, key=bluntend_pair
#         self._state_fingerprint = None
#
#     def state_fingerprint(self):
#         if self._state_fingerprint is None:
#             hashes = frozenset(hash((cmplx.strands_fingerprint(),
#                                      cmplx.hybridization_fingerprint(),
#                                      cmplx.stacking_fingerprint()))
#                                for cmplx in self.complexes)
#             self._state_fingerprint = hashes
#         return self._state_fingerprint
#
#     def reset_state_fingerprint(self):
#         self._state_fingerprint = None
