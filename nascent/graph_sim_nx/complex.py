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

# pylint: disable=C0103,W0142

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
import numpy as np
# Debug imports:
import pdb
import inspect
import traceback

# Relative imports
from .connected_multigraph import ConnectedMultiGraph
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .nx_utils import draw_graph_and_save
from .constants import (
    PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION,
    DIRECTION_SYMMETRIC, DIRECTION_UPSTREAM, DIRECTION_DOWNSTREAM,
)
from .debug import printd, pprintd, do_print

# Module-level constants and variables:
MAX_ICID_RADIUS = 10
MAX_COMPLEX_SYMMETRY = 10
# Sequential id generators: get next sequential id with next(make_helix_sequential_id)
make_sequential_id = sequential_number_generator()
helix_sequential_id_generator = sequential_number_generator()
# supercomplex_sequential_id_gen = sequential_number_generator()
loopid_sequential_number_generator = sequential_number_generator()





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
    Currently domain-level representation: Each node is a domain.
    ## TODO: Would it make sense to have this be a regular graph of DomainEnds instead of a MultiGraph of domains?
    ##       (see discussion below).

    ## TODO: (Optimization) Instead of creating a new complex whenever two strands hybridize, consider having a
    ##       "delegation" system where one strand's complex simply delegates responsibility to another strand's complex.

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

    Note: We are currently relying on Complex being a *directed* MultiDiGraph
    for e.g. stacking_edges() where we rely on stacking interactions being
    *from* the 5' upstream domain *to* the 3' downstream domain and guaranteed to not be the other way around.
    This is only needed because complex is a domain-level graph and would not be needed for a DomainEnd-level graph.
    (I'm also maintaining the independent self.stacked_pairs which has the same content..)
    TODO: Consolidate use of graph vs stacked_pairs/etc for storing interaction information!

    It might be a lot simpler to have a DomainEnd graph as the basis...

    Discussion: How to generate proper state fingerprints?
    """
    def __init__(self, data=None, strands=None, origin="o", reaction_deque_size=1000): #, reaction_spec_pair=None):
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
        # domain.in_complex_identifier. It will not be 100% certain/accurate, but it should be close enough.
        # Counters in python 3.3+ have nice unary operators. E.g. ```+mycounter``` to get all elements with positive
        # count, and ```-mycounter``` to get all negative, by adding mycounter to a new empty Counter.
        self.domain_species_counter = Counter()
        # We are already calculating a strand-species count, might as well make it permanent running counter:
        self.strand_species_counter = Counter()
        self.domain_strand_specie_counter = Counter()   # Entries are (strand.name, domain.name)

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

        ### Complex energy: ###
        # Includes: Volume energies (for every bi-molecular reaction), Shape/Loop energy and of course,
        # stacking- and hybridization energies. For each, we have both dH and dS.
        # self.energies_dHdS = {'volume': [0, 0], 'shape': [0, 0],
        #                       'stacking': [0, 0], 'hybridization': [0, 0]}
        # Enthalpies in self.energies_dHdS[0], entropies in self.energies_dHdS[1]
        # self.energies_dHdS = [{'volume': 0.0, 'shape': 0.0, 'stacking': 0.0, 'hybridization': 0.0},
        #                       {'volume': 0.0, 'shape': 0.0, 'stacking': 0.0, 'hybridization': 0.0}]
        # Ordered by contribution then enthalpy/entropy:
        self.energies_dHdS = {'volume': [0., 0.], 'shape': [0., 0.], 'stacking': [0., 0.], 'hybridization': [0., 0.]}
        # Sum with:[sum(vals) for vals in zip(*cmplx.energies_dHdS.values())]

        self.energy_total_dHdS = [0, 0]
        # energies_dHdS = np.array([(0.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 0.0)],
        #     dtype=[('volume', float), ('shape', float), ('stacking', float), ('hybridization', float)])
        # Edit: numpy's "structured array" is almost useless. Use ndarray if you want speed.
        # energies_dHdS = np.ndarray((4,2))  # energies_dHdS[<contribution_idx>, <dH or dS>]

        ## Per-interaction energies:
        # self.volume_energies = {}  # strand: (dH=0, dS) # Edit: Just asserting that dS_vol * n_strands = dS_vol_total
        # self.volume_energy = # better to just use energies_dHdS['volume'] ?
        ## TODO: Consolidate whether you are using the complex graph to store interactions or per-interaction-type dicts
        ## Note: You have a global domainend system graph which you could also use to store energies...
        self.hybridization_energies = {} # frozenset({dom1, dom2}): (dH, dS)   # Same as self.hybridized_pairs
        self.stacking_energies = {} # frozenset((h1end3p, h2end5p), (h1end3p, h2end5p)): (dH, dS) # self.stacked_pairs
        self.loop_energies = {} # loopid: (dH, dS)
        self.energy_contributions = {'hybridization': self.hybridization_energies,
                                     'stacking': self.stacking_energies,
                                     'shape': self.loop_energies,
                                     # 'volume': self.volume_energies,
                                    }
        # Note: Above energy contributions keys are instances, not fingerprints. Do not use in reaction graph.

        ## TODO: Use numpy arrays (or recarrays) for storing energies to make calculations easier.


        ## TODO: Is the above sufficiently unique? Consider replacing with a DomainEnd approach below:
        # self.stacked_domain_end_pairs = set()
        # self.hybridized_domain_ends = set()  # Could probably just be hybridized_domains
        # self.backbone_connected_domain_end_pairs = set()

        # State fingerprint is used for cache lookup
        # To reduce recalculation overhead, it is separated into three parts:
        # state = strands + hybridizations + stacking fingerprints.
        self._state_fingerprint = None
        self._strands_fingerprint = None
        self._hybridization_fingerprint = None
        self._merged_complexes_historic_fingerprints = []
        self._stacking_fingerprint = None
        self._nodes_have_been_labelled = None

        ## DEBUG attributes: ##
        self._historic_strands = []
        self._historic_fingerprints = [0]
        # edit: deque are not good for arbitrary access, using a list
        self.reaction_deque = deque(maxlen=reaction_deque_size) # deque with... (reaction_pair, reaction_attr)?
        #self.reaction_deque = [0]    # deque(maxlen=reaction_deque_size)
        # self.reaction_deque_size = reaction_deque_size
        self.reaction_invocation_count = Counter()  # maps reaction_key => Number of times that reaction has occured.
        self.reaction_throttle_cache = {} # For per-complex throttling
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


        ### Keeping track of loops: ###
        # Dict with loops. Keyed by LoopID ? Maybe a set instead?
        # No, loop is a dict (mutable). We don't have frozendict in Python, but there is the MappingProxyType obj class.
        # Having loops as dicts indexed by an ID is nearly identical to having loops being objects in Python.
        # In Python, objects can be a little slow because they are not just data containers but also have methods etc.
        # In e.g. Julia, having a separate "Loop" type would be faster than Dict. (Julia is ALL about types).
        self.loops = {}  # "loops_by_id; LoopID => loop dict
        self.loop_by_loopid = self.loops # alias.
        # What is a loop? It is primarily made up of InterfaceNodes, because that's the only thing that makes sense.
        # For each loop, we have a dict with entries (suggested):
        #     loopid: uuid for this loop.
        #     path: list of ifnodes, formed by reaction between the first and last element, but considered cyclic!
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
        self.loopid_by_hash = {}   # loop hash => loop_id. See calculate_loop_hash function.
        # Should we use a set or list? We probably have random deletion of loops, so set is arguably better.
        # Delegatee: A person designated to act for or represent another or others.
        # Delegator: A person who is delegating or has delegated responsibility to someone else.
        self.loop_delegations = {}     # loop0id => {loop1id, loop2id}   (delegator => delegatees)
        self.loop_delegations_rev = {} # loop1id => {set of loop0s split by loop1} (delegatee => delegators)

        ## IFNODES (which makes up a path): ##
        self.ifnode_by_hash = {}  # ifnode-state-hash => ifnode-instance.
        # ifnode_by_hash is reset with complex and should be updated at every change to the complex..
        #


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

    @property
    def volume_energy(self):
        # return (self.energies_dHdS[0]['volume'], self.energies_dHdS[1]['volume'])
        return self.energies_dHdS['volume']

    @property
    def volume_entropy(self):
        # return self.energies_dHdS[1]['volume']
        return self.energies_dHdS['volume'][1]

    @volume_entropy.setter
    def volume_entropy(self, entropy):
        self.energies_dHdS['volume'][1] = entropy


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

    # def domains_gen(self):
    #     return (domain for strand in self.strands for domain in strand.domains)



    # @state_should_change
    def add_strand(self, strand, update_graph=False):
        """ We keep track of strands for use with fingerprinting, etc. """
        # printd("%r: Adding strand %r..." % (self, strand))
        self.strands.add(strand)
        self.strands_by_name[strand.name].add(strand)
        self.strand_species_counter[strand.name] += 1
        self.domain_strand_specie_counter.update(d.domain_strand_specie for d in strand.domains)
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
        # self.reset_state_fingerprint()


    # @state_should_change
    def remove_strand(self, strand, update_graph=False, update_edge_pairs=True):
        """ We keep track of strands for use with fingerprinting, etc. """
        # printd("%r: Removing strand %r..." % (self, strand))
        all_removed_hybridization_pairs, all_removed_stacking_pairs = set(), set()
        self.strands.remove(strand)
        self.strands_by_name[strand.name].remove(strand)
        self.strand_species_counter[strand.name] -= 1
        self.domain_strand_specie_counter.subtract(d.domain_strand_specie for d in strand.domains)
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
            self.domain_strand_specie_counter.update(d.domain_strand_specie for d in strand.domains)
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
    def remove_strands(self, strands, update_graph=False, update_edge_pairs=True, update_energies=True):
        """ Strands must be a set. """
        # printd("%r: Removing strands %s..." % (self, strands))
        # # self.history.append("remove_strands: Removing strands %s (update_graph=%s)" % (strands, update_graph))
        # all_removed_hybridization_pairs, all_removed_stacking_pairs = set(), set()
        removed_hybridization_energy, removed_stacking_energy = {}, {}
        removed_volume_energy, removed_loops, removed_loop_energy = {}, {}, {}
        ## TODO: Loop removal not implemented
        for strand in strands:
            self.strands_by_name[strand.name].remove(strand)
            self.strand_species_counter[strand.name] -= 1
            self.domain_strand_specie_counter.subtract(d.domain_strand_specie for d in strand.domains)
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
                    # all_removed_hybridization_pairs |= obsolete_hybridization_pairs
                    # all_removed_stacking_pairs |= obsolete_stacking_pairs
                    for pair in obsolete_hybridization_pairs:
                        removed_hybridization_energy = self.hybridization_energies[pair]
                        del self.hybridization_energies[pair]
                    for pair in obsolete_stacking_pairs:
                        removed_stacking_energy = self.stacking_energies[pair]
                        del self.stacking_energies[pair]
            if update_graph:
                # strand is also a (linear) graph of domains:
                self.remove_nodes_from(strand)
                # self.ends5p3p_graph.remove_nodes_from(strand.ends5p3p_graph)
                # self.strand_graph.remove_node(strand)
            # removed_volume_energy[strand] = self.volume_energies[strand]
            # del self.volume_energies[strand]
        # Edit: A complex cannot know the volume energy in general, but it can guess using:
        self.volume_entropy = self.volume_entropy*((len(self.strands)-len(strands)-1)/(len(self.strands)-1))
        if not isinstance(strands, set):
            strands = set(strands)
        self.strands -= strands
        self._historic_strands.append(sorted(self.strands, key=lambda s: (s.name, s.suid)))
        self.N_strand_changes += 1

        ## TODO: Remove loops (and calculate energy)
        # return all_removed_hybridization_pairs, all_removed_stacking_pairs
        return removed_hybridization_energy, removed_stacking_energy, removed_loop_energy, removed_loops


    # @state_should_change
    def add_hybridization_edge(self, domain_pair, hyb_energy):
        """
        Add a hybridization edge to the graph.
        If hyb_energy is provided, also add an entry to hybridization_energies.
        """
        # self.history.append("add_hybridization_edge: domain_pair = %s" % (domain_pair,))
        if isinstance(domain_pair, frozenset):
            domain1, domain2 = tuple(domain_pair)
        else:
            domain1, domain2 = domain_pair
            domain_pair = frozenset(domain_pair)
        self.add_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION, direction=DIRECTION_SYMMETRIC)
        if self.is_directed():
            self.add_edge(domain2, domain1, key=HYBRIDIZATION_INTERACTION, direction=DIRECTION_SYMMETRIC)
        self.hybridized_pairs.add(domain_pair)
        self.hybridization_energies[domain_pair] = hyb_energy
        ## TODO: This could be used to calculate fingerprints faster, but probably not worth it:
        # hybridized_domain_name_pairs = frozenset(frozenset(d.name for d in pair) for pair in self.hybridized_pairs)
        # if len(hybridized_domain_name_pairs) == len(self.hybridized_pairs):
        #     # Complex has a unique combination of hybridization edges with domain names as nodes:
        #     self._hybridization_fingerprint = hash(hybridized_domain_name_pairs)
        self._hybridization_fingerprint = None
        self._state_fingerprint = None

    # @state_should_change
    def remove_hybridization_edge(self, domain_pair):
        # self.history.append("add_hybridization_edge: domain_pair = %s" % (domain_pair,))
        if isinstance(domain_pair, frozenset):
            domain1, domain2 = tuple(domain_pair)
        else:
            domain1, domain2 = domain_pair
            domain_pair = frozenset(domain_pair)
        self.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        if self.is_directed():
            self.remove_edge(domain2, domain1, key=HYBRIDIZATION_INTERACTION)
        self.hybridized_pairs.remove(domain_pair)
        del self.hybridization_energies[domain_pair]
        self._hybridization_fingerprint = None
        self._state_fingerprint = None


    # @state_should_change
    def add_stacking_edge(self, stacking_pair, stacking_energy):
        """
        Stacking pair must be tuple ((h1end3p, h2end5p), (h2end3p, h1end5p))
        or frozenset((h1end3p, h2end5p), (h2end3p, h1end5p)).
        """
        # self.history.append("add_stacking_edge: stacking_pair = %s" % (stacking_pair,))
        if isinstance(stacking_pair, frozenset):
            (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(stacking_pair)
        else:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
            stacking_pair = frozenset(stacking_pair)
        self.add_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION, direction=DIRECTION_DOWNSTREAM)
        self.add_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION, direction=DIRECTION_DOWNSTREAM)
        if self.is_directed():
            self.add_edge(h1end5p.domain, h1end3p.domain, key=STACKING_INTERACTION, direction=DIRECTION_UPSTREAM)
            self.add_edge(h2end5p.domain, h2end3p.domain, key=STACKING_INTERACTION, direction=DIRECTION_UPSTREAM)
        self.stacked_pairs.add((h1end3p.domain, h1end5p.domain))
        self.stacked_pairs.add((h2end3p.domain, h2end5p.domain))
        self.stacking_energies[stacking_pair] = stacking_energy
        self._stacking_fingerprint = None
        self._state_fingerprint = None


    # @state_should_change
    def remove_stacking_edge(self, stacking_pair):
        """ Remove stacking edge from complex graph and stacked_pairs set.
        Stacking pair must be tuple ((h1end3p, h2end5p), (h2end3p, h1end5p)). """
        # self.history.append("remove_stacking_edge: stacking_pair = %s" % (stacking_pair,))
        if isinstance(stacking_pair, frozenset):
            (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(stacking_pair)
        else:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
            stacking_pair = frozenset(stacking_pair)
        self.remove_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION)
        self.remove_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION)
        if self.is_directed():
            self.remove_edge(h1end5p.domain, h1end3p.domain, key=STACKING_INTERACTION)
            self.remove_edge(h2end5p.domain, h2end3p.domain, key=STACKING_INTERACTION)
        self.stacked_pairs.remove((h1end3p.domain, h1end5p.domain))
        self.stacked_pairs.remove((h2end3p.domain, h2end5p.domain))
        del self.stacking_energies[stacking_pair]
        self._stacking_fingerprint = None
        self._state_fingerprint = None


    def check_complex_energy_change(self, dH_hyb, dS_hyb, dH_stack, dS_stack, dS_shape, dS_volume,
                                    reaction_attr, reacted_pair, volume_entropy):
        """ Update complex energy sub-totals and energy_total_dHdS total. """
        # Note: reaction_attr.is_intra is *always* true for dehybridize/unstack reactions;
        # use result['case'] to determine if the volume energy of the complex is changed.
        # reacted_spec_pair will only occour in self._statedependent_dH_dS when dehybridization_rate_constant
        # has been called. Also: Is this really state dependent? Can't we just let de-hybridization and unstacking
        # be independent of complex state?
        #dH, dS = self._statedependent_dH_dS[reacted_spec_pair]
        # Made this more simple by inverting values once depending on is_forming, when calculating dH_hyb etc.

        # if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        #     cmplx.energies_dHdS[0]['hybridization'] += dH_hyb
        #     cmplx.energies_dHdS[1]['hybridization'] += dS_hyb
        # elif reaction_attr.reaction_type is STACKING_INTERACTION:
        #     cmplx.energies_dHdS[0]['stacking'] += dH_stack
        #     cmplx.energies_dHdS[1]['stacking'] += dS_stack
        # else:
        #     raise ValueError("Un-supported reaction type %s" % reaction_attr.reaction_type)
        # # Cases:
        # #   Formation Case 0/1:   IntRA-STRAND/COMPLEX hybridization.
        # #   Formation Case 2/3/4: IntER-complex hybridization between two complexes/strands.
        # #   Breaking  Case 0/1:   Complex is still intact.
        # #   Breaking  Case 2/3/4: Complex is split in two.
        # if reaction_result['case'] <= 1:
        #     # IntRA-complex reaction. A single complex both before and after reaction.
        #     cmplx.energies_dHdS[1]['shape'] += dS_shape
        # else:
        #     # IntER-strand/complex reaction; Two strands/complexes coming together or splitting up.
        #     cmplx.energies_dHdS[1]['volume'] += dS_volume
        subtotals_before = {k: tuple(v) for k, v in self.energies_dHdS.items()}  # {'volume': (dH, dS)}
        total_before = tuple(self.energy_total_dHdS)  # (dH, dS)
        if dH_hyb:
            assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
            # We group energies_dHdS first by enthalpy/entropy, then by contribution type
            # because then it's easier to summarize total enthalpy/entropy.
            # self.energies_dHdS[0]['hybridization'] += dH_hyb
            # self.energies_dHdS[1]['hybridization'] += dS_hyb
            self.energies_dHdS['hybridization'][0] += dH_hyb
            self.energies_dHdS['hybridization'][1] += dS_hyb
            if reacted_pair:
                if reaction_attr.is_forming:
                    self.hybridization_energies[reacted_pair] = (dH_hyb, dS_hyb)
                else:
                    del self.hybridization_energies[reacted_pair]
        if dH_stack:
            assert reaction_attr.reaction_type is STACKING_INTERACTION
            self.energies_dHdS['stacking'][0] += dH_stack
            self.energies_dHdS['stacking'][1] += dS_stack
            if reacted_pair:
                if reaction_attr.is_forming:
                    self.stacking_energies[reacted_pair] = (dH_stack, dS_stack)
                else:
                    del self.stacking_energies[reacted_pair]
        if dS_shape:
            assert reaction_attr.is_intra is True
            self.energies_dHdS['shape'][1] += dS_shape
            # Need to update contributions from all loops... better to just keep them separate...
        if dS_volume:
            assert reaction_attr.is_intra is False or reaction_attr.is_forming is False
            self.energies_dHdS['volume'][1] += dS_volume
            assert (int(self.energies_dHdS['volume'][1] / dS_volume * (1 if reaction_attr.is_forming else -1))
                    == len(self.strands)-1)

        ## TODO: Use a more appropriate rounding, not just integer. Using tuple to indicate immutable.
        # self.energy_total_dHdS = tuple([int(sum(d.values())) for d in self.energies_dHdS])
        self.energy_total_dHdS = tuple([int(sum(vals)) for vals in zip(*self.energies_dHdS.values())]) # (dH, dS)


    def recalculate_complex_energy(self, volume_entropy, dS_shape=None):
        """
        Re-calculate complex energy based on all individual energy contributions
        (hybridizations, stackings, loops and volume entropies)
        param :volume_entropy: is the entropy gained when *releasing* two components (should be positive).
        """
        # energy_contributions['hybridization'] = {domain_pair: [dH, dS], ...}
        for contrib_key, entries in self.energy_contributions.items():
            if contrib_key == 'shape':
                # we currently do not track loops (TODO). Instead, we calculate shape entropy using reaction dS_shape:
                if dS_shape:
                    self.energies_dHdS[contrib_key][1] += dS_shape
            else:
                self.energies_dHdS[contrib_key] = [sum(vals) for vals in zip(*entries.values())] if entries else [0, 0]
        # energies_dHdS = {'hybridization': [dH, dS], 'stacking': [dH, dS], ...}
        self.volume_entropy = -volume_entropy * (len(self.strands) - 1)
        self.energy_total_dHdS = [int(sum(vals)) for vals in zip(*self.energies_dHdS.values())]


    def check_complex(self, volume_entropy):
        """
        Check that complex is correct. Mostly for debugging.
        param :volume_entropy: is the entropy gained when *releasing* two components (should be positive).
        """
        # energy_contributions['hybridization'] = {domain_pair: [dH, dS], ...}
        for contrib_key, entries in self.energy_contributions.items():
            # e.g. contribution = 'hybridization', dHdS_vals = {frozenset((d1, d2)): (dH, dS)}
            if contrib_key == 'shape': # we currently do not track loops (TODO)
                continue
            assert all(np.isclose(self.energies_dHdS[contrib_key], ([sum(vals) for vals in zip(*entries.values())] if entries else [0, 0])))
        # energies_dHdS = {'hybridization': [dH, dS], 'stacking': [dH, dS], ...}
        # Check volume entropy:
        assert np.isclose(self.volume_entropy, -volume_entropy*(len(self.strands)-1))
        assert all(np.isclose(self.energy_total_dHdS, [int(sum(vals)) for vals in zip(*self.energies_dHdS.values())]))


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
        if self._state_fingerprint is None:
            # Check that all domains in complex have also been reset:
            # assert all(d._specie_state_fingerprint is None for d in self.domains())
            old_domain_fingerprints = {domain: domain._specie_state_fingerprint
                                       for domain in self.domains()
                                       if domain._specie_state_fingerprint is not None}
            if not self._nodes_have_been_labelled:
                self.label_complex_nodes()
                if old_domain_fingerprints:
                    for domain, old_fingerprint in old_domain_fingerprints.items():
                        # If complex state fingerprint is None, then domains cannot have valid state fingerprint.
                        # Question is: Have the old, obsolete domain state fingerprint been used for anything?
                        print("WARNING: domain %s state fingerprint" % (domain, ),
                              "had not been reset prior to calling Complex.state_fingerprint()!")
                        if domain.state_fingerprint() != old_fingerprint:
                            print(" -- AND NEW STATE FINGERPRINT DIFFERS FROM OLD: %s vs %s" %
                                  (domain._specie_state_fingerprint, old_fingerprint))
                        else:
                            print(" -- redundant label_complex_nodes: new fingerprint same as old.")
                        from inspect import currentframe, getframeinfo
                        frameinfo = getframeinfo(currentframe().f_back)
                        print(" -- called from %s:%s" % (frameinfo.filename, frameinfo.lineno))
            else:
                # Check that the labelling is invariant:
                old_icids = {d: d.icid for d in self.domains()}
                self.label_complex_nodes()
                new_icids = {d: d.icid for d in self.domains()}
                assert old_icids == new_icids
            ## State fingerprint must not include anything unique to this complex instance such as str(self)
            # self._state_fingerprint = hash((
            #     self.strands_fingerprint(),
            #     self.hybridization_fingerprint(),
            #     self.stacking_fingerprint()
            #     )) % 100000  # TODO: Remove modulus when done debugging.
            # The above state-fingerprint method was based on a "minimize required calculations" strategy.
            # It works really well as long as there are no symmetric parts in the complex and I could basically
            # just label each edge by {source.name, target.name} then there would not be any identical edges.
            # Obviously, this breaks down really fast, and we have to assign a unique in-complex-identifier
            # to each node/domain. The new node-labelling scheme should be fast for the case above where all
            # (h, {source.name, target.name}) or (b, (source.name, target.name)) edges are unique.
            # However, it must be fully re-calculated for any kind of reaction.
            # And, when using icid labelling, we MUST recalculate state fingerprint from scratch;
            # we cannot re-use any of the old partial hybridization/stacking fingerprints,
            # because the icid for many nodes could have changed.
            # We have to re-calculate domain icids after every reaction, because it can be hard to know if
            # a reaction makes the complex symmetrix.
            # The only case where we may not need to do this is when there is only a single copy of each
            # domain/strand in the complex. (This might actually be a pretty common case).
            # In that case, we can indeed make optimizations by separating the state fingerprint
            # into different unique parts, where a stacking reaction will not affect the other partial fingerprint.

            # Each edges is a tuple of (node-pair, interaction, direction):
            edges = frozenset((frozenset(((d1.strand.name, d1.name, d1.icid), (d2.strand.name, d2.name, d2.icid))),
                               ekey, eattr['direction']) for d1, d2, ekey, eattr in self.edges(keys=True, data=True))
            self._state_fingerprint = hash(edges) % 100000  # TODO: Remove modulus when done debugging.

            # print("\nRe-calculated state fingerprint for %r: %s" % (self, self._state_fingerprint))
            #pdb.set_trace()
            # historic fingerprints are now added only through assert_state_change:
            # self._historic_fingerprints.append((self._state_fingerprint,))  # TODO: Remove when done debugging.
            # self.history.append("state_fingerprint: Calculated fingerprint: %r" % (self._state_fingerprint,))
        return self._state_fingerprint


    # def adjust_icid_radius_or_use_instance(self, in_complex_id_counts=None, n_tries=4):
    #     """
    #     TODO: Consider a better way to do this.
    #     """
    #     if in_complex_id_counts is None:
    #         # icid of 0 means "this domain is the only of its kind, no need to calculate an icid."
    #         # We can have as many of these as we'd like.
    #         nonzero_icids = [icid for icid in [d.in_complex_identifier() for d in self.nodes()] if icid != 0]
    #         in_complex_id_counts = Counter(nonzero_icids)
    #     while n_tries > 0:
    #         # if 0 in in_complex_id_counts:
    #         #     del in_complex_id_counts[0]
    #         self.icid_radius *= 2 # Double the range (up to n_tries=4 times: 2, 4, 8, 16 times original)
    #         for domain in self.nodes():
    #             domain.state_change_reset()
    #         nonzero_icids = [icid for icid in [d.in_complex_identifier() for d in self.nodes()] if icid != 0]
    #         in_complex_id_counts = Counter(nonzero_icids)
    #         n_tries -= 1
    #         if all(count < 2 for icid, count in in_complex_id_counts.items()):
    #             print("adjust_icid_radius_or_use_instance: Unique icid found at icid_radius %s for %r" %
    #                   (self.icid_radius, self))
    #             pdb.set_trace()
    #             break
    #     else:
    #         # Increasing icid range did not help; fall back to using domain instances...
    #         print("Complex: Increasing icid range did not help; falling back to using domain instances...")
    #         # If Complex.icid_use_instance has been set to True, fall back to using domain instances as
    #         # in_complex_identifier (icid) for *all* domains. If icid_use_instance is not True but is
    #         # boolean True, it is assumed to be a set of domains for which to use domain instance as icid.
    #         self.icid_use_instance = True
    #         for domain in self.nodes():
    #             domain.state_change_reset()
    #         #pdb.set_trace()


    def label_complex_nodes(self):
        """
        Make sure all nodes in complex have a unique label.
        Currently, this is done for Domain nodes (not DomainEnds).
        Note that a node label is only state-specific when combined with the parent complex' state_fingerprint!
        """
        ### TODOS:
        ## TODO: Assert that backbone up/downstream are properly affecting neighboring hashrings
        ##       I'm concerned that the edge directionality of the directed MultiDiGraph complex domain_graph
        ##       can interferre.
        ## TODO: There should be a difference between backbone connecting upstream or downstream.
        ##       I've added UP/DOWN STREAM enumerations to constants module.

        self.symmetry_fold = 1
        non_unique_nodes = self.domains()   # list if rebuilding, or set if updating in-place.
        for node in non_unique_nodes:
            node.hash_ring.clear()
            node.icid = node.name  # or use domain_strand_specie = (strand.name, domain.name) ?
            node._sym_break_label = 0
            node._asymmetric_unit_label = 0 # (0,)
            node.symmetric_icids.clear()   # Can be used to find equivalent nodes from one asymmetric unit to another.
        node_icids = {node: node.icid for node in non_unique_nodes}
        node_icid_count = Counter(node_icids.values())
        non_unique_nodes = [node for node in non_unique_nodes if node_icid_count[node.icid] > 1]
        n_previous_nonunique_nodes = len(non_unique_nodes) + 1
        total_loops = 0

        while non_unique_nodes and self.symmetry_fold <= MAX_COMPLEX_SYMMETRY:
            ring_number = 0
            loops_for_this_symmetry_fold = 0
            n_previous_nonunique_nodes += 1 # Make sure we give it a try
            while non_unique_nodes and (len(non_unique_nodes) < n_previous_nonunique_nodes): # or ring_number < MAX_ICID_RADIUS):
                # ring_number =+ 1
                total_loops += 1
                loops_for_this_symmetry_fold += 1
                n_previous_nonunique_nodes = len(non_unique_nodes)
                for node in non_unique_nodes:  # append new hash_ring to each node
                    # use neightbor_rings, which only gives unique neighbors? Nah.

                    # node.hash_ring[ring_number] = hash(frozenset(
                    #     (neighbor.icid, eattr['interaction']) for neighbor, eattr in self.adj[node].items()))

                    # Alternatively, keep track for debugging:
                    node.hash_ring = {
                        neighbor: (ekey,
                                   # No interaction key for hybridization edges; we use edge_keys
                                   # for interactions in Complex and system domain_graph (MultiDiGraph).
                                   # eattr['interaction'],
                                   neighbor.icid, neighbor._sym_break_label, neighbor._asymmetric_unit_label)
                        for neighbor, edges_by_key in self.adj[node].items() for ekey, eattr in edges_by_key.items()
                        # for neighbor, eattr in self.adj[node].items()  # Edit: We have a MultiDiGraph with keys.
                    }
                print("Remaning non-unique nodes: %s" % non_unique_nodes) # , ", ".join(
                print("\nNodes hash_ring: (symmetry fold %s, %s loops for this symmetry, %s loops total)" %
                      (self.symmetry_fold, loops_for_this_symmetry_fold, total_loops))
                pprint({node: node.hash_ring for node in non_unique_nodes})
                for node in non_unique_nodes: # Then update all nodes' icid:
                    # node.icid = node.hash_ring[-1]   # if not keeping track
                    # Make sure to add current node.icid to next icid when hashing:
                    node.icid = hash((node.icid, frozenset(node.hash_ring.values())))  # if keeping track for debugging
                    node.icid %= 100000  ## TODO: Remove modulus when done debugging!
                    hash_ring_values = list(node.hash_ring.values())
                    # Asym label is always a tuple of (length_of_next_asym_tuple, next_asym_tuple, ...)
                    # Asym unit label propagates from the symmetry-breaking labelled node outwards.
                    # asym_labels = [asym_label for ekey, icid, sym_break_label, asym_label in hash_ring_values]
                    asym_labels = [val[3] for val in hash_ring_values]
                    # asym_labels.append(node._asymmetric_unit_label) # Include the node itself in the comparison.
                    # Or, do:
                    max_neighboring_asym_label = max(asym_labels) #, key=lambda asym_label: (len(asym_label), asym_label))
                    if max_neighboring_asym_label > node._asymmetric_unit_label:
                        node._asymmetric_unit_label = max_neighboring_asym_label
                print("Node icid (in-complex-identifiers):")
                pprint({node: node.icid for node in non_unique_nodes})
                print("Node icid symmetry breaking labels:")
                pprint({node: node._sym_break_label for node in non_unique_nodes})
                print("Node._asymmetric_unit_label:")
                pprint({node: node._asymmetric_unit_label for node in non_unique_nodes})

                # Then update the list of non_unique nodes:
                node_icids = {node: node.icid for node in non_unique_nodes}
                node_icid_count = Counter(node_icids.values())
                non_unique_nodes = [node for node in non_unique_nodes if node_icid_count[node.icid] > 1]

            # Are there cases where non_unique_nodes does not decrease, but will still eventually produce uniquely-labelled nodes?
            # If that's the case, we need to relax the len(non_unique_nodes) < n_previous_nonunique_nodes criteria in the while loop.
            # Perhaps or instead of and: (ring_number < max_icid_radius or len(non_unique_nodes) < n_previous_nonunique_nodes) ?
            if non_unique_nodes:
                # We still have non_unique nodes:
                # Need to have a consistent and invariant labelling method.
                # To do this, start by labelling a single of the non_unique_nodes.
                # The selection of this node should be the same for all complexes of the same state.
                sym_break_label, node_uuid, selected_node = min(
                    ((node.name, node._asymmetric_unit_label, node.icid, self.symmetry_fold), node.uuid, node)
                    for node in non_unique_nodes)
                # Increase complex symmetry:
                self.symmetry_fold += 1
                # Also, it would be nice to be able to label the different asymmetric units.
                # Add a "asymetric unit" label to each node of each of the symmetric parts of complex:
                # Using domain._asymmetric_unit_label to propagate this
                selected_node._sym_break_label = sym_break_label
                selected_node._asymmetric_unit_label = self.symmetry_fold*2 # sym_break_label
                print("Increasing symmetry fold to %s by labelling node %s with symmetry-breaking label %s" %
                      (self.symmetry_fold, selected_node, sym_break_label))
                # pdb.set_trace()
                for node in non_unique_nodes:
                    # Add current icid to node.symmetric_icids so you can later find symmetric domains:
                    node.symmetric_icids.append(node.icid)

                # What is the best way to "increment" the asymmetric unit label?
                # Could also be byte string or similar or octal or hex values.
                # Or, instead of tuples, just do regular addition.
                # But, for now I want to retain more info for debugging.
            # end inner "while non_unique_nodes are being reduced" while loop.
        # end outer "try to break symmetry" while loop
        self._nodes_have_been_labelled = True

        if self.symmetry_fold >= 10:
            print("WARNING: Unusually high symmetry fold %s encountered!" % self.symmetry_fold,
                  "In-complex-identifier for reactants used in reaction spec pair (state hashes) may be ambiguous!")
        else:
            if do_print:
                print("All nodes have been uniquely labelled in complex %s with symmetry fold %s in %s loops.\n" %
                      (self, self.symmetry_fold, total_loops))


    def assert_state_change(self, reaction_spec_pair, reaction_attr, reset=True,
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
        print("DEBUG: reset_state_fingerprint called with %s\nFROM:" % (locals(),))
        traceback.print_stack(limit=4)
        if reset: ## TODO: Determine if we should always reset here
            ## We could reset state fingerprints (strands, hybridizations, stackings);
            ## we can use reaction_attr.reaction_type to know whether to reset hyb or stack,
            ## but we would need to probe result['case'] in order to know whether to reset
            ## strands/backbone fingerprint.
            self.reset_state_fingerprint()
            self.rebuild_ifnode_by_hash_index()
            self.rebuild_loopid_by_hash_index()
        try:
            assert reaction_spec_pair != self.reaction_deque[-1] != None
        except IndexError:
            pass
        state_fingerprint = self.state_fingerprint()
        assert state_fingerprint != self._historic_fingerprints[-1] != None
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

        # Update ifnode_by_hash cache (doing it here for now, should be done lazily eventually...)
        self.rebuild_ifnode_by_hash_index()

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
        self.reaction_deque.append((reaction_spec_pair, reaction_attr))
        #     return None
        return state_fingerprint, None


    def get_all_fingerprints(self):
        return (self._state_fingerprint, self._strands_fingerprint, self._stacking_fingerprint)


    def reset_state_fingerprint(self, reset_strands=True, reset_hybridizations=True, reset_stacking=True,
                                reset_domains=True, reset_loops=True):
        """ Properly reset complex state fingerprint (and also all domains, by default). """
        ## TODO: Make sure reset_stacking is properly set when required!
        # self.history.append("reset_state_fingerprint: Unsetting fingerprints: %r" % (locals(),))
        # if self._state_fingerprint is None:
        #     from inspect import currentframe, getframeinfo
        #     frameinfo = getframeinfo(currentframe().f_back)
        #     print("Possible excessive reset_state_fingerprint(%s, %s, %s, %s) called from %s:%s" %
        #           (reset_strands, reset_hybridizations, reset_stacking, reset_domains,
        #            frameinfo.filename, frameinfo.lineno))
        #     frameinfo = getframeinfo(currentframe().f_back.f_back)
        #     print(" -- Previous frameinfo file and line: %s:%s" % (frameinfo.filename, frameinfo.lineno))
        print("DEBUG: reset_state_fingerprint called with %s\nFROM:" % (locals(),))
        traceback.print_stack(limit=4)
        self._state_fingerprint = None
        if reset_strands:
            self._strands_fingerprint = None
        if reset_hybridizations:
            self._hybridization_fingerprint = None
        if reset_stacking:
            self._stacking_fingerprint = None
        self._nodes_have_been_labelled = False
        if reset_domains:
            for domain in self.domains():
                domain.state_change_reset()
        if reset_loops or True:
            self.ifnode_by_hash = None
            self.loopid_by_hash = None
        # Make sure to update Complex's ifnode_by_hash and loopid_by_hash indexes:
        # self.rebuild_ifnode_by_hash_index()
        # self.rebuild_loopid_by_hash_index()



    #### LOOP TRACKING LOGIC: ####

    def get_ifnode_by_hash(self, ifnode_statespec):
        """ Get ifnode instance by it's state-specific hash (fingerprint). """
        # Q: Is it really worth bothering with caching of this? We are only going to use it once, I think.
        try:
            return self.ifnode_by_hash[ifnode_statespec]
        except KeyError:
            print("ERROR: Could not get ifnode by ifnode_statespec %s" % ifnode_statespec)
            print(" -- resetting self.ifnode_by_hash")
            pdb.set_trace()
            # Update the map, check that no ifnode gives the same hash.
            ifnodes = [end.ifnode for domain in self.domains() for end in (domain.end5p, domain.end3p)]
            top_ifnodes = {ifnode for ifnode in ifnodes if ifnode.delegatee is None}
            self.ifnode_by_hash = {ifnode.state_fingerprint(): ifnode for ifnode in top_ifnodes}
            try:
                return self.ifnode_by_hash[ifnode_statespec]
            except KeyError:
                print("ERROR: Still could not get ifnode by ifnode_statespec %s" % ifnode_statespec)
                print(" -- invoking reset_state_fingerprint and then resetting self.ifnode_by_hash")
                self.reset_state_fingerprint()
                ifnodes = [end.ifnode for domain in self.domains() for end in (domain.end5p, domain.end3p)]
                top_ifnodes = {ifnode for ifnode in ifnodes if ifnode.delegatee is None}
                self.ifnode_by_hash = {ifnode.state_fingerprint(): ifnode for ifnode in top_ifnodes}
                return self.ifnode_by_hash[ifnode_statespec]

    def calculate_loop_path_spec(self, loop_path):
        """
        Return the list:
            [ifnode.state_fingerprint() for ifnode in loop_path]
        Using this method will make it easier to cache the ifnode values at a later point
        (although arguably the cached value should be stored in the ifnode - but premature optimization).
        Also, calling this will ensure that the ifnode_by_hash index hash been built when ifnode hashes are used.
        """
        if self.ifnode_by_hash is None:
            self.rebuild_ifnode_by_hash_index()
        return [ifnode.state_fingerprint() for ifnode in loop_path]

    def calculate_loop_hash(self, loop_path_spec):
        """
        Calculate a state-dependent hash for the given loop spec (list of ifnode hashes).
        The calculation algorithm will "cyclirize" the list (connecting the start and end nodes)
        and then use all edges between nodes as the hash:
            hash(frozenset(frozenset((src, tgt)) for src, tgt in edgelist))
        An alternative would be to simply calculate the hash based on the set of all ifnodes:
            hash(frozenset(frozenset(loop_path_spec))
        and hope that is sufficiently unique to distinguish the loop by its hash.
        Is there any reason why it wouldn't be?
        """
        if self.loopid_by_hash is None:
            # If we are ever creating a loop_hash, we want to make sure Complex.loopid_by_hash is also available.
            self.rebuild_loopid_by_hash_index()
        edgelist = list(zip(loop_path_spec[:-1], loop_path_spec[1:]))
        edgelist.append((loop_path_spec[-1], loop_path_spec[0]))
        edges = frozenset(frozenset((src, tgt)) for src, tgt in edgelist)
        return hash(edges)

    def calculate_loop_hash_slow(self, loop_path, loop=None, loopid=None):
        """
        Calculate a state-dependent hash for the given loop.
        Complex state fingerprint and domain icid labelling must be up-to-date,
        but does not otherwise rely on any cached values.
        Use calculate_loop_hash(loop_path_spec) for a fast version of this.
        """
        if loop_path is None:
            if loop is None:
                assert loopid is not None
                loop = self.loop_by_loopid[loopid]
            loop_path = loop['path']
        edgelist = list(zip(loop_path[:-1], loop_path[1:]))
        edgelist.append((loop_path[-1], loop_path[0]))
        # Converting each tuple to frozenset because order doesn't matter.
        edges = frozenset(frozenset((src.state_fingerprint(), tgt.state_fingerprint()))
                          for src, tgt in edgelist)
        return hash(edges)


    def rebuild_ifnode_by_hash_index(self):
        """ Rebuild ifnode_by_hash dict. """
        # Go over all ifnodes (for all domain ends):
        cmplx_ifnodes = [end.ifnode for d in self.domains() for end in (d.end5p, d.end3p)
                         if end.ifnode.delegatee is None]
        self.ifnode_by_hash = {ifnode.state_fingerprint(): ifnode for ifnode in cmplx_ifnodes}


    def recreate_loop_ifnodes_from_spec(self, ifnodes_specs):
        """
        Recreate a loop path with ifnodes (instances) from a list of ifnodes fingerprints.
        """
        ## TODO: Use a Blist to store path-nodes (you will be doing lots of arbitrary inserts/deletions)
        if self.ifnode_by_hash is None:
            self.rebuild_ifnode_by_hash_index()
        # assert all(spec in self.ifnode_by_hash for spec in ifnodes_specs)
        return [self.ifnode_by_hash[spec] for spec in ifnodes_specs]


    def rebuild_loopid_by_hash_index(self, update_ifnodes=None):
        """ Update all loop hashes. Will also update ifnodes if update_ifnodes is True (default). """
        if update_ifnodes or self.ifnode_by_hash is None:
            self.rebuild_ifnode_by_hash_index()
        self.loopid_by_hash = {}
        # pdb.set_trace()
        # old_loop_hash_by_id = {loopid: loop_hash for loop_hash, loopid in self.loopid_by_hash.items()}
        for loopid, loop in self.loop_by_loopid.items():
            new_path_spec = [ifnode.state_fingerprint() for ifnode in loop['path']]
            new_loop_hash = self.calculate_loop_hash(new_path_spec)
            loop['path_spec'] = new_path_spec
            loop['loop_hash'] = new_loop_hash
            self.loopid_by_hash[new_loop_hash] = loopid
            # del self.loopid_by_hash[old_loop_hash_by_id[loopid]] # Checking that the old one was present
        # self.loopid_by_hash = loopid_by_hash # overwrite the old map


    def effectuate_loop_changes(self, loop_effects, is_forming):
        """
        loop_effects is a dict describing what the effects of
        a new loop being formed or an existing loop being broken.
        it should contain a 'changed_loops_by_hash' entry which is a dict with:
        loop_hash:

        Note: The loop hashes applies to the state *before* the loop is formed.

        The timing of this needs to be exactly right. Make sure you keep track of when and where
        hashes (state specs) are re-calculated for ifnodes and loops:
        ifnode state fingerprint is reset and re-calculated:
            IfNode.state_fingerprint() is called by:
                Complex.calculate_loop_hash(loop_path) similarly calls ifnode.state_fingerprint() for all src,tgt.
                Complex.rebuild_loopid_by_hash_index() calls ifnode.state_fingerprint(), which simply calls
            ifnode.state_fingerprint() simply calls:
                ifnode.domain_end.state_fingerprint() for ifnode in IfNode.delegated_edges.keys()
            DomainEnd.state_fingerprint() calls:
                (DomainEnd.domain.state_fingerprint(), DomainEnd.end, DomainEnd.stack_partner is not None)
            Domain.state_fingerprint() uses:
                Domain._specie_state_fingerprint - which is recalculated if None as:
                Domain._specie_state_fingerprint = (dspecie, self.partner is not None, c_state, in_complex_identifier)
            Domain._specie_state_fingerprint reset by (and only by):
                Domain.state_change_reset()
            Domain.state_change_reset is currently called from multiple places:
                Complex.reset_state_fingerprint()
                ReactionMgr.hybridize_and_process (but not stack_and_process?) - debugging...
                ComponentMgr.join/break_complex_at() - if Domain.strand.complex is None
                ReactionMgr.post_reaction_processing - but only for free strands (!)
            Complex.reset_state_fingerprint is called from:
                ComponentMgr.join/break_complex_at()
                Complex.assert_state_change()
                # Complex.get_ifnode_by_hash -- but just for debugging
                # ReactionMgr.hybridize_and_process (but not stack_and_process?) - but only temporarily, for debugging...

        Now, we need to be able to update Complex loops based on the abstract hash-based
        loop_hash and path_spec values loop_effects.
        These hashes are all based on the Complex state *before* performing the reaction.
        Thus, we need to either effectuate_loop_changes(loop_effects) or at least update loop_effects
        before calling Complex.reset_state_fingerprint(), i.e.:
            in ComponentMgr.join/break_complex_at, or
            in ReactionMgr.post_reaction_processing() before calling Complex.assert_state_change.

        I think I will try to effectuate_loop_changes() it in ComponentMgr.join/break_complex_at
            Edit: Is a little tedious since we do not have reaction_spec_pair
            ComponentMgr generally doesn't know much about reactions...
            We DO calculate dHdS for stacking and hybridization in ComponentMgr.hybridize/stack,
            but they are state independent and does not rely on reaction_spec_pair.

        Moved effectuate_loop_changes invocation to ReactionMgr.post_reaction_processing

        Note: Tools to generate call graphs:
            pycallgraph:
            pyan: https://github.com/dafyddcrosby/pyan, https://github.com/davidfraser/pyan
            yapppi: https://bitbucket.org/sumerc/yappi/
        Profiling blog posts:
            TalkPython podcast #28, http://blog.thehumangeo.com/2015/07/28/profiling-in-python/
            PyCharm professional has a very nice profiler with call-graph:
            http://blog.jetbrains.com/pycharm/2015/05/pycharm-4-5-eap-build-141-988-introducing-python-profiler/

        """
        # pdb.set_trace()
        if is_forming:
            # 1a. Add new loop:
            new_loop = loop_effects['new_loop']
            new_loop_id = next(loopid_sequential_number_generator)
            new_loop_hash = new_loop['loop_hash']
            self.loops[new_loop_id] = new_loop
            self.loopid_by_hash[new_loop_hash] = new_loop_id

            # Make sure ifnode_by_hash is up-to-date:
            # TODO: Check if this is needed or even a good idea. Shouldn't the index already be up-to-date?
            # Yes, it is definitely a bad idea to re-build ifnode_by_hash index.
            # Even if most ifnode state hashes haven't changed, the ifnode delegation may have changed,
            # and we use THAT when re-building the index.
            # self.rebuild_ifnode_by_hash_index()

            # Update loop path to use this Complex's ifnodes (instances):
            new_loop['path'] = self.recreate_loop_ifnodes_from_spec(new_loop['path_spec'])

            # Testing for changes in ifnode delegation:
            assert new_loop['path'][0].top_delegate() == new_loop['path'][-1].top_delegate()
            new_loop['path'][0] = new_loop['path'][0].top_delegate()
            new_loop['path'].pop() # Remove the last ifnode (which is also the first after reaction)
            # Make sure path consists of top_delegates only:
            assert [ifnode.top_delegate() for ifnode in new_loop['path']] == new_loop['path']
            if max(Counter(ifnode.top_delegate() for ifnode in new_loop['path']).values()) > 1:
                print("\n\nWARNING: new_loop['path'] has duplicate top_delegate ifnodes:")
                print(new_loop['path'])
                print([ifnode.top_delegate() for ifnode in new_loop['path']])
                print(Counter(ifnode.top_delegate() for ifnode in new_loop['path']).most_common(3))
                pdb.set_trace()
                # Isn't it natural that delegation will have changed? The loop effects were calculated before
                # the reaction was done, but delegation has now been updated.
                # Either (a) Take changes in delegation into account when calculating the path,
                # or (b) update the path to account for delegation changes.
                # (a) is probably hard since we can't know in advance which node is selected as top_delegate.
                # We don't consider the edge between reactants when calculating activity, so we don't need to modify.
                # In fact, we could just say that the loop path is just always


            # Make sure the loop hash haven't changed:
            # new_loop_hash_from_path = self.calculate_loop_hash(new_loop['path'])
            # new_loop_hash_from_path_spec = self.calculate_loop_hash(None, new_loop['path_spec'])
            # assert new_loop['loop_hash'] == self.calculate_loop_hash(new_loop['path'])
            # Edit: you cannot re-calculate the loop hash from the path ifnodes and expect to get the same hash
            # if you have performed a reaction because then the ifnode's delegation scheme will have changed
            # and that is used when calling src/tgt.state_fingerprint().
            # If you want to rely on a cached hash, you have to use self.loops[loopid]['loop_hash']
            assert new_loop['loop_hash'] == self.calculate_loop_hash(new_loop['path_spec'])
            for ifnode in new_loop['path']:
                self.loopids_by_interface[ifnode].add(new_loop_id)
        else:
            # 1b: Remove existing loop
            del_loop_hash = loop_effects['del_loop_hash']
            del_loop_id = self.loopid_by_hash[del_loop_hash]
            if del_loop_id == loop_effects['del_loop_id']:
                print("del_loop_id == loop_effects['del_loop_id']: %s == %s" %
                      (del_loop_id, loop_effects['del_loop_id']))
            del_loop = self.loops[del_loop_id]
            for ifnode in del_loop['path']:
                self.loopids_by_interface[ifnode].remove(del_loop_id)
            assert not any(del_loop_id in loopids_set for loopids_set in self.loopids_by_interface.values())
            del self.loops[del_loop_id]
            del self.loopid_by_hash[del_loop_hash]

            # Update ifnodes hashes: - No. See above reason for why not.
            # self.rebuild_ifnode_by_hash_index()

        # 2. Update existing loops:
        for loop_hash, loop_info in loop_effects['changed_loops_by_hash'].items():
            if loop_hash is None:
                # None is used to indicate a new loop. but only when used as loop_id, not hash.
                print("\n\nThis shouldn't happen!!\n")
                pdb.set_trace()
            # Update loop ifnodes:
            loop_id = self.loopid_by_hash[loop_hash]
            old_loop_info = self.loop_by_loopid[loop_id]
            old_loop_nodes = set(old_loop_info['path'])
            new_loop_nodes = set(loop_info['path'])
            # Remove old ifnode entries from self.loopids_by_interface and add new:
            for ifnode in (old_loop_nodes - new_loop_nodes):
                self.loopids_by_interface[ifnode].remove(loop_id)
            for ifnode in (new_loop_nodes - old_loop_nodes):
                self.loopids_by_interface[ifnode].add(loop_id)

            assert all(ifnode in self.loopids_by_interface and loop_id in self.loopids_by_interface[ifnode]
                       for ifnode in new_loop_nodes)
            old_loop_info.update(loop_info)

        for loop in self.loops.values():
            top_delegate_path = [ifnode.top_delegate() for ifnode in loop['path']]
            assert top_delegate_path == loop['path']
            assert max(Counter(top_delegate_path).values()) == 1

        # Update loop hashes (state specs):
        # Note: We need the old loop state specs (fingerprints, hashes) when updating existing loops above
        # self.rebuild_loopid_by_hash_index(update_ifnodes=False)

        # TODO: Do we need to update loops when ifnode delegation changes?
        # E.g. if stacking, then four ifnodes are represented as one ifnode.
        # All four ifnodes will have the same hash, because ifnode.state_fingerprint uses all delegates.
        # So if a path includes the edge {src, tgt} and the edge is collapsed by e.g. stacking to a single ifnode,
        # then that is just {ifnode}, which should be OK.
        # What if we have a path that includes ifnode2 and ifnode2 is expanded to {ifnode2, ifnode3}?
        # The system-level InterfaceGraph only includes edges between the top delegates.
        # InterfaceGraph.delegate() will move the delegator's edges to the top delegatee.
        # So new paths DOES NOT include
        # But either way, the all affected paths should be updated to reflect the new ifnode top_delegate path, right?



    # def register_new_loop(self, loop_path_spec, loop_activity, loop_path=None, replacing_loop_spec=None):
    #     """ Register a new loop. """
    #     if loop_path_spec is None:
    #         loop_path_spec = [ifnode.state_fingerprint() for ifnode in loop_path]
    #     else:
    #         loop_path = self.recreate_loop_ifnodes_from_spec(loop_path_spec)
    #     # recreate_loop_ifnodes_from_spec will update self.ifnode_by_hash
    #     # You can then do: ifnodes = [self.ifnode_by_hash[spec] for spec in ifnodes_specs]
    #     loop_hash = self.calculate_loop_hash(loop_path)
    #     # To get a loop dict from a loop hash, first get the id with loopid_by_hash,
    #     # then use self.loops[loopid] (loops_by_id) to get loop info dict.
    #     assert loop_hash not in self.loopid_by_hash
    #     # Get new id for loop: (It would be a little easier if loops were object instances)
    #     loopid = next(loopid_sequential_number_generator)  # Is currently a uuid
    #     # What's in a loop (dict/object)?
    #     loop = {
    #         'loopid': loopid,
    #         'path': loop_path,
    #         'path_spec': loop_path_spec,
    #         'node_set': set(loop_path),  # Not sure this is worth it...
    #         'original_path': tuple(loop_path), # For later reference..
    #         'activity': loop_activity,  # a1 in the calculations, not "a"?
    #         'original_activity': loop_activity,  # For later reference..
    #         'original_loop_hash': loop_hash,  # The current loop hash is in self.loopid_by_hash.
    #         # For loop delegation/overrides see self.loop_delegations(_rev).
    #     }
    #     assert loopid not in self.loop_by_loopid
    #     self.loop_by_loopid[loopid] = loop
    #
    #
    # def update_changed_loop(self, loop0_hash, replacement_loop, refresh_cache=True):
    #     """ Update an existing loop to use a new, shorter path: """
    #     if refresh_cache:
    #         self.update_ifnodes_and_loops_specs()
    #     loop0id = self.loopid_by_hash[loop0_hash]
    #     loop0 = self.loop_by_loopid[loop0id]
    #     if 'loop_hash' in loop0: # This should usually just be a lookup entry in self.loopid_by_hash
    #         try:
    #             del self.loopid_by_hash[loop0['loop_hash']]
    #         except KeyError:
    #             print("Could not delete old entry in loopid_by_hash using old loop_hash %s" % loop0['loop_hash'])
    #     loop0.update(replacement_loop)  # Or you could just overwrite the old entry in loop_by_loopid.




    ## OLD FINGERPRINT METHODS: ##

    def strands_species_count(self):
        """ Count the number of strand species. Used as part of complex state finger-printing. """
        strand_species_counts = {}
        domain_species_counts = {}
        for strand in self.strands:
            if strand.name not in strand_species_counts:
                strand_species_counts[strand.name] = 1
            else:
                strand_species_counts[strand.name] += 1
            for domain in strand.domains:
                if domain.name not in domain_species_counts:
                    domain_species_counts[domain.name] = 1
                else:
                    domain_species_counts[domain.name] += 1
        # Remove entries from strand_species_counter with zero count:
        depleted = [sname for sname, count in self.strand_species_counter.items() if count < 1]
        for sname in depleted:
            del self.strand_species_counter[sname]

        assert strand_species_counts == self.strand_species_counter

        # Remove entries from domain_species_counter with zero count:
        depleted = [name for name, count in self.domain_species_counter.items() if count < 1]
        for name in depleted:
            del self.domain_species_counter[name]
        assert domain_species_counts == self.domain_species_counter
        #return strand_species_counts
        # Used for hashing, so return a hashable frozenset(((specie1, count), ...))
        return frozenset(strand_species_counts.items())

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


    def hybridization_edges(self, direction=None):
        """
        How to get a list of domain hybridization?
        * self.edges() will generate a list of all connections: backbone, hybridization and stacking.

        """
        # you can loop over all edges and filter by type:
        if direction is None:
            edges = [(d1, d2) for d1, d2, key, interaction in self.edges(keys=True, data='interaction')
                     if key == HYBRIDIZATION_INTERACTION or interaction == HYBRIDIZATION_INTERACTION]
        else:
            edges = [(d1, d2) for d1, d2, key, eattr in self.edges(keys=True, data=True)
                     if key == HYBRIDIZATION_INTERACTION and eattr['direction'] == direction]
        # or maybe it is easier to keep a dict with hybridization connections:
        assert set(frozenset(tup) for tup in edges) == self.hybridized_pairs
        ## TODO: Remove check and only use one or the other
        ## (I might not need self.hybridized_pairs with new domain icid labelling scheme..
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


    def stacking_edges(self, direction=DIRECTION_DOWNSTREAM):
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
        if direction is None:
            edges = [(d1, d2) for d1, d2, key, interaction in self.edges(keys=True, data='interaction')
                     if key == STACKING_INTERACTION or interaction == STACKING_INTERACTION]
        else:
            edges = [(d1, d2) for d1, d2, key, eattr in self.edges(keys=True, data=True)
                     if key == STACKING_INTERACTION and eattr['direction'] == direction]
        # Only if Complex is a DiGraph:
        assert set(edges) == self.stacked_pairs # and len(set(edges)) == len(self.stacked_pairs)
        # For now, I have to keep a dict with hybridization connections:
        return edges
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




# class SuperComplex(ConnectedMultiGraph):
#     """
#     Each node is a complex. For complexes that are making temporary stacking interactions.
#     Since stacking interactions are unlikely, the idea was to join them by forming a "SuperComplex",
#     instead of performing a full join_complex_at followed by a costly break_complex_at cycle for a transient interaction.
#     However, the added complexity currently isn't worth the expected gain (premature optimization).
#     Perhaps I will do this later, or perhaps I will find a way to make break_complex_at cheaper,
#     e.g. by better caching.
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
