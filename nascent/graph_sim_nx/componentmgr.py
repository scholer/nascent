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

# pylint: disable=C0103,W0142,W0212,R0902

"""

Module for managing the whole system.

I think there might be a lot of code involved in managing the system (graphs and structure),
and having all that in the simulator makes it a bit hard to read.

Splitting the "system state" code out to a separate module allows the simulator to be focused
on just the "stochastic simulation" part and be mostly agnostic on the system details.

Question: Who manages "reaction propensity" etc?
 - I'm thinking the system manager takes care of everything except the simulation steps / temperature control.


The "management system" is composed of three classes:
* Graph Manager - Takes care of system-level graphs
* Component Manager - Takes care of "components", i.e. domains, strands, complexes - as well as hybridize/dehybridize.
* Reaction Manager - Everything related to reactions, i.e. energies, c_j, propensity functions, etc.



## This manager will take care of: ##
System state:
 * Strands,
 * Domains,
 * Complexes
 * System graphs
 * Structural elements
Reactions: (this could be a separate object, but for now system state and reactions are integrated)
 * Possible hybridization reactions
 * Energies, energy model
 * Hybridization and dehybridization rates
 * Propensity functions
 *

"""

from __future__ import absolute_import, print_function, division
import os
#import random
from collections import defaultdict
#import math
#from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
from networkx.algorithms.components import connected_components, connected_component_subgraphs
# import numpy as np
import pdb

# from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
#, R # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K
# from .constants import R, N_AVOGADRO, AVOGADRO_VOLUME_NM3
from .constants import HYBRIDIZATION_INTERACTION, PHOSPHATEBACKBONE_INTERACTION, STACKING_INTERACTION
from .constants import HELIX_XOVER_DIST, HELIX_STACKING_DIST, HELIX_WIDTH
from .complex import Complex
from .graph_manager import GraphManager
from .debug import printd, pprintd
 # Enthalpies in units of R*K, entropies in units of R = 1.987 cal/mol/K
from nascent.energymodels.biopython import DNA_NN4_R
from .utils import (sequential_number_generator, sequential_uuid_gen)






class ComponentMgr(GraphManager):
    """
    System-manager class to manage the system state and state changes
    (but it does not control *which* state changes are induced).

    The GraphManager super-class provides everything related to structural analysis:
    loops and intra-complex activity calculation, etc.
    """

    def __init__(self, strands, params, domain_pairs=None):
        # super(ComponentMgr, self).__init__(strands=strands)
        GraphManager.__init__(self, strands=strands)
        self.params = params
        self.stacking_joins_complexes = params.get('stacking_joins_complexes', True)
        self.enable_helix_bending = params.get('enable_helix_bending', False)
        # super-complexes: two or more complexes held together by stacking interactions.
        # Edit: I've made a proper SuperComplex class. Edit2: I've dropped using supercomplexes for now.
        # self.supercomplexes = set()
        # complexes are assigned sequential complex-unique IDs upon instantiation, no need to keep track of order here.
        # Since we remove obsolete (degraded) complexes in arbitrary order, we don't want a regular list.
        self.complexes = set()
        self.removed_complexes = [] # But it might be interesting to keep track of deletion order.
        self.strands = strands
        self.strands_by_name = defaultdict(list)
        for strand in strands:
            self.strands_by_name[strand.name].append(strand)
        print("Strands in self.strands_by_name:")
        print("\n".join("- %10s: %s species" % (sname, len(strands))
                        for sname, strands in self.strands_by_name.items()))
        self.domains_list = [domain for strand in strands for domain in strand.domains]
        self.domains = set(self.domains_list)  # doesn't change

        # Stats - counts
        self.N_domains = len(self.domains)
        self.N_strands = len(self.strands)

        self.N_domains_hybridized = sum(1 for domain in self.domains_list if domain.partner is not None)
        # self.N_strands_hybridized = sum(1 for oligo in self.strands if oligo.is_hybridized())
        self.N_stacked_ends = 0 # Running count (for performance and system checks)
        self.N_stacked_ends = sum(1 for domain in self.domains_list for end in (domain.end3p, domain.end5p)
                                  if end.stack_partner is not None)
        self.N_stacked_ends = self.n_stacked_ends()

        self.domains_by_name = defaultdict(list)
        self.unhybridized_domains_by_name = defaultdict(set)
        self.hybridized_domains_by_name = defaultdict(set)

        for d in self.domains:
            self.domains_by_name[d.name].append(d)
            if d.partner is None:
                self.unhybridized_domains_by_name[d.name].add(d)
            else:
                self.hybridized_domains_by_name[d.name].add(d)
        print("Domains in self.domains_by_name:")
        print("\n".join("- %10s: %s species" % (dname, len(domains))
                        for dname, domains in self.domains_by_name.items()))
        if domain_pairs is None:
            # mapping: dom_a -> dom_A, dom_A -> dom_a
            # TODO: This could perhaps be a list, if you want to have different types of domains interacting,
            # E.g. dom_a could be perfect match for dom_A, while dom_ax has 1 mismatch:
            # domain_pairs[dom_A.name] = [dom_a.name, dom_ax.name]
            # Or it could be a set of sets: {{da, dA}, {dA, dax}} and then generate partners by:
            # partners_species = set(chain(pair for pair in domain_pairs if dA in pair)) - {dA}
            # However, might as well only do this once and save the list!
            # Also, if you change domain_pairs mapping, remember to adjust domain_dHdS cache as well.
            domain_pairs = {d.name: [d.name.lower()] if d.name == d.name.upper() else [d.name.upper()]
                            for d in self.domains_list}
            # remove pairs without partner:
            domain_pairs = {d1name: d2names for d1name, d2names in domain_pairs.items()
                            if d2names[0] in self.domains_by_name}
        # allow self-complementarity?
        assert not any(k == v for k, v in domain_pairs.items())
        self.domain_pairs = domain_pairs
        self.hyb_dehyb_file = None
        if params.get('save_hybdehyb_to_file'):
            fn = os.path.join(params.get('working_directory', '.'), "hyb_dehyb.py") \
                if params['save_hybdehyb_to_file'] is True else params['save_hybdehyb_to_file']
            self.hyb_dehyb_file = open(fn, 'w')


    def n_hybridizable_domains(self):
        return sum(1 for domain in self.domains_list if domain.name in self.domain_pairs)

    def n_hybridized_domains(self):
        """ Count the number of hybridized domains. """
        count = sum(1 for domain in self.domains_list if domain.partner is not None)
        if not count % 2 == 0:
            print("Weird - n_hybridized_domains counts to %s (should be an even number)" % count)
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.domains_list
                                                   if domain.partner is not None))
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.hybridized_domains_by_name))
        assert self.N_domains_hybridized == count
        return count

    def n_stacked_ends(self):
        count = sum(1 for domain in self.domains_list for end in (domain.end3p, domain.end5p)
                    if end.stack_partner is not None)
        n_stacking_edges = sum(1 for source, target, key in self.ends5p3p_graph.edges(keys=True)
                               if key == STACKING_INTERACTION)
        assert count == n_stacking_edges*2
        assert count == self.N_stacked_ends
        return count

    def n_partially_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for strand in self.strands if strand.is_hybridized())

    def n_fully_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for strand in self.strands if strand.is_fully_hybridized())



    def hybridize(self, domain1, domain2, join_complex=True):
        """
        Splitting out logic in preparation for Julia implementation.
        :is_forming: and :is_intra: is only used for checking.
        returns
            result dict with:
            changed_complexes, new_complexes, obsolete_complexes, free_strands, case
        Where changed_complexes includes new_complexes but not obsolete_complexes.
        Edit: changed_complexes DOES NOT include new_complexes.
        Cases:
            0: Intra-strand hybridization.
            1: Intra-complex hybridization; no change.
            2: Inter-complex hybridization; merge two complexes into one.
            3: Complex+strand hybridization; add strand to complex.
            4: Inter-strand hybridization; create new complex
        """
        # printd("%s.hybridize(%s, %s) invoked..." % (type(self).__name__, domain1, domain2))
        if self.hyb_dehyb_file:
            print("domain1, domain2 = (domains_by_duid[%s], domains_by_duid[%s])" % (domain1.duid, domain2.duid),
                  file=self.hyb_dehyb_file)
            print("assert domain1.domain_strand_specie == %s and domain2.domain_strand_specie == %s" %
                  (domain1.domain_strand_specie, domain2.domain_strand_specie), file=self.hyb_dehyb_file)
            print("sysmgr.hybridize(domain1, domain2)", file=self.hyb_dehyb_file)
        assert domain1 != domain2
        assert domain1.partner is None
        assert domain2.partner is None

        #dset = frozenset((domain1, domain2))
        #sset = frozenset((domain1.strand, domain2.strand))
        domain1.partner = domain2
        domain2.partner = domain1
        strand1 = domain1.strand
        strand2 = domain2.strand

        # Update system-level graphs:
        edge_kwargs = {"interaction": HYBRIDIZATION_INTERACTION,
                       #"len": 6, # With of ds helix ~ 2 nm ~ 6 bp
                       #"weight": 2,
                       #"key": HYBRIDIZATION_INTERACTION
                       'dist_ee_sq': HELIX_WIDTH**2,
                       'dist_ee_nm': HELIX_WIDTH,
                       'stiffness': 1,
                      }
        #key = (domain1.universal_name, domain2.universal_name, HYBRIDIZATION_INTERACTION)
        s_edge_key = (frozenset((domain1.universal_name, domain2.universal_name)), HYBRIDIZATION_INTERACTION)
        # printd("Adding strand_graph edge (%s, %s, key=%s)" % (strand1, strand2, s_edge_key))
        edge_key = HYBRIDIZATION_INTERACTION
        self.strand_graph.add_edge(strand1, strand2, key=s_edge_key, interaction=edge_key)
        self.domain_graph.add_edge(domain1, domain2, key=edge_key, attr_dict=edge_kwargs)
        self.ends5p3p_graph.add_edge(domain1.end5p, domain2.end3p, key=edge_key, attr_dict=edge_kwargs)
        self.ends5p3p_graph.add_edge(domain2.end5p, domain1.end3p, key=edge_key, attr_dict=edge_kwargs)
        # Hybridization => merge two DomainEnds into a single node by delegating representation of one to the other.
        end1_delegatee = self.interface_graph.merge(domain1.end5p.ifnode, domain2.end3p.ifnode)
        end2_delegatee = self.interface_graph.merge(domain1.end3p.ifnode, domain2.end5p.ifnode)
        # You could just use interface_graph[end1_delegate][end2_delegate] instead of looping...
        # Remember, interface_graph.adj reflects the current *representation*. After merging/delegating ifnodes,
        # the edges from delegator is no longer available. Either update edge before merge or only update delegatee.
        self.interface_graph[end1_delegatee][end2_delegatee]['dist_ee_nm'] = domain1.ds_dist_ee_nm
        self.interface_graph[end1_delegatee][end2_delegatee]['dist_ee_sq'] = domain1.ds_dist_ee_sq
        self.interface_graph[end1_delegatee][end2_delegatee]['stiffness'] = 1

        # Update edge dist_ee_nm, dist_ee_sq and stiffness attrs for edge between d.end5p and d.end3p (for both d):
        for d in (domain1, domain2):
            if self.ends5p3p_graph.is_multigraph():
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['dist_ee_nm'] = d.ds_dist_ee_nm
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['dist_ee_sq'] = d.ds_dist_ee_sq
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['stiffness'] = 1
            else:
                self.ends5p3p_graph[d.end5p][d.end3p]['dist_ee_nm'] = d.ds_dist_ee_nm  # domain ds end-to-end dist
                self.ends5p3p_graph[d.end5p][d.end3p]['dist_ee_sq'] = d.ds_dist_ee_sq  # cached squared result, nm2
                self.ends5p3p_graph[d.end5p][d.end3p]['stiffness'] = 1  # ds helix has stiffness 1
            # Edit: only updating edges from delegatee not delegator edge.
            # self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['dist_ee_nm'] = d.ds_dist_ee_nm
            # self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['dist_ee_sq'] = d.ds_dist_ee_sq
            # self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['stiffness'] = 1

        if join_complex:
            result = self.join_complex_at(domain_pair=(domain1, domain2), edge_kwargs=edge_kwargs)
        else:
            c1 = strand1.complex
            result = {'changed_complexes': [c1] if c1 is not None else None,
                      'new_complexes': None, 'obsolete_complexes': None,
                      'free_strands': None, 'case': -1}
            if c1 is not None:
                # c1.history.append("hybridize: Adding edge (%r, %r, %r)." % (domain1, domain2, edge_key))
                # I'm experimenting having Complex being a DiGraph...
                c1.add_hybridization_edge((domain1, domain2))
                # c1.add_edge(domain1, domain2, key=edge_key, attr_dict=edge_kwargs)
                # c1.add_edge(domain2, domain1, key=edge_key, attr_dict=edge_kwargs)
                # c1._hybridization_fingerprint = None
                # c1._state_fingerprint = None


        assert strand1.complex == strand2.complex != None

        if self.hyb_dehyb_file:
            print("print('- hybridize complete.')", file=self.hyb_dehyb_file)

        self.N_domains_hybridized += 2

        return result


    def dehybridize(self, domain1, domain2, break_complex=True):
        """
        Dehybridize domain2 from domain1.
        Returns a result dict with keys:
            changed_complexes, new_complexes, obsolete_complexes, free_strands, case
        Three cases:
        ## Case 0: No complexes; Intra-strand (e.g. hairpin) de-hybridization.
        ## Case 1/(a) Two smaller complexes - must create a new complex for detached domain:
        ## Case 2/(b) One complex and one unhybridized strand - no need to do much further
        ## Case 3/(c) Two unhybridized strands
        """
        # printd("%s.dehybridize(%s, %s) invoked..." % (type(self).__name__, domain1, domain2))
        if self.hyb_dehyb_file:
            print("domain1, domain2 = (domains_by_duid[%s], domains_by_duid[%s])" % (domain1.duid, domain2.duid),
                  file=self.hyb_dehyb_file)
            print("assert domain1.domain_strand_specie = %s and domain1.domain_strand_specie = %s" %
                  (domain1.domain_strand_specie, domain2.domain_strand_specie), file=self.hyb_dehyb_file)
            print("sysmgr.dehybridize(domain1, domain2)", file=self.hyb_dehyb_file)
        assert domain1 != domain2
        assert domain1.partner == domain2 != None
        assert domain2.partner == domain1 != None
        assert domain1.end3p.hyb_partner is domain2.end5p and domain2.end5p.hyb_partner is domain1.end3p
        assert domain2.end3p.hyb_partner is domain1.end5p and domain1.end5p.hyb_partner is domain2.end3p
        ## Edit/new: Domains CANNOT dehybridize if they are stacked
        assert all(end.stack_partner is None for end in
                   (domain1.end3p, domain1.end5p, domain2.end3p, domain2.end5p))

        # dset = frozenset((domain1, domain2))
        #sset = frozenset(domain1.strand, domain2.strand)
        domain1.partner = None
        domain2.partner = None

        strand1 = domain1.strand
        strand2 = domain2.strand
        c = strand1.complex
        c.domains_distances = {}    # Reset distances.
        assert c == strand2.complex

        # Update system-level graphs:
        s_edge_key = (frozenset((domain1.universal_name, domain2.universal_name)), HYBRIDIZATION_INTERACTION)
        # printd("%s: Removing strand_graph edge (%s, %s, key=%s)" % (type(self).__name__, strand1, strand2, s_edge_key))
        self.strand_graph.remove_edge(strand1, strand2, key=s_edge_key)
        edge_key = HYBRIDIZATION_INTERACTION
        self.domain_graph.remove_edge(domain1, domain2, key=edge_key)
        self.ends5p3p_graph.remove_edge(domain1.end5p, domain2.end3p, key=edge_key)
        self.ends5p3p_graph.remove_edge(domain2.end5p, domain1.end3p, key=edge_key)
        # Dehybridization => Undo the effect of current InterfaceNode delegation:
        self.interface_graph.split(domain1.end5p.ifnode, domain2.end3p.ifnode)
        self.interface_graph.split(domain1.end3p.ifnode, domain2.end5p.ifnode)

        # Update edge dist_ee_nm, dist_ee_sq and stiffness attrs for edge between d.end5p and d.end3p (for both d):
        for d in (domain1, domain2):
            if self.ends5p3p_graph.is_multigraph():
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['dist_ee_nm'] = d.ss_dist_ee_nm
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['dist_ee_sq'] = d.ss_dist_ee_sq
                self.ends5p3p_graph[d.end5p][d.end3p][PHOSPHATEBACKBONE_INTERACTION]['stiffness'] = 0
                # ss backbone has zero stiffness
            else:
                self.ends5p3p_graph[d.end5p][d.end3p]['dist_ee_nm'] = d.ss_dist_ee_nm
                self.ends5p3p_graph[d.end5p][d.end3p]['dist_ee_sq'] = d.ss_dist_ee_sq
                self.ends5p3p_graph[d.end5p][d.end3p]['stiffness'] = 0  # ss backbone has zero stiffness
            self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['dist_ee_nm'] = d.ss_dist_ee_nm
            self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['dist_ee_sq'] = d.ss_dist_ee_sq
            self.interface_graph[d.end5p.ifnode][d.end3p.ifnode]['stiffness'] = 0

        ## If domain is stacked, break the stacking interaction before breaking complex:
        ## Edit/new: Domains CANNOT dehybridize if they are stacked
        assert all(endp.stack_partner is None for endp in
                   (domain1.end3p, domain2.end3p, domain2.end5p, domain1.end5p))
        # unstacking_results = {}
        # for d, h2d in ((domain1, domain2), (domain2, domain1)):
        #     h1end3p = d.end3p
        #     h2end5p = h2d.end5p  # Is "unpartnered" above, so cannot use end.partner
        #     #             h1end3p         h1end5p
        #     # Helix 1   ----------3' : 5'----------
        #     # Helix 2   ----------5' : 3'----------
        #     #             h2end5p         h2end3p
        #     if h1end3p.stack_partner is not None:
        #         ## Variable unpacking:
        #         h1end5p = h1end3p.stack_partner
        #         h2end3p = h2end5p.stack_partner
        #         assert h2end3p == h1end5p.hyb_partner  # Other stacking end is still intact.
        #         # Edit: When I do domain.partner = None above, it also sets domain.end5p/3p = None! (property)
        #         # assert h2end5p is h1end3p.hyb_partner
        #         # With break_complex=False, unstack() will update system-graphs and complex domain graph,
        #         # but not check if complex is broken apart.
        #         # But wait - we SHOULD check if complex breaks apart. That can EASILY happen.
        #         ## TODO: FIX THIS!!
        #         ## Hypothesis: We can ONLY break a hybridization that is NOT stacked?
        #         ## - Must give correct result for the "1 vs 2 vs N split" example case.
        #         ## Can domains dehybridize if stacked? If not, why?
        #         ## - We currently only consider stacking of hybridized domains (duplexes), NOT single-stranded domains.
        #         ## - If we cannot stack domains that are not hybridized, then we cannot dehybridize domains that are
        #         ##   stacked (the reverse rx). All reactions must be reversible,
        #         ##   otherwise we would get a net cyclic transport in our reaction network.
        #         unstacking_results[frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))] = \
        #             self.unstack(h1end3p, h2end5p, h2end3p, h1end5p, break_complex=False)
        #         # printd("Unstack result for domain %r:" % (d,))
        #         # pprintd(unstacking_results[-1][1])
        #     # else:
        #         # printd("Domain %r 3p end is not stacked...")

        if break_complex:
            result = self.break_complex_at(domain_pair=(domain1, domain2))
            # printd("Result for breaking complex %s between %r and %r:" % (c, domain1, domain2))
            # pprintd(result)
        else:
            # Manually remove complex edges (but don't break it into two):
            c1 = strand1.complex
            result = {'changed_complexes': [c1],
                      'new_complexes': None, 'obsolete_complexes': None,
                      'free_strands': None, 'case': -1}
            if c1 is not None:
                c1.history.append("dehybridize: Removing edge (%r, %r, %r)." %
                                  (domain1, domain2, HYBRIDIZATION_INTERACTION))
                c1.remove_hybridization_edge((domain1, domain2))
                # c1.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
                # c1.hybridized_pairs.remove(frozenset((domain1, domain2)))
                # c1._hybridization_fingerprint = None
                # c1._state_fingerprint = None


        # result['unstacking_results'] = unstacking_results
        # We return unstacking_results as a dict of: {<ends four-tuple>: result), ....}
        # for use in ReactionMgr.hybridize_and_process which invokes ReactionMgr.update_possible_stacking_reactions(...)

        if self.hyb_dehyb_file:
            print("print('- dehybridize complete.')", file=self.hyb_dehyb_file)

        self.N_domains_hybridized -= 2

        return result


    def stack(self, h1end3p, h2end5p, h2end3p, h1end5p, join_complex=True):
        """
        Form a stacking interaction.
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p

        Note: Domains on the same helix may or may not be also connected by their phosphate backbone.
        E.g. you could have a hinge, where one helix is backbone-connected and the other one not.
        This is probably the most common case, e.g. in N-way junctions.
        """
        if self.hyb_dehyb_file:
            print("h1end3p, h2end5p, h2end3p, h1end5p =",
                  "getattr(domains_by_duid[%s], 'end%s')," % (h1end3p.domain.duid, h1end3p.end),
                  "getattr(domains_by_duid[%s], 'end%s')," % (h1end5p.domain.duid, h1end5p.end),
                  "getattr(domains_by_duid[%s], 'end%s')," % (h2end3p.domain.duid, h2end3p.end),
                  "getattr(domains_by_duid[%s], 'end%s') " % (h2end5p.domain.duid, h2end5p.end),
                  file=self.hyb_dehyb_file)
            print("assert h1end3p.domain.domain_strand_specie == ", h1end3p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h1end5p.domain.domain_strand_specie == ", h1end5p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h2end3p.domain.domain_strand_specie == ", h2end3p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h2end5p.domain.domain_strand_specie == ", h2end5p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("sysmgr.stack(h1end3p, h2end5p, h2end3p, h1end5p, join_complex=%s)" % join_complex,
                  file=self.hyb_dehyb_file)
        stacking_tuple = ((h1end3p, h2end5p), (h2end3p, h1end5p))
        stacking_pair = frozenset(stacking_tuple)
        ## Variable unpacking:
        h1strand1 = h1end3p.domain.strand
        h1strand2 = h1end5p.domain.strand
        h2strand1 = h2end3p.domain.strand
        h2strand2 = h2end5p.domain.strand
        c1 = h1strand1.complex
        c2 = h1strand2.complex
        if c1:
            # Edit: We can have stacking between complexes
            # assert c1 == h2strand1.complex
            # assert c2 == h2strand2.complex
            pass

        ## Assertions:
        assert h1end3p.stack_partner is None
        assert h1end5p.stack_partner is None
        assert h2end3p.stack_partner is None
        assert h2end5p.stack_partner is None

        ## Update domain end attributes:
        h1end3p.stack_partner, h1end5p.stack_partner = h1end5p, h1end3p
        h2end3p.stack_partner, h2end5p.stack_partner = h2end5p, h2end3p

        stack_string = "%s%s/%s%s" % (h1end3p.base, h1end5p.base, h2end5p.base, h2end3p.base)
        # Stacking string is:  "CA/GT" - i.e. {h1end3p}{h1end5p}/{h2end5p}{h2end3p} as above
        #          h1end3p     h1end5p
        # h1 5'---ATGCATG|C : A|TGCATGC---3'
        # h2 3'---TACGTAC|G : T|ACGTACG---5'
        #          h2end5p ^^^ h2end3p
        if stack_string not in DNA_NN4_R:
            stack_string = stack_string[::-1]
            assert stack_string in DNA_NN4_R

        h1end3p.stack_string = h1end5p.stack_string = h2end3p.stack_string = h2end5p.stack_string = stack_string

        ## Update system-level graphs:
        edge_kwargs = {"interaction": STACKING_INTERACTION,
                       ## You can have separate attrs, length_nm and length_bp, but there should be just one "len":
                       #"len": 0.5, # With of ds helix ~ 2 nm ~ 6 bp
                       #"weight": 2,
                       #"key": HYBRIDIZATION_INTERACTION
                       "dist_ee_nm": HELIX_STACKING_DIST,
                       "dist_ee_sq": HELIX_STACKING_DIST**2,
                       "stiffness": 1,
                      }
        # h1_edge_key = (frozenset((h1end3p.instance_name, h1end5p.instance_name)), STACKING_INTERACTION)
        # h2_edge_key = (frozenset((h2end3p.instance_name, h2end5p.instance_name)), STACKING_INTERACTION)
        # Edit: domain-level stacking interactions have directionality, use tuple:
        h1_edge_key = (STACKING_INTERACTION, (h1end3p.instance_name, h1end5p.instance_name))
        h2_edge_key = (STACKING_INTERACTION, (h2end3p.instance_name, h2end5p.instance_name))

        # printd("Adding strand_graph edge (%s, %s, key=%s)" % (h1strand1, h1strand2, h1_edge_key))
        self.strand_graph.add_edge(h1strand1, h1strand2, key=h1_edge_key, interaction=STACKING_INTERACTION)
        # printd("Adding strand_graph edge (%s, %s, key=%s)" % (h2strand1, h2strand2, h2_edge_key))
        self.strand_graph.add_edge(h2strand1, h2strand2, key=h2_edge_key, interaction=STACKING_INTERACTION)

        self.domain_graph.add_edge(h1end3p.domain, h1end5p.domain, key=h1_edge_key, attr_dict=edge_kwargs)
        self.domain_graph.add_edge(h2end3p.domain, h2end5p.domain, key=h2_edge_key, attr_dict=edge_kwargs)
        self.ends5p3p_graph.add_edge(h1end3p, h1end5p, key=STACKING_INTERACTION, attr_dict=edge_kwargs)
        self.ends5p3p_graph.add_edge(h2end3p, h2end5p, key=STACKING_INTERACTION, attr_dict=edge_kwargs)
        # Stacking => merge two DomainEnds into a single node by delegating representation of one to the other.
        delegatee = self.interface_graph.merge(h1end3p.ifnode.top_delegate(), h2end3p.ifnode.top_delegate())

        if join_complex and self.stacking_joins_complexes:
            result = self.join_complex_at(stacking_pair=stacking_tuple, edge_kwargs=edge_kwargs)
        else:
            assert False  # Making sure this doesn't happen during testing...
            result = {'changed_complexes': [c1],
                      'new_complexes': None, 'obsolete_complexes': None,
                      'free_strands': None, 'case': -1}
            if c1 is not None:
                c1.history.append(("stack: Manually adding stacking edges (%r, %r) (%r, %r) "
                                   "but not performing merge (if relevant).") %
                                  (h1end3p.domain, h1end5p.domain, h2end3p.domain, h2end5p.domain))
                c1.add_stacking_edge(stacking_tuple)
                # c1.add_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION, attr_dict=edge_kwargs)
                # c1.add_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION, attr_dict=edge_kwargs)
                # c1._stacking_fingerprint = None  # Reset hybridization fingerprint
                # c1._state_fingerprint = None

        self.N_stacked_ends += 4

        if self.hyb_dehyb_file:
            print("print('- stacking complete.')", file=self.hyb_dehyb_file)

        return result


    def unstack(self, h1end3p, h2end5p, h2end3p, h1end5p, break_complex=True):
        """
        Break a stacking interaction.
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p
        """
        if self.hyb_dehyb_file:
            print("h1end3p, h2end5p, h2end3p, h1end5p =",
                  "getattr(domains_by_duid[%s], 'end%s')," % (h1end3p.domain.duid, h1end3p.end),
                  "getattr(domains_by_duid[%s], 'end%s')," % (h1end5p.domain.duid, h1end5p.end),
                  "getattr(domains_by_duid[%s], 'end%s')," % (h2end3p.domain.duid, h2end3p.end),
                  "getattr(domains_by_duid[%s], 'end%s') " % (h2end5p.domain.duid, h2end5p.end),
                  file=self.hyb_dehyb_file)
            print("assert h1end3p.domain.domain_strand_specie ==", h1end3p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h1end5p.domain.domain_strand_specie ==", h1end5p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h2end3p.domain.domain_strand_specie ==", h2end3p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("assert h2end5p.domain.domain_strand_specie ==", h2end5p.domain.domain_strand_specie,
                  file=self.hyb_dehyb_file)
            print("sysmgr.unstack(h1end3p, h2end5p, h2end3p, h1end5p, break_complex=%s)" % break_complex,
                  file=self.hyb_dehyb_file)

        stacking_tuple = ((h1end3p, h2end5p), (h2end3p, h1end5p))
        stacking_pair = frozenset(stacking_tuple)

        ## Variable unpacking:
        h1strand1 = h1end3p.domain.strand
        h1strand2 = h1end5p.domain.strand
        h2strand1 = h2end3p.domain.strand
        h2strand2 = h2end5p.domain.strand
        c1 = h1strand1.complex
        c2 = h1strand2.complex
        if c1:
            assert c1 == h2strand1.complex
            assert c2 == h2strand2.complex

        ## Assertions:
        assert h1end3p.stack_partner is h1end5p
        assert h1end5p.stack_partner is h1end3p
        assert h2end3p.stack_partner is h2end5p
        assert h2end5p.stack_partner is h2end3p

        ## Update domain end attributes:
        h1end3p.stack_partner, h1end5p.stack_partner = None, None
        h2end3p.stack_partner, h2end5p.stack_partner = None, None
        h1end3p.stack_string = h1end5p.stack_string = h2end3p.stack_string = h2end5p.stack_string = None

        ## Update system-level graphs:
        # h1_edge_key = (frozenset((h1end3p.instance_name, h1end5p.instance_name)), STACKING_INTERACTION)
        # h2_edge_key = (frozenset((h2end3p.instance_name, h2end5p.instance_name)), STACKING_INTERACTION)
        # Edit: stacking keys have directionality. Use tuple.
        h1_edge_key = (STACKING_INTERACTION, (h1end3p.instance_name, h1end5p.instance_name))
        h2_edge_key = (STACKING_INTERACTION, (h2end3p.instance_name, h2end5p.instance_name))
        # printd("Removing strand_graph edge (%s, %s, key=%s)" % (h1strand1, h1strand2, h1_edge_key))
        self.strand_graph.remove_edge(h1strand1, h1strand2, key=h1_edge_key)
        # printd("Removing strand_graph edge (%s, %s, key=%s)" % (h2strand1, h2strand2, h2_edge_key))
        self.strand_graph.remove_edge(h2strand1, h2strand2, key=h2_edge_key)

        self.domain_graph.remove_edge(h1end3p.domain, h1end5p.domain, key=h1_edge_key)
        self.domain_graph.remove_edge(h2end3p.domain, h2end5p.domain, key=h2_edge_key)
        self.ends5p3p_graph.remove_edge(h1end3p, h1end5p, key=STACKING_INTERACTION)
        self.ends5p3p_graph.remove_edge(h2end3p, h2end5p, key=STACKING_INTERACTION)
        #             h1end3p         h1end5p
        # Helix 1   ----------3' : 5'----------
        # Helix 2   ----------5' : 3'----------
        #             h2end5p         h2end3p
        # Stacking => merge two DomainEnds into a single node by delegating representation of one to the other.
        # Just because h1end3p.ifnode.delegatee is not None, does not mean that h2end5p *is* the delegatee.
        # When the duplexes are stacked, both may well not be None.
        end1_delegatee = h1end3p.ifnode if h2end5p.ifnode.delegatee is h1end3p.ifnode else h2end5p.ifnode
        end2_delegatee = h2end3p.ifnode if h1end5p.ifnode.delegatee is h2end3p.ifnode else h1end5p.ifnode
        # Stacking should be the top delegation, there shouldn't be further layers:
        assert (end1_delegatee.delegatee is None and end2_delegatee.delegatee is end1_delegatee) \
            or (end2_delegatee.delegatee is None and end1_delegatee.delegatee is end2_delegatee)
        # Cannot use top_delegate, that would just get the same for both, have to find the delegate for each duplex.
        delegatee = self.interface_graph.split(end1_delegatee, end2_delegatee)

        if break_complex and self.stacking_joins_complexes:
            result = self.break_complex_at(stacking_pair=stacking_tuple)
        else:
            # Manually update complex domain graph, but do not process/break the complex into two (if applicable).
            # assert False # checking that this does not happen during testing..
            ## Changed: Even if we don't explicitly break the complex (e.g. for stacking/unstacking),
            ##    we should still return a results dict including a 'changed_complexes' entry,
            ##    and maybe use case = -1 to indicate that no breaking was attempted.
            result = {'changed_complexes': [c1],
                      'new_complexes': None, 'obsolete_complexes': None,
                      'free_strands': None, 'case': -1}
            if c1 is not None:
                # c1.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
                # c1._hybridization_fingerprint = None
                # old_fingerprints1 = c1.get_all_fingerprints()
                # c1.reset_state_fingerprint()
                # old_fingerprints2 = c1.get_all_fingerprints()
                c1.history.append(("unstack: Manually removing stacking edges (%r, %r) (%r, %r) "
                                   "but not breaking complex (if relevant)...") %
                                  (h1end3p.domain, h1end5p.domain, h2end3p.domain, h2end5p.domain))
                c1.remove_stacking_edge(stacking_tuple)
                # c1.remove_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION)
                # c1.remove_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION)
                # c1.reset_state_fingerprint()
                # c1._stacking_fingerprint = None  # Reset hybridization fingerprint
                # c1._state_fingerprint = None
                # for end in h1end3p, h2end5p, h2end3p, h1end5p:
                #     end.domain.state_change_reset(reset_complex=False)
                # new_fingerprints = c1.get_all_fingerprints()
                # assert new_fingerprints != old_fingerprints2

        self.N_stacked_ends -= 4

        if self.hyb_dehyb_file:
            print("print('- un-stacking complete.')", file=self.hyb_dehyb_file)

        return result


    def join_complex_at(self, domain_pair=None, stacking_pair=None, edge_kwargs=None):
        """
        Join one or more complexes, by forming one or more edges.
        Return result dict with:
        result = {'changed_complexes': None,
                  'new_complexes': None,
                  'obsolete_complexes': None,
                  'free_strands': None,
                  'case': None}
        'case' is any of:
            Case -1: No complex merge/break was attempted,
                     used e.g. for stacking where we don't always want to look for inter-complex reactions.
            Case 0 : intRA-STRAND hybridization.
            Case 1 : IntRA-complex hybridization
            Case 2 : IntER-complex hybridization between two complexes. Merge the two complexs:
            Case 3a: domain2/strand2 is not in a complex; use c1
            Case 3b: domain1/strand1 is not in a complex; use c2
            Case 4 : Neither strands are in existing complex; create new complex
        Thus:
            Case 0-1: Intra/uni-molecular reactions; no change in volume energy.
        """
        if domain_pair is None and stacking_pair is None:
            raise ValueError("Must provide either domain_pair or stacking_pair.")
        if domain_pair is not None and stacking_pair is not None:
            raise NotImplementedError("Providing both domain_pair and stacking_pair is not yet implemented!")

        if domain_pair:
            domain1, domain2 = domain_pair
            # NOTE: Complex graph edge attributes doesn't matter very much!
            if edge_kwargs is None:
                edge_kwargs = {"interaction": HYBRIDIZATION_INTERACTION,
                               #"len": 6, # With of ds helix ~ 2 nm ~ 6 bp
                               #"weight": 2,
                               'dist_ee_nm': HELIX_WIDTH,
                               'dist_ee_sq': HELIX_WIDTH**2,
                               'stiffness': 1
                              }
            domain_edge_kwargs = edge_kwargs
            #h1end3p = h1end5p = (h2end3p, h2end5p) = ()
        else:
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
            if edge_kwargs is None:
                edge_kwargs = {"interaction": STACKING_INTERACTION,
                               #"len": 6, # With of ds helix ~ 2 nm ~ 6 bp
                               #"weight": 2,
                               'dist_ee_nm': HELIX_STACKING_DIST,
                               'dist_ee_sq': HELIX_STACKING_DIST**2,
                               'stiffness': 1
                              }
            stacking_edge_kwargs = edge_kwargs
            domain1, domain2 = h1end3p.domain, h2end3p.domain
        ## Variable unpacking:
        strand1 = domain1.strand
        strand2 = domain2.strand
        c1 = strand1.complex
        c2 = strand2.complex

        # changed_complexes, new_complexes, obsolete_complexes, free_strands = None, None, None, []
        result = {'changed_complexes': None,
                  'new_complexes': None,
                  'obsolete_complexes': None,
                  'free_strands': None,
                  'case': None}

        if strand1 == strand2:
            # If forming an intra-strand connection, we don't need to check if we make or merge any Complexes
            # print("\nhybridize case 0: intra-strand hybridization/stacking.")
            # if domain_pair:
            #     print(domain_pair)
            # if stacking_pair:
            #     print(stacking_pair)
            result['case'] = 0
            assert c1 == c2
            if c1 is None:
                # No complex to update, just return:
                return result
            else:
                c_major = c1
                result['changed_complexes'] = [c1]
        ## Else: Analyse what happened to the complex
        elif c1 and c2:
            ## Both domains are in a complex.
            if c1 == c2:
                ## Case (1): Intra-complex hybridization
                # printd("hybridize case 1: intra-complex hybridization.")
                result['case'] = 1
                assert strand1 in c2.strands and strand2 in c1.strands
                # c1.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION) # done below
                c_major = c1
                result['changed_complexes'] = [c_major]
                c_major.history.append("join_complex_at: Intra-complex hybridization.")
            else:
                ## Case (2): Inter-complex hybridization between two complexes. Merge the two complexs:
                # printd("hybridize case 2: inter-complex hybridization.")
                result['case'] = 2
                c_major, c_minor = (c1, c2) if (len(c1.nodes()) >= len(c2.nodes())) else (c2, c1)
                c_major.history.append("join_complex_at: Inter-complex hybridization, merging with c_minor = %r" %
                                       (c_minor,))
                # Import nodes and edges to major complex:
                # add_strands updates: strand.complex, c_major.strands, c_major.strands_by_name
                c_major.add_strands(c_minor.strands, update_graph=False)
                # We use the c_minor graph - rather than strands ^ - to update c_major graph:
                c_major.add_nodes_from(c_minor.nodes(data=True))
                c_major.add_edges_from(c_minor.edges(keys=True, data=True))
                c_major.N_strand_changes += c_minor.N_strand_changes # Currently not keeping track of strand changes.
                c_major.history.append(c_minor.history)
                c_major.hybridized_pairs |= c_minor.hybridized_pairs
                c_major.stacked_pairs |= c_minor.stacked_pairs
                # Complex.add_strands takes care of updating counters:
                #c_major.domain_species_counter += c_minor.domain_species_counter
                #c_major.strand_species_counter += c_minor.strand_species_counter
                ## Delete the minor complex:
                c_minor.strands.clear()     # Clear strands
                c_minor.strands_by_name.clear()
                c_minor.domain_species_counter.clear()
                c_minor.strand_species_counter.clear()
                c_minor.hybridized_pairs.clear()
                c_minor.stacked_pairs.clear()
                c_minor.node.clear()        # Clear nodes
                c_minor.adj.clear()         # Clear edges (graph.adj attribute contains edge data for nx graphs)
                c_minor.history.append("This complex is now merged with %r" % (c_major,))
                result['obsolete_complexes'] = [c_minor]
                result['changed_complexes'] = [c_major]
        elif c1:
            ## Case 3a: domain2/strand2 is not in a complex; use c1
            # printd("hybridize case 3a: strand hybridizing to complex.")
            result['case'] = 3
            c_major = c1
            c_major.history.append("join_complex_at: Strand+complex hybridization.")
            c_major.add_strand(strand2, update_graph=True)
            c_major._strands_fingerprint = None
            result['changed_complexes'] = [c_major]
        elif c2:
            ## Case 3b: domain1/strand1 is not in a complex; use c2
            # printd("hybridize case 3b: strand hybridizing to complex.")
            result['case'] = 3
            c_major = c2
            c_major.history.append("join_complex_at: Strand+complex hybridization.")
            c_major.add_strand(strand1, update_graph=True)
            c_major._strands_fingerprint = None
            result['changed_complexes'] = [c_major]
        else:
            ## Case 4: Neither strands are in existing complex; create new complex
            result['case'] = 4
            # printd("hybridize case 4: inter-strand hybridization (forming a new complex).")
            new_complex = Complex(strands=[strand1, strand2])
            new_complex.history.append("join_complex_at: Complex created from strands %r, %r." % (strand1, strand2))
            # new_complex.strands |= {strand1, strand2}
            # strand1.complex = strand2.complex = new_complex
            # new_complex.add_strand(strand1)
            # new_complex.add_strand(strand2)
            c_major = new_complex
            result['new_complexes'] = [new_complex]

        # Create the hybridization connection in the major complex graph:
        if domain_pair:
            c_major.history.append("join_complex_at: Adding edge (%r, %r, %r)." %
                                   (domain1, domain2, HYBRIDIZATION_INTERACTION))
            c_major.add_hybridization_edge((domain1, domain2))
        if stacking_pair:
            c_major.history.append("join_complex_at: Adding edges (%r, %r, %r) and (%r, %r, %r)." %
                                   (h1end3p.domain, h1end5p.domain, STACKING_INTERACTION,
                                    h2end3p.domain, h2end5p.domain, STACKING_INTERACTION))
            c_major.add_stacking_edge(stacking_pair)

        # Reset complex: (maybe not needed here)
        c_major.domain_distances = {} # Reset distances
        c_major.reset_state_fingerprint()

        return result


    def break_complex_at(self, domain_pair=None, stacking_pair=None, interaction=None):
        """
        Break complex at one or more edges and determine how it splits up.
        Note: This only deals with the complex (although it uses system graphs for analysis).
        Make sure to break edges in system graphs before calling this method to determine and
        process complex breakage.

        Case enumeration, conforming to cases returned by join_complex_at:
            Case -1:  Join/break_complex was not attempted (e.g. if check is disabled).
            CASE  0:  NO COMPLEX PRESENT, intra-strand reaction.
            Case  1:  The two strands are still connected: No need to do anything further
            Case  2:  Two smaller complexes - must create a new complex for detached domain:
            Case  3:  One complex and one unhybridized strand - no need to do much further
            Case  4:  Two unhybridized strands
        Cases 0-1 has no change in volume entropy; cases 2-4 has 1 unit of volume entropy change.

        Old case enumeration:
            ## Case 0: No complexes; Intra-strand (e.g. hairpin) de-hybridization.
            ## Case 1/(a) Two smaller complexes - must create a new complex for detached domain:
            ## Case 2/(b) One complex and one unhybridized strand - no need to do much further
            ## Case 3/(c) Two unhybridized strands
        """
        if domain_pair is None and stacking_pair is None:
            raise ValueError("Must provide either domain_pair or stacking_pair.")
        if domain_pair is not None and stacking_pair is not None:
            raise NotImplementedError("Providing both domain_pair and stacking_pair is not yet implemented!")
        if stacking_pair:
            (h1end3p, h2end5p), (h2end3p, h1end5p) = stacking_pair
        if domain_pair:
            domain1, domain2 = domain_pair
            if interaction is None:
                interaction = HYBRIDIZATION_INTERACTION
            #h1end3p = h1end5p = (h2end3p, h2end5p) = ()
        else:
            domain1, domain2 = h1end3p.domain, h2end3p.domain
            if interaction is None:
                interaction = STACKING_INTERACTION

        edge_key = interaction

        strand1 = domain1.strand
        strand2 = domain2.strand
        c = strand1.complex
        assert c == strand2.complex

        result = {'changed_complexes': None,
                  'new_complexes': None,
                  'obsolete_complexes': None,
                  'free_strands': None,
                  'case': None,
                  'is_forming': False,
                 }

        ## CASE 0: NO COMPLEX PRESENT
        if strand1 == strand2 and c is None:
            # The domains are on the same strand and we don't have any complex to update
            # Case 0?
            result['case'] = 0
            result['free_strands'] = [strand1]
            return result

        assert c is not None

        ## First, BREAK CONNECTION:
        if domain_pair is not None:
            ## Update complex graph, breaking the d1-d2 hybridization edge:
            # c.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
            # c.hybridized_pairs.remove(frozenset(domain_pair))
            c.remove_hybridization_edge(domain_pair)
            c.history.append("break_complex_at: Removed edge (%r, %r, %r) and unsetting _hybridization_fingerprint..." %
                             (domain1, domain2, HYBRIDIZATION_INTERACTION))
            #c.strand_graph.remove_edge(strand1, strand2, s_edge_key)  ## complex strand_graph is obsolete
            c._hybridization_fingerprint = None  # Reset hybridization fingerprint
        if stacking_pair is not None:
            ## Update complex graph, breaking (h1end3p, h1end5p) and (h2end3p, h2end5p) stacking edges:
            # c.remove_edge(h1end3p.domain, h1end5p.domain, key=STACKING_INTERACTION)
            # c.remove_edge(h2end3p.domain, h2end5p.domain, key=STACKING_INTERACTION)
            # # Stacked_pairs have direction (h1nd3p, h1end5p) Unless you want it to be the same stacking pair as here?
            # c.stacked_pairs.remove((h1end3p.domain, h1end5p.domain))
            # c.stacked_pairs.remove((h1end3p.domain, h2end5p.domain))
            c.remove_stacking_edge(stacking_pair)
            c.history.append("break_complex_at: Removed edges (%r, %r, %r) and (%r, %r, %r) and unsetting _stacking_fingerprint..." %
                             (h1end3p.domain, h1end5p.domain, STACKING_INTERACTION,
                              h2end3p.domain, h2end5p.domain, STACKING_INTERACTION))
            c._stacking_fingerprint = None  # Reset hybridization fingerprint
        c.history.append(" - break_complex_at: also unsetting _state_fingerprint (was: %s)" % (c._state_fingerprint,))
        c._state_fingerprint = None
        c.domain_distances = {}  # Reset distances:

        ## Then, SEE WHAT HAS HAPPENED TO THE COMPLEX:

        ## Determine the connected component for strand 1:
        dom1_cc_oligos = nx.node_connected_component(self.strand_graph, strand1)
        # Could also use nx.connected_components_subgraphs(c)

        ## Case 1: The two strands are still connected: No need to do anything further
        if strand2 in dom1_cc_oligos:
            ## The two strands are still connected: No need to do anything further
            c.history.append("( - break_complex_at: complex is still intact.)")
            result['case'] = 1
            result['changed_complexes'] = [c]
            # printd("Dehybridize case 0: Complex still intact.")
            return result

        #### The two strands are no longer connected: ####
        c._strands_fingerprint = None
        c.reset_state_fingerprint()

        ## Need to split up. Three cases:
        ## Case 2(a) Two smaller complexes - must create a new complex for detached domain:
        ## Case 3(b) One complex and one unhybridized strand - no need to do much further
        ## Case 4(c) Two unhybridized strands

        dom2_cc_oligos = nx.node_connected_component(self.strand_graph, strand2)
        assert strand2 not in dom1_cc_oligos
        assert strand1 not in dom2_cc_oligos
        dom1_cc_size = len(dom1_cc_oligos)  # cc = connected component
        dom2_cc_size = len(dom2_cc_oligos)
        # printd("dom1_cc_size=%s, dom2_cc_size=%s, len(c.strands)=%s" % (dom1_cc_size, dom2_cc_size, len(c.strands)))

        assert len(c.strands) == dom1_cc_size + dom2_cc_size

        if dom2_cc_size > 1 and dom1_cc_size > 1:
            ## Case 2(a) Two smaller complexes - must create a new complex for detached domain:
            # Determine which of the complex fragments is the major and which is the minor:
            ## TODO: I don't really need the domain-graph, the strand-level should suffice.
            ## (although it is a good check to have while debuggin...)
            # Make sure NOT to make a copy of graph attributes (nodes/edges/etc)
            # connected_component_subgraphs is not implemented for directed graphs!
            # However, you can just do it manually using system graphs instead of complex domain graph:
            # cc_subgraphs = list(connected_component_subgraphs(c, copy=False))
            # if len(cc_subgraphs) != 2:
            #     print("Unexpected length %s of connected_component_subgraphs(c):" % len(cc_subgraphs))
            #     pprint(cc_subgraphs)
            cc1_nodes = nx.node_connected_component(self.domain_graph, domain1)
            cc2_nodes = nx.node_connected_component(self.domain_graph, domain2)
            #graph_minor, graph_major = sorted(cc_subgraphs, key=lambda g: len(g.nodes()))
            if len(cc1_nodes) < len(cc2_nodes):
                graph_minor = self.domain_graph.subgraph(cc1_nodes)
                # Remember to add the departing domain.strand to the new_complex_oligos list:
                new_complex_oligos = set(dom1_cc_oligos)
            else:
                graph_minor = self.domain_graph.subgraph(cc2_nodes)
                new_complex_oligos = set(dom2_cc_oligos)
            # Use graph_minor to initialize; then all other is obsolete
            # c.strands -= new_complex_oligos
            # c.remove_nodes_from(graph_minor.nodes())
            c.history.append("break_complex_at: Breaking off minor complex...")
            removed_hybridization_pairs, removed_stacking_pairs = \
                c.remove_strands(new_complex_oligos, update_graph=True, update_edge_pairs=True)
            c_new = Complex(data=graph_minor, strands=new_complex_oligos)
            # Remember to add the hybridized_pairs and stacked_pairs removed from the major complex.
            c_new.hybridized_pairs |= removed_hybridization_pairs
            c_new.stacked_pairs |= removed_stacking_pairs
            # printd("Dehybridize case (a) - De-hybridization caused splitting into two complexes:")
            # printd(" - New complex: %s, nodes = %s" % (c_new, c_new.nodes()))
            # printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
            # printd(" - graph_minor nodes:", graph_minor.nodes())
            c.history.append(" - break_complex_at: minor complex %r with %s strands broken off." % (c_new, new_complex_oligos))
            # Case 2: Two smaller complexes - must create a new complex for detached domain:
            result['case'] = 2
            result['changed_complexes'] = [c]
            result['new_complexes'] = [c_new]
        elif dom2_cc_size > 1 or dom1_cc_size > 1:
            ## Case 3(b) one complex and one unhybridized strand - no need to do much further
            # Which-ever complex has more than 1 strands is the major complex:
            domain_minor = domain1 if dom1_cc_size == 1 else domain2
            c.history.append("break_complex_at: Breaking off strand...")
            c.remove_strand(domain_minor.strand, update_graph=True)
            # printd("Dehybridize case (b) - De-hybridization caused a free strand to split away:")
            # printd(" - Free strand: %s, nodes = %s" % (domain_minor.strand, domain_minor.strand.nodes()))
            # printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
            result['case'] = 3
            result['changed_complexes'] = [c]
            result['free_strands'] = [domain_minor.strand]
        else:
            ## Case 4(c) Two unhybridized strands
            result['case'] = 4
            result['obsolete_complexes'] = [c]
            result['free_strands'] = [domain1.strand, domain2.strand]
            c.history.append("break_complex_at: Breaking now-obsolete complex into two free strands...")
            c.remove_strands((strand1, strand2), update_graph=True)
            # c.remove_nodes_from(strand1) # iter(nx.Graph) yields nodes.
            # c.remove_nodes_from(strand2)
            assert c.strands == set()
            ## This sometimes fails:
            if not all(len(strandset) == 0 for strandset in c.strands_by_name.values()):
                print(" FAIL: all(len(strandset) == 0 for strandset in c.strands_by_name.values())")
                pprint(c.strands_by_name)
                pprint(c.nodes())
                pprint(c.edges())
            assert all(len(strandset) == 0 for strandset in c.strands_by_name.values())
            assert len(c.nodes()) == 0
            strand1.complex, strand2.complex = None, None
            # printd("Dehybridize case (c) - De-hybridization caused complex to split into two free strands:")
            # printd(" - Free strands 1: %s, nodes1 = %s" % (domain1.strand, domain1.strand.nodes()))
            # printd(" - Free strands 1: %s, nodes1 = %s" % (domain2.strand, domain2.strand.nodes()))
            # printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))

        if domain_pair:
            assert domain1.partner is None
            assert domain2.partner is None
        if stacking_pair:
            assert all(end.stack_partner is None for end in (h1end3p, h2end5p, h2end3p, h1end5p))

        return result
