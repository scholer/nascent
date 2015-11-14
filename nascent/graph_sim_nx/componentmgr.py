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

import os
#import random
from collections import defaultdict, namedtuple
ReactionAttrs = namedtuple('ReactionAttrs', ['is_hybridizing', 'is_intra'])
#import math
#from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
from networkx.algorithms.components import connected_components, connected_component_subgraphs
import numpy as np
# import pdb

# from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
#, R # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K
# from .constants import R, N_AVOGADRO, AVOGADRO_VOLUME_NM3
from .constants import HYBRIDIZATION_INTERACTION, PHOSPHATEBACKBONE_INTERACTION, STACKING_INTERACTION
from .complex import Complex
from .graph_manager import GraphManager
from .debug import printd, pprintd
from nascent.energymodels.biopython import DNA_NN4, energy_tables_in_units_of_R
from .utils import (sequential_number_generator, sequential_uuid_gen)

supercomplex_sequential_id = sequential_number_generator()





class ComponentMgr(GraphManager):
    """
    System-manager class to manage the system state and state changes
    (but it does not control *which* state changes are induced).

    The GraphManager super-class provides everything related to structural analysis:
    loops and intra-complex activity calculation, etc.
    """

    def __init__(self, strands, params, domain_pairs=None):
        GraphManager.__init__(self, strands=strands)
        self.params = params
        # super-complexes: two or more complexes held together by stacking interactions.
        # supercomplex id (scid) => {'complexes': {...}, 'stacking_edges'}
        self.supercomplexes = defaultdict(set)
        # complexes are assigned sequential complex-unique IDs upon instantiation, no need to keep track of order here.
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
        self.N_strands_hybridized = sum(1 for oligo in self.strands if oligo.is_hybridized())

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
            domain_pairs = {d.name: d.name.lower() if d.name == d.name.upper() else d.name.upper()
                            for d in self.domains_list}
            # remove pairs without partner:
            domain_pairs = {d1name: d2name for d1name, d2name in domain_pairs.items()
                            if d2name in self.domains_by_name}
        # allow self-complementarity?
        assert not any(k == v for k, v in domain_pairs.items())
        self.domain_pairs = domain_pairs
        self.hyb_dehyb_file = open(os.path.join(params.get('working_directory', '.'), "hyb_dehyb.py"), 'w')


    def hybridize(self, domain1, domain2):
        """
        Splitting out logic in preparation for Julia implementation.
        returns
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        Where changed_complexes includes new_complexes but not obsolete_complexes.
        Edit: changed_complexes DOES NOT include new_complexes.
        """
        printd("%s.hybridize(%s, %s) invoked..." % (type(self).__name__, domain1, domain2))
        print("domain1, domain2 = (domains_by_duid[%s], domains_by_duid[%s])" % (domain1.duid, domain2.duid),
              file=self.hyb_dehyb_file)
        print("assert domain1.domain_strand_specie = %s and domain2.domain_strand_specie = %s" %
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
        c1 = strand1.complex
        c2 = strand2.complex

        # Update system-level graphs:
        edge_kwargs = {"interaction": HYBRIDIZATION_INTERACTION,
                       "len": 6, # With of ds helix ~ 2 nm ~ 6 bp
                       "weight": 2,
                       #"key": HYBRIDIZATION_INTERACTION
                      }
        #key = (domain1.universal_name, domain2.universal_name, HYBRIDIZATION_INTERACTION)
        s_edge_key = (frozenset((domain1.universal_name, domain2.universal_name)), HYBRIDIZATION_INTERACTION)
        printd("Adding strand_graph edge (%s, %s, key=%s)" % (strand1, strand2, s_edge_key))
        self.strand_graph.add_edge(strand1, strand2, key=s_edge_key, interaction=HYBRIDIZATION_INTERACTION)
        self.domain_graph.add_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION, **edge_kwargs)
        self.ends5p3p_graph.add_edge(domain1.end5p, domain2.end3p, key=HYBRIDIZATION_INTERACTION, **edge_kwargs)
        self.ends5p3p_graph.add_edge(domain2.end5p, domain1.end3p, key=HYBRIDIZATION_INTERACTION, **edge_kwargs)

        # changed_complexes, new_complexes, obsolete_complexes, free_strands = None, None, None, []
        result = {'changed_complexes': None,
                  'new_complexes': None,
                  'obsolete_complexes': None,
                  'free_strands': None,
                  'case': None}

        if strand1 == strand2:
            # If forming an intra-strand connection, no need to make or merge any Complexes
            print("hybridize case 0: intra-strand hybridization.")
            result['case'] = 0
            if strand1.complex:
                assert strand1.complex == strand2.complex
                strand1.complex.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION)
                result['changed_complexes'] = [c1]
            return result

        ## Update complex:
        if c1 and c2:
            ## Both domains are in a complex.
            if c1 == c2:
                ## Case (1): Intra-complex hybridization
                printd("hybridize case 1: intra-complex hybridization.")
                result['case'] = 1
                assert strand1 in c2.strands and strand2 in c1.strands
                # c1.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION) # done below
                c_major = c1
                result['changed_complexes'] = [c_major]
            else:
                ## Case (2): Inter-complex hybridization between two complexes. Merge the two complexs:
                printd("hybridize case 2: inter-complex hybridization.")
                result['case'] = 2
                c_major, c_minor = (c1, c2) if (len(c1.nodes()) >= len(c2.nodes())) else (c2, c1)
                # Import nodes and edges to major complex:
                # add_strands updates: strand.complex, c_major.strands, c_major.strands_by_name
                c_major.add_strands(c_minor.strands, update_graph=False)
                # We use the c_minor graph - rather than strands ^ - to update c_major graph:
                c_major.add_nodes_from(c_minor.nodes(data=True))
                c_major.add_edges_from(c_minor.edges(keys=True, data=True))
                c_major.N_strand_changes += c_minor.N_strand_changes # Currently not keeping track of strand changes.
                ## Delete the minor complex:
                c_minor.strands.clear()     # Clear strands
                c_minor.strands_by_name.clear()
                c_minor.node.clear()        # Clear nodes
                c_minor.adj.clear()         # Clear edges (graph.adj attribute contains edge data for nx graphs)
                result['obsolete_complexes'] = [c_minor]
                result['changed_complexes'] = [c_major]
        elif c1:
            ## Case 3a: domain2/strand2 is not in a complex; use c1
            printd("hybridize case 3a: strand hybridizing to complex.")
            result['case'] = 3
            c1.add_strand(strand2, update_graph=True)
            c_major = c1
            result['changed_complexes'] = [c_major]
        elif c2:
            ## Case 3b: domain1/strand1 is not in a complex; use c2
            printd("hybridize case 3b: strand hybridizing to complex.")
            result['case'] = 3
            c2.add_strand(strand1, update_graph=True)
            c_major = c2
            result['changed_complexes'] = [c_major]
        else:
            ## Case 4: Neither strands are in existing complex; create new complex
            result['case'] = 4
            printd("hybridize case 4: inter-strand hybridization (forming a new complex).")
            new_complex = Complex(strands=[strand1, strand2])
            # new_complex.strands |= {strand1, strand2}
            # strand1.complex = strand2.complex = new_complex
            # new_complex.add_strand(strand1)
            # new_complex.add_strand(strand2)
            c_major = new_complex
            result['new_complexes'] = [new_complex]

        # Create the hybridization connection in the major complex graph:
        c_major.add_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION, **edge_kwargs)
        c_major.domain_distances = {} # Reset distances

        assert strand1.complex == strand2.complex != None

        print("print('- hybridize complete.')", file=self.hyb_dehyb_file)
        return result


    def dehybridize(self, domain1, domain2):
        """
        Dehybridize domain2 from domain1.
        Returns:
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        """
        printd("%s.dehybridize(%s, %s) invoked..." % (type(self).__name__, domain1, domain2))
        print("domain1, domain2 = (domains_by_duid[%s], domains_by_duid[%s])" % (domain1.duid, domain2.duid),
              file=self.hyb_dehyb_file)
        print("assert domain1.domain_strand_specie = %s and domain1.domain_strand_specie = %s" %
              (domain1.domain_strand_specie, domain2.domain_strand_specie), file=self.hyb_dehyb_file)
        print("sysmgr.dehybridize(domain1, domain2)", file=self.hyb_dehyb_file)
        assert domain1 != domain2
        assert domain1.partner == domain2 != None
        assert domain2.partner == domain1 != None

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
        printd("%s: Removing strand_graph edge (%s, %s, key=%s)" % (type(self).__name__, strand1, strand2, s_edge_key))
        self.strand_graph.remove_edge(strand1, strand2, key=s_edge_key)
        self.domain_graph.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        self.ends5p3p_graph.remove_edge(domain1.end5p, domain2.end3p, key=HYBRIDIZATION_INTERACTION)
        self.ends5p3p_graph.remove_edge(domain2.end5p, domain1.end3p, key=HYBRIDIZATION_INTERACTION)

        if strand1 == strand2 and c is None:
            # The domains are on the same strand and we don't have any complex to update
            return None, None, None, [strand1]
        assert c is not None

        result = {'changed_complexes': None,
                  'new_complexes': None,
                  'obsolete_complexes': None,
                  'free_strands': None,
                  'case': None}

        ## Update complex graph, breaking the d1-d2 hybridization edge:
        c.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        #c.strand_graph.remove_edge(strand1, strand2, s_edge_key)  ## complex strand_graph is obsolete

        c.domain_distances = {}  # Reset distances:
        c._hybridization_fingerprint = None  # Reset hybridization fingerprint

        ## Determine the connected component for strand 1:
        dom1_cc_oligos = nx.node_connected_component(self.strand_graph, strand1)
        # Could also use nx.connected_components_subgraphs(c)

        if strand2 in dom1_cc_oligos:
            ## The two strands are still connected: No need to do anything further
            result['case'] = 0
            result['changed_complexes'] = [c]
            printd("Dehybridize case 0: Complex still intact.")
            return result


        #### The two strands are no longer connected: ####
        c._strands_fingerprint = None

        ## Need to split up. Three cases:
        ## Case (a) Two smaller complexes - must create a new complex for detached domain:
        ## Case (b) One complex and one unhybridized strand - no need to do much further
        ## Case (c) Two unhybridized strands

        dom2_cc_oligos = nx.node_connected_component(self.strand_graph, strand2)
        assert strand2 not in dom1_cc_oligos
        assert strand1 not in dom2_cc_oligos
        dom1_cc_size = len(dom1_cc_oligos)  # cc = connected component
        dom2_cc_size = len(dom2_cc_oligos)
        # printd("dom1_cc_size=%s, dom2_cc_size=%s, len(c.strands)=%s" % (dom1_cc_size, dom2_cc_size, len(c.strands)))

        assert len(c.strands) == dom1_cc_size + dom2_cc_size

        if dom2_cc_size > 1 and dom1_cc_size > 1:
            # Case (a) Two smaller complexes - must create a new complex for detached domain:
            # Determine which of the complex fragments is the major and which is the minor:
            ## TODO: I don't really need the domain-graph, the strand-level should suffice.
            ## (although it is a good check to have while debuggin...)
            # Make sure NOT to make a copy of graph attributes (nodes/edges/etc)
            cc_subgraphs = list(connected_component_subgraphs(c, copy=False))
            if len(cc_subgraphs) != 2:
                print("Unexpected length %s of connected_component_subgraphs(c):" % len(cc_subgraphs))
                pprint(cc_subgraphs)
            graph_minor, graph_major = sorted(cc_subgraphs, key=lambda g: len(g.nodes()))
            # Remember to add the departing domain.strand to the new_complex_oligos list:
            new_complex_oligos = set(dom2_cc_oligos if domain1 in graph_major.nodes() else dom1_cc_oligos)
            # Use graph_minor to initialize; then all other is obsolete
            # c.strands -= new_complex_oligos
            # c.remove_nodes_from(graph_minor.nodes())
            c.remove_strands(new_complex_oligos, update_graph=True)
            c_new = Complex(data=graph_minor, strands=new_complex_oligos)
            printd("Dehybridize case (a) - De-hybridization caused splitting into two complexes:")
            printd(" - New complex: %s, nodes = %s" % (c_new, c_new.nodes()))
            printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
            printd(" - graph_minor nodes:", graph_minor.nodes())
            printd(" - graph_major nodes:", graph_major.nodes())
            result['case'] = 1
            result['changed_complexes'] = [c]
            result['new_complexes'] = [c_new]
            #changed_complexes.append(c_new)
        elif dom2_cc_size > 1 or dom1_cc_size > 1:
            # Case (b) one complex and one unhybridized strand - no need to do much further
            # Which-ever complex has more than 1 strands is the major complex:
            domain_minor = domain1 if dom1_cc_size == 1 else domain2
            c.remove_strand(domain_minor.strand, update_graph=True)
            printd("Dehybridize case (b) - De-hybridization caused a free strand to split away:")
            printd(" - Free strand: %s, nodes = %s" % (domain_minor.strand, domain_minor.strand.nodes()))
            printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
            result['case'] = 2
            result['changed_complexes'] = [c]
            result['free_strands'] = [domain_minor.strand]
        else:
            # Case (c) Two unhybridized strands
            result['case'] = 3
            result['obsolete_complexes'] = [c]
            result['free_strands'] = [domain1.strand, domain2.strand]
            c.remove_strands({strand1, strand2})
            c.remove_nodes_from(strand1) # iter(nx.Graph) yields nodes.
            c.remove_nodes_from(strand2)
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
            printd("Dehybridize case (c) - De-hybridization caused complex to split into two free strands:")
            printd(" - Free strands 1: %s, nodes1 = %s" % (domain1.strand, domain1.strand.nodes()))
            printd(" - Free strands 1: %s, nodes1 = %s" % (domain2.strand, domain2.strand.nodes()))
            printd(" - Old complex: %s, nodes = %s" % (c, c.nodes()))

        assert domain1.partner is None
        assert domain2.partner is None
        print("print('- dehybridize complete.')", file=self.hyb_dehyb_file)
        return result


    def n_hybridized_domains(self):
        """ Count the number of hybridized domains. """
        count = sum(1 for domain in self.domains_list if domain.partner is not None)
        if not count % 2 == 0:
            print("Weird - n_hybridized_domains counts to %s (should be an even number)" % count)
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.domains_list
                                                   if domain.partner is not None))
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.hybridized_domains_by_name))
        return count


    def n_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for oligo in self.strands if oligo.is_hybridized())


    def stack(self, h1end3p, h1end5p, h2end3p, h2end5p):
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

        ## Variable unpacking:
        h1domain1 = h1end3p.domain
        h1domain2 = h1end5p.domain
        h2domain1 = h2end3p.domain
        h2domain2 = h2end5p.domain
        h1strand1 = h1domain1.strand
        h1strand2 = h1domain2.strand
        h2strand1 = h2domain1.strand
        h2strand2 = h2domain2.strand
        c1 = h1strand1.complex
        c2 = h1strand2.complex
        if c1:
            assert c1 == h2strand1.complex
            assert c2 == h2strand2.complex

        ## Assertions:
        assert h1end3p.stack_partner is None
        assert h1end5p.stack_partner is None
        assert h2end3p.stack_partner is None
        assert h2end5p.stack_partner is None

        ## Update domain end attributes:
        h1end3p.stack_partner, h1end5p.stack_partner = h1end5p, h1end3p
        h2end3p.stack_partner, h2end5p.stack_partner = h2end5p, h2end3p

        stack_string = "%s%s/%s%s" % (h1end3p.base, h1end5p.base, h2end3p.base, h2end5p.base)
        if stack_string not in DNA_NN4:
            stack_string = stack_string[::-1]
            assert stack_string in DNA_NN4

        h1end3p.stack_string = h1end5p.stack_string = h2end3p.stack_string = h2end5p.stack_string = stack_string

        ## Update system-level graphs:
        edge_kwargs = {"interaction": STACKING_INTERACTION,
                       ## TODO: Decide whether "len" is in bp or nm.
                       ## You could have separate attrs, length_nm and length_bp, but there should be just one "len".
                       "len": 0.5, # With of ds helix ~ 2 nm ~ 6 bp
                       "weight": 2,
                       #"key": HYBRIDIZATION_INTERACTION
                      }
        h1_edge_key = (frozenset((h1end3p.instance_name, h1end5p.instance_name)), STACKING_INTERACTION)
        h2_edge_key = (frozenset((h2end3p.instance_name, h2end5p.instance_name)), STACKING_INTERACTION)
        printd("Adding strand_graph edge (%s, %s, key=%s)" % (h1strand1, h1strand2, h1_edge_key))
        self.strand_graph.add_edge(h1strand1, h1strand2, key=h1_edge_key, interaction=STACKING_INTERACTION)
        printd("Adding strand_graph edge (%s, %s, key=%s)" % (h2strand1, h2strand2, h2_edge_key))
        self.strand_graph.add_edge(h2strand1, h2strand2, key=h2_edge_key, interaction=STACKING_INTERACTION)

        self.domain_graph.add_edge(h1domain1, h1domain2, key=STACKING_INTERACTION, **edge_kwargs)
        self.domain_graph.add_edge(h2domain1, h2domain2, key=STACKING_INTERACTION, **edge_kwargs)
        self.ends5p3p_graph.add_edge(h1end3p, h1end5p, key=STACKING_INTERACTION, **edge_kwargs)
        self.ends5p3p_graph.add_edge(h2end3p, h2end5p, key=STACKING_INTERACTION, **edge_kwargs)

        # changed_complexes, new_complexes, obsolete_complexes, free_strands = None, None, None, []
        result = {'changed_complexes': None,
                  #'new_complexes': None,
                  #'obsolete_complexes': None,
                  #'free_strands': None,
                  #'case': None
                 }
        if c1 is not None and c2 is not None:
            # We may or may not want to update the complex reactions if stacking caused significant change...
            # This is also relevant for interactions (stacking and hybridization) between two complexes.
            result['changed_complexes'] = [c1] if c1 is c2 else [c1, c2]

        return result


    def unstack(self, h1end3p, h1end5p, h2end3p, h2end5p):
        """
        Break a stacking interaction.
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p
        """

        ## Variable unpacking:
        h1domain1 = h1end3p.domain
        h1domain2 = h1end5p.domain
        h2domain1 = h2end3p.domain
        h2domain2 = h2end5p.domain
        h1strand1 = h1domain1.strand
        h1strand2 = h1domain2.strand
        h2strand1 = h2domain1.strand
        h2strand2 = h2domain2.strand
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
        h1_edge_key = (frozenset((h1end3p.instance_name, h1end5p.instance_name)), STACKING_INTERACTION)
        h2_edge_key = (frozenset((h2end3p.instance_name, h2end5p.instance_name)), STACKING_INTERACTION)
        printd("Removing strand_graph edge (%s, %s, key=%s)" % (h1strand1, h1strand2, h1_edge_key))
        self.strand_graph.remove_edge(h1strand1, h1strand2, key=h1_edge_key, interaction=STACKING_INTERACTION)
        printd("Removing strand_graph edge (%s, %s, key=%s)" % (h2strand1, h2strand2, h2_edge_key))
        self.strand_graph.remove_edge(h2strand1, h2strand2, key=h2_edge_key, interaction=STACKING_INTERACTION)

        self.domain_graph.remove_edge(h1domain1, h1domain2, key=STACKING_INTERACTION)
        self.domain_graph.remove_edge(h2domain1, h2domain2, key=STACKING_INTERACTION)
        self.ends5p3p_graph.remove_edge(h1end3p, h1end5p, key=STACKING_INTERACTION)
        self.ends5p3p_graph.remove_edge(h2end3p, h2end5p, key=STACKING_INTERACTION)

        # changed_complexes, new_complexes, obsolete_complexes, free_strands = None, None, None, []
        result = {'changed_complexes': None}
        if c1 is not None and c2 is not None:
            # We may or may not want to update the complex reactions if stacking caused significant change...
            # This is also relevant for interactions (stacking and hybridization) between two complexes.
            result['changed_complexes'] = [c1] if c1 is c2 else [c1, c2]

        return result
