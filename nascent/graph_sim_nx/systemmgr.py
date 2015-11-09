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

# pylint: disable=C0103,W0142,W0212

"""

Module for managing the whole system.

I think there might be a lot of code involved in managing the system (graphs and structure),
and having all that in the simulator makes it a bit hard to read.

Splitting the "system state" code out to a separate module allows the simulator to be focused
on just the "stochastic simulation" part and be mostly agnostic on the system details.

Question: Who manages "reaction propensity" etc?
 - I'm thinking the system manager takes care of everything except the simulation steps / temperature control.

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
 * Propencity functions
 *

"""

import os
#import random
from collections import defaultdict, namedtuple
ReactionAttrs = namedtuple('ReactionAttrs', ['is_hybridizing', 'is_intra'])
#import math
from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
from networkx.algorithms.components import connected_components, connected_component_subgraphs
import numpy as np
import pdb

from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
from .constants import R, N_AVOGADRO, AVOGADRO_VOLUME_NM3 #, R # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K
from .constants import PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION
from .complex import Complex
from .structure_analyzer import StructureAnalyzer
from .nx_utils import draw_graph_and_save







class SystemMgr(StructureAnalyzer):
    """
    System-manager class to manage the system state and state changes
    (but it does not control *which* state changes are induced).

    The StructureAnalyzer super-class provides everything related to structural analysis:
    loops and intra-complex activity calculation, etc.
    """

    def __init__(self, volume, strands, params, domain_pairs=None):
        StructureAnalyzer.__init__(self, strands=strands)
        self.params = params
        self.temperature = params.get('temperature', 300)
        self.volume = volume or params.get('volume')
        self.fnprefix = ""
        # Base single-molecule activity of two free single reactants in volume.
        self.specific_bimolecular_activity = 1/self.volume/N_AVOGADRO # × M
        self.include_steric_repulsion = False # True
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
        # Keep track of domain state depletions.
        # If the (grouped) domain species count occilates between 0 and 1, then the grouped approach might be
        # disadvantageous. Instead, propensities should be evaluated for each domain instance.
        # Or you can have a mixed model where hybridized/duplexed domains and domains on free strands are grouped
        # while unhybridized domains within a complex are considered individually.
        # Or you could consider intra-complex reactions individually while group other reactions.
        self.N_state_depletions = 0
        self.N_state_repletions = 0

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
        #self.upper_case_domains = [d for d in self.domains if d == d.upper()]
        #self.lower_case_domains = [d for d in self.domains if d == d.lower()]

        print("Simulation system manager initiated at V=%s with %s strands spanning %s domains." \
              % (self.volume, len(self.strands), len(self.domains)))


        ## Symbol nomenclature:
        ## Sᵢ, S₁, S₂ - domain species (domain name). Domains with the same sequence are of the same specie.
        ## Fᵢ, F₁, F₂ - domain state fingerprint - domains of the same specie and in the same strand/complex state.
        ## Dᵢ, D₁, D₂ - domain instances, individual domain molecules.

        ### Caches: ###
        self.cache = {} # defaultdict(dict)
        # Standard enthalpy and entropy of hybridization,
        # indexed as [frozenset((d1.name, d2.name))][0 for enthalpy, 1 for entropy]
        # - Note: index with frozenset((a,b)) or just cache[a, b] = cache[b, a] = value? # d[1,2] same as d[(1,2)]
        # --> Creating a set() or frozenset() takes about 10x longer than to make tuple.
        # --> dict assignment with frozenset is 0.4/0.5 us vs 0.17/0.05 for the "double tuple entry" (python/pypy),
        #       so if you have the memory for it, go for the faster tuples which takes 2x memory.
        # Note: Whenever I use "name", I'm refering to the name of the specie - domain or strand specie.
        self.cache['domain_hybridization_energy'] = self.domain_dHdS = {} # indexed simply as Sᵢ or {S₁, S₂} ?

        ## Relative activity - used to moderate the activity of domains based on their accessibility in a complex.
        # Note: Not sure whether to use this or use a "loop energy" type of activation energy.
        # Relative activities, indexed as:
        #  - [({domain1-specie, domain2-specie}, complex-state-fingerprint)] for intra-complex reactions
        #  - [(domain1-specie, domain1-complex-state-fingerprint)] for inter-complex reactions.
        # relative activity is 1 for free strands
        # For inter-complex reactions (including free strands), the combined relative activity is simply the product:
        #   rel_activity = dom1_rel_activity * dom2_rel_activity

        self.cache['intracomplex_activity'] = {}    # indexed by {F₁, F₂}
        self.cache['stochastic_rate_constant'] = {} # indexed by {F₁, F₂}
        # For surface/steric/electrostatic reduction in k_on.
        # TODO: Check surface-hybridization kinetics papers whether reduction in k_on yields same reduction in k_off.
        #       I.e. if the surface reduces Tm, hybridization energy or equivalent.
        #       For electrostatic and possibly steric, I would expect a reuction in Tm but not for surfaces,
        #       so it is a question whether these can in fact be treated equally.
        self.cache['domain_accessibility'] = {}     # indexed by F₁
        # Note: F₁, F₂ is domain-state-fingerprints. It is supposed to be invariant for domain instances
        # of the same domain specie in the same complex state. E.g. two free "domain A" should be in the same state.
        # Alternative name: dom_specie_spec

        #self.relative_activity_cache = defaultdict(dict)

        ## State-dependent hybridization energy cache
        # I think this was for hybridization that required e.g. bending or zippering...
        self._statedependent_dH_dS = {}  # indexed as {Fᵢ}

        self.possible_hybridization_reactions = {}  # {d1, d2} => c_j
        self.reaction_attrs = {}                    # {d1, d2} => (is_hybridizing, is_intra)
        # When we are not grouping reactions by domain state species, then propensity functions are just the same
        # as propensity constants... (X₁==X₂==X₁₂==1)
        self.propensity_functions = self.possible_hybridization_reactions
        if strands:
            self.init_possible_reactions()



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




    def init_possible_reactions(self):
        """
        Reactions have:
            reaction_spec => propensity_constant c_j, state_change_vector v_j
        However, I do:
            ({domain1, domain2}, is_hybridizing, is_intracomplex) => propensity_constant.
        Where domspec is a domain-state fingerprint, and
            {domspec1, domspec2} = domspec_pair = {F₁, F₂}
            F₁, F₂ are symbols for state-dependent "finger-prints" aka "state-species" for the domains:
                hash((dspecie, c_state, self.in_complex_identifier()))

        We need to find possible reactions. Possible ways to do that:
         1. Use self.domain_pairs - domain species names.
         2. Go over all domains.

        Although I expect all domains to be in the same start at the onset of the simulation,
        I figure it is best to go through them all, determining their state-specie,
        just to make sure I have it right.

        I am grouping domains into sub-populations according to their state.
        I do this to reduce the memory required for caching and increasing the chance of cache hits.

        I am basically just using a "domain subspecie state hash" (domspec) to denote a set of
        particular domains all in the same state.
        E.g. a domspec 928457 could be denote any/all "domain1 on strandA that are free (not in a complex)"
        To get a list of all actual domains matching a certain domain subpopulation domspec, use:
            self.domain_state_subspecies[domspec] - will you all domains matching <domspec>
        Maybe a better name would be domain_subpopulation_by_state?
        """
        print("\nInitializing possible reactions (UNGROUPED)...")
        Rxs = self.possible_hybridization_reactions
        for d1 in self.domains:
            # d1_domspec = d1.domain_state_fingerprint()
            if d1.partner is not None:
                is_hyb, is_intra = False, True
                # If d1 is hybridized, then there is only one possible reaction channel: dehybridizing
                d2 = d1.partner
                # d2_domspec = d2.domain_state_fingerprint()
                domain_pair = frozenset((d1, d2))
                #reaction_spec = (domspec_pair, is_hyb, is_intra)
                Rxs[domain_pair] = self.calculate_c_j(d1, d2, is_hybridizing=is_hyb, is_intra=is_intra)
                self.reaction_attrs = (is_hyb, is_intra)
                # self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                # self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
            else:
                # find all possible reaction channels for domain d1
                # self.domain_pairs[dname] = [list of complementary domain names]
                # all possible hybridization partners/candidates:
                is_hyb = True
                if d1.name not in self.domain_pairs:
                    # No domain partners for this domain
                    continue
                for d2 in (d for cname in self.domain_pairs[d1.name]
                           for d in self.unhybridized_domains_by_name[cname]):
                    assert d2.partner is None
                    # d2_domspec = d2.domain_state_fingerprint()
                    domain_pair = frozenset((d1, d2))
                    is_intra = (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex) or \
                               (d1.strand == d2.strand)  # intra-complex OR intra-strand reaction
                    #reaction_spec = (domspec_pair, is_hyb, is_intra)
                    if domain_pair not in Rxs:
                        # R_j = (c_j, v_j) - propensity constant for reaction j
                        Rxs[domain_pair] = self.calculate_c_j(d1, d2, is_hybridizing=True, is_intra=is_intra)
                    # self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                    # self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
        print(len(self.possible_hybridization_reactions), "possible hybridization reactions initialized.")


    def update_possible_reactions(self, changed_domains, d1, d2, is_hybridizing, reaction_doms_specs=None):
        """
        Maybe it is better to just re-calculate a_j for all changed domains,
        and then filter out reactions with a_j = 0,
        rather than bothering with keeping track of depleted domain state species?

        Wow, this is getting really complex.
        I'm almost contemplating going back to the old approach where we do not try to group and "count"
        domain states, but instead we just consider them all independent.
        Pros for simulating without grouping domains by domain state:
        - We don't have to keep track of possible_hybridization_reactions, etc.
        - Each domain instance has its own propensity function; (just beware of (d1, d2) and (d2, d1) duplicates).
        - At least for small copy numbers, grouping and keeping track of domain state count is likely inefficient.
        - Because of the large number of possible states, even for large copy numbers, the state count might occilate
            between 0 and 1. That is particularly unfavorable for the "grouped" strategy.
        - We can still use "domain state fingerprints" for caching.
        Cons:
        - If we have large copy numbers, a_j becomes very large. Well... It becomes the number of domains.
            Although, even if we have 1000 domain instances, it still only takes 0.1 ms to find a suitable a_j.
        - If we have 10 domain "a" copies and 10 complementary domain "A" copies, then there are 10*10=100 combinations.
            If
        - Grouping has other advantages, e.g. caching.

        Instead, I guess we could have

        # Set vs list:
        # - Lists are slightly faster to iterate over;
        # - Sets are significantly faster for lookups and removing arbitrary elements (about 80 times faster)

        """
        print(("\nupdate_possible_reactions invoked with d1=%s, d2=%s, is_hybridizing=%s, reaction_doms_specs=%s, "
               "changed_domains:") % (d1, d2, is_hybridizing, reaction_doms_specs))
        pprint(changed_domains)


        n_possible_start = len(self.possible_hybridization_reactions)
        # n_species_start = len(self.domain_state_subspecies)

        # Add changed domains to new species
        depleted_domspecs = set()
        changed_domspecs = set()
        new_domspecs = set()
        # old_d1d2_doms_specs = frozenset((d1._specie_state_fingerprint, d2._specie_state_fingerprint))
        updated_reactions = set()

        if is_hybridizing:
            self.unhybridized_domains_by_name[d1.name].remove(d1)
            self.unhybridized_domains_by_name[d2.name].remove(d2)
            self.hybridized_domains_by_name[d1.name].add(d1)
            self.hybridized_domains_by_name[d2.name].add(d2)
        else:
            self.hybridized_domains_by_name[d1.name].remove(d1)
            self.hybridized_domains_by_name[d2.name].remove(d2)
            self.unhybridized_domains_by_name[d1.name].add(d1)
            self.unhybridized_domains_by_name[d2.name].add(d2)

        # If d1, d2 is hybridizing, then we need to eliminate reactions involving d1, d2 from
        # possible_reactions.
        if is_hybridizing:
            obsolete_reactions = [domain_pair for domain_pair in self.possible_hybridization_reactions
                                  if d1 in domain_pair or d2 in domain_pair]
            for domain_pair in obsolete_reactions:
                del self.possible_hybridization_reactions[domain_pair]
                del self.reaction_attrs[domain_pair]
            # You can keep track of domain hybridization reactions, but then
            # please note that the reaction can have become obsolete by the other domain being hybridized.
            # It is thus only a set of possible reactions that can be deleted.
            # Add de-hybridization reaction for d1, d2: (will be skipped when processing all changed domains)
            domain_pair = frozenset((d1, d2))
            self.possible_hybridization_reactions[domain_pair] = \
                self.calculate_c_j(d1, d2, is_hybridizing=False, is_intra=True)
            self.reaction_attrs[domain_pair] = ReactionAttrs(is_hybridizing=False, is_intra=True)
        else:
            # d1 and d2 have just been de-hybridized;
            # Find new possible partners for d1 and d2:
            # d1 and d2 are processed together with the other changed domains
            pass


        for domain in changed_domains:
            # IMPORTANT: changed_domains must not contain any
            old_domspec = domain._specie_state_fingerprint
            # For d1, d2, old_domspec yields the new domspec instead!
            # print("Re-setting and re-calculating domain_state_fingerprint for %s - old is: %s"
            #       % (domain, domain._specie_state_fingerprint))
            domain.state_change_reset()
            new_domspec = domain.domain_state_fingerprint()
            if new_domspec == old_domspec:
                print("\nWeird: new_domspec == old_domspec for changed domain %s: %s == %s" %
                      (domain, new_domspec, old_domspec))

            if domain.partner is not None:
                # Currently, we don't re-calculate hybridized domains that are already hybridized (other than d1, d2)
                continue
            else:
                ## No partner -- consider hybridization reactions with all other unhybridized complementary domains.
                # for d2 in [d for cname in self.domain_pairs[d1.name]
                #            for d in self.unhybridized_domains_by_name[cname]]:
                #pdb.set_trace()
                if domain.name in self.domain_pairs:
                    for cname in self.domain_pairs[domain.name]:
                        for domain2 in self.unhybridized_domains_by_name[cname]:
                            # Remove old reaction propensity:
                            # (TODO: Consider only updating/removing intra-complex reactions)
                            domain_pair = frozenset((domain, domain2))
                            is_intra = (domain.strand.complex is not None and \
                                        domain.strand.complex == domain2.strand.complex) \
                                        or (domain.strand == domain2.strand)  # intra-complex OR intra-strand reaction
                            self.possible_hybridization_reactions[domain_pair] = \
                                self.calculate_c_j(domain, domain2, is_hybridizing=True, is_intra=is_intra)
                            self.reaction_attrs[domain_pair] = ReactionAttrs(is_hybridizing=True, is_intra=is_intra)
                            if domain_pair not in updated_reactions:
                                updated_reactions.add(domain_pair)


    def print_reaction_stats(self):
        """ Print information on reaction pathways. """
        print("Reaction pathways:")
        for reaction_spec, c_j,  in self.possible_hybridization_reactions.items():
            #domspec1, domspec2 = tuple(reaction_spec[0])
            domain1, domain2 = tuple(reaction_spec)
            is_hybridizing = domain1.partner is None
            print("%s %s %s" % (domain1, "hybridize" if is_hybridizing else "de-hybrid", domain2),
                  (": %0.03e x 1 x 1 = %03e" % (
                      c_j, self.propensity_functions[reaction_spec])
                  ) if is_hybridizing else (":   %0.03e x   1   = %03e" % (
                      c_j, self.propensity_functions[reaction_spec]))
                 )


    def map_reaction_specs_by_domspec(self):
        """
        Generate an equivalent to self.hybridization_reactions_by_domspec
        from self.possible_hybridization_reactions, i.e. a map of
            domspec => {set of reaction_specs involving this domspec}
        """
        reactions_by_domspec = defaultdict(set)
        for domain_pair in self.possible_hybridization_reactions:
            for domspec in domain_pair:
                reactions_by_domspec[domspec].add(domain_pair)
        return reactions_by_domspec


    def get_reactions_by_domspec(self, domspec):
        """ Filter possible_hybridization_reactions, returning only reactions involving :domspec:. """
        return {reaction_spec: v
                for reaction_spec, v in self.possible_hybridization_reactions.items()
                if domspec in reaction_spec[0]}


    def calculate_c_j(self, d1, d2, is_hybridizing, is_intra=None):
        """
        :is_intra: whether the reaction is intra-molecular; that is, either intra-complex or intra-strand.
        Calculate propensity constant, c_j, for hybridization or dehybridization of d1 and d2.
        Edit: Changed name from "propensity constant" to just "c_j", to avoid confusion with propensity *function*, a_j.
        Alternative names:
        * specific/molecular/instance/individual reaction frequency
        * specific: specific to a single pair of molecules from S₁, S₂.
        * molecular/instance/individual: as opposed to the total propensity for the entire species population, a_j.
        * "frequency" because of unit s⁻¹.
        ** Could also be "propensity", but could be confused with propensity *function*, a_j.
        * "stochastic rate constant"  ⁽¹,²⁾
        * "specific probability rate constant" (³)
        Regarding "probability" vs "frequency":
        * c_j has units of s⁻¹, so probability (which is unitless) wouldn't be correct.
        * The product "c_j∙dt" gives the probability that a selected pair of instances
            will undergo reaction R_j in timespan dt.

        For bi-molecular reactions, the propensity constant is k_on rate constant times concentration:
            c = k_on * [domain]   # c is propensity, not concentration.
        For uni-molecular reactions, the propensity constant is simply the k_off rate constant:
            c = k_off

        Differences between rate constant, k_j, and c_j:
        * c_j is stochastic, k_j is *mass-action* rate constant.⁽¹⁾
        * c_j = k_j / (N_A ∙ volume)ᴺ⁻¹   # N is the number of reactants, e.g. 2 for bimolecular reactions.
        * v_j = k_j [S₁] [S₁]   # reaction rate, has unit of M/s. (Note: not state change vector ν_j, nu)
        * a_j = c_j x₁ x₂       # propensity function, has unit of s⁻¹

        Question: Where would it be appropriate to include effects of intra-complex reactions:
        * Effects related to effective volume / effective/relative activity should be included
            when calculating c_j. c_j is volume dependent, k_j is independent on volume.
        * Effects related to hybridization energy should be included in k_j,
            since k_j *is* dependent on hybridization energy. This could include:
            * Bending strain energy (decreasing k_on) - although bending could also be interpreted as influencing
                the spatial probability distribution function and thus be kind of volumetric in nature.
            * Zippering strain (increasing k_off),
            * Electrostatic interactions.
        Not sure about steric interactions. That is probably volumetric, affecting spatial probability
        districution function, to be exact. The same could be said about electrostatic...
        I guess if an interaction affects k_off, then it is energetic in nature, while if
        it only affects k_on, then it is volumetric in nature (spatial pdf).
        For a first approximation, maybe only consider the effects that can be interpreted in terms
        of effective concentration / relative activity. Then you can assume constant rate constants,
        k_on/k_off (except for stacking affecting k_off), but otherwise have all effect apply only
        to c_j and never k_j.
        These effect could even be isolated to a single factor, relative_activity(d1, d2).

        Question: What is the best way to capture volume pdf effects?

        Question: What should relative activity entail?
            Should rel act be 1 for two free strands - or 1/volume? -- Thus basically the stochastic activity..
            Maybe start by implementing it here and not worry about "relative" vs actual activity.

        TODO: Account for INTRA-complex reactions.
            Question: Should this be done here or in hybridization_rate_constant(d1, d2) ?
        TODO: Account for steric interferrence in INTER-complex reactions.

        Refs:
            [1]: http://people.uleth.ca/~roussel/C4000biochemnet/slides/14CME.pdf
            [2]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1303974/ "stochastic rate constant cμ"
            [3]: http://www.async.ece.utah.edu/BioBook/lec4.pdf
        """
        if is_hybridizing:
            # if d1.strand.complex == d2.strand.complex != None:
            # Edit: We might have a single strand interacting with itself
            if is_intra:
                assert d1.strand.complex == d2.strand.complex
                assert d1.strand.complex is not None or d1.strand == d2.strand
                # Intra-complex reaction:
                try:
                    stochastic_activity = self.intracomplex_activity(d1, d2)
                except nx.NetworkXNoPath as e:
                    c = d1.strand.complex
                    print(("ERROR: No path between d1 %s and d2 %s "
                           "which are both in complex %s. ") % (d1, d2, c))
                    print("complex.nodes():")
                    pprint(c.nodes())
                    print("complex.edges():")
                    pprint(c.edges())
                    workdir = self.params.get("working_directory", os.path.expanduser("~"))
                    draw_graph_and_save(c, outputfn=os.path.join(workdir, "complex_graph.png"))
                    draw_graph_and_save(self.strand_graph, outputfn=os.path.join(workdir, "strand_graph.png"))
                    draw_graph_and_save(self.domain_graph, outputfn=os.path.join(workdir, "domain_graph.png"))
                    draw_graph_and_save(self.ends5p3p_graph, outputfn=os.path.join(workdir, "ends5p3p_graph.png"))
                    raise e
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
            if self.include_steric_repulsion:
                steric_activity_factor = np.prod([1 if d.strand.complex is None else self.steric_activity_factor(d)
                                                  for d in (d1, d2)])
                stochastic_activity *= steric_activity_factor

            k_on = self.hybridization_rate_constant(d1, d2)
            # hybridization rates of free strands can be approximated as being constant
            # (somewhere between 1e5 /M/s and 1e6 /M/s). This is because the activation energy is fairly negligible.
            # Hybridization rate constant, k, is in unit of /M/s = L/mol/s.
            # Activity is always in units of M - so resultant
            c_j = k_on * stochastic_activity # should be e.g. 0.1 /s

            ## Should we include stacking interactions in k_on? (That would increase k_on...)
            # a) Yes: We include them in k_off, so they must also always be included in k_on (because...)
            # b1) No: Stacking and unstacking is separate from hybridization/dehybridization:
            #           (unhybridized, unstacked) ←→ (hybridized, unstacked) ←→ (hybridized, stacked)
            #           To de-hybridize, a domain must first un-stack.
            #           AND WE ONLY CONSIDER DE-HYBRIDIZATION REACTION FOR UN-STACKED DOMAINS
            #           (No.. because stacking reactions might be much faster than hybridization)
            # b2) No: We have a separate reaction for stacking/unstacking.
            #           If the domain is stacked, that is included in k_off, but the domain can spontaneously unstack.
            # c) From a thermodynamic/kinetic perspective: Not sure.
            #           On one hand, I wouldn't expect stacking to increase hybridization rate.
            #           On second thought:
            # d) Thermodynamically: If you include stacking interactions in k_off, you CANNOT also include them in k_on,
            #           as this would effectively include it twice:
            #           ∆G/RT = ln(K) = ln(k_on/k_off) = ln(k_on) - ln(k_off)
            #                 = ln(k_on0*k_stacking) - ln(k_off0/k_stacking) = ln(k_on) - ln(k_off) + 2*ln(k_stacking)
        else:
            # De-hybridization reaction;
            # k_off depends on ΔG°  (at least when we are at T < Tm at molar concentrations where ΔG° < 0)
            k_off = self.dehybridization_rate_constant(d1, d2)
            # For uni-molecular reactions, c_j = k_j
            c_j = k_off
        return c_j #, is_hybridizing

    def steric_activity_factor(self, d):
        """
        """
        # TODO: Implement steric_activity_factor
        return 1

    def intracomplex_activity(self, d1, d2):
        """
        Return the activity for hybridization of two domains within a complex.
        TODO: Add steric/electrostatic/surface repulsion (individual activity contribution for d1 and d2,
              so that activity = intercomplex_activity*d1_activity*d2_activity)
        """
        assert d1.strand.complex == d2.strand.complex != None
        cache_key = frozenset([d.domain_state_fingerprint() for d in (d1, d2)])
        cache = self.cache['intracomplex_activity']
        if cache_key in cache:
            return cache[cache_key]
        #activity = d1.strand.complex.intracomplex_activity(d1, d2)
        activity = super().intracomplex_activity(d1, d2)
        cache[cache_key] = activity
        return activity


    def hybridization_rate_constant(self, d1, d2):
        #T = self.temperature
        # At T >> Tm (ΔG° ~ 0), this depends on ΔG°.
        # Also depends on:
        # * salt / ionic strength
        # * length of d1, d2 (and the full length of their strands)
        # But for now, just return default rate constant of 1e5 /s/M
        # Another question: What about intra-complex hybridization?
        #   Should this be accounted for here at the hybridization rate constant,
        #   or at the propensity constant?
        return 1e5


    def dehybridization_rate_constant(self, d1, d2):
        """
        Calculate dehybridization rate constant:
            K = exp(-ΔG°/RT) = exp(-(ΔH°-TΔS°)/RT) = exp(ΔS°/R-ΔH°/R/T) = exp(ΔS°-ΔH°/T) # ΔH°, ΔS° in units of R, R/K
            k_off = k_on/K,
                  = k_on * exp(+ΔG°/RT)
        If T < Tm at molar concentrations (ΔG° < 0), hybridization rate, k_on, is approximately constant at 1e5-1e6.
        If T > Tm at molar concentrations (ΔG° > 0), hybridization rate goes down.
        In practice, we don't have to worry too much about.
        """
        T = self.temperature
        # dH, dS, dG at standard conditions
        dH, dS = self.dHdS_from_state_cache(d1, d2) # Units of R, R/K
        #dG = dH - T*dS # No need to perform extra operation, just use dH and dS:

        # Consideration: Do we really need to have a state-specific energy cache?
        #   a) Yes: Because we might have bending energies, etc which have complex state dependence
        #   b) No:  It should be pleanty fast to add additional energy contributions on the fly.

        ## Add additional energy contributions:
        ## - Stacking energies
        ## - Other? Bending?

        stack_dH, stack_dS = self.stacking_energy(d1, d2)
        # dH += stack_dH
        # dS += stack_dS

        K = exp(dS - dH/T)
        k_on = 1e5  # Is only valid for if K > 1
        k_off = k_on/K
        print("Hybridization energy for domain %s with sequence %s" % (d1, d1.sequence))
        print("  dS=%0.03g, dH=%.03g, T=%s, dS-dH/T=%.03g K=%0.02g, k_off=%0.02g" % (dS, dH, T, dS - dH/T, K, k_off))
        return k_off

    def stacking_energy(self, d1, d2):
        """ Calculate stacking energy for d1 and d2. """
        avg_stack_dH, avg_stack_dS = -4000, -10
        # TODO: Add base-specific stacking interactions
        stack_dH, stack_dS = 0, 0
        for d in (d1, d2):
            if d.end5p.stack_partner:
                stack_dH += avg_stack_dH
                stack_dS += avg_stack_dS
            if d.end3p.stack_partner:
                stack_dH += avg_stack_dH
                stack_dS += avg_stack_dS
        return stack_dH, stack_dS

    def dHdS_from_state_cache(self, d1, d2):
        """
        Cache hybridization energy
        Returns
            dH, dS
        where dH is in units of gas constant R, and dS in units of R/K
        """
        doms_state_hash = frozenset(d.domain_state_fingerprint() for d in (d1, d2))
        if doms_state_hash not in self._statedependent_dH_dS:
            if (d1.name, d2.name) not in self.domain_dHdS:
                dH, dS = hybridization_dH_dS(d1.sequence, c_seq=d2.sequence[::-1])
                print("Energy (dH, dS) for %s hybridizing to %s: %0.4g kcal/mol, %0.4g cal/mol/K" % (d1, d2, dH, dS))
                print("     ", d1.sequence, "\n     ", d2.sequence[::-1])
                # Convert to units of R, R/K:
                dH, dS = dH*1000/R, dS/R
                # print(" - In units of R: %0.4g R, %0.4g R/K" % (dH, dS))
                self.domain_dHdS[(d1.name, d2.name)] = dH, dS
            else:
                dH, dS = self.domain_dHdS[(d1.name, d2.name)]
            ## Also add for individual domains, just for good measure.
            for d in (d1, d2):
                if d1.name not in self.domain_dHdS:
                    self.domain_dHdS[d.name] = hybridization_dH_dS(d1.sequence)
            # TODO: make adjustments, e.g. due to:
            # * dangling ends
            # * stacking
            # * mechanical bending of the helix
            # * electrostatic repulsion
            self._statedependent_dH_dS[doms_state_hash] = (dH, dS)
        return self._statedependent_dH_dS[doms_state_hash]


    def get_relative_activity(self, d1, d2):
        # Historical domain reaction rate constants
        # [T][{(domain1-specie, complex-state), (domain2-specie, complex-state)}]
        # [T][{(domain1-specie, complex-state), (domain2-specie, complex-state)}]
        # Question: Is it really needed to save all permutations of free and complexed rate constants?
        # At least when the two domains are in separate complexes, it seems excessive.
        # But is is probably nice to have for when the two domain species occupy the same complex
        # (where a detailed calculation actually matter)
        # maybe just index by:
        #  - For intra-complex reactions:
        #    relative_intra-complex-reactivity for d1<->d2 =
        #        [({domain1-specie, domain2-specie}, complex-state-fingerprint)]
        #  - For inter-complex reactions, you extract the relative reactivity of each:
        #     d1_relative_activity = [(domain1-specie, domain1-complex-state-fingerprint)]
        #     d2_relative_activity = [(domain2-specie, domain2-complex-state)]
        #  - If the domain strand is not in a complex, relative activity is set to 1.
        #      [{domain1-specie, complex-statedomain2-specie}, complex-fingerprint, intra-complex-reactivity)]
        if d1.strand.complex is None:
            d1_rel_act = 1
        if d2.strand.complex is None:
            d2_rel_act = 1
        if d1.strand.complex == d2.strand.complex != None:
            rel_act = 1
        return rel_act


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


    def hybridize(self, domain1, domain2):
        """
        Splitting out logic in preparation for Julia implementation.
        returns
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        Where changed_complexes includes new_complexes but not obsolete_complexes.
        Edit: changed_complexes DOES NOT include new_complexes.
        """
        print("SystemMgr.hybridize(%s, %s)" % (domain1, domain2))
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
        print("Adding strand_graph edge (%s, %s, key=%s)" % (strand1, strand2, s_edge_key))
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
                result['case'] = 1
                assert strand1 in c2.strands and strand2 in c1.strands
                # c1.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION) # done below
                c_major = c1
                result['changed_complexes'] = [c_major]
            else:
                ## Case (2): Inter-complex hybridization between two complexes. Merge the two complexs:
                result['case'] = 2
                c_major, c_minor = (c1, c2) if (len(c1.nodes()) >= len(c2.nodes())) else (c2, c1)
                for strand in c_minor.strands:
                    strand.complex = c_major
                c_major.strands |= c_minor.strands
                # c_major.N_strand_changes += c_minor.N_strand_changes # Currently not keeping track of strand changes.
                ## Delete the minor complex:
                c_major.add_nodes_from(c_minor.nodes(data=True))
                c_major.add_edges_from(c_minor.edges(keys=True, data=True))
                c_minor.strands = set()
                c_minor.adj = {}        # For networkx, the graph.adj attribute is the main edge data store
                result['obsolete_complexes'] = [c_minor]
                result['changed_complexes'] = [c_major]
        elif c1:
            ## Case 3a: domain2/strand2 is not in a complex; use c1
            result['case'] = 3
            c1.add_strand(strand2, update_graph=True)
            c_major = c1
            result['changed_complexes'] = [c_major]
        elif c2:
            ## Case 3b: domain1/strand1 is not in a complex; use c2
            result['case'] = 3
            c2.add_strand(strand1, update_graph=True)
            c_major = c2
            result['changed_complexes'] = [c_major]
        else:
            ## Case 4: Neither strands are in existing complex; create new complex
            result['case'] = 4
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

        return result


    def dehybridize(self, domain1, domain2):
        """
        Dehybridize domain2 from domain1.
        Returns:
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        """

        print("SystemMgr.dehybridize(%s, %s)" % (domain1, domain2))
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
        print("Removing strand_graph edge (%s, %s, key=%s)" % (strand1, strand2, s_edge_key))
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
        # print("dom1_cc_size=%s, dom2_cc_size=%s, len(c.strands)=%s" % (dom1_cc_size, dom2_cc_size, len(c.strands)))

        assert len(c.strands) == dom1_cc_size + dom2_cc_size

        if dom2_cc_size > 1 and dom1_cc_size > 1:
            # Case (a) Two smaller complexes - must create a new complex for detached domain:
            # Determine which of the complex fragments is the major and which is the minor:
            ## TODO: I don't really need the domain-graph, the strand-level should suffice.
            ## (although it is a good check to have while debuggin...)
            cc_subgraphs = list(connected_component_subgraphs(c))
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
            print("Case (a) - De-hybridization caused splitting into two complexes:")
            print(" - New complex: %s, nodes = %s" % (c_new, c_new.nodes()))
            print(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
            print(" - graph_minor nodes:", graph_minor.nodes())
            print(" - graph_major nodes:", graph_major.nodes())
            result['case'] = 1
            result['changed_complexes'] = [c]
            result['new_complexes'] = [c_new]
            #changed_complexes.append(c_new)
        elif dom2_cc_size > 1 or dom1_cc_size > 1:
            # Case (b) one complex and one unhybridized strand - no need to do much further
            # Which-ever complex has more than 1 strands is the major complex:
            domain_minor = domain1 if dom1_cc_size == 1 else domain2
            c.remove_strand(domain_minor.strand, update_graph=True)
            print("Case (b) - De-hybridization caused a free strand to split away:")
            print(" - Free strand: %s, nodes = %s" % (domain_minor.strand, domain_minor.strand.nodes()))
            print(" - Old complex: %s, nodes = %s" % (c, c.nodes()))
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
            print("Case (c) - De-hybridization caused complex to split into two free strands:")
            print(" - Free strands 1: %s, nodes1 = %s" % (domain1.strand, domain1.strand.nodes()))
            print(" - Free strands 1: %s, nodes1 = %s" % (domain2.strand, domain2.strand.nodes()))
            print(" - Old complex: %s, nodes = %s" % (c, c.nodes()))

        assert domain1.partner is None
        assert domain2.partner is None
        return result


    def react_and_process(self, domain_pair, is_hybridizing, is_intra="Not implemented"):
        """
        Will select a random pair of domain instances from the domain species pair
        for hybridization reaction, or a random duplex in case of dehybridization reactions.
        Reaction specie consists of:
            ({domspec1, domspec2}, is_hybridizing, is_intracomplex)
        """
        print("react_and_process invoked with args: is_hybridizing=%s, is_intra=%s, domain_pair:" %
              (is_hybridizing, is_intra))
        pprint(domain_pair)

        d1, d2 = tuple(domain_pair)
        if is_hybridizing is None:
            is_hybridizing = d1.partner is None
        if is_hybridizing:
            assert d1.partner is None and d2.partner is None
        else:
            assert d1.partner == d2 and d2.partner == d1


        if is_hybridizing:
            print("Hybridizing domain %s %s and %s %s..." %
                  (repr(d1), d1._specie_state_fingerprint, repr(d2), d2._specie_state_fingerprint))
            hyb_result = self.hybridize(d1, d2)
            print("Completed hybridization of domain %s and %s..." % (d1, d2))
        else:
            print("De-hybridizing domain %s %s and %s %s..." %
                  (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
            hyb_result = self.dehybridize(d1, d2)
            print("Completed de-hybridization of domain %s and %s..." % (d1, d2))


        # 4: Update/re-calculate possible_hybridization_reactions and propensity_functions
        # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
        # - possible_hybridization_reactions
        # - propensity_functions
        # Note: If evaluating whether an object is boolean False, the steps include:
        # - Does it have a __len__ attribute? - Yes? Return bool(len(obj))
        # Whereas evaluating whether "obj is None" is about 10 times faster.
        # Fixed: Excluding domains that have no partners -- these will never undergo hybridization reaction.

        print("hyb_result: (is_hybridizing: %s)" % is_hybridizing)
        pprint(hyb_result)
        if hyb_result['free_strands']:
            print("free strands domains:")
            pprint([s.domains for s in hyb_result['free_strands']])

        changed_domains = []
        if hyb_result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in hyb_result['changed_complexes'] for domain in cmplx.nodes()
                                if domain.partner is None and domain.name in self.domain_pairs]
            print("Changed complexes domains:")
            pprint(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if hyb_result['new_complexes']:
            self.complexes |= set(hyb_result['new_complexes'])
            print("Adding new complexes %s to sysmgr.complexes:" % hyb_result['new_complexes'])
            pprint(self.complexes)
            new_cmplx_domains = [domain for cmplx in hyb_result['new_complexes'] for domain in cmplx.nodes()
                                 if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += new_cmplx_domains
            print("New complexes domains:")
            pprint(new_cmplx_domains)
        if hyb_result['free_strands']:
            free_st_domains = [domain for strand in hyb_result['free_strands'] for domain in strand.domains
                               if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += free_st_domains
        if hyb_result['obsolete_complexes']:
            self.complexes -= set(hyb_result['obsolete_complexes'])
            self.removed_complexes += hyb_result['obsolete_complexes']
            print("Removing obsolete complexes %s from sysmgr.complexes:" % hyb_result['obsolete_complexes'])
            pprint(self.complexes)
            print("sysmgr.removed_complexes:")
            pprint(self.removed_complexes)

        if is_hybridizing:
            # d1, d2 are partnered and not included in changed_domains (which currently excludes hybridized domains)
            changed_domains += [d1, d2]

        if len(changed_domains) != len(set(changed_domains)):
            print("WARNING: changed_domains contains duplicates!! THIS SHOULD NOT HAPPEN!")
            print("changed_domains:")
            pprint(changed_domains)
            print("hyb_result['changed_complexes']:")
            pprint(hyb_result['changed_complexes'])
            print("changed_complexes domains: (using complex.nodes())")
            pprint([domain for cmplx in hyb_result['changed_complexes'] for domain in cmplx.nodes()])
            print("changed_complexes domains: (using complex.strands)")
            pprint([domain for cmplx in hyb_result['changed_complexes']
                    for s in cmplx.strands for domain in s.domains])
            print("free_strands:")
            pprint(hyb_result['free_strands'])
            print("free_strands domains:")
            pprint([d for s in hyb_result['free_strands'] for d in s.domains])

            print("( %s reaction between %s (%s) and %s (%s) )" %
                  ("Hybridization" if is_hybridizing else "De-hybridization", d1, repr(d1), d2, repr(d2)))
            print("-----------------")

        self.update_possible_reactions(changed_domains, d1, d2, is_hybridizing, domain_pair)
        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        if d1.strand.complex:
            d1.strand.complex.reset_state_fingerprint()
        if d2.strand.complex:
            d2.strand.complex.reset_state_fingerprint()

        print("changed_domains _specie_state_fingerprint:")
        pprint([(d, d._specie_state_fingerprint) for d in changed_domains])
        print("changed_domains specie_state_fingerprint():")
        pprint([(d, d.domain_state_fingerprint()) for d in changed_domains])

        assert set(self.propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        # Return the selected domains that were actually hybridized/dehybridized
        return d1, d2
