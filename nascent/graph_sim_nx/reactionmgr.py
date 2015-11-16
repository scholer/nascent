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
# Using namedtuple is quite slow in CPython, but faster than dicts in PyPy.
# (Attribute accessing namedtuple.a is faster than dict['a']; not sure about instantiation.)
ReactionAttrsOld = namedtuple('ReactionAttrsOld', ['is_forming', 'is_intra'])
ReactionAttrs = namedtuple('ReactionAttrs', ['reaction_type', 'is_forming', 'is_intra'])
#import math
from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
# from networkx.algorithms.components import connected_components, connected_component_subgraphs
import numpy as np
import pdb

from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS, energy_tables_in_units_of_R
DNA_NN4_R = energy_tables_in_units_of_R['DNA_NN4']
from .constants import R, N_AVOGADRO # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K:
# from .constants import AVOGADRO_VOLUME_NM3
from .constants import PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION
from .constants import REACTION_NAMES
# from .complex import Complex
from .componentmgr import ComponentMgr
from .nx_utils import draw_graph_and_save
from .debug import printd, pprintd
from .domain import Domain, DomainEnd





class ReactionMgr(ComponentMgr):
    """
    System-manager class to manage the system state and state changes
    (but it does not control *which* state changes are induced).

    The GraphManager super-class provides everything related to structural analysis:
    loops and intra-complex activity calculation, etc.
    """

    def __init__(self, volume, strands, params, domain_pairs=None):
        ComponentMgr.__init__(self, strands=strands, params=params, domain_pairs=None)
        self.params = params
        self.temperature = params.get('temperature', 300)
        self.volume = volume or params.get('volume')
        self.stacking_rate_constant = 1e5
        # Base single-molecule activity of two free single reactants in volume.
        self.specific_bimolecular_activity = 1/self.volume/N_AVOGADRO # × M
        self.include_steric_repulsion = False # True

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
        self.reaction_attrs = {}                    # {d1, d2} => (is_forming, is_intra)
        # TODO: Change to generic (reaction_type, is_forming, is_intra)
        self.reaction_pairs_by_domain = defaultdict(set) # domain => {d1, d2} domain_pair of valid reactions.

        self.stacking_interactions = {}         # {(end5p, end3p), (end5p, end3p)}
        self.possible_stacking_reactions = {}   # {(end5p, end3p), (end5p, end3p)} => c_j
        ## What stacking interactions to keep track of?
        ## end5p:end3p
        ## Only ends of hybridized domains (duplexes)
        ## Do we include ALL possible stacking interactions?
        ## or only to nearby ends?
        ## Strategies:
        ## a) Calculate equivalently to domain hybridization?
        ## b) (For intra-complex stacking): Go over all stacking interactions and
        ##      shuffle stacking state according to partition function.
        ## Will single-stranded neighbours produce steric interferrence? Will double-stranded?
        ## I need stacking state ends finger-prints. perhaps just FE = (FD, "5p"/"3p")
        ## where FD is domain state fingerprint.

        # When we are not grouping reactions by domain state species, then propensity functions are just the same
        # as propensity constants... (X₁==X₂==X₁₂==1)
        self.propensity_functions = self.possible_hybridization_reactions
        self.stacking_propensity_functions = self.possible_stacking_reactions

        if strands is not None:
            self.init_possible_reactions()

        self.invoked_reactions_file = None
        if params.get('save_invoked_reactions_to_file'):
            fn = os.path.join(params.get('working_directory', '.'), "invoked_reactions.py") \
                if params['save_invoked_reactions_to_file'] is True else params['save_invoked_reactions_to_file']
            self.invoked_reactions_file = open(fn, 'w')


    def init_possible_reactions(self):
        """
        Reactions have:
            reaction_spec => propensity_constant c_j, state_change_vector v_j
        However, I do:
            ({domain1, domain2}, is_forming, is_intracomplex) => propensity_constant.
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
        printd("\nInitializing possible reactions (UNGROUPED)...")
        Rxs = self.possible_hybridization_reactions
        for d1 in self.domains:
            # d1_domspec = d1.state_fingerprint()
            if d1.partner is not None:
                is_forming, is_intra = False, True
                # If d1 is hybridized, then there is only one possible reaction channel: dehybridizing
                d2 = d1.partner
                # d2_domspec = d2.state_fingerprint()
                domain_pair = frozenset((d1, d2))
                #reaction_spec = (domspec_pair, is_forming, is_intra)
                Rxs[domain_pair] = self.calculate_c_j(d1, d2, is_forming=is_forming, is_intra=is_intra)
                self.reaction_attrs[domain_pair] = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION,
                                                                 is_forming=is_forming, is_intra=is_intra)
                self.reaction_pairs_by_domain[d1].add(domain_pair)
                self.reaction_pairs_by_domain[d2].add(domain_pair)
                # self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                # self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
            else:
                # find all possible reaction channels for domain d1
                # self.domain_pairs[dname] = [list of complementary domain names]
                # all possible hybridization partners/candidates:
                is_forming = True
                if d1.name not in self.domain_pairs:
                    # No domain partners for this domain
                    continue
                for d2 in (d for cname in self.domain_pairs[d1.name]
                           for d in self.unhybridized_domains_by_name[cname]):
                    assert d2.partner is None
                    # d2_domspec = d2.state_fingerprint()
                    domain_pair = frozenset((d1, d2))
                    is_intra = (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex) or \
                               (d1.strand == d2.strand)  # intra-complex OR intra-strand reaction
                    #reaction_spec = (domspec_pair, is_forming, is_intra)
                    if domain_pair not in Rxs:
                        # R_j = (c_j, v_j) - propensity constant for reaction j
                        Rxs[domain_pair] = self.calculate_c_j(d1, d2, is_forming=True, is_intra=is_intra)
                        self.reaction_attrs[domain_pair] = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION,
                                                                         is_forming=is_forming, is_intra=is_intra)
                        self.reaction_pairs_by_domain[d1].add(domain_pair)
                        self.reaction_pairs_by_domain[d2].add(domain_pair)
                    # self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                    # self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
        print(len(self.possible_hybridization_reactions), "possible hybridization reactions initialized.")


    def update_possible_reactions(self, changed_domains, reacted_pair, reaction_attr, reaction_domspec_pair=None):
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
        d1, d2 = tuple(reacted_pair)
        is_forming = reaction_attr.is_forming
        ## TODO: Consolidate this and init_possible_reactions into a single method.
        printd(("\nupdate_possible_reactions invoked with d1=%s, d2=%s, is_forming=%s, reaction_domspec_pair=%s, "
                "changed_domains:") % (d1, d2, is_forming, reaction_domspec_pair))
        pprintd(changed_domains)

        # n_possible_start = len(self.possible_hybridization_reactions)
        # n_species_start = len(self.domain_state_subspecies)

        # old_d1d2_doms_specs = frozenset((d1._specie_state_fingerprint, d2._specie_state_fingerprint))
        updated_reactions = set()
        old_domspecs = {domain: domain._specie_state_fingerprint
                        for domain in changed_domains}

        if is_forming:
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
        if is_forming:
            obsolete_reactions = [domain_pair for domain_pair in self.possible_hybridization_reactions
                                  if d1 in domain_pair or d2 in domain_pair]
            expected_reactions = self.reaction_pairs_by_domain[d1] | self.reaction_pairs_by_domain[d2]
            # the tracked "expected" can easily include reactions that have been removed.
            # For instance: We add {A#1, a#1} to the set. Then A#1 is hybridized to a#2.
            # We delete reaction_pairs_by_domain[A#1] and [a#2], but a#1 is still there.
            # Thus, reaction_pairs_by_domain is only a "suggested" set of reactions that may or may not still be valid.
            if len(set(obsolete_reactions) - expected_reactions) > 0:
                print("\nlen(set(obsolete_reactions) - expected_reactions) > 0:")
                print("set(obsolete_reactions):")
                pprint(set(obsolete_reactions))
                print("expected_reactions:")
                pprint(expected_reactions)
                print("set(obsolete_reactions)-expected_reactions:")
                pprint(set(obsolete_reactions)-expected_reactions)
                print("expected_reactions-set(obsolete_reactions):")
                pprint(expected_reactions-set(obsolete_reactions))
            else:
                printd("\n obsolete_reactions MATCHES expected_reactions (%s elements)"
                       % len(expected_reactions))

            self.reaction_pairs_by_domain[d1].clear()
            self.reaction_pairs_by_domain[d2].clear()
            for domain_pair in obsolete_reactions:
                del self.possible_hybridization_reactions[domain_pair]
                del self.reaction_attrs[domain_pair]
            # You can keep track of domain hybridization reactions, but then
            # please note that the reaction can have become obsolete by the other domain being hybridized.
            # It is thus only a set of possible reactions that can be deleted.
            # Add de-hybridization reaction for d1, d2: (will be skipped when processing all changed domains)
            domain_pair = frozenset((d1, d2))
            # Calling calculate_c_j or any other can produce a call to domain.state_fingerprint()
            self.possible_hybridization_reactions[domain_pair] = \
                self.calculate_c_j(d1, d2, is_forming=False, is_intra=True)
            self.reaction_attrs[domain_pair] = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION,
                                                             is_forming=False, is_intra=True)
            self.reaction_pairs_by_domain[d1].add(domain_pair)
            self.reaction_pairs_by_domain[d2].add(domain_pair)
        else:
            # d1 and d2 have just been de-hybridized;
            # Find new possible partners for d1 and d2:
            # d1 and d2 are processed together with the other changed domains
            pass


        for domain in changed_domains:
            # IMPORTANT: changed_domains must not contain any
            old_domspec = old_domspecs[domain]
            # For d1, d2, old_domspec yields the new domspec instead!
            # printd("Re-setting and re-calculating state_fingerprint for %s - old is: %s"
            #       % (domain, domain._specie_state_fingerprint))
            domain.state_change_reset()
            new_domspec = domain.state_fingerprint()
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
                                self.calculate_c_j(domain, domain2, is_forming=True, is_intra=is_intra)
                            self.reaction_attrs[domain_pair] = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION,
                                                                             is_forming=True, is_intra=is_intra)
                            self.reaction_pairs_by_domain[domain].add(domain_pair)
                            self.reaction_pairs_by_domain[domain2].add(domain_pair)
                            if domain_pair not in updated_reactions:
                                updated_reactions.add(domain_pair)


    def print_reaction_stats(self):
        """ Print information on reaction pathways. """
        printd("Reaction pathways (ungrouped reactions):")
        #keyfun = lambda kv: sorted([d.instance_name for d in kv[0]])
        # We sort the domains to make the stats easier to read:
        # Inner sorted sorts the domain pair (d1, d2 vs d2, d1), outer sorts the reactions by domain
        for reaction_spec_str, c_j, reaction_spec in sorted((sorted([repr(d) for d in k]), v, k) \
            for k, v in self.possible_hybridization_reactions.items()):
            #domspec1, domspec2 = tuple(reaction_spec[0])
            d1, d2 = tuple(reaction_spec)
            domain1, domain2 = reaction_spec_str
            is_forming = d1.partner is None
            is_intra = (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex) or \
                       (d1.strand == d2.strand)  # intra-complex OR intra-strand reaction
            reaction_attr = self.reaction_attrs[reaction_spec]
            reaction_str = REACTION_NAMES[reaction_attr.is_forming][reaction_attr.reaction_type]
            # Use %10s or %18r for domain str/repr
            printd("%18s %9s (%s) %18s" % (domain1,
                                           "hybridize" if is_forming else "de-hybrid",
                                           "intra" if is_intra else "inter", domain2),
                   (": %0.03e x 1 x 1 = %03e" % (
                       c_j, self.propensity_functions[reaction_spec])
                   ) if is_forming else (": %0.03e   x   1 = %03e" % (
                       c_j, self.propensity_functions[reaction_spec])),
                   "" if reaction_attr.is_forming == is_forming and reaction_attr.is_intra == is_intra
                   else (" - EXPECTED: %s" % reaction_attr)
                  )
        # Inner sorted sorts the domain pair (d1, d2 vs d2, d1), outer sorts the reactions by domain
        for reaction_spec_str, c_j, reaction_spec in sorted((sorted([repr(d) for d in k]), v, k) \
            for k, v in self.possible_stacking_reactions.items()):
            #domspec1, domspec2 = tuple(reaction_spec[0])
            d1, d2 = tuple(reaction_spec)
            domain1, domain2 = reaction_spec_str
            reaction_attr = self.reaction_attrs[reaction_spec]
            is_forming = reaction_attr.is_forming
            is_intra = reaction_attr.is_intra
            #reaction_str = ("" if reaction_attr.is_forming else "de-") + reaction_attr.reaction_type
            reaction_str = REACTION_NAMES[reaction_attr.is_forming][reaction_attr.reaction_type]
            # Use %10s or %18r for domain str/repr
            printd("%44s %9s (%s) %44s" % (domain1, reaction_str,
                                           "intra" if is_intra else "inter", domain2),
                   (": %0.03e x 1 x 1 = %03e" % (
                       c_j, self.stacking_propensity_functions[reaction_spec])
                   ) if is_forming else (": %0.03e   x   1 = %03e" % (
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


    def calculate_c_j(self, d1, d2, is_forming, is_intra=None):
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
        if is_forming:
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
        return c_j #, is_forming


    def steric_activity_factor(self, d):
        """
        """
        # TODO: Implement steric_activity_factor
        return 1


    def intracomplex_activity(self, elem1, elem2, reaction_type=HYBRIDIZATION_INTERACTION):
        """
        Return the activity for hybridization of two domains within a complex.
        TODO: Add steric/electrostatic/surface repulsion (individual activity contribution for d1 and d2,
              so that activity = intercomplex_activity*d1_activity*d2_activity)
        """
        if isinstance(elem1, Domain):
            d1, d2 = elem1, elem2
        else:
            d1, d2 = elem1.domain, elem2.domain
        assert d1.strand.complex == d2.strand.complex != None
        cache_key = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
        cache = self.cache['intracomplex_activity']
        if cache_key in cache:
            return cache[cache_key]
        #activity = d1.strand.complex.intracomplex_activity(d1, d2)
        activity = super().intracomplex_activity(elem1, elem2)
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
        ## TODO: Make hybridization k_on depend on domain length.
        ## This ensures that if we break a domain in two, the equilibrium remais the same:
        ## K1 = k1_on/k1_off = exp(ΔG1°/RT)
        ## K2 = k2_on/k2_off = exp((ΔG2a°+ΔG2a°)/RT)
        ## K2 = k2_on/k2_off * k2_on/k2_off = exp(ΔG2a°/RT) * exp(ΔG2b°/RT)
        ## But it is not an AND combination, it's an OR: We only need either 2A or 2B to be hybridized.
        ##                  /-- k1b on/off ->  (B:b, C, c)  -- k2c on/off --\
        ## (B, C) + (b, c) <                                                 >  (B:b, C:c)
        ##                  \-- k1c on/off ->  (C:c, B, b)  -- k2b on/off --/
        ## That combination gives us:
        ## K2 = ([B_h, C_u] + [B_u, C_h]) / ([B_u]*[C_u]) = k1b_on/k1b_off + k1c_on/k1c_off
        ## K1 == K2
        ## k1_on = k1_off*exp(dS1-dH1/T), k2_on = k2_off*exp(dS2-dH2/T)
        ## Essentially, the off rate for case 2 (unhybridize B:b AND C:c) is the same as
        ## the off rate for case 1 (unhybridize A:a).
        ## However, since we have two reaction pathways (and the same number of each starting domain),
        ## if k_on was the same, we would have twice as much reacting in case (2) in any time interval.
        return 5e4*min((len(d1), len(d2)))      # A 20 nt duplex will have k_on = 1e6 /M/s


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

        stack_dH, stack_dS = self.duplex_stacking_energy(d1, d2)
        # dH += stack_dH
        # dS += stack_dS

        K = exp(dS - dH/T)
        k_on = self.hybridization_rate_constant(d1, d2)  # 1e5  # Is only valid for if K > 1
        k_off = k_on/K
        printd("Hybridization energy for domain %s with sequence %s" % (d1, d1.sequence))
        printd("  dS=%0.03g, dH=%.03g, T=%s, dS-dH/T=%.03g K=%0.02g, k_off=%0.02g" % (dS, dH, T, dS - dH/T, K, k_off))
        return k_off

    def duplex_stacking_energy(self, d1, d2):
        """
        Calculate current stacking energy dH, dS (in units of R, R/K)
        for duplex of d1 and d2.
        """
        # TODO: Cache this.
        # avg_stack_dH, avg_stack_dS = -4000, -10
        # TODO: Check base-specific stacking interactions
        stack_dH, stack_dS = 0, 0
        # for d in (d1, d2):
            # if d.end5p.stack_partner:
            #     stack_dH += avg_stack_dH
            #     stack_dS += avg_stack_dS
            # if d.end3p.stack_partner:
            #     stack_dH += avg_stack_dH
            #     stack_dS += avg_stack_dS
        # We only need one end; all ends that are part of a stacking interaction has a copy of the base stacking string.
        if d1.end3p.stack_partner is not None:
            dH, dS = DNA_NN4_R[d1.end3p.stack_string]
            stack_dH, stack_dS = stack_dH + dH, stack_dS + dS
        if d2.end3p.stack_partner is not None:
            dH, dS = DNA_NN4_R[d2.end3p.stack_string]
            stack_dH, stack_dS = stack_dH + dH, stack_dS + dS
        return stack_dH, stack_dS

    def dHdS_from_state_cache(self, d1, d2):
        """
        Cache hybridization energy
        Returns
            dH, dS
        where dH is in units of gas constant R, and dS in units of R/K
        """
        doms_state_hash = frozenset(d.state_fingerprint() for d in (d1, d2))
        if doms_state_hash not in self._statedependent_dH_dS:
            if (d1.name, d2.name) not in self.domain_dHdS:
                dH, dS = hybridization_dH_dS(d1.sequence, c_seq=d2.sequence[::-1])
                print("Energy (dH, dS) for %s hybridizing to %s: %0.4g kcal/mol, %0.4g cal/mol/K" % (d1, d2, dH, dS))
                print("     ", d1.sequence, "\n     ", d2.sequence[::-1])
                # Convert to units of R, R/K:
                dH, dS = dH*1000/R, dS/R
                # printd(" - In units of R: %0.4g R, %0.4g R/K" % (dH, dS))
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


    def react_and_process(self, domain_pair, reaction_attr):
        """
        Will select a random pair of domain instances from the domain species pair
        for hybridization reaction, or a random duplex in case of dehybridization reactions.
        Reaction specie consists of:
            ({domspec1, domspec2}, is_forming, is_intracomplex)
        """
        print("\nreact_and_process invoked with args: domain_pair = %s, reaction_attr = %s" %
              (domain_pair, reaction_attr))
        printd("domain domspecs/fingerprints (before reaction):", [d.state_fingerprint() for d in domain_pair])
        if self.invoked_reactions_file:
            print("domain_pair = frozenset((domains_by_duid[%s], domains_by_duid[%s]))" %
                  tuple([d.duid for d in domain_pair]),
                  file=self.invoked_reactions_file)
            print("sysmgr.react_and_process(domain_pair, %s)" % (reaction_attr, ),
                  file=self.invoked_reactions_file)

        d1, d2 = tuple(domain_pair)
        is_forming = reaction_attr.is_forming
        is_intra = reaction_attr.is_intra
        #is_forming = reaction_attr.is_forming
        #is_intra = reaction_attr.is_intra
        if is_forming:
            assert d1.partner is None and d2.partner is None
        else:
            assert d1.partner == d2 and d2.partner == d1

        if is_forming:
            printd("Hybridizing domain %s %s and %s %s..." %
                   (repr(d1), d1._specie_state_fingerprint, repr(d2), d2._specie_state_fingerprint))
            result = self.hybridize(d1, d2)
            printd("Completed hybridization of domain %s and %s..." % (d1, d2))
        else:
            printd("De-hybridizing domain %s %s and %s %s..." %
                   (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
            result = self.dehybridize(d1, d2)
            printd("Completed de-hybridization of domain %s and %s..." % (d1, d2))

        # 4: Update/re-calculate possible_hybridization_reactions and propensity_functions
        # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
        # - possible_hybridization_reactions
        # - propensity_functions
        # Note: If evaluating whether an object is boolean False, the steps include:
        # - Does it have a __len__ attribute? - Yes? Return bool(len(obj))
        # Whereas evaluating whether "obj is None" is about 10 times faster.
        # Fixed: Excluding domains that have no partners -- these will never undergo hybridization reaction.

        printd("result: (is_forming: %s)" % is_forming)
        pprintd(result)
        if result['free_strands']:
            printd("free strands domains:")
            pprintd([s.domains for s in result['free_strands']])

        changed_domains = []
        if result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()
                                if domain.partner is None and domain.name in self.domain_pairs]
            printd("Changed complexes domains:")
            pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['new_complexes']:
            self.complexes |= set(result['new_complexes'])
            printd("Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            pprintd(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()
                                 if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += new_cmplx_domains
            printd("New complexes domains:")
            pprintd(new_cmplx_domains)
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains
                               if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += free_st_domains
        if result['obsolete_complexes']:
            self.complexes -= set(result['obsolete_complexes'])
            self.removed_complexes += result['obsolete_complexes']
            printd("Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            pprintd(self.complexes)
            printd("sysmgr.removed_complexes:")
            pprintd(self.removed_complexes)

        if is_forming:
            # d1, d2 are partnered and not included in changed_domains (which currently excludes hybridized domains)
            changed_domains += [d1, d2]

        if len(changed_domains) != len(set(changed_domains)):
            print("WARNING: changed_domains contains duplicates!! THIS SHOULD NOT HAPPEN!")
            print("changed_domains:")
            pprint(changed_domains)
            print("result['changed_complexes']:")
            pprint(result['changed_complexes'])
            print("changed_complexes domains: (using complex.nodes())")
            pprint([domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()])
            print("changed_complexes domains: (using complex.strands)")
            pprint([domain for cmplx in result['changed_complexes']
                    for s in cmplx.strands for domain in s.domains])
            print("free_strands:")
            pprint(result['free_strands'])
            print("free_strands domains:")
            pprint([d for s in result['free_strands'] for d in s.domains])

            print("( %s reaction between %s (%s) and %s (%s) )" %
                  ("Hybridization" if is_forming else "De-hybridization", d1, repr(d1), d2, repr(d2)))
            print("-----------------")


        ## Update reactions for changed domains:
        self.update_possible_reactions(changed_domains, reacted_pair=domain_pair, reaction_attr=reaction_attr)
        # Only provide reacted_pair
        self.update_possible_stacking_reactions(for_domains=changed_domains)


        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        if d1.strand.complex:
            d1.strand.complex.reset_state_fingerprint()
        if d2.strand.complex:
            d2.strand.complex.reset_state_fingerprint()

        printd("changed_domains _specie_state_fingerprint:")
        pprintd([(d, d._specie_state_fingerprint) for d in changed_domains])
        printd("changed_domains specie_state_fingerprint():")
        pprintd([(d, d.state_fingerprint()) for d in changed_domains])

        assert set(self.propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        if self.invoked_reactions_file:
            print("print('- react_and_process complete. (execfile)')", file=self.invoked_reactions_file)
        # Return the selected domains that were actually hybridized/dehybridized
        return (d1, d2), result


    ################################
    ###    STACKING REACTIONS    ###
    ################################

    def update_possible_stacking_reactions(self, for_domains=None, reacted_pair=None,
                                           reaction_attr=None):
        """
        The stacking equivalent to init_possible_reactions.
        If changed_domains is None, consider ALL domains.
        (This makes it equivalent to init_possible_stacking_reactions)
        :reacted_pair: If provided, should be a frozenset of {(h1end3p, h2end5p), (h2end3p, h1end5p)}.
        Only provided reacted_pair when processing a recent stacking/unstacking reaction.
        DO NOT provide reacted_pair when processing a hybridization/dehybridization reaction.
        (The hybridized/dehybridized domains are listed in "for_domains" and processed accordingly.)
        """
        if for_domains is None:
            for_domains = self.domains
        print("Updating possible reactions for domains:")
        pprint(for_domains)
        print("Reacted pair:", reacted_pair)
        print("Reaction attr:", reaction_attr)

        ## NOTE: The reaction could have been either stacking/unstacking OR HYBRIDIZATION/DEHYBRIDIZATION.

        ## First, determine reactions that are obsolete by stacking/un-stacking:
        if reacted_pair is not None:
            # reacted_pair MUST be duplex ends for a recent stacking/unstacking reaction; it should not be
            # given when updating after hybridization/dehybridization reactions.
            if reaction_attr is None:
                reaction_attr = self.reaction_attrs[reacted_pair]
            #duplexend1, duplexend2 = reacted_pair  # for stacking/unstacking reaction
            #duplexend1, duplexend2 = reacted_pair  # for stacking/unstacking reaction
            if reaction_attr.is_forming:
                assert reaction_attr.reaction_type == STACKING_INTERACTION
                # Note: We don't process for if reaction_attr.reaction_type == HYBRIDIZATION_INTERACTION:
                # This case, where we have formed a new duplex, is processed under "for h1end3p in stacked_ends" below.
                obsolete_reactions = [other_pair for other_pair in self.possible_stacking_reactions
                                      if len(reacted_pair & other_pair) > 0]
                                      # set.intersection(other_set) is slower than '&' for small sets:
                for reaction_pair in obsolete_reactions:
                    del self.possible_hybridization_reactions[reaction_pair]
                    del self.reaction_attrs[reaction_pair]
            else:
                # Stacking interaction between duplex ends was broken and is available for new rections.
                # This is treated together with the changed domains below:
                pass

        # Only consider hybridized domains:
        hybridized_domains = [domain for domain in for_domains if domain.partner is not None]
        stacked_ends = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is not None]
        unstacked_ends3p = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is None]
        #unstacked_ends5p = [domain.end5p for domain in hybridized_domains if domain.end5p.stack_partner is None]
        # for domain in hybridized_domains:
        #     if domain.end3p.stack_partner is None:
        updated_reactions = set()
        for h1end3p in stacked_ends:
            h1end5p = h1end3p.stack_partner
            h2end5p = h1end3p.hyb_partner
            h2end3p = h2end5p.stack_partner
            ## Assertions:
            assert h2end3p == h1end5p.hyb_partner
            # Multi-graph NetworkX interface is graph[node1][node2][key] => edge_attrs
            assert STACKING_INTERACTION in self.ends5p3p_graph[h1end5p][h1end3p]
            assert STACKING_INTERACTION in self.ends5p3p_graph[h1end3p][h1end5p] # Should be a given.
            assert STACKING_INTERACTION in self.ends5p3p_graph[h2end5p][h2end3p]
            assert STACKING_INTERACTION in self.ends5p3p_graph[h2end3p][h2end5p]

            stacking_pair = frozenset((h1end3p, h1end5p), (h2end3p, h2end5p))
            if stacking_pair not in self.possible_stacking_reactions:
                self.possible_stacking_reactions[stacking_pair] = \
                    self.calculate_stacking_c_j(h1end3p, h1end5p, h2end3p, h2end5p,
                                                is_forming=False, is_intra=True)
                self.reaction_attrs[stacking_pair] = ReactionAttrs(reaction_type=STACKING_INTERACTION,
                                                                   is_forming=False, is_intra=True)
        for h1end3p in unstacked_ends3p:
            h2end5p = h1end3p.hyb_partner
            h1s1 = h1end3p.domain.strand
            h1c1 = h1s1.complex
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            for h2end3p in unstacked_ends3p:
                if h2end3p is h1end3p:
                    continue
                h1end5p = h2end3p.hyb_partner
                ## Assertions:
                assert all(end is not None for end in [h1end5p, h1end3p, h2end3p, h2end5p])
                if h1end3p in self.ends5p3p_graph[h1end5p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h1end5p][h1end3p]
                if h1end5p in self.ends5p3p_graph[h1end3p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h1end3p][h1end5p]
                if h2end3p in self.ends5p3p_graph[h2end5p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h2end5p][h2end3p]
                if h2end5p in self.ends5p3p_graph[h2end3p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h2end3p][h2end5p]

                stacking_pair = frozenset(((h1end3p, h1end5p), (h2end3p, h2end5p)))
                if stacking_pair not in self.possible_stacking_reactions:
                    h1s2 = h1end5p.domain.strand
                    h1c2 = h1s2.complex
                    is_intra = False
                    if h1s1 is h1s2:
                        is_intra = 'strand'
                    elif h1c1 is not None and h1c2 is not None:
                        if h1c1 is h1c2:
                            is_intra = 'complex'
                        elif h1c1.supercomplex is not None and h1c1.supercomplex is h1c2.supercomplex:
                            is_intra = 'superc'
                    self.possible_stacking_reactions[stacking_pair] = \
                        self.calculate_stacking_c_j(h1end3p, h1end5p, h2end3p, h2end5p,
                                                    is_forming=True, is_intra=is_intra)
                    self.reaction_attrs[stacking_pair] = ReactionAttrs(reaction_type=STACKING_INTERACTION,
                                                                        is_forming=True, is_intra=True)



    def stack_and_process(self, stacking_pair):
        """
        The stacking equivalent to react_and_process.
        Perform a stacking reaction, determine which domains have changed,
        and forward that info to reaction updater method.
        stacking_pair = frozenset((h1end3p, h1end5p), (h2end3p, h2end5p))
        Returns stacking_pair, result
        """
        ## TODO: Consolidate stacking and hybridization
        reaction_attr = self.reaction_attrs[stacking_pair]
        reaction_str = ("" if reaction_attr.is_forming else "de-") + reaction_attr.reaction_type

        (h1end3p, h1end5p), (h2end3p, h2end5p) = tuple(stacking_pair)
        d1 = h1end3p.domain
        d2 = h2end3p.domain
        # Assertions:
        assert h2end5p == h1end3p.hyb_partner
        assert h2end3p == h1end5p.hyb_partner
        if reaction_attr.is_forming:
            # assert d1.partner is None and d2.partner is None
            assert h1end5p.stack_partner is None and  h1end3p.stack_partner is None
            assert h2end3p.stack_partner is None and  h2end5p.stack_partner is None
        else:
            assert h1end5p == h1end3p.stack_partner
            assert h2end3p == h2end5p.stack_partner
            # assert d1.partner == d2 and d2.partner == d1

        printd("%s h1end3p, h1end5p, h2end3p, h2end5p - %s %s and %s %s..." %
               (reaction_str, h1end3p, h1end5p, h2end3p, h2end5p))
        if reaction_attr.is_forming:
            result = self.stack(h1end3p, h1end5p, h2end3p, h2end5p)
        else:
            result = self.unstack(h1end3p, h1end5p, h2end3p, h2end5p)
        printd("Completed %s of h1end3p, h1end5p, h2end3p, h2end5p - %s %s and %s %s" %
               (reaction_str, h1end3p, h1end5p, h2end3p, h2end5p))

        printd("%s result:" % reaction_str)
        pprintd(result)

        changed_domains = []
        if result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()
                                if domain.partner is None and domain.name in self.domain_pairs]
            printd("Changed complexes domains:")
            pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains
                               if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += free_st_domains



        self.update_possible_stacking_reactions(changed_domains, stacking_pair, reaction_attr)
        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        if h1end3p.domain.strand.complex:
            h1end3p.domain.strand.complex.reset_state_fingerprint()
        if h1end3p.domain.strand.complex:
            h1end3p.domain.strand.complex.reset_state_fingerprint()

        return stacking_pair, result


    def unstacking_rate_constant(self, h1end3p, h2end3p, T=None):
        """
        Rate for breaking stacking interaction between duplex ends
        (giving the domain 3p end of involved duplexes).
        Note: stacking_rate_constant is a constant attribute, not a function.
        """
        if T is None:
            T = self.temperature
        stack_dH, stack_dS = DNA_NN4_R[h1end3p.stack_string]
        h1end5p = h2end3p.hyb_partner
        h2end5p = h1end3p.hyb_partner
        assert h1end3p.stack_string == h2end3p.stack_string == h1end5p.stack_string == h2end5p.stack_string
        ## TODO: Remove stack_string equality assertion
        # K = exp(stack_dS - stack_dH/T)
        # k_off = self.stacking_rate_constant/K
        k_off = self.stacking_rate_constant * exp(stack_dH/T - stack_dS)
        return k_off


    def calculate_stacking_c_j(self, h1end3p, h1end5p, h2end3p, h2end5p, is_forming, is_intra,
                               reaction_type=STACKING_INTERACTION, reaction_pair_fingerprint=None):
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
        # TODO: Add caching of stacking c_j
        d1 = h1end3p.domain
        d2 = h2end3p.domain
        assert h1end3p.hyb_partner == h2end5p
        assert h1end3p == h2end5p.hyb_partner
        assert h2end3p.hyb_partner == h1end5p
        assert h2end3p == h1end5p.hyb_partner
        if is_forming:
            # if d1.strand.complex == d2.strand.complex != None:
            # Edit: We might have a single strand interacting with itself
            if is_intra:
                assert d1.strand.complex == d2.strand.complex
                assert d1.strand.complex is not None or d1.strand == d2.strand
                # Intra-complex reaction:
                # This should really be the averaae of dist(h1end3p, h1end5p), dist(h2end3p, h2end5p).
                # However, except for the case where the helices are already backbone connected,
                # we can approximate this by just one distance, dist(h1end3p, h2end3p):
                if h1end3p.pb_downstream == h1end5p and h2end3p.pb_downstream == h2end5p:
                    stochastic_activity = 1
                stochastic_activity = self.intracomplex_activity(h1end3p, h2end3p, reaction_type=STACKING_INTERACTION)
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
            if self.include_steric_repulsion:
                steric_activity_factor = ((1 if d1.strand.complex is None else self.steric_activity_factor(d1)) *
                                          (1 if d2.strand.complex is None else self.steric_activity_factor(d2)))
                stochastic_activity *= steric_activity_factor

            k_on = self.stacking_rate_constant # (h1end3p, h2end3p)    # is fairly constant as long as ΔG° < 0
            c_j = k_on * stochastic_activity
        else:
            # De-hybridization reaction;
            k_off = self.unstacking_rate_constant(h1end3p, h2end3p)  # k_off depends on ΔG°
            c_j = k_off # For uni-molecular reactions, c_j = k_j
        return c_j #, is_forming



    def check_system(self):
        is_good = True

        ### Check self.(un)hybridized_domains_by_name lookup tables: ###

        unhybridized_domains_by_name = defaultdict(set)
        hybridized_domains_by_name = defaultdict(set)
        for name in self.unhybridized_domains_by_name:
            _ = unhybridized_domains_by_name[name]
        for name in self.hybridized_domains_by_name:
            _ = hybridized_domains_by_name[name]
        for name, domains in self.domains_by_name.items():
            for d in domains:
                if d.partner is None:
                    unhybridized_domains_by_name[name].add(d)
                else:
                    hybridized_domains_by_name[name].add(d)
        if self.unhybridized_domains_by_name != unhybridized_domains_by_name:
            print("\nSystem-check: self.unhybridized_domains_by_name != unhybridized_domains_by_name")
            print("self.unhybridized_domains_by_name:")
            pprint(self.unhybridized_domains_by_name)
            print("unhybridized_domains_by_name (real):")
            pprint(unhybridized_domains_by_name)
            is_good = False
        if self.hybridized_domains_by_name != hybridized_domains_by_name:
            print("\nSystem-check: self.hybridized_domains_by_name != hybridized_domains_by_name")
            print("self.hybridized_domains_by_name:")
            pprint(self.hybridized_domains_by_name)
            print("hybridized_domains_by_name (real):")
            pprint(hybridized_domains_by_name)
            is_good = False

        ### Check complexes: ###
        for cmplx in self.complexes:
            # Check complex strands:
            strands_by_name = defaultdict(set)
            for name in cmplx.strands_by_name:
                _ = strands_by_name[name]
            for strand in cmplx.strands:
                strands_by_name[strand.name].add(strand)
            if cmplx.strands_by_name != strands_by_name:
                print("\nSystem-check: Problem with complex %r" % cmplx)
                print("System-check: cmplx.strands_by_name != strands_by_name")
                print("cmplx.strands_by_name:")
                pprint(cmplx.strands_by_name)
                print("strands_by_name (from cmplx.strands):")
                pprint(strands_by_name)
                is_good = False

            # Check complex domains:
            strand_domains = {d for strand in cmplx.strands for d in strand.domains}
            if strand_domains != set(cmplx.nodes()):
                is_good = False
                print("\nSystem-check: Problem with complex %r" % cmplx)
                print("Domains from complex.strands:")
                pprint(strand_domains)
                print("complex.nodes():")
                pprint(cmplx.nodes())


        ### Check self.reaction_pairs_by_domain lookup tables: ###


        ### Check self.possible_reactions: ###


        ### Check self.propensity functions (grouped sysmgr): ###

        if not is_good:
            pdb.set_trace()
