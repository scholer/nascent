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
from collections import defaultdict, namedtuple, Counter, deque
from itertools import zip_longest
# Using namedtuple is quite slow in CPython, but faster than dicts in PyPy.
# (Attribute accessing namedtuple.a is faster than dict['a']; not sure about instantiation.)
#import math
from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
# from networkx.algorithms.components import connected_components, connected_component_subgraphs
# import numpy as np
import pdb
import random

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
ReactionAttrs = namedtuple('ReactionAttrs', ['reaction_type', 'is_forming', 'is_intra'])





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
        # Base single-molecule activity of two free single reactants in volume.
        self.specific_bimolecular_activity = 1/self.volume/N_AVOGADRO # × M
        self.include_steric_repulsion = False # True

        ## Reaction settings/parameters:
        # When a complex is changed, update hybridization reaction for all pairs,
        # including hybridization to other complexes:
        self.reaction_update_pair_against_all = params.get("reaction_update_pair_against_all", True)

        ## Hybridization reaction parameters:
        self.hybridization_rate_constant_per_nt = params.get("hybridization_rate_constant_per_nt", 5e4)
        self.reaction_dehybridization_include_stacking_energy = params.get(
            "reaction_dehybridization_include_stacking_energy", True)

        ## Stacking reaction parameters:
        self.stacking_rate_constant = params.get('stacking_rate_constant', 1e5)
        self.enable_intercomplex_stacking = params.get("enable_intercomplex_stacking", False)


        ## Reaction throttle:  (doesn't work, yet)
        self.reaction_throttle = params.get("reaction_throttle", False)
        self.reaction_throttle_fun = params.get("reaction_throttle_fun")
        self.reaction_throttle_offset = params.get("reaction_throttle_offset", 0)
        self.reaction_throttle_rolling_fraction = params.get("reaction_throttle_rolling_fraction", False)
        # Use relative cache instead of calculating new based on absolute Nric
        self.reaction_throttle_use_cache = params.get("reaction_throttle_use_cache", self.reaction_throttle)
        # Nric = normalized_reaction_invocation_count is calculated as:
        # sysmgr.reaction_invocation_count[reaction_spec]/sum(len(sysmgr.domains_by_name[d.name] for d in (d1, d2)))
        if self.reaction_throttle:
            if self.reaction_throttle_fun is None:
                self.reaction_throttle_fun = lambda Nric: exp(-Nric/10)
            elif isinstance(self.reaction_throttle_fun, str):
                self.reaction_throttle_fun = eval(self.reaction_throttle_fun)
        self.reaction_throttle_cache = defaultdict(lambda : 1.0)


        ## Stats:
        # Every time a complex changes, record the new fingerprint here:
        self.complex_state_encounters = Counter()  # keyed by <c state>, value incremented at every occurrence.
        # How many times has a reaction been fired, by (reaction_pair, reaction_attr)
        self.reaction_invocation_count = Counter() # keyed by <pair>, value incremented for every reaction invocation
        # Keep track of the last N reactions:
        self.reaction_invocation_deque = deque(maxlen=self.params.get("reaction_throttle_rolling_window_size", 40))
        # This is slightly more light-weight than the "complex states at time t" collected by StatsManager.


        ## Symbol nomenclature:
        ## Sᵢ, S₁, S₂ - domain species (domain name). Domains with the same sequence are of the same specie.
        ## Fᵢ, F₁, F₂ - domain state fingerprint - domains of the same specie and in the same strand/complex state.
        ## Dᵢ, D₁, D₂ - domain instances, individual domain molecules.

        ### Caches: ###
        self.cache = {} # defaultdict(dict)
        # Standard enthalpy and entropy of hybridization, in units of R or R/T,
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

        # Activities and c_j caches, indexed by {d₁.F, d₂.F} or {(h1e3p.F, h2e5p.F), (h2e3p.F, h1e5p.F)}
        self.cache['intracomplex_activity'] = {}
        self.cache['stochastic_rate_constant'] = {}
        self.cache['stochastic_rate_constant_attrs'] = defaultdict(list)
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
        # as propensity constants: a_j == X₁ X₂ c_j == c_j, since  X₁==X₂==X₁₂==1 for ungrouped reactions.
        self.propensity_functions = self.possible_hybridization_reactions
        self.stacking_propensity_functions = self.possible_stacking_reactions

        if strands is not None:
            self.init_possible_reactions()

        ## Invoked reactions file, for re-playing the simulation:
        self.invoked_reactions_file = None
        if params.get('save_invoked_reactions_to_file'):
            fn = os.path.join(params.get('working_directory', '.'), "invoked_reactions.py") \
                if params['save_invoked_reactions_to_file'] is True else params['save_invoked_reactions_to_file']
            self.invoked_reactions_file = open(fn, 'w')

        print("Simulation system manager initiated at V=%s with %s strands spanning %s domains." \
              % (self.volume, len(self.strands), len(self.domains)))


    def reset_temperature(self, T):
        """ Set system temperature to T and and reset specified caches. """
        self.temperature = T
        if self.params.get("reaction_throttle_reset_on_temperature_change"):
            self.reaction_invocation_count.clear()
        self.cache['intracomplex_activity'].clear()
        self.cache['stochastic_rate_constant'].clear()


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
        # printd("\nInitializing possible reactions (UNGROUPED)...")
        self.update_possible_hybridization_reactions(changed_domains=None, recalculate_hybridized=True)
        self.update_possible_stacking_reactions(changed_domains=None, recalculate_stacked=True)
        print(len(self.possible_hybridization_reactions), "possible hybridization reactions initialized.")


    def update_possible_hybridization_reactions(self, changed_domains, reacted_pair=None,
                                                reaction_attr=None, reacted_domspec_pair=None,
                                                pair_against_all=True, recalculate_hybridized=None,
                                                recalculate_these=None):
        """
        First, process the reacted_pair: if pair was formed (hybridized), remove each domain
        from other possible reactions involving either of the two domains.

        Second, process changed_domains. If :changed_domains: is None, it defaults to ALL domains.
        We only consider un-hybridized domains, unless :recalculate_hybridized: is set to True.

        Arguments:
        :changed_domains:   A list/set of domains that were affected by a recent reaction.
            If None, re-calculate *all* domains.
        :reacted_pair:      The domains that were just reacted, prompting this re-calculation.
        :reaction_attr:     ReactionAttr namedtuple for the reaction that was just triggered.
        :reacted_domspec_pair: Should be frozenset(d for d in reacted_pair), can be provided if you already have it.

        :pair_against_all: If True, consider/recalculate reactions against all other un-hybridized domains.
            Otherwise, only consider/recalculate reactions on within the same complex.
        :recalculate_hybridized: If True, re-calculate all changed domains, even hybridized domains,
            which are by default not re-calculated.
        :recalculate_these: (list/set/tuple) If given, domains in this set are always re-calculated, regardless of
            whether the domain is hybridized or not.

        Variable nomenclature:
            d1, d2: The domains that were just reacted.
            domain, domain2: Loop variables for re-calculating hybridization reactions for changed domain-domain pairs.
        """
        if recalculate_these is None:
            recalculate_these = set()
        if changed_domains is None:
            changed_domains = self.domains
            if recalculate_hybridized is None:
                recalculate_hybridized = True
        # if reacted_domspec_pair is None:
        #     reacted_domspec_pair = frozenset((d.state_fingerprint() for d in reacted_pair))
        ## TODO: Consolidate this and init_possible_reactions into a single method.
        # printd(("\nupdate_possible_hybridization_reactions invoked with reacted_pair=%s, reaction_attr=%s, reacted_domspec_pair=%s, "
        #         "changed_domains:") % (reacted_pair, reaction_attr, reacted_domspec_pair))
        # pprintd(changed_domains)

        # n_possible_start = len(self.possible_hybridization_reactions)
        # n_species_start = len(self.domain_state_subspecies)

        # old_d1d2_doms_specs = frozenset((d1._specie_state_fingerprint, d2._specie_state_fingerprint))
        updated_reactions = set()
        old_domspecs = {domain: domain._specie_state_fingerprint
                        for domain in changed_domains}

        ## PROCESS THE REACTED PAIR:
        if reacted_pair:
            d1, d2 = tuple(reacted_pair)
            is_forming = reaction_attr.is_forming

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
                    # printd("\n obsolete_reactions MATCHES expected_reactions (%s elements)" % len(expected_reactions))
                    pass

                self.reaction_pairs_by_domain[d1].clear()
                self.reaction_pairs_by_domain[d2].clear()
                for domain_pair in obsolete_reactions:
                    del self.possible_hybridization_reactions[domain_pair]
                    del self.reaction_attrs[domain_pair]
                # You can keep track of domain hybridization reactions, but then
                # please note that the reaction can have become obsolete by the other domain being hybridized.
                # It is thus only a set of possible reactions that can be deleted.
                # Add de-hybridization reaction for d1, d2: (will be skipped when processing all changed domains)
                domain_pair = reacted_pair # frozenset((d1, d2))
                domain_spec_pair = frozenset((d1.state_fingerprint(), d2.state_fingerprint()))
                # Calling calculate_hybridization_c_j or any other can produce a call to domain.state_fingerprint()
                self.reaction_attrs[domain_pair] = ra = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION,
                                                                      is_forming=False, is_intra=True)
                # Use calculate_c_j(elem1, elem2, reaction_attr, reaction_spec_pair) to enable throttling and caching.
                self.possible_hybridization_reactions[domain_pair] = (
                    #self.calculate_hybridization_c_j(d1, d2, is_forming=False, is_intra=True)
                    self.calculate_c_j(d1, d2, reaction_attr=ra, reaction_spec_pair=domain_spec_pair))
                self.reaction_pairs_by_domain[d1].add(domain_pair)
                self.reaction_pairs_by_domain[d2].add(domain_pair)
            else:
                # d1 and d2 have just been de-hybridized;
                # Find new possible partners for d1 and d2:
                # d1 and d2 are processed together with the other changed domains
                pass
            del d1
            del d2

        ## PROCESS CHANGED DOMAINS (all domains in all affected complexes):
        for domain in changed_domains:
            # IMPORTANT: changed_domains must not contain any... what?
            # TODO: Checking old vs new domspec is only needed if reacted_pair is not None
            old_domspec = old_domspecs[domain]
            # For d1, d2, old_domspec yields the new domspec instead!
            # printd("Re-setting and re-calculating state_fingerprint for %s - old is: %s"
            #       % (domain, domain._specie_state_fingerprint))
            domain.state_change_reset()
            new_domspec = domain.state_fingerprint()
            if new_domspec == old_domspec and reacted_pair is not None:
                c = domain.strand.complex
                print("\nWeird: new_domspec == old_domspec for changed domain %s: %s == %s" %
                      (domain, new_domspec, old_domspec))
                print(" - complex._state_fingerprint =", c._state_fingerprint)
                print(" - complex.strands_fingerprint() =", c.strands_fingerprint())
                print(" - complex.hybridization_fingerprint() =", c.hybridization_fingerprint())
                print(" - complex.stacking_fingerprint() =", c.stacking_fingerprint())
                print(" - complex.reset_state_fingerprint() ...")
                print(" - complex.state_fingerprint() =", c.state_fingerprint())

            if domain.partner is None and domain.name in self.domain_pairs:
                ## No partner (but is has potential partners at the specie level specified in domain_pairs)
                ## consider hybridization reactions with all other unhybridized complementary domains.
                ## Or maybe just re-calculate for complementary domains in changed domains?
                # for d2 in [d for cname in self.domain_pairs[d1.name]
                #            for d in self.unhybridized_domains_by_name[cname]]:
                #pdb.set_trace()
                is_forming = True
                for cname in self.domain_pairs[domain.name]:
                    for domain2 in self.unhybridized_domains_by_name[cname]:
                        # Remove old reaction propensity:
                        # TODO: Consider only updating/removing intra-complex hybridization reactions: Inter-complex reactions are only affected if steric/electrostatic repulsion is included.
                        if (pair_against_all or domain2 in changed_domains) and domain is not domain2:
                            domain_pair = frozenset((domain, domain2))
                            if domain_pair in updated_reactions:
                                continue
                            domain_spec_pair = frozenset((domain.state_fingerprint(), domain2.state_fingerprint()))
                            # is_intra: intra-complex OR intra-strand reaction:
                            is_intra = (domain.strand.complex is not None and \
                                        domain.strand.complex == domain2.strand.complex) \
                                        or (domain.strand == domain2.strand)
                            self.reaction_attrs[domain_pair] = ra = ReactionAttrs(
                                reaction_type=HYBRIDIZATION_INTERACTION, is_forming=True, is_intra=is_intra)
                            self.possible_hybridization_reactions[domain_pair] = (
                                #self.calculate_hybridization_c_j(domain, domain2, is_forming=True, is_intra=is_intra)
                                self.calculate_c_j(domain, domain2, reaction_attr=ra,
                                                   reaction_spec_pair=domain_spec_pair))
                            self.reaction_pairs_by_domain[domain].add(domain_pair)
                            self.reaction_pairs_by_domain[domain2].add(domain_pair)
                            #if domain_pair not in updated_reactions:
                            updated_reactions.add(domain_pair)
            elif domain.partner is not None and (recalculate_hybridized or domain in recalculate_these):
                # Currently, we don't re-calculate hybridized domains that are already hybridized (other than d1, d2)
                # - unless explicitly told to via recalculate_hybridized argument.
                # If we have a particular pair of domains that we would like to re-calculate (e.g. because of recent
                # stacking reaction), then these can be specified with recalculate_these
                # d2_domspec = d2.state_fingerprint()
                domain2 = domain.partner
                is_forming, is_intra = False, True
                domain_pair = frozenset((domain, domain2))
                domain_spec_pair = frozenset((domain.state_fingerprint(), domain2.state_fingerprint()))
                if domain_pair in updated_reactions:
                    continue
                #reaction_spec = (domspec_pair, is_forming, is_intra)
                self.reaction_attrs[domain_pair] = ra = ReactionAttrs(
                    reaction_type=HYBRIDIZATION_INTERACTION, is_forming=is_forming, is_intra=is_intra)
                self.possible_hybridization_reactions[domain_pair] = (
                    #self.calculate_hybridization_c_j(d1, d2, is_forming=is_forming, is_intra=is_intra)
                    self.calculate_c_j(domain, domain2, reaction_attr=ra, reaction_spec_pair=domain_spec_pair))
                self.reaction_pairs_by_domain[domain].add(domain_pair)
                self.reaction_pairs_by_domain[domain2].add(domain_pair)
                updated_reactions.add(domain_pair)


    def print_reaction_stats(self):
        """ Print information on reaction pathways. """
        print("Reaction pathways (ungrouped reactions):")
        print("\n".join(self.reaction_stats_strs()))


    def reaction_stats_strs(self, ):
        reaction_stats_strs = []
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
            stat = ("%18s %9s %s %18s" % (domain1,
                                          "hybrid" if is_forming else "de-hyb",
                                          "intra" if is_intra else "inter", domain2),
                    (": %0.03e x 1 x 1 = %03e" % (
                        c_j, self.propensity_functions[reaction_spec])
                    ) if is_forming else (": %0.03e   x   1 = %03e" % (
                        c_j, self.propensity_functions[reaction_spec])),
                    "" if reaction_attr.is_forming == is_forming and reaction_attr.is_intra == is_intra
                    else (" - EXPECTED: %s" % reaction_attr)
                   )
            reaction_stats_strs.append("".join(stat))
        # Inner sorted sorts the domain pair (d1, d2 vs d2, d1), outer sorts the reactions by domain
        for reaction_spec_str, c_j, reaction_spec in sorted((sorted([repr(d) for d in k]), v, k) \
            for k, v in self.possible_stacking_reactions.items()):
            #domspec1, domspec2 = tuple(reaction_spec[0])
            d1, d2 = tuple(reaction_spec)
            domain1, domain2 = reaction_spec_str
            reaction_attr = self.reaction_attrs[reaction_spec]
            is_forming = reaction_attr.is_forming
            #reaction_str = ("" if reaction_attr.is_forming else "de-") + reaction_attr.reaction_type
            reaction_str = REACTION_NAMES[reaction_attr.is_forming][reaction_attr.reaction_type]
            # Use %10s or %18r for domain str/repr
            stat = ("%42s %9s %s %42s" % (domain1, reaction_str[:8], #reaction_attr.is_intra, domain2),
                                          "intra" if reaction_attr.is_intra else "inter", domain2),
                    (": %0.03e x 1 x 1 = %03e" % (c_j, self.stacking_propensity_functions[reaction_spec]))
                    if reaction_attr.is_forming else
                    (": %0.03e   x   1 = %03e" % (c_j, self.stacking_propensity_functions[reaction_spec]))
                   )
            reaction_stats_strs.append("".join(stat))
        return reaction_stats_strs


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


    def calculate_c_j(self, elem1, elem2, reaction_attr, reaction_spec_pair=None):
        """
        General method for calculation of c_j stochastic rate constant.
        Unlike the specialized calculate_hybridization_c_j and calculate_stacking_c_j,
        this method also takes care of caching and throttling.
        Cached values does *not* include throttle factor.
        :elem1:, :elem2: Either two domains (if hybridization reaction) or
                         two two-tuples of (h-end3p, h-end5p) pairs (if stacking interaction).
        :reaction_attr: A ReactionAttr namedtuple with (reaction_type, is_forming, is_intra) values.
        :reaction_spec_pair:    A frozenset of state fingerprints for each instance in elem1, elem2.

        """
        assert reaction_spec_pair != frozenset((elem1, elem2))
        assert reaction_attr.reaction_type != HYBRIDIZATION_INTERACTION or \
            reaction_attr.is_forming or reaction_attr.is_intra
        reaction_type, is_forming, is_intra = reaction_attr
        if reaction_spec_pair is not None and reaction_spec_pair in self.cache['stochastic_rate_constant']:
            # TOOD: Make sure that this includes stacking state!
            # TODO: We actually have two kinds of species fingerprints: One that allows symmetry and one that don't.
            c_j, ra = self.cache['stochastic_rate_constant'][reaction_spec_pair]
            assert ra == reaction_attr
        else:
            if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                if reaction_spec_pair is None:
                    reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
                c_j = self.calculate_hybridization_c_j(elem1, elem2, is_forming=is_forming, is_intra=is_intra,
                                                       reaction_spec_pair=reaction_spec_pair)
            elif reaction_attr.reaction_type is STACKING_INTERACTION:
                if reaction_spec_pair is None:
                    reaction_spec_pair = frozenset(((elem1[0].state_fingerprint(), elem1[0].state_fingerprint()),
                                                    (elem2[0].state_fingerprint(), elem2[0].state_fingerprint())))
                c_j = self.calculate_stacking_c_j(elem1[0], elem1[1], elem2[0], elem2[1],
                                                  is_forming=is_forming, is_intra=is_intra,
                                                  reaction_spec_pair=reaction_spec_pair)
            else:
                raise ValueError("Unknown reaction type %r for reaction_attr %s between %s and %s" %
                                 (reaction_attr.reaction_type, reaction_attr, elem1, elem2))
            self.cache['stochastic_rate_constant'][reaction_spec_pair] = (c_j, reaction_attr)
            # raise NotImplementedError("Not fully implemented yet...")
        # if self.reaction_throttle:
        #     if self.reaction_throttle_use_cache:
        #         throttle_factor = self.reaction_throttle_cache[reaction_spec_pair]
        #     else:
        #         throttle_factor = self.calculate_throttle_factor(elem1, elem2, reaction_type, reaction_spec_pair)
        #         if throttle_factor < 1 or (False or random.random() < max((0.1, throttle_factor))):
        #             print(("Throttling {t:.02e} for {inter_intra} {de}-{reaction_type} reaction: {elem1} <-> {elem2}"
        #                    "\n(spec: {spec} [{ric} invocations])").format(
        #                 t=throttle_factor, inter_intra="intra" if reaction_attr.is_intra else "inter",
        #                 de="  " if reaction_attr.is_forming else "de", reaction_type=reaction_attr.reaction_type,
        #                 elem1=elem1, elem2=elem2, ric=self.reaction_invocation_count[reaction_spec_pair],
        #                 spec=set(reaction_spec_pair)))
        #             # print(" - reaction_spec_pair:", set(reaction_spec_pair))
        #     c_j *= throttle_factor
        # else:
        #     print("NOT throttling reaction:", elem1, elem2, reaction_attr)
        return c_j


    def calculate_throttle_factor(self, elem1, elem2, reaction_type, reaction_spec_pair):
        """
        :reaction_spec_pair: Must be a state-species pair; not instances. Symmetric degeneracy is acceptable.
        """
        assert reaction_spec_pair is not None and reaction_type is not None
        if reaction_type is HYBRIDIZATION_INTERACTION:
            N_dom_copies = sum(len(self.domains_by_name[d.name]) for d in (elem1, elem2))
        elif reaction_type is STACKING_INTERACTION:
            N_dom_copies = sum(len(self.domains_by_name[end.domain.name])
                                   for endpair in (elem1, elem2) for end in endpair)
        Nric = self.reaction_invocation_count[reaction_spec_pair]/N_dom_copies
        # pprint(locals())
        #pprint(self.reaction_invocation_count)
        #print(Nric)
        #pdb.set_trace()
        if Nric < self.reaction_throttle_offset:
            return 1
        if self.reaction_throttle_rolling_fraction:
            frac = self.reaction_invocation_deque.count(reaction_spec_pair)/len(self.reaction_invocation_deque)
            if frac < self.reaction_throttle_rolling_fraction:
                return 1
        return self.reaction_throttle_fun(Nric-self.reaction_throttle_offset)


    def calculate_hybridization_c_j(self, d1, d2, is_forming, is_intra, reaction_type=HYBRIDIZATION_INTERACTION,
                                    reaction_spec_pair=None):
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
        ## Done: Make c_j a cached value [done through calculate_c_j generic method]

        assert is_forming or is_intra   # If we are not forming, we are de-hybridizing; is_intra should be True

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
                #
                # steric_activity_factor = np.prod([1 if d.strand.complex is None else self.steric_activity_factor(d)
                #                                   for d in (d1, d2)])
                steric_activity_factor = ((1 if d1.strand.complex is None else self.steric_activity_factor(d1)) *
                                          (1 if d2.strand.complex is None else self.steric_activity_factor(d2)))
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
            # TODO: Ensure that k_off depends on end stacking interactions!
            k_off = self.dehybridization_rate_constant(d1, d2)
            # For uni-molecular reactions, c_j = k_j
            c_j = k_off
        return c_j #, is_forming


    def steric_activity_factor(self, d):
        """
        """
        # TODO: Implement steric_activity_factor
        return 1


    def intracomplex_activity(self, elem1, elem2, reaction_type=HYBRIDIZATION_INTERACTION, reaction_spec_pair=None):
        """
        Return the activity for hybridization of two domains within a complex.
        :elem1:, :elem2: are either two domains, or two pairs of duplex domain ends:
        For stacking activity, use
            :elem1: = (h1end3p, h2end5p)  and  :elem2: = (h2end3p, h1end5p)   [new stacking pair grouping?]

        TODO: Add steric/electrostatic/surface repulsion (individual activity contribution for d1 and d2,
              so that activity = intercomplex_activity*d1_activity*d2_activity)
        """
        if isinstance(elem1, tuple):
            #assert reaction_type == STACKING_INTERACTION
            if reaction_spec_pair is None:
                reaction_spec_pair = frozenset(((elem1[0].state_fingerprint(), elem1[1].state_fingerprint()),
                                                (elem2[0].state_fingerprint(), elem2[1].state_fingerprint())))
            else: # TODO: Remove assertion when done debugging.
                assert reaction_spec_pair == frozenset(((elem1[0].state_fingerprint(), elem1[1].state_fingerprint()),
                                                        (elem2[0].state_fingerprint(), elem2[1].state_fingerprint())))
            d1, d2 = elem1[0].domain, elem2[0].domain
        else:
            if reaction_spec_pair is None:
                reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
            else: # TODO: Remove assertion when done debugging.
                assert reaction_spec_pair == frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
            if isinstance(elem1, Domain):
                d1, d2 = elem1, elem2
                assert reaction_type == HYBRIDIZATION_INTERACTION
            elif isinstance(elem1, DomainEnd):
                d1, d2 = elem1.domain, elem2.domain
        assert d1.strand.complex == d2.strand.complex != None
        # cache = self.cache['intracomplex_activity']
        if reaction_spec_pair in self.cache['intracomplex_activity']:
            return self.cache['intracomplex_activity'][reaction_spec_pair]
        #activity = d1.strand.complex.intracomplex_activity(d1, d2)
        if isinstance(elem1, tuple): # stacking domain-ends
            # activity = (super().intracomplex_activity(elem1[0], elem2[0]) +
            #             super().intracomplex_activity(elem1[1], elem2[1]))/2
            activity = super().intracomplex_activity(elem1, elem2)   # h1end3p, h2end3p
        else: # hybridizing domains
            activity = super().intracomplex_activity(elem1, elem2)
        self.cache['intracomplex_activity'][reaction_spec_pair] = activity
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
        return self.hybridization_rate_constant_per_nt*min((len(d1), len(d2)))      # A 20 nt duplex will have k_on = 1e6 /M/s


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

        if self.reaction_dehybridization_include_stacking_energy:
            stack_dH, stack_dS = self.duplex_stacking_energy(d1, d2)
            if any((stack_dH, stack_dS)):
                print("\nAdding stack_dH, stack_dS = %s, %s to dehybridization energy (%s and %s)" %
                      (stack_dH, stack_dS, d1, d2))
            dH += stack_dH
            dS += stack_dS

        K = exp(dS - dH/T)
        k_on = self.hybridization_rate_constant(d1, d2)  # 1e5  # Is only valid for if K > 1
        k_off = k_on/K
        # printd("Hybridization energy for domain %s with sequence %s" % (d1, d1.sequence))
        # printd("  dS=%0.03g, dH=%.03g, T=%s, dS-dH/T=%.03g K=%0.02g, k_off=%0.02g" % (dS, dH, T, dS - dH/T, K, k_off))
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


    def react_and_process(self, domain_pair, reaction_attr, reaction_spec_pair=None):
        """
        Will select a random pair of domain instances from the domain species pair
        for hybridization reaction, or a random duplex in case of dehybridization reactions.
        Reaction specie consists of:
            ({domspec1, domspec2}, is_forming, is_intracomplex)
        """
        ## TODO: Merge this with stack_and_process
        # printd("\nreact_and_process invoked with args: domain_pair = %s, reaction_attr = %s" % (domain_pair, reaction_attr))
        # printd("domain domspecs/fingerprints (before reaction):", [d.state_fingerprint() for d in domain_pair])
        if self.invoked_reactions_file:
            print("domain_pair = frozenset((domains_by_duid[%s], domains_by_duid[%s]))" %
                  tuple([d.duid for d in domain_pair]),
                  file=self.invoked_reactions_file)
            print("sysmgr.react_and_process(domain_pair, %s)" % (reaction_attr, ),
                  file=self.invoked_reactions_file)


        d1, d2 = tuple(domain_pair)
        if reaction_spec_pair is None:
            ## TODO: Add a global reaction_spec_pair by reaction_pair cache (making a frozenset is a bit slow).
            reaction_spec_pair = frozenset((d1.state_fingerprint(), d2.state_fingerprint()))
        else: ## TODO: Remove assertion when done debugging.
            assert reaction_spec_pair == frozenset((d1.state_fingerprint(), d2.state_fingerprint()))

        self.reaction_invocation_count[(reaction_spec_pair, reaction_attr)] += 1
        if self.reaction_throttle_use_cache:
            self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)] = \
                0.9995 * self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)]


        assert reaction_attr.reaction_type == HYBRIDIZATION_INTERACTION
        if reaction_attr.is_intra:
            assert d1.strand.complex == d2.strand.complex
            if d1.strand.complex is None:
                assert d1.strand == d2.strand
        else:
            assert (d1.strand.complex != d2.strand.complex) or (d1.strand.complex is None or d2.strand.complex is None)

        if reaction_attr.is_forming:
            # printd("Hybridizing domain %s %s and %s %s..." % (repr(d1), d1._specie_state_fingerprint, repr(d2), d2._specie_state_fingerprint))
            assert d1.partner is None and d2.partner is None
            result = self.hybridize(d1, d2)
            # printd("Completed hybridization of domain %s and %s..." % (d1, d2))
        else:
            # printd("De-hybridizing domain %s %s and %s %s..." % (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
            assert d1.partner == d2 and d2.partner == d1
            result = self.dehybridize(d1, d2)
            # printd("Completed de-hybridization of domain %s and %s..." % (d1, d2))

        # 4: Update/re-calculate possible_hybridization_reactions and propensity_functions
        # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
        # - possible_hybridization_reactions
        # - propensity_functions
        # Note: If evaluating whether an object is boolean False, the steps include:
        # - Does it have a __len__ attribute? - Yes? Return bool(len(obj))
        # Whereas evaluating whether "obj is None" is about 10 times faster.
        # Fixed: Excluding domains that have no partners -- these will never undergo hybridization reaction.

        # printd("result: (is_forming: %s)" % is_forming)
        # pprintd(result)
        # if result['free_strands']:
            # printd("free strands domains:")
            # pprintd([s.domains for s in result['free_strands']])

        changed_domains = []
        if result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()
                                # if domain.partner is None and domain.name in self.domain_pairs
                                # edit: we need to have all, not just unhybridized - for stacking
                               ]
            # printd("Changed complexes domains:")
            # pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['new_complexes']:
            self.complexes |= set(result['new_complexes'])
            # printd("Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            # pprintd(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()
                                 # if domain.partner is None and domain.name in self.domain_pairs
                                 # edit: we need to have all, not just unhybridized - for stacking
                                ]
            changed_domains += new_cmplx_domains
            # printd("New complexes domains:")
            # pprintd(new_cmplx_domains)
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains
                               #if domain.partner is None and domain.name in self.domain_pairs
                               # edit: we need to have all, not just unhybridized - for stacking
                              ]
            changed_domains += free_st_domains
        if result['obsolete_complexes']:
            self.complexes -= set(result['obsolete_complexes'])
            self.removed_complexes += result['obsolete_complexes']
            # printd("Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            # pprintd(self.complexes)
            # printd("sysmgr.removed_complexes:")
            # pprintd(self.removed_complexes)

        # if reaction_attr.is_forming:
        #     # d1, d2 are partnered and not included in changed_domains (which currently excludes hybridized domains)
        #     # edit: we now include *all* domains in changed_domains, hybridized and un-hybridized.
        #     changed_domains += [d1, d2]

        if len(changed_domains) != len(set(changed_domains)):
            print("\nWARNING: changed_domains contains duplicates!! THIS SHOULD NOT HAPPEN!\n")
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
                  ("Hybridization" if reaction_attr.is_forming else "De-hybridization", d1, repr(d1), d2, repr(d2)))
            print("-----------------")


        ## Update reactions for changed domains:
        self.update_possible_hybridization_reactions(changed_domains, reacted_pair=domain_pair, reaction_attr=reaction_attr)
        # Only provide reacted_pairs if a reaction really has occured:
        # HEY, IF UNHYBRIDIZE() SETS d1.partner = None, then that currently also sets d1.end3p.hyb_partner = None, etc!
        # Then, if update_possible_stacking_reactions tries to "assert h2end5p == h1end3p.hyb_partner", it will fail.
        # Actually, it won't fail until it tries to stack the ends that have just been de-hybridized.
        # if 'unstacking_results' in result:
        #     reacted_stacking_pairs = result['unstacking_results']
        #     # unstacking_results is a dict of: {(h1end3p, h2end5p, h2end3p, h1end5p): result), ....}
        # else:
        #     reacted_stacking_pairs = None
        # Using dehybridized_duplex_ends is much better than trying to pretend that a (stacking) reaction has occoured:
        dehybridized_ends = [d1.end5p, d1.end3p, d2.end5p, d2.end3p] if not reaction_attr.is_forming else None
        self.update_possible_stacking_reactions(changed_domains=changed_domains,
                                                dehybridized_ends=dehybridized_ends)
        # Note: We still use 'unstacking_results' for forwarding these side-effect unstacking reactions to dispatcher.

        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        if d1.strand.complex:
            d1.strand.complex.reset_state_fingerprint()
        if d2.strand.complex not in (d1.strand.complex, None):
            d2.strand.complex.reset_state_fingerprint()

        # printd("changed_domains _specie_state_fingerprint:")
        # pprintd([(d, d._specie_state_fingerprint) for d in changed_domains])
        # printd("changed_domains specie_state_fingerprint():")
        # pprintd([(d, d.state_fingerprint()) for d in changed_domains])

        assert set(self.propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        if self.invoked_reactions_file:
            print("print('- react_and_process complete. (execfile)')", file=self.invoked_reactions_file)
        # Return the selected domains that were actually hybridized/dehybridized
        return (d1, d2), result


    ################################
    ###    STACKING REACTIONS    ###
    ################################

    def update_possible_stacking_reactions(self, changed_domains=None, reacted_pairs=None, reaction_attrs=None,
                                           dehybridized_ends=None, recalculate_stacked=False,
                                           reaction_endspec_pair=None):
        """
        The stacking equivalent to update_possible_hybridization_reactions.
        If changed_domains is None, consider ALL domains.
        (This makes it equivalent to init_possible_stacking_reactions)
        :reacted_pair: If provided, should be a frozenset of {(h1end3p, h2end5p), (h2end3p, h1end5p)}.
        Only provide reacted_pair when processing a recent stacking/unstacking reaction.
        DO NOT provide reacted_pair when processing a hybridization/dehybridization reaction.
        (The hybridized/dehybridized domains are listed in "changed_domains" and processed accordingly.)

        :dehybridized_duplex_ends: is a list of ends [d1end3p, d2end5p, ...]
            for domains that have just been dehybridized (and thus cannot undergo stacking reactions).
        """
        # printd("Updating possible stacking reactions for domains: (None=All)")
        # pprintd(changed_domains)
        # printd("Reacted pairs:", reacted_pairs)
        # printd("Reaction attrs:", reaction_attrs)
        if changed_domains is None:
            # Update possible stacking reactions for *all* domains:
            changed_domains = list(self.domains)

        ## NOTE: The reaction could have been either stacking/unstacking OR HYBRIDIZATION/DEHYBRIDIZATION.

        ## First, determine reactions that are obsolete by stacking/un-stacking:
        if reacted_pairs is not None:
            # reacted_pair MUST be duplex ends for a recent stacking/unstacking reaction; it should not be
            # given when updating after hybridization/dehybridization reactions, except if a de-hybridization
            # resulted in broken stacking interactions:
            if reaction_attrs is None:
                reaction_attrs = [self.reaction_attrs[reacted_pair] for reacted_pair in reacted_pairs]
            newly_formed_stacks = [(reacted_pair, reaction_attr)
                                   for reacted_pair, reaction_attr in zip_longest(reacted_pairs, reaction_attrs)
                                   if reaction_attr.is_forming]
            # We are only considering newly formed stacking pairs; newly broken stacking pairs
            # should be in changed_domains and treated together with the other changed domains:
            # Stacking interaction between duplex ends was broken and is available for new rections.
            # This is treated together with the changed domains below:
            for reacted_pair, reaction_attr in newly_formed_stacks:
                #duplexend1, duplexend2 = reacted_pair  # for stacking/unstacking reaction
                #duplexend1, duplexend2 = reacted_pair  # for stacking/unstacking reaction
                assert reaction_attr.reaction_type == STACKING_INTERACTION
                # Note: We don't process for if reaction_attr.reaction_type == HYBRIDIZATION_INTERACTION:
                # This case, where we have formed a new duplex, is processed under "for h1end3p in stacked_ends" below.
                obsolete_reactions = [other_pair for other_pair in self.possible_stacking_reactions
                                      if len(reacted_pair & other_pair) > 0]
                # "if len(reacted_pair & other_pair) > 0" checks if the two stacking_pairs have any ends in common.
                # set.intersection(other_set) is slower than '&' for small sets:
                for reaction_pair in obsolete_reactions:
                    # printd("Removing from possible_stacking_reactions obsolete stacking_pair:", reaction_pair)
                    del self.possible_stacking_reactions[reaction_pair]
                    del self.reaction_attrs[reaction_pair]

            ## TODO: Check to make sure that all reacted_pairs domains are in changed_domains
            # {(h1end3p, h1end5p), (h2end3p, h2nd5p)}
            for end in (end for reacted_pair in reacted_pairs for tup in reacted_pair for end in tup):
                if not end.domain in changed_domains:
                    # printd(" - update_possible_stacking_reactions: Adding %s end's domain to changed_domains." % end)
                    #changed_domains.add(end.domain)
                    changed_domains.append(end.domain)

        if dehybridized_ends is not None:
            obsolete_reactions = [stacking_pair for stacking_pair in self.possible_stacking_reactions
                                  if any(end in ends_tup for end in dehybridized_ends
                                         for ends_tup in stacking_pair)]
            # printd("Removing reactions obsolete by dehybridized_duplex_ends:", dehybridized_ends)
            # pprintd(obsolete_reactions)
            for stacking_pair in obsolete_reactions:
                del self.possible_stacking_reactions[stacking_pair]


        # Only consider hybridized domains:
        hybridized_domains = [domain for domain in changed_domains if domain.partner is not None]
        stacked_ends = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is not None]
        # printd("stacked_ends:")
        # pprintd(stacked_ends)
        unstacked_ends3p = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is None]
        if self.enable_intercomplex_stacking:
            # Match against *all* unstacked domains (not just intra-complex)
            all_unstacked_ends3p = [domain.end3p for domain in self.domains
                                    if domain.partner is not None and domain.end3p.stack_partner is None]
        else:
            # Only consider stacking within changed complexes:
            all_unstacked_ends3p = unstacked_ends3p
        #unstacked_ends5p = [domain.end5p for domain in hybridized_domains if domain.end5p.stack_partner is None]
        # for domain in hybridized_domains:
        #     if domain.end3p.stack_partner is None:
        updated_reactions = set()
        # stacked_ends is a list of stacked DomainEnd3p on hybridized domains in changed_domains:
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

            if h1end3p.stack_partner == h1end3p.pb_downstream and h2end3p.stack_partner == h2end3p.pb_downstream:
                # Performance hack for multiple domains on the same strand on the same duplex. Don't unstack.
                continue

            stacking_pair = frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
            # if stacking_pair in updated_reactions:
            #     continue
            # For currently-stacked ends, we only need to update if its not in possible_stacking_reactions;
            # Once a stacking_pair is in possible_stacking_reactions, it doesn't change except when unstacked.
            # What if stacking_pair is already in possible_stacking_reactions, but the complex has changed
            # enough to make the first c_j obsolete? Uhm, see above...
            # And for the case looking at forming stacking reaction: the complex has changed,
            # so we *MUST* to re-calculate c_j.
            # stackspec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
            #                             (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
            # if stacking_pair in self.reaction_attrs:
            #     if stackspec_pair == self.reaction_attrs[stacking_pair].stackspec_pair:
            #         # No change since last:
            #         continue
            if stacking_pair not in self.possible_stacking_reactions or \
                (recalculate_stacked and stacking_pair not in updated_reactions):
                ## TODO: Reaction throttle should probably be inserted here. Edit: Use calculate_c_j() to add throttle.
                reaction_attr = ReactionAttrs(reaction_type=STACKING_INTERACTION,
                                              is_forming=False, is_intra=True) #, reactionspec_pair=stackspec_pair)
                self.possible_stacking_reactions[stacking_pair] = \
                    self.calculate_stacking_c_j(h1end3p, h2end5p, h2end3p, h1end5p,
                                                is_forming=False, is_intra=True)
                self.reaction_attrs[stacking_pair] = reaction_attr
                updated_reactions.add(stacking_pair)
        # unstacked_ends is a list of unstacked DomainEnd3p on hybridized domains in changed_domains:
        for h1end3p in unstacked_ends3p:
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            #
            for h2end3p in all_unstacked_ends3p:
                if h2end3p is h1end3p:
                    continue # Ends cannot stack to them self

                # unpack variables:
                h2end5p = h1end3p.hyb_partner
                assert h1end3p.domain.partner == h2end5p.domain
                h1s1 = h1end3p.domain.strand
                h1c1 = h1s1.complex
                h1end5p = h2end3p.hyb_partner
                h1s2 = h1end5p.domain.strand
                h1c2 = h1s2.complex
                if not self.enable_intercomplex_stacking and (h1c1 != h1c2 if h1c1 is not None else h1s1 != h1s2):
                    # Interaction between different complexes (or strands):
                    continue
                if h2end3p.hyb_partner == h1end3p.stacked_upstream() and not self.enable_helix_bending:
                    # Don't allow duplex to hybridize to itself (opposite ends).
                    #                 h1end5p    h1end3p
                    # This end --> 5'--------------------3' <-- cannot hybridize to this.
                    # Helix 2   -> 3'--------------------5' <-
                    #                h2end3p    h2end5p
                    # printd("update_possible_stacking_reactions: Ignoring stacking between %r and %r" % (h1end3p, h2end3p))
                    continue
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

                stacking_pair = frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
                if stacking_pair not in updated_reactions: # self.possible_stacking_reactions:
                    is_intra = False
                    if h1c1 is not None and h1c1 is h1c2:
                        is_intra = 'complex'
                    elif h1s1 is h1s2:
                        is_intra = 'strand'
                    # print("  %s -: :- %s \n  %s -: :- %s  " % (h1end3p, h1end5p, h2end5p, h2end3p))
                    c_j = self.calculate_stacking_c_j(h1end3p, h2end5p, h2end3p, h1end5p,
                                                      is_forming=True, is_intra=is_intra)
                    self.possible_stacking_reactions[stacking_pair] = c_j
                    self.reaction_attrs[stacking_pair] = ReactionAttrs(reaction_type=STACKING_INTERACTION,
                                                                       is_forming=True, is_intra=is_intra)
                    updated_reactions.add(stacking_pair)
        # pprint(self.reaction_attrs)
        # pprint(locals())
        # pdb.set_trace()



    def stack_and_process(self, stacking_pair, reaction_spec_pair=None):
        """
        The stacking equivalent to react_and_process.
        Perform a stacking reaction, determine which domains have changed,
        and forward that info to reaction updater method.
        stacking_pair = frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
        Ends annotation:
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p
        Returns stacking_pair, result
        """
        ## TODO: Consolidate stacking and hybridization

        # printd("stack_and_process invoked with stacking_pair:")
        # pprintd(stacking_pair)
        (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(stacking_pair)
        if self.invoked_reactions_file:
            print(("stacking_pair = frozenset(("
                   "(getattr(domains_by_duid[%s], 'end%s'), getattr(domains_by_duid[%s], 'end%s'))"
                   "(getattr(domains_by_duid[%s], 'end%s'), getattr(domains_by_duid[%s], 'end%s'))))") %
                  (h1end3p.domain.duid, h1end3p.end,
                   h1end5p.domain.duid, h1end5p.end,
                   h2end3p.domain.duid, h2end3p.end,
                   h2end5p.domain.duid, h2end5p.end),
                  file=self.invoked_reactions_file)
            print("sysmgr.stack_and_process(stacking_pair)", file=self.invoked_reactions_file)

        reaction_attr = self.reaction_attrs[stacking_pair]
        if reaction_spec_pair is None:
            reaction_spec_pair = frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
        else:
            assert reaction_spec_pair == frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
        reaction_str = ("" if reaction_attr.is_forming else "de-") + reaction_attr.reaction_type

        self.reaction_invocation_count[(reaction_spec_pair, reaction_attr)] += 1
        # if self.reaction_throttle_use_cache:
        #     self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)] = \
        #         0.95 * self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)]

        # d1 = h1end3p.domain
        # d2 = h2end3p.domain
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
        # printd("stack_and_process: %s h1end3p, h2end5p, h2end3p, h1end5p - %s %s and %s %s..." % (reaction_str, h1end3p, h2end5p, h2end3p, h1end5p))
        assert reaction_attr.reaction_type == STACKING_INTERACTION
        if reaction_attr.is_intra:
            assert len(set(e.domain.strand.complex for e in (h1end3p, h2end5p, h2end3p, h1end5p))) == 1
            if h1end3p.domain.strand.complex is None:
                assert len(set(e.domain.strand for e in (h1end3p, h2end5p, h2end3p, h1end5p))) == 1

        if reaction_attr.is_forming:
            assert all(e.stack_partner is None for e in (h1end3p, h2end5p, h2end3p, h1end5p))
            result = self.stack(h1end3p, h2end5p, h2end3p, h1end5p)
        else:
            assert all(e.stack_partner is not None for e in (h1end3p, h2end5p, h2end3p, h1end5p))
            result = self.unstack(h1end3p, h2end5p, h2end3p, h1end5p)
        # printd("stack_and_process: Completed %s of h1end3p, h2end5p, h2end3p, h1end5p - %s %s and %s %s" % (reaction_str, h1end3p, h2end5p, h2end3p, h1end5p))

        # printd("stack_and_process: %s result:" % reaction_str)
        # pprintd(result)

        changed_domains = []
        if result['changed_complexes']:
            # Edit: For stacking, we want duplexed domains, regardless of domain_pairs!!
            # No longer filtering with "if domain.partner is None and domain.name in self.domain_pairs"
            # Instead, do appropriate filtering in update_possible_hybridization_reactions.
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()]
            # printd("stack_and_process: Changed complexes domains:")
            # pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains]
            changed_domains += free_st_domains
        if result['new_complexes']:
            self.complexes |= set(result['new_complexes'])
            # printd("stack_and_process: Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            # pprintd(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()]
            changed_domains += new_cmplx_domains
            # printd("stack_and_process: New complexes domains:")
            # pprintd(new_cmplx_domains)
        if result['obsolete_complexes']:
            self.complexes -= set(result['obsolete_complexes'])
            self.removed_complexes += result['obsolete_complexes']
            # printd("stack_and_process: Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            # pprintd(self.complexes)
            # printd("stack_and_process: sysmgr.removed_complexes:")
            # pprintd(self.removed_complexes)


        assert all(isinstance(domain, Domain) for domain in changed_domains)

        # Update hybridization reactions:
        self.update_possible_hybridization_reactions(
            changed_domains, recalculate_these={e.domain for e in (h1end3p, h2end5p, h2end3p, h1end5p)})
        # Update stacking reactions:
        self.update_possible_stacking_reactions(changed_domains, reacted_pairs=(stacking_pair,),
                                                reaction_attrs=(reaction_attr,))
        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this to hybridization/dehybridization methods and apply conditionally.
        if h1end3p.domain.strand.complex:
            h1end3p.domain.strand.complex.reset_state_fingerprint()
        if h2end3p.domain.strand.complex:
            h2end3p.domain.strand.complex.reset_state_fingerprint()

        return stacking_pair, result


    def unstacking_rate_constant(self, h1end3p, h2end3p, T=None):
        """
        Rate for breaking stacking interaction between duplex ends
        (giving the domain 3p end of involved duplexes).
        Note: stacking_rate_constant is a constant attribute, not a function.
        """
        if T is None:
            T = self.temperature
        h1end5p = h2end3p.hyb_partner
        h2end5p = h1end3p.hyb_partner
        ## TODO: Remove assertions
        assert h1end3p.stack_string == h2end3p.stack_string
        assert h1end3p.stack_string == h2end3p.stack_string == h1end5p.stack_string == h2end5p.stack_string
        assert h1end3p.stack_string is not None
        stack_dH, stack_dS = DNA_NN4_R[h1end3p.stack_string] # Why does this have None as key?
        ## TODO: Remove stack_string equality assertion
        # K = exp(stack_dS - stack_dH/T)
        # k_off = self.stacking_rate_constant/K
        k_off = self.stacking_rate_constant * exp(stack_dH/T - stack_dS)
        return k_off


    def calculate_stacking_c_j(self, h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra,
                               reaction_type=STACKING_INTERACTION, reaction_spec_pair=None):
        """
        Calculate stochastic rate constant c_j for forming/breaking a stacking interaction.
        Ends annotation:
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p

        Note: Domains on the same helix may or may not be also connected by their phosphate backbone.
        E.g. you could have a hinge, where one helix is backbone-connected and the other one not.
        This is probably the most common case, e.g. in N-way junctions.
        Nomenclature:
        :reaction_spec_pair: (For state species - previously "reaction_pair_fingerprint")
         - For stacking interactions, this is called "stackspec_pair", for domains its called "domspec_pair".
         - The general case its called "reaction_spec_pair".

        Q: Where do we apply throttling? Here? Or under the "undate_*_reactions" methods?
        A: Maybe you could have a "calculate_c_j" wrapper which takes care of caching and throttling?
        """
        # TODO: Add caching of stacking c_j
        # Caching is mostly needed for intracomplex reactions, and we already cache intracomplex_activity results!
        if reaction_spec_pair is None:
            reaction_spec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                            (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        else:
            ## TODO: Remove excessive assertions
            assert reaction_spec_pair == frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                                    (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        # Comment out during testing:
        if reaction_spec_pair in self.cache['stochastic_rate_constant']:
            return self.cache['stochastic_rate_constant'][reaction_spec_pair]

        # Comment out in production:
        # Reset complex (and domain) state fingerprints and re-check:
        # complexes = {c for c in (end.domain.strand.complex for end in (h1end3p, h2end5p, h2end3p, h1end5p))
        #              if c is not None}
        # for c in complexes:
        #     c.reset_state_fingerprint()
        # assert reaction_spec_pair == frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
        #                                         (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        d1 = h1end3p.domain
        d2 = h2end3p.domain
        # printd("Re-setting complex state and re-calculating state fingerprint:")
        # d1.state_change_reset()
        # d2.state_change_reset()
        # d1.strand.complex.reset_state_fingerprint(reset_domains=True)
        # d2.strand.complex.reset_state_fingerprint(reset_domains=True)
        # reaction_spec_pair_fresh = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
        #                                       (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        # try:
        #     assert reaction_spec_pair == reaction_spec_pair_fresh
        # except AssertionError as e:
        #     print("Got different reaction_spec_pair fingerprints after resetting domains/complexes!")
        #     raise e
        # if d1.strand.complex is not None:
        #     d1.strand.complex.reset_state_fingerprint()
        # if d2.strand.complex is not None:
        #     d2.strand.complex.reset_state_fingerprint()
        assert h1end3p.hyb_partner == h2end5p
        assert h1end3p == h2end5p.hyb_partner
        assert h2end3p.hyb_partner == h1end5p
        assert h2end3p == h1end5p.hyb_partner
        # printd("calculate_stacking_c_j invoked with h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra:")
        # printd(", ".join(str(p) for p in (h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra)))
        if is_forming:
            # if d1.strand.complex == d2.strand.complex != None:
            # Edit: We might have a single strand interacting with itself
            if is_intra:
                assert d1.strand.complex == d2.strand.complex
                assert d1.strand.complex is not None or d1.strand == d2.strand
                # Intra-complex reaction:
                # This should really be the average of dist(h1end3p, h1end5p), dist(h2end3p, h2end5p).
                # However, except for the case where the helices are already backbone connected,
                # we can approximate this by just one distance, dist(h1end3p, h2end3p):
                if h1end3p.pb_downstream == h1end5p and h2end3p.pb_downstream == h2end5p:
                    stochastic_activity = 1
                # (h1end3p, h2end5p), (h2end3p, h1end5p)
                #stochastic_activity = self.intracomplex_activity(h1end3p, h2end3p, reaction_type=STACKING_INTERACTION)
                stochastic_activity = self.intracomplex_activity((h1end3p, h2end5p), (h2end3p, h1end5p),
                                                                 reaction_type=STACKING_INTERACTION)
                # printd("reaction_spec_pair:")
                # pprintd(reaction_spec_pair)
                # printd(" - stacking stochastic_activity (intra reaction) = %s" % stochastic_activity)
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
                # printd(" - stacking stochastic_activity (inter reaction) = %s" % stochastic_activity)
            if self.include_steric_repulsion:
                steric_activity_factor = ((1 if d1.strand.complex is None else self.steric_activity_factor(d1)) *
                                          (1 if d2.strand.complex is None else self.steric_activity_factor(d2)))
                stochastic_activity *= steric_activity_factor

            k_on = self.stacking_rate_constant # (h1end3p, h2end3p)    # is fairly constant as long as ΔG° < 0
            c_j = k_on * stochastic_activity
            # printd(" - stacking stochastic rate constant c_j = %s * %s = %s" % (k_on, stochastic_activity, c_j))
        else:
            # De-hybridization reaction;
            k_off = self.unstacking_rate_constant(h1end3p, h2end3p)  # k_off depends on ΔG°
            c_j = k_off # For uni-molecular reactions, c_j = k_j
            # printd(" - un-stacking stochastic rate constant c_j = k_off = %s" % (c_j,))

        reaction_attrdict = {'ends': (h1end3p, h2end5p, h2end3p, h1end5p),
                             'is_forming': is_forming,
                             'is_intra': is_intra,
                             'reaction_type': reaction_type,
                             'c_j': c_j,
                             ('k_%s' % ("on" if is_forming else "off")): k_on if is_forming else k_off,
                             'stack_seq (all - None if forming)': tuple(end.stack_string for end in
                                                                        (h1end3p, h2end5p, h2end3p, h1end5p)),
                             'steric_activity_factor': steric_activity_factor
                                                       if (is_forming and self.include_steric_repulsion) else None,
                             'stochastic_activity': stochastic_activity if is_forming else "1 (unimolecular reac.)",
                             'domains (str, frozen):': tuple(repr(end.domain)
                                                             for end in (h1end3p, h2end5p, h2end3p, h1end5p)),
                             'domains (now):': tuple(end.domain for end in (h1end3p, h2end5p, h2end3p, h1end5p)),
                             'complexes stacking fingerprint': "c stack fp", # Key too long to pretty print...
                             #"c stack fp": tuple(c.stacking_fingerprint() for c in complexes),
                             #"c.stacking_edges": tuple(repr(c.stacking_edges()) for c in complexes),
                            }

        if reaction_spec_pair in self.cache['stochastic_rate_constant']:
            try:
                c_j_cache = self.cache['stochastic_rate_constant'][reaction_spec_pair]
                #assert np.isclose(c_j, self.cache['stochastic_rate_constant'][reaction_spec_pair])
                atol, rtol = 1e-3, 1e-3
                #assert abs(c_j_cache - c_j) < 1e-2 and abs(c_j_cache - c_j)/c_j < 1e-5
                assert abs(c_j_cache - c_j) < (atol + rtol*abs(c_j)) # avoid zero-division
            except AssertionError as e:
                print(("\n\nERROR: Newly calculated c_j %0.03e for %r reaction_spec_pair %s "
                       "is different from old value %0.03e") %
                      (c_j, reaction_type, reaction_spec_pair,
                       self.cache['stochastic_rate_constant'][reaction_spec_pair]))
                print("Re-setting complex state and re-calculating state fingerprint:")
                d1.state_change_reset()
                d2.state_change_reset()
                reaction_spec_pair_fresh = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                                      (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
                print("Fresh state fingerprint: %s" % (reaction_spec_pair_fresh,))
                print(" - same as old reaction_spec_pair?", reaction_spec_pair_fresh == reaction_spec_pair)
                print("stochastic_rate_constant_attrs:")
                pprint(reaction_attrdict)
                print("Previous stochastic_rate_constant_attrs:")
                pprint(self.cache['stochastic_rate_constant_attrs'][reaction_spec_pair])
                print("\nd1.print_history:")
                d1.print_history()
                print("\nd2.print_history:")
                d2.print_history()
                c1, c2 = d1.strand.complex, d1.strand.complex
                if c1 or c2:
                    if c1:
                        print("\nc1.print_history:")
                        c1.print_history()
                    if c2 == c1:
                        print("\nc1 == c2")
                    else:
                        print("\nc2.print_history:")
                        c2.print_history()
                raise e
        else:
            self.cache['stochastic_rate_constant'][reaction_spec_pair] = c_j
            self.cache['stochastic_rate_constant_attrs'][reaction_spec_pair].append(reaction_attrdict)

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
