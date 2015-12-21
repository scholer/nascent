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
 * Propensity functions
 *

"""

import os
#import random
from collections import defaultdict
#import math
from math import exp #, log as ln
#from datetime import datetime
from pprint import pprint
import networkx as nx
from networkx.algorithms.components import connected_components, connected_component_subgraphs
# import numpy as np
import pdb

from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
from .constants import R, N_AVOGADRO, AVOGADRO_VOLUME_NM3 #, R # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K
from .constants import PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION
from .complex import Complex
from .nx_utils import draw_graph_and_save
from .reactionmgr import ReactionMgr






class ReactionMgrGrouped(ReactionMgr):
    """
    This sub-class of ReactionMgr will do about the same as reactionmgr (same interfaces, etc),
    but it will group the domains by state specie.

    This enables the "molecule count" dependent propensity functions used by Gillepsie DM:
        a_j = c_j * x₁ * x₂     (for bi-molecular hybridization reactions), or
        a_j = c_j * x₁₂         (for uni-molecular dehybridization reactions).

    If state-specie grouping is not used, we can assume a_j = c_j for all reactions,
    but we have to iterate over *all* combinations of domain instance pairs,
    instead of iterating over combinations of domain state species.

    Grouping is theoretically faster for simulations with high copy-numbers and a relatively
    small distribution of actual states.
    However, for low copy-numbers or systems that "spread out" over a large state-permutation space,
    the benefit is minor and does not offset the cost of managing domain state species reactions.

    ## Symbol nomenclature:
    ## Sᵢ, S₁, S₂ - domain super-species (domain name). Domains with the same sequence are of the same specie.
    ## Fᵢ, F₁, F₂ - domain state fingerprint - domains of the same super-specie and in the same strand/complex state.
    ## Dᵢ, D₁, D₂ - domain instances, individual domain molecules.

    """

    def __init__(self, volume, strands, params, domain_pairs=None, **kwargs):
        self.fnprefix = "grouped"

        # Up-to-date list of hybridizations that are currently possible:
        # contains [dom_specie_spec, c_j, is_forming] tuples
        # {(domain1, cstate), (domain2, cstate)} => c_j, is_forming
        self.possible_hybridization_reactions = {}  # Rj, j=0..M - however, I index by reaction_spec
        self.depleted_hybridization_reactions = {}
        # Quick lookup with {F₁: [{F₁, F₂}, ...]}; equivalent to
        # {domspec for k in possible_hybridization_reactions if domspec in k}
        self.hybridization_reactions_by_domspec = defaultdict(set)

        # Keep track of domain state depletions.
        # If the (grouped) domain species count occilates between 0 and 1, then the grouped approach might be
        # disadvantageous. Instead, propensities should be evaluated for each domain instance.
        # Or you can have a mixed model where hybridized/duplexed domains and domains on free strands are grouped
        # while unhybridized domains within a complex are considered individually.
        # Or you could consider intra-complex reactions individually while group other reactions.
        self.N_state_depletions = 0
        self.N_state_repletions = 0


        # This class has some attributes with same name as the non-grouped class, but
        # with different content. This includes:
        # - possible_hybridization_reactions[reaction_spec] => c_j  ('is_forming' is now part of reaction_spec)
        # - depleted_hybridization_reactions[reaction_spec] => c_j
        # - hybridization_reactions_by_domspec[domspec]     => reaction_spec
        # - hybridization_propensity_functions[reaction_spec]             => a_j
        # In general, where the un-grouped super-class uses pairs of {domain1, domain2} instances,
        # this class uses pairs of {domspec1, domspec2} domain species.
        # ("domspec" usually refer to a domain state "fingerprint" that is identical for
        #  all domains in a particular state and different for domains in different states.)
        # Thus, a reaction_spec is no longer: ({domain1, domain2}, is_hyb, is_intra)
        # but rather ({domspec1, domspec2}, is_hyb, is_intra)

        super().__init__(volume, strands, params, domain_pairs=None)

        ## Domain state species related attributes: ##

        # mapping with the different domains in different states
        # indexed as: Fᵢ => list of domains in this conformation
        # For simple complexes, Fᵢ is just (domain-strand-specie, complex-state-fingerprint)
        # However, if the complex has multiple copies of that strand, the above is insufficient to specifying
        # which strand-domain within the complex we are talking about.
        # Using Fᵢ = domain.state_fingerprint() will give us a proper hash.
        self.domain_state_subspecies = defaultdict(set)
        for domain in self.domains_list:
            self.domain_state_subspecies[domain.state_fingerprint()].add(domain)

        # Propensity functions:  aj(x)
        # {(domain1, cstate), (domain2, cstate)} => a
        self.hybridization_propensity_functions = {} # Indexed by indexed by ({domain1, domain2}, is_hyb, is_intra)
        if strands:
            self.init_all_propensity_functions()



    def init_possible_reactions(self):
        """
        Reactions have:
            reaction_spec => propensity_constant c_j, state_change_vector v_j
        However, I do:
            ({domspec1, domspec2}, is_forming, is_intracomplex) => propensity_constant.
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
        print("Initializing possible reactions (GROUPED)...")
        Rxs = {}
        for d1 in self.domains:
            d1_domspec = d1.state_fingerprint()
            if d1.partner is not None:
                is_hyb, is_intra = False, True
                # If d1 is hybridized, then there is only one possible reaction channel: dehybridizing
                d2 = d1.partner
                d2_domspec = d2.state_fingerprint()
                domspec_pair = frozenset((d1_domspec, d2_domspec))
                reaction_spec = (domspec_pair, is_hyb, is_intra)
                if reaction_spec not in Rxs:
                    # tuple values are (propensity constant c_j, is_forming)
                    Rxs[reaction_spec] = self.calculate_hybridization_c_j(d1, d2, is_forming=is_hyb, is_intra=is_intra)
                self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
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
                    d2_domspec = d2.state_fingerprint()
                    domspec_pair = frozenset((d1_domspec, d2_domspec))
                    is_intra = (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex) or \
                               (d1.strand == d2.strand)  # intra-complex OR intra-strand reaction
                    reaction_spec = (domspec_pair, is_hyb, is_intra)
                    if reaction_spec not in Rxs:
                        # R_j = (c_j, v_j) - propensity constant for reaction j
                        Rxs[reaction_spec] = self.calculate_hybridization_c_j(d1, d2, is_forming=True, is_intra=is_intra)
                    self.hybridization_reactions_by_domspec[d1_domspec].add(reaction_spec)
                    self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
        self.possible_hybridization_reactions = Rxs
        print(len(self.possible_hybridization_reactions), "possible hybridization reactions initialized.")


    def update_possible_hybridization_reactions(self, changed_domains, reacted_pair, reaction_attr, reaction_domspec_pair=None):
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
        """
        print(("\nupdate_possible_hybridization_reactions invoked with d1=%s, d2=%s, is_forming=%s, reaction_domspec_pair=%s, "
               "changed_domains:") % (d1, d2, is_forming, reaction_domspec_pair))
        pprint(changed_domains)

        ### INITIAL ASSERTIONS. FAIL FAST AND FAIL HARD. ###
        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        generated_reactions_by_domspec_map = self.map_reaction_specs_by_domspec()

        if self.hybridization_reactions_by_domspec != generated_reactions_by_domspec_map:
            print("\n!! self.hybridization_reactions_by_domspec != self.map_reaction_specs_by_domspec()")
            print("self.hybridization_reactions_by_domspec:")
            pprint(self.hybridization_reactions_by_domspec)
            print("self.map_reaction_doms_specs_by_domspecs():")
            pprint(generated_reactions_by_domspec_map)
            pdb.set_trace()

        n_possible_start = len(self.possible_hybridization_reactions)
        n_propensity_start = len(self.hybridization_propensity_functions)
        n_species_start = len(self.domain_state_subspecies)

        # Add changed domains to new species
        depleted_domspecs = set()
        changed_domspecs = set()
        new_domspecs = set()
        old_d1d2_doms_specs = frozenset((d1._specie_state_fingerprint, d2._specie_state_fingerprint))

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


        for domain in changed_domains:
            # print("update_possible_hybridization_reactions: processing changed domain:", domain)
            # IMPORTANT: changed_domains must not contain any
            # Set vs list:
            # - Lists are slightly faster to iterate over;
            # - Sets are significantly faster for lookups and removing arbitrary elements (about 80 times faster)
            old_domspec = domain._specie_state_fingerprint
            # For d1, d2, old_domspec yields the new domspec instead!
            try:
                # This is really expensive for lists so domain_species should be a set
                print("\nRemoving domain %s from old domain_state_subspecies group %s" %
                      (domain, old_domspec))
                self.domain_state_subspecies[old_domspec].remove(domain)
                print(" - it now has %s elements." % len(self.domain_state_subspecies[old_domspec]))
            except KeyError as e:
                print("\ndomain %s (%s) not in self.domain_state_subspecies[%s] (%s elements):" % (
                    domain, repr(domain), old_domspec, len(self.domain_state_subspecies[old_domspec])))
                pprint(self.domain_state_subspecies[old_domspec])
                print("\nself.domain_state_subspecies:")
                pprint(self.domain_state_subspecies)
                for domspec_i, spec_domain_set in self.domain_state_subspecies.items():
                    if domain in spec_domain_set:
                        print(" - domain %s was found in self.domain_state_subspecies[%s]:" % (domain, domspec_i),
                              self.domain_state_subspecies[domspec_i])
                if domain not in (d1, d2) or True:
                    raise e
            if len(self.domain_state_subspecies[old_domspec]) == 0:
                depleted_domspecs.add(old_domspec)
                del self.domain_state_subspecies[old_domspec]
            else:
                changed_domspecs.add(old_domspec)  # change in count of this domspec
            # print("Re-setting and re-calculating state_fingerprint for %s - old is: %s"
            #       % (domain, domain._specie_state_fingerprint))
            domain.state_change_reset()
            new_domspec = domain.state_fingerprint()
            if new_domspec == old_domspec:
                print("\nWeird: new_domspec == old_domspec for changed domain %s: %s == %s" %
                      (domain, new_domspec, old_domspec))
            if new_domspec not in self.domain_state_subspecies:
                # perhaps use self.hybridization_reactions_by_domspec ??
                new_domspecs.add(new_domspec)
            else:
                changed_domspecs.add(new_domspec)
            print("\nAdding domain %s to domain_state_subspecies under new domspec: %s" % \
                  (domain, new_domspec))
            self.domain_state_subspecies[new_domspec].add(domain)
            print(" - it now has %s elements." % len(self.domain_state_subspecies[new_domspec]))

        ## Do we need to update hybridized domains (duplexes)?
        ## - As long as we don't have bending that affects k_off rates: No.
        ## But then... how to keep track of hybridized domains??
        ## Things like stacking/unstacking will affect their k_off rate constants, so the count depends on
        ## the state. Unless I do it completely different and just say that
        ## the domspec "fingerprint" for a hybridized domain does not depend on the complex state at all,
        ## but only depends on whether the domain is stacked or not.
        ## E.g. fingerprint = (domain_species_name, (end5p.stack_partner is not None, end3p.stack_partner is not None)

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())


        ## Consideration: Do we really need to check for depleated states all the time?? Maybe just do it
        ## once every 10th run or so.

        ## Depleted domain states:
        for domspec in depleted_domspecs:
            # This could be replaced by using hybridization_reactions_by_domspec
            del self.hybridization_reactions_by_domspec[domspec]
            # reaction_spec[0] = domspec_pair is frozendict((domspec1, domspec2))
            depleted_reactions = {reaction_spec: v for reaction_spec, v in self.possible_hybridization_reactions.items()
                                  if domspec in reaction_spec[0]}
            for reaction_spec, v in depleted_reactions.items():
                domspec_pair = reaction_spec[0]
                try:
                    print("Deleting reaction_spec %s from possible_hybridization_reactions..." % (reaction_spec,))
                    del self.possible_hybridization_reactions[reaction_spec]
                    print("Deleting reaction_spec %s from hybridization_propensity_functions..." % (reaction_spec,))
                    del self.hybridization_propensity_functions[reaction_spec]
                    domspec2 = next(F for F in domspec_pair if F != domspec)
                    print("Removing reaction_spec from hybridization_reactions_by_domspec[domspec1/2]...")
                    #self.hybridization_reactions_by_domspec[domspec].remove(reaction_spec) # We have already deleted.
                    self.hybridization_reactions_by_domspec[domspec2].remove(reaction_spec)
                except KeyError as e:
                    print("Unexpected KeyError:", e)
                    print("Starting pdb...")
                    pdb.set_trace()
                self.depleted_hybridization_reactions[reaction_spec] = v
                self.N_state_depletions += 1

        ## TODO: Check that domspecs in hybridization_reactions_by_domspec matches possible_hybridization_reactions
        ##       and hybridization_propensity_functions.
        ## TODO: Consider consolidating possible_hybridization_reactions and hybridization_propensity_functions to a single dict:
        ##          {{F1, F2}: is_forming, c_j, a_j}
        ##          However, that would make it harder to (a) do expansions, or (b) know if a_j has been calculated.

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())


        # Keep track of which reaction paths we have already updates so we won't re-calculate them.
        # E.g. if we have re-calculated reaction (c hybridize C), don't re-calculate (C hybridize c).
        updated_reactions = set()

        ## New domain states:
        for domspec in new_domspecs:
            repleted_reactions = {reaction_spec: v for reaction_spec, v in self.depleted_hybridization_reactions.items()
                                  if domspec in reaction_spec[0] and all(
                                      len(self.domain_state_subspecies[F]) > 0 for F in reaction_spec[0])}
            for reaction_spec, c_j in repleted_reactions.items():
                del self.depleted_hybridization_reactions[reaction_spec]
                self.possible_hybridization_reactions[reaction_spec] = c_j
                self.N_state_repletions += 1
                for ds in reaction_spec[0]:
                    self.hybridization_reactions_by_domspec[ds].add(reaction_spec)
            # When do we add domain with domspec to domain_state_subspecies ?
            ## TODO: DO NOT OVERWRITE OUTER 'd1' and 'd2' VARIABLES !!
            domain1 = next(iter(self.domain_state_subspecies[domspec]))
            ## TODO: Include hybridization state in domspec. Then you have that as well strand name and domain name,
            ##          meaning you don't have to find a concrete domain instance.
            if domain1.partner is not None:
                # is_hyb not same as is_forming (specifies what just happened to d1, d2))
                is_hyb, is_intra = False, True
                domain2 = domain1.partner
                d2_domspec = domain2.state_fingerprint()
                domspec_pair = frozenset((domspec, d2_domspec))
                reaction_spec = (domspec_pair, is_hyb, is_intra)
                self.hybridization_reactions_by_domspec[domspec].add(reaction_spec)
                self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
                if reaction_spec not in self.possible_hybridization_reactions:
                    # tuple values are (propensity constant c_j, is_forming)
                    # domain1 has a partner (domain2), so this must be a de-hybridization reaction:
                    print(("de-hybrid reaction_spec %s for new domspec %s is not in possible_hybridization_reactions, "
                           "calculating c_j...") % (reaction_spec, domspec))
                    self.possible_hybridization_reactions[reaction_spec] = \
                        self.calculate_hybridization_c_j(domain1, domain2, is_forming=False, is_intra=is_intra)
                elif reaction_spec not in updated_reactions:
                    # Should not happen - it is supposedly a "new domain state".
                    # Edit: Can happen when c and C is changed: (C hybridize c) AND (c hybridize C)
                    print(("Weird: reaction_spec %s for hybridized domains domain1, domain2 = (%s, %s) is "
                           "already in possible_reactions - that should not happen.") % (reaction_spec, domain1, domain2))
                else:
                    print(" GOOD: reaction_spec %s for new domspec %s is already in possible_hybridization_reactions." \
                          % (reaction_spec, domspec))
                if reaction_spec not in updated_reactions:
                    self.recalculate_propensity_function(reaction_spec)
                    updated_reactions.add(reaction_spec)
            else:
                ## No partner -- consider hybridization reactions with all other unhybridized complementary domains.
                # for domain2 in [d for cname in self.domain_pairs[domain1.name]
                #            for d in self.unhybridized_domains_by_name[cname]]:
                #pdb.set_trace()
                if domain1.name in self.domain_pairs:
                    is_hyb = True
                    try:
                        for cname in self.domain_pairs[domain1.name]:
                            for domain2 in self.unhybridized_domains_by_name[cname]:
                                d2_domspec = domain2.state_fingerprint()
                                domspec_pair = frozenset((domspec, d2_domspec))
                                is_intra = (domain1.strand.complex is not None and \
                                            domain1.strand.complex == domain2.strand.complex) or \
                                           (domain1.strand == domain2.strand)  # intra-complex OR intra-strand reaction
                                reaction_spec = (domspec_pair, is_hyb, is_intra)
                                if reaction_spec not in self.possible_hybridization_reactions:
                                    # R_j = (c_j, v_j) - propensity constant for reaction j
                                    print(("hybridizing reaction_spec %s for new domspec %s is not in possible"
                                           "_hybridization_reactions, calculating c_j...") % (reaction_spec, domspec))
                                    self.possible_hybridization_reactions[reaction_spec] = \
                                        self.calculate_hybridization_c_j(domain1, domain2, is_forming=True, is_intra=is_intra)
                                else:
                                    print(" GOOD: reaction_spec %s for new domspec %s is already in possible_hybridization_reactions." \
                                          % (reaction_spec, domspec))
                                self.hybridization_reactions_by_domspec[domspec].add(reaction_spec)
                                self.hybridization_reactions_by_domspec[d2_domspec].add(reaction_spec)
                                if reaction_spec not in updated_reactions:
                                    self.recalculate_propensity_function(reaction_spec)
                                    updated_reactions.add(reaction_spec)
                    except KeyError as e:
                        print("KeyError:", e)
                        print("domain1, domain1.name:", domain1, ",", domain1.name)
                        print("self.domain_pairs:")
                        pprint(self.domain_pairs)
                        print("self.domain_pairs[domain1.name]:")
                        pprint(self.domain_pairs[domain1.name])
                        print("self.unhybridized_domains_by_name:")
                        pprint(self.unhybridized_domains_by_name)
                        print("All domains that should be unhybridized (partner is None):")
                        real_unhybridized_domains = [d for d in self.domains if d.partner is None]
                        pprint(real_unhybridized_domains)
                        print(" - %s elements" % len(real_unhybridized_domains))
                        print("Difference (symmetric):")
                        print(set(self.unhybridized_domains_by_name.keys()) ^
                              set(d.name for d in real_unhybridized_domains))
                        print("domain2, domain2.name:", domain2, ",", domain2.name)
                        pdb.set_trace()

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        ## Changed domains states (i.e. only a change in number of domains in that state)
        ## Just update hybridization_propensity_functions for those domspec
        for domspec in changed_domspecs: #| new_domspecs:
            ## TODO: Make sure we do not calculate the full product "matrix", but only one of the symmetric halfs.
            # print(("Re-evaluating propensity for all hybridization reactions involving domspec "
            #        "%s in changed_domspecs | new_domspecs..." % (domspec, )))
            for reaction_spec in self.hybridization_reactions_by_domspec[domspec]:
                assert reaction_spec in self.possible_hybridization_reactions
                print("Re-calculating propensity for reaction_spec %s ..." % (reaction_spec,))
                if not isinstance(reaction_spec, tuple) and not isinstance(reaction_spec[0], frozenset):
                    pdb.set_trace()
                if reaction_spec not in updated_reactions:
                    self.recalculate_propensity_function(reaction_spec)
                    updated_reactions.add(reaction_spec)

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        n_possible_end = len(self.possible_hybridization_reactions)
        n_propensity_end = len(self.hybridization_propensity_functions)
        n_species_end = len(self.domain_state_subspecies)
        print_debug_info = False

        if n_possible_end != n_propensity_end:
            print("\n!! n_possible_end != n_propensity_end: %s != %s" % (n_possible_end, n_propensity_end))
            print_debug_info = True

        generated_reactions_by_domspec_map = self.map_reaction_specs_by_domspec()
        if self.hybridization_reactions_by_domspec != generated_reactions_by_domspec_map:
            print("\n!! self.hybridization_reactions_by_domspec != self.map_reaction_specs_by_domspec()")
            print("self.hybridization_reactions_by_domspec:")
            pprint(self.hybridization_reactions_by_domspec)
            print("self.map_reaction_specs_by_domspec():")
            pprint(generated_reactions_by_domspec_map)
            print_debug_info = True

        if print_debug_info:
            print("--- FURTHER DEBUG INFO: ---")
            print("n_possible_start, n_propensity_start:", n_possible_start, n_propensity_start)
            print("n_species_start, n_species_end:", n_species_start, n_species_end)
            print("changed_domains:", changed_domains)
            print("domain1, domain2:", domain1, domain2)
            print("reaction_domspec_pair:", reaction_domspec_pair)
            print("is_forming:", is_forming)
            print("self.possible_hybridization_reactions: (%s elements)" % len(self.possible_hybridization_reactions))
            pprint(list(self.possible_hybridization_reactions))
            print("self.hybridization_propensity_functions: (%s elements)" % len(self.hybridization_propensity_functions))
            pprint(list(self.hybridization_propensity_functions))
            print("Shared keys: (%s elements)" %
                  len(set(self.possible_hybridization_reactions.keys()) & set(self.hybridization_propensity_functions.keys())))
            pprint(list(set(self.possible_hybridization_reactions.keys()) & set(self.hybridization_propensity_functions.keys())))
            print("Differing keys: (symmetric difference, %s elements)" %
                  len(set(self.possible_hybridization_reactions.keys()) ^ set(self.hybridization_propensity_functions.keys())))
            pprint(list(set(self.possible_hybridization_reactions.keys()) ^ set(self.hybridization_propensity_functions.keys())))

            print("self.depleted_hybridization_reactions: (%s elements)" % len(self.depleted_hybridization_reactions))
            pprint(list(self.depleted_hybridization_reactions.keys()))

            print("self.domain_state_subspecies.keys(): (%s elements)" % len(self.domain_state_subspecies))
            pprint({k: len(v) for k, v in self.domain_state_subspecies.items()})
            #pprint(list(self.domain_state_subspecies.keys()))

            print("new_domspecs: (%s elements)" % len(new_domspecs))
            pprint(new_domspecs)
            print("changed_domspecs: (%s elements)" % len(changed_domspecs))
            pprint(changed_domspecs)
            print("depleted_domspecs: (%s elements)" % len(depleted_domspecs))
            pprint(depleted_domspecs)

            print("self.hybridization_reactions_by_domspec:")
            pprint(self.hybridization_reactions_by_domspec)
            pdb.set_trace()

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        self.print_reaction_stats()


    def print_reaction_stats(self):
        """ Print information on reaction pathways. """
        print("Reaction pathways: (grouped reactions)")
        for reaction_spec, c_j in self.possible_hybridization_reactions.items():
            domspec1, domspec2 = tuple(reaction_spec[0])
            is_forming = reaction_spec[1]
            print("%20s %s (%s) %20s" % (domspec1[0],
                                     "hybridize" if is_forming else "de-hybrid",
                                     "intra" if reaction_spec[2] else "inter",
                                     domspec2[0]),
                  (": %0.03e x%02s x%02s = %03e" % (
                      c_j, len(self.domain_state_subspecies[domspec1]), len(self.domain_state_subspecies[domspec2]),
                      self.hybridization_propensity_functions[reaction_spec])
                  ) if is_forming else (":   %0.03e x%02s   = %03e" % (
                      c_j, len(self.domain_state_subspecies[domspec1]), self.hybridization_propensity_functions[reaction_spec]))
                 )


    def map_reaction_specs_by_domspec(self):
        """
        Generate an equivalent to self.hybridization_reactions_by_domspec
        from self.possible_hybridization_reactions, i.e. a map of
            domspec => {set of reaction_specs involving this domspec}
        """
        reactions_by_domspec = defaultdict(set)
        for reaction_spec in self.possible_hybridization_reactions:
            for domspec in reaction_spec[0]:
                reactions_by_domspec[domspec].add(reaction_spec)
        return reactions_by_domspec


    def get_reactions_by_domspec(self, domspec):
        """ Filter possible_hybridization_reactions, returning only reactions involving :domspec:. """
        return {reaction_spec: v
                for reaction_spec, v in self.possible_hybridization_reactions.items()
                if domspec in reaction_spec[0]}



    def recalculate_propensity_function(self, reaction_spec):
        """
        Recalculate propensity function for a single reaction_spec =
            (frozenset((F1, F2)), is_forming, is_intra)
        Assumes that c_j, is_forming is already available in
        self.possible_hybridization_reactions.
        """
        try:
            c_j = self.possible_hybridization_reactions[reaction_spec]
        except KeyError as e:
            print("KeyError %s for self.possible_hybridization_reactions[%s]" % (e, reaction_spec,))
            print("self.possible_hybridization_reactions: (%s elements)" % len(self.possible_hybridization_reactions))
            pprint(list(self.possible_hybridization_reactions))
            print("self.hybridization_propensity_functions: (%s elements)" % len(self.hybridization_propensity_functions))
            pprint(list(self.hybridization_propensity_functions))
            print("self.domain_state_subspecies.keys(): (%s elements)" % len(self.domain_state_subspecies))
            pprint({k: len(v) for k, v in self.domain_state_subspecies.items()})
            print("self.hybridization_reactions_by_domspec:")
            pprint(self.hybridization_reactions_by_domspec)

            pdb.set_trace()
        is_forming = reaction_spec[1]
        if is_forming:
            # a_j = c_j * x₁ * x₂
            self.hybridization_propensity_functions[reaction_spec] = (c_j *
                #np.prod([len(self.domain_state_subspecies[ds]) for ds in reaction_spec[0]])
                (len(self.domain_state_subspecies[reaction_spec[0][0]]) *
                 len(self.domain_state_subspecies[reaction_spec[0][1]])))
        else:
            # a_j = c_j * x₃ ,      x₃ is number of duplexes
            #self.hybridization_propensity_functions[domspec_pair] = c_j * len(self.domain_state_subspecies[domspec_pair[0]])
            self.hybridization_propensity_functions[reaction_spec] = \
                c_j * len(self.domain_state_subspecies[next(iter(reaction_spec[0]))])




    def init_all_propensity_functions(self):
        """
        Reaction propensity specifies how likely a certain reaction of one or more species is to occour.
        The product "propensity × dt" denotes the probability that a reaction will occour in infinitesimal time dt.
        This is just the reaction rate times the number of reactant species present (concentration):
            propensity = k_on [domain1] [domain2]   # if hybridizing
            propensity = k_off * [duplex]           # if dehybridizing
        However, since our simulation operates stochastically at the molecular level, we prefer to express reaction
        propencities in term of discrete molecule count, as:
            a_j = c_j * X₁ * X₂         # for hybridization reactions
            a_j = c_j * X₁              # for dehybridization reactions

        Where c_j is the "stochastic rate constant" (which is essentially k_j * stochastic_activity), and
        X₁, X₂ are the number of domains of species X₁, X₂ or duplexes of species X₁.

        Gillespie calls it "propensity function", since propensity is a function of state (concentration/species counts)
        and uses the symbol a_j to indicate the propensity of reaction j. This should not be confused with "activity".
        Gillespie uses a_0 to indicate sum(a_j for a_j in a) - I just use "a_sum".

        The propencify function sum, a_0 or a_sum, is simply: a_sum = sum(a)

        Calculation of initial propencities is done here as:
            a = [c_j * (prod(len(domain_subspecies[ds]) for ds in doms_spec)  # equivalnt to: c_j * X₁ * X₂, or c_j * X₁
                        if is_forming else len(self.domain_state_subspecies[reaction_spec[0][0]]))
                 for doms_spec, c_j, is_forming in possible_reactions_lst]
        Where doms_spec is a set of Fᵢ = (domain-specie, cstate) for each domain in the reaction, used to specify
        all domains of a particular species (e.g. "domain A") in a particular state (e.g. free or in any complex state).
        The "prod(len(domain_subspecies[ds]) for ds in doms_spec)" expression is equivalent to c_j * X₁ * X₂,
        for hybridization reactions, while "len(self.domain_state_subspecies[reaction_spec[0][0]]))" is equivalent to c_j * X₁
        for dehybridzation reactions.

        For uni-molecular reactions, c_j is simply the reaction rate constant, k_j.
        For bi-molecular reactions, c_j is k_j *times* the concentration of a single reactant molecule.
        E.g. for the bi-molecular reaction j between two species, S₁ and S₂, with population/count x₁ and x₂
        and rate constant k_j, we calculate c_j as:
            c_j = k_j / N_AVOGADRO / volume    # c_j should be in units of s⁻¹
        At a constant concentration, as we increase the simulation volume, we decrease the concentration of each
        individual molecule (instance), so c_j decreases. But since the number of specie molecules/instances
        increases by the same factor, the product c_j * x₁ is constant, and the propensity a_j = c_j * x₁ * x₂
        actually *increases*. This is intuitive: If we have a larger reaction volume at a constant concentration,
        we expect more reactions to occour within a given timespan Δt. This is even simpler to see for uni-molecular
        reactions (Sᵢ -> ...), where c_j = k_j and a_j increases linearly with xᵢ, which again is linearly proportional
        to the volume. Meaning that in a given timespan Δt, the larger the volume, the more Sᵢ will react.

        We typically save the "unit-normalized" c_j instead of k_j, because:
            (1) it prevents repeatedly doing the floating-point calculation c_j = k_j / N_AVOGADRO / volume, and
            (2) we don't need to worry about uni- vs bi- vs tri-molecular etc, we just need to multiply c_j with
                the reactant species population count, x₁ (x₂, x₃, etc) to get the propensity function a_j.
        """
        print("Initializing propensity functions...")
        # domspec_pair (plural) is: frozenset({(domain1, cstate), (domain2, cstate)})
        a = {reaction_spec: (c_j * #(np.prod([len(self.domain_state_subspecies[ds]) for ds in reaction_spec[0]])
                             len(self.domain_state_subspecies[reaction_spec[0][0]]) *
                             len(self.domain_state_subspecies[reaction_spec[0][1]]))
                                   #if is_forming else len(self.domain_state_subspecies[reaction_spec[0][0]]))
                                   # reaction_spec[1] = is_forming
                            if reaction_spec[1] else
                            len(self.domain_state_subspecies[next(iter(reaction_spec[0]))])
             for reaction_spec, c_j in self.possible_hybridization_reactions.items()}
        self.hybridization_propensity_functions = a
        print(len(self.hybridization_propensity_functions), "propensity functions initialized.")




    def hybridize_and_process(self, domspec_pair, reaction_attr):
        """
        Will select a random pair of domain instances from the domain species pair
        for hybridization reaction, or a random duplex in case of dehybridization reactions.
        Reaction specie consists of:
            ({domspec1, domspec2}, is_forming, is_intracomplex)
        """
        print("\nreact_and_process invoked with args: domspec_pair = %s, reaction_attr = %s" %
              (domspec_pair, reaction_attr))
        if not all(domspec in self.domain_state_subspecies for domspec in domspec_pair):
            print("domspec_pair:", domspec_pair)
            print("Not all domspec in domspec_pair are in self.domain_state_subspecies:")
            pprint(self.domain_state_subspecies)
        assert all(domspec in self.domain_state_subspecies for domspec in domspec_pair)
        is_forming = reaction_attr.is_forming
        ## HERE WE SELECT (formerly REMOVE/POP) arbitrary domains FROM self.domain_state_subspecies sets ##
        if is_forming:
            # For each domspec in domspec_pair, select a domain instance belonging to that specie population.
            # Hmm... how do we tell whether the reaction is supposed to be intra-complex hybridization
            # vs between two domains in two separate complexes? I need to have a flag for that as well.
            # Argh, another good reason NOT to have domain-specie based reactions but just have a reaction
            # path for *each* combination of domain instace pairs.
            # Obviously, you could still use domspecs for caching.
            d1, d2 = [next(iter(species_list))  # species_list.pop()
                      for species_list in
                      [self.domain_state_subspecies[dom_spec] for dom_spec in domspec_pair]]
        else:
            # Duplex de-hybridization reaction:
            d1 = next(iter(self.domain_state_subspecies[next(iter(domspec_pair))]))
            d2 = d1.partner
            assert d2 is not None

        domain_species_keys_before = {k: len(v) for k, v in self.domain_state_subspecies.items()}
        if is_forming:
            print("Hybridizing domain %s %s and %s %s..." %
                  (repr(d1), d1._specie_state_fingerprint, repr(d2), d2._specie_state_fingerprint))
            # changed_complexes, new_complexes, obsolete_complexes, free_strands = self.hybridize(d1, d2)
            result = self.hybridize(d1, d2)
            print("Completed hybridization of domain %s and %s..." % (d1, d2))
            # print("Completed hybridization of domain %s (%s) and %s (%s)..." %
            #       (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
        else:
            print("De-hybridizing domain %s %s and %s %s..." %
                  (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
            #changed_complexes, new_complexes, obsolete_complexes, free_strands = self.dehybridize(d1, d2)
            result = self.dehybridize(d1, d2)
            print("Completed de-hybridization of domain %s and %s..." % (d1, d2))
            # print("Completed de-hybridization of domain %s (%s) and %s (%s)..." %
            #       (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
        #domain_species_keys_after = set(self.domain_state_subspecies.keys())
        domain_species_keys_after = {k: len(v) for k, v in self.domain_state_subspecies.items()}
        if domain_species_keys_after != domain_species_keys_before:
            print("domain_species_keys_after != domain_species_keys_before (%s vs %s elements)" %
                  (len(domain_species_keys_after), len(domain_species_keys_before)))
            print("domain_species_keys_before:")
            pprint(domain_species_keys_before)
            print("domain_species_keys_after:")
            pprint(domain_species_keys_after)
        # else:
        #     print("domain_species_keys_after == domain_species_keys_before (%s vs %s elements)" %
        #           (len(domain_species_keys_after), len(domain_species_keys_before)))
        #     print("domain_species_keys_after hybridization/dehybridization:")
        #     pprint(domain_species_keys_after)



        # 4: Update/re-calculate possible_hybridization_reactions and hybridization_propensity_functions
        # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
        # - possible_hybridization_reactions
        # - hybridization_propensity_functions
        # Note: If evaluating whether an object is boolean False, the steps include:
        # - Does it have a __len__ attribute? - Yes? Return bool(len(obj))
        # Whereas evaluating whether "obj is None" is about 10 times faster.
        # Fixed: Excluding domains that have no partners -- these will never undergo hybridization reaction.

        print("result: (is_forming: %s)" % is_forming)
        pprint(result)
        if result['free_strands']:
            print("free strands domains:")
            pprint([s.domains for s in result['free_strands']])

        changed_domains = []
        if result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()
                                if domain.partner is None and domain.name in self.domain_pairs]
            print("Changed complexes domains:")
            pprint(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['new_complexes']:
            self.complexes |= set(result['new_complexes'])
            print("Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            pprint(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()
                                 if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += new_cmplx_domains
            print("New complexes domains:")
            pprint(new_cmplx_domains)
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains
                               if domain.partner is None and domain.name in self.domain_pairs]
            changed_domains += free_st_domains
        if result['obsolete_complexes']:
            self.complexes -= set(result['obsolete_complexes'])
            self.removed_complexes += result['obsolete_complexes']
            print("Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            pprint(self.complexes)
            print("sysmgr.removed_complexes:")
            pprint(self.removed_complexes)

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

        self.update_possible_hybridization_reactions(changed_domains, reacted_pair=(d1, d2), is_forming=is_forming,
                                       reaction_domspec_pair=domspec_pair)
        # DEBUGGING: Resetting complex fingerprint.
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        if d1.strand.complex:
            d1.strand.complex.reset_state_fingerprint()
        if d2.strand.complex:
            d2.strand.complex.reset_state_fingerprint()

        print("changed_domains _specie_state_fingerprint:")
        pprint([(d, d._specie_state_fingerprint) for d in changed_domains])
        print("changed_domains specie_state_fingerprint():")
        pprint([(d, d.state_fingerprint()) for d in changed_domains])

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        # Return the selected domains that were actually hybridized/dehybridized
        return (d1, d2), result
