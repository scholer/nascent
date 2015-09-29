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

# pylint: disable=C0103

"""

Module for

"""

import os
import random
from collections import defaultdict
import math
from math import exp, log as ln
from datetime import datetime

import networkx as nx
import numpy as np


from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
from .thermodynamic_utils import thermodynamic_meltingcurve
from .simulator import Simulator


# Module-level constants and variables
N_AVOGADRO = 6.022e23   # /mol
R = 1.987  # universal gas constant in cal/mol/K
VERBOSE = 0



def complex_sizes_hist(complexes):
    hist = {}
    for c in complexes:
        N = len(c.Strands)
        if N not in hist:
            hist[N] = 0
        hist[N] += 1
    return hist




class DM_Simulator(Simulator):
    """

    I generally use a "domain state subspecies" hash to identify subpopulations in the same state.
    This "domspec" hash is composed of:
        - domain-strand-specie name
        - complex state hash
        - in-complex-domain-identifier  - in case the complex has multiple copies of the domain-strand

    I use this "domspec" hash for many things, including
        - grouping domains in the same state
        - cache keys

    Obviously, much of the data regards a combination of two domains.
    To group these combinations, I use a frozenset of the two domspecs:
        doms_specs = frozenset({domspec1, domspec2})

    Primary data structures:
        - possible_hybridization_reactions[doms_specs] = propencity_constant cj, is_hybridizing
            Contains up-to-date hybridization and dehybridization reactions
            dict, keyed by doms_specs = {domain1-state, domain2-state}
            This is "R" and "R_j" in Gillespie's formulas
        - propencity_functions[doms_specs] = actual propencity for reaction between domain1 and domain2
            with domain-states doms_specs = {domain1-state, domain2-state}
            This is "a" or "a(x)" in Gillespie's formulas.


    Overview of the different caches:
    To avoid doing a lot of re-calculation, I generally save all results.
    If the result involves only one domain, I store it as domspec for that domain.
    If the result involves two domains, I store it under the doms_specs = {domain1-domspec, domain2-domspec}


    1-domain caches, state-dependent (domspec):
        - domain activity for inter-complex reactions.
            If the domains are not in part of the same complex, I can usually just take a "quick look around"
            and estimate steric interactions.

    2-domain caches, state-dependent (doms_specs):
        - domain activities for intra-complex reactions.
            If the two domains are part of the same complex, estimating their mutual activity is a lot harder.
            I generally see if they are mechanically, physically separated, whether they are close or far,
            and whether the domains between them form a rigid or flexible separation.
        - reaction constants k_on and k_off at a particular temperature.
            reaction constants depend on temperature, so this cache must be reset
            whenever the temperature changes.
            Reaction constants also depends on ionic strength, the length of the interacting strands,
            and to a lesser extend GC content stabilizing meta-stable intermediates.
        - _statedependent_dH_dS: Usually we consider dH and dS to be constant.
            However, in some states, e.g. for a duplex where the 5p of one strand and 3p of the other
            is being pulled perpendicularly to the helical axis, thus "zippering" them apart (or zipped when forming).
            This mechanical stress will change dH and dS of hybridization, favoring the dehybridized state.
            _statedependent_dH_dS cache should only be used for intra-complex reactions, i.e.
            duplexes within the same complex.

    2-domain caches using domain-specie *name* and not state:
        - domain_dH_dS: Standard dH and dS for hybridization between a domain of species "da"
            with another interacting domain of type "dA", both in their "free solution state".
            I considered using sequences instead of domain names, but using names seems more consistent.


    """


    def __init__(self, **kwargs):
        #super().__init__(**kwargs)
        Simulator.__init__(self, **kwargs)

        # Relative activities, indexed as:
        #  - [({domain1-specie, domain2-specie}, complex-state-fingerprint)] for intra-complex reactions
        #  - [(domain1-specie, domain1-complex-state-fingerprint)] for inter-complex reactions.
        # relative activity is 1 for free strands
        self.relative_activity_cache = defaultdict(dict)


        # mapping with the different domains in different states
        # indexed as: (domain-strand-specie, complex-state-fingerprint) => list of domains in this conformation
        # However, again we have to problem of fully specifying exactly WHAT strand-domain within the
        # complex we are talking about, if the complex has multiple copies of that strand.
        #
        self.domain_state_subspecies = defaultdict(dict)
        for domain in self.domains_list:
            self.domain_state_subspecies[domain.domain_state_fingerprint()] = domain


        # Up-to-date list of hybridizations that are currently possible:
        # contains [dom_specie_spec, c_j, is_hybridizing] tuples
        # {(domain1, cstate), (domain2, cstate)} => c_j, is_hybridizing
        self.possible_hybridization_reactions = {}  # Rj, j=0..M - however, I index by doms_specs

        # Propencity functions:  aj(x)
        # {(domain1, cstate), (domain2, cstate)} => a
        self.propencity_functions = {}

        # Futher caches:
        self._statedependent_dH_dS = {}


    def init_possible_reactions(self):
        """
        Reactions have:
            doms_specs => propencity_constant c_j, state_change_vector v_j
            {(domain1-specie, cstate), (domain2-specie, cstate)} => c_j, is_hybridizing

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
        Rxs = {}
        for d1 in self.domains:
            d1_statespecie = d1.domain_state_fingerprint()
            if d1.partner:
                # If d1 is hybridized, then there is only one possible reaction channel: dehybridizing
                d2 = d1.partner
                d2_statespecie = d2.domain_state_fingerprint()
                doms_specs = frozenset((d1_statespecie, d2_statespecie))
                if doms_specs not in Rxs:
                    # tuple values are (propencity constant c_j, is_hybridizing)
                    Rxs[doms_specs] = self.calculate_propencity_constant(d1, d2, is_hybridizing=False)
            else:
                # find all possible reaction channels for domain d1
                # self.domain_pairs[dname] = [list of complementary domain names]
                # all possible hybridization partners/candidates:
                for d2 in (d for cname in self.domain_pairs[d1.name] for d in self.domains_by_name[cname]):
                    d2_statespecie = d2.domain_state_fingerprint()
                    doms_specs = frozenset((d1_statespecie, d2_statespecie))
                    if doms_specs not in Rxs:
                        # R_j = (c_j, v_j) - propencity constant for reaction j
                        Rxs[doms_specs] = self.calculate_propencity_constant(d1, d2, is_hybridizing=True)
        self.possible_hybridization_reactions = R


    def init_all_propencities(self):
        """
        a = [c_j * prod(len(domain_subspecies(ds)) for ds in doms_spec)
             for doms_spec, c_j, is_hybridizing in possible_reactions_lst]  # doms_spec is the (domain-specie, cstate) above.
        # to get the sum:
        a_sum = sum(a)
        """
        # doms_specs (plural) is: frozenset({(domain1, cstate), (domain2, cstate)})
        # propensity = k_on [domain1] [domain2] if hybridizing else k_off * [duplex]
        # TODO: Check for factor 2 on [duplex] vs [domain1-hybridized] + [domain2-hybridized]
        #a = {doms_specs: c_j * (np.prod([len(self.domain_state_subspecies[ds]) for ds in doms_specs])
        #                        if id_hybridizing else len(self.domain_state_subspecies[doms_specs[0]]))
        # Edit: either using product (for [domain1]*[domain2]) or sum (for [domain1]+[domain2])
        a = {doms_specs: c_j * \
             (np.prod if is_hybridizing else sum)([len(self.domain_state_subspecies[ds]) for ds in doms_specs])
             # TODO: propencities must also include the relative activity of d1 against d2
             # (if this is not already included in c_j for that reaction)
                # reduces to prod([N_domain1, N_domain2]) or sum([N_domain1, N_domain2])
             for doms_specs, (c_j, is_hybridizing) in self.possible_hybridization_reactions.items()}  # doms_spec is the (domain-specie, cstate) above.
        self.propencity_functions = a


    def calculate_propencity_constant(self, d1, d2, is_hybridizing):
        volume = self.volume
        if is_hybridizing:
            k_on = self.hybridization_rate_constant(d1, d2)
            # hybridization rate constant is in unit of /M/s = L/mol/s.
            # We need it in number of molecules, so divide by NA
            c = k_on/volume/N_AVOGADRO # should be e.g. 0.1 /s
        else:
            k_off = self.dehybridization_rate_constant(d1, d2)
            # k_off depends on ΔG°  (at least when we are at T < Tm at molar concentrations where ΔG° < 0)
            c = k_off
        return c, is_hybridizing


    def hybridization_rate_constant(self, d1, d2):
        #T = self.temperature
        # At T >> Tm (ΔG° ~ 0), this depends on ΔG°.
        # Also depends on:
        # * salt / ionic strength
        # * length of d1, d2 (and the full length of their strands)
        # But for now, just return default rate constant of 1e5 /s/M
        return 1e5


    def dehybridization_rate_constant(self, d1, d2):
        """
        Calculate dehybridization rate constant:
            k_off = k_on/K,   K = exp(-ΔG°/RT)
                  = k_on * exp(+ΔG°/RT)
        """
        T = self.temperature
        # dH, dS, dG at standard conditions
        dH, dS = self.dHdS_from_state_cache(d1, d2) # Units of R, R/K
        #dG = dH - T*dS # No need to perform extra operation, just use dH and dS:
        k_off = 1e5 * exp(dH/T - dS)
        return 1e5


    def dHdS_from_state_cache(self, d1, d2):
        """
        Cache hybridization energy
        Returns
            dH, dS
        where dH is in units of gas constant R, and dS in units of R/K
        """
        doms_state_hash = frozenset(d.domain_state_fingerprint() for d in (d1, d2))
        if doms_state_hash not in self._statedependent_dH_dS:
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
        if d1.Complex == d2.Complex != None:
            rel_act = 1
        if d1.Complex is None:
            d1_rel_act = 1



    def simulate(self, T=None, n_steps_max=100000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        if T is None:
            T = self.temperature
        else:
            self.temperature = T

        n_done = 0
        while n_done < n_steps_max:

            # Step 1. (I expect that propencity_functions is up-to-date)
            reactions, propencity_functions = zip(list(self.propencity_functions.items()))
            # reactions: list of possible reactions
            a = propencity_functions # propencity_functions[j]: propencity for reaction[j]
            a0_sum = sum(propencity_functions)

            # Step 2: Generate values for τ and j:  - easy.
            r1, r2 = random.random(), random.random()
            dt = ln(1/r1)/a0_sum
            # find j:
            breaking_point = r2*a0_sum
            j = 0 # python 0-based index: a[0] = j_1
            sum_j = 0
            while sum_j < breaking_point:
                sum_j += a[j]
                j += 1
            reaction = reactions[j]

            # 3. Effect the next reaction by replacing t ← t + τ and x̄ ← x̄ + νj.
            self.simulation_time += dt
            # dom_spec, c_j, is_hybridizing = possible_reactions_lst[j]
            # Edit, is indexed by doms_spec
            # doms_specs, (c_j, is_hybridizing) in self.possible_hybridization_reactions.items()
            c_j, is_hybridizing =



            # Update and re-calculate:
            # - self.possible_hybridization_reactions
            # - self.propencity_functions


            ## FINISH OF AND LOOP TO STEP (1)
            n_done += 1
            self.N_steps += 1
            if n_done % 10000 == 0:
                print("Simulated %s of %s steps at T=%s K (%0.0f C). %s state changes with %s selections in %s total steps." % \
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))



