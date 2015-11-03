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

import sys
import random
#from collections import defaultdict
from collections import deque
#import math
from math import log as ln # exp
#from datetime import datetime
import time
if sys.platform == "win32":
    timer = time.clock
else:
    timer = time.time


#import networkx as nx
#import numpy as np


#from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
#from .thermodynamic_utils import thermodynamic_meltingcurve
from .simulator import Simulator
# N_AVOGADRO in /mol,
# Universal Gas constant in cal/mol/K
#from .constants import N_AVOGADRO, R

# Module-level constants and variables
VERBOSE = 0




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
        - propensity_functions[doms_specs] = actual propencity for reaction between domain1 and domain2
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
        self.timings = {}  # Performance profiling of the simulator (real-world or cpu time)
        self.system_stats['tau_deque'] = deque(maxlen=10)
        self.print_post_step_fmt = ("\r{self.N_steps: 5} "
                                    "[{self.sim_system_time:0.02e}"
                                    " {stats[tau]:0.02e}"
                                    " {stats[tau_mean]:0.02e}] "
                                    "[{timings[step_time]:0.02e}] "
                                    "[{timings[step_time]:0.02e}] "
                                    )


    def simulate(self, T=None, n_steps_max=100000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        self.timings['simulation_start_time'] = timer()
        self.timings['step_start_time'] = self.timings['step_end_time'] = self.timings['simulation_start_time']

        sysmgr = self.systemmgr
        if T is None:
            T = sysmgr.temperature
        else:
            sysmgr.temperature = T

        n_done = 0

        while n_done < n_steps_max:

            # sysmgr.possible_hybridization_reactions is dict with  {doms_specs: (c_j, is_hybridizing)}
            # sysmgr.propensity_functions is dict with              {doms_specs: a_j}

            # Step 1. (I expect that propensity_functions is up-to-date)
            reactions, propensity_functions = zip(sysmgr.propensity_functions.items())
            # reactions: list of possible reactions
            a = propensity_functions # propensity_functions[j]: propencity for reaction[j]
            a0_sum = sum(propensity_functions)

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
            # Propensity functions are only used to select which reaction path to fire.
            # Now that we have that, we just have to select an domain pair with doms_specs
            reaction_doms_specs = reactions[j]  # is a doms_specs {F₁, F₂}

            # 3. Effect the next reaction by updating time and counts (replacing t ← t + τ and x̄ ← x̄ + νj).
            # 3a: Update system time.
            # 3b: Hybridize/dehybridize domains and Update graphs/stats/etc

            # 3a: Update system time:
            self.sim_system_time += dt

            # 3b: Hybridize/dehybridize:
            c_j, is_hybridizing = sysmgr.possible_hybridization_reactions[reaction_doms_specs]
            d1, d2 = [species_list.pop() for species_list in
                      [sysmgr.domain_state_subspecies[dom_spec] for dom_spec in reaction_doms_specs]]
            if is_hybridizing:
                changed_complexes, new_complexes, obsolete_complexes, free_strands = sysmgr.hybridize(d1, d2)
            else:
                changed_complexes, new_complexes, obsolete_complexes, free_strands = sysmgr.dehybridize(d1, d2)

            # 4: Update/re-calculate possible_hybridization_reactions and propensity_functions
            # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
            # - possible_hybridization_reactions
            # - propensity_functions
            changed_domains = [domain for cmplx in changed_complexes for domain in cmplx.nodes()
                               if not domain.partner]+\
                              [domain for strand in free_strands for domain in strand.domains if not domain.partner]
            if is_hybridizing:
                # d1, d2 are partnered and not included in changed_domains (which filtered out hybridized domains)
                changed_domains.extend([d1, d2])
            sysmgr.update_possible_reactions(changed_domains, d1, d2, is_hybridizing)


            ## FINISH OF AND LOOP TO STEP (1)
            n_done += 1
            self.N_steps += 1
            self.system_stats['tau'] = dt
            self.system_stats['tau_deque'].append(dt)
            self.system_stats['tau_mean'] = sum(self.system_stats['tau_deque'])/len(self.system_stats['tau_deque'])
            self.timings['start_time'], self.timings['end_time'] = self.timings['end_time'], timer()
            self.timings['step_time'] = self.timings['start_time'] - self.timings['end_time']

            print(self.print_post_step_fmt.format(
                self=self, stats=self.system_stats, timings=self.timings, sysmgr=self.systemmgr),
                  end="")
            if n_done % 10000 == 0:
                print(("Simulated %s of %s steps at T=%s K (%0.0f C). "+
                       "%s state changes with %s selections in %s total steps.") %
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))




