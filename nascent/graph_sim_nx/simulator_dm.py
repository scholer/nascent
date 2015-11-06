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

# pylint: disable=C0103,W0212

"""

Module for

"""

import sys
import random
#from collections import defaultdict
from collections import deque
#import math
from math import log as ln # exp
from datetime import datetime
from pprint import pprint
import pdb
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
from .constants import (N_AVOGADRO, R,
                        HYBRIDIZATION_INTERACTION,
                        PHOSPHATEBACKBONE_INTERACTION,
                        STACKING_INTERACTION)

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
        - possible_hybridization_reactions[doms_specs] = propensity_constant cj, is_hybridizing
            Contains up-to-date hybridization and dehybridization reactions
            dict, keyed by doms_specs = {domain1-state, domain2-state}
            This is "R" and "R_j" in Gillespie's formulas
        - propensity_functions[doms_specs] = actual propensity for reaction between domain1 and domain2
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


    def __init__(self, params, **kwargs):
        #super().__init__(**kwargs)
        Simulator.__init__(self, params=params, **kwargs)
        self.time_per_T = params.get('time_per_T', 10) # seconds
        self.timings = {}  # Performance profiling of the simulator (real-world or cpu time)
        self.system_stats['tau_deque'] = deque(maxlen=10)
        self.print_post_step_fmt = ("\r{self.N_steps: 5} "
                                    "[{self.sim_system_time:0.02e}"
                                    " {stats[tau]:0.02e}"
                                    " {stats[tau_mean]:0.02e}] "
                                    "[{timings[step_time]:0.02e}] "
                                    #"[{timings[step_time]:0.02e}] "
                                   )
        # Must specify: time, tau, T, nodes
        self.state_change_hybridization_template = {'change_type': 1,
                                                    'forming': 1,
                                                    'interaction': HYBRIDIZATION_INTERACTION,
                                                    'multi': False}
        self.state_change_dehybridization_template = {'change_type': 1,
                                                      'forming': 0,
                                                      'interaction': HYBRIDIZATION_INTERACTION,
                                                      'multi': False}


    def simulate(self, T=None, n_steps_max=100000, simulation_time=None, systime_max=None):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        self.timings['simulation_start_time'] = timer()
        self.timings['step_start_time'] = self.timings['step_end_time'] = self.timings['simulation_start_time']
        if systime_max is None:
            systime_max = self.sim_system_time + simulation_time

        sysmgr = self.systemmgr
        if T is None:
            T = sysmgr.temperature
        else:
            sysmgr.temperature = T

        n_done = 0

        print("sysmgr.domain_pairs:")
        pprint(sysmgr.domain_pairs)

        while n_done < n_steps_max:

            print("\n\n------ Step %03s -----------\n" % n_done)

            # sysmgr.possible_hybridization_reactions is dict with  {doms_specs: (c_j, is_hybridizing)}
            # sysmgr.propensity_functions is dict with              {doms_specs: a_j}

            # Step 1. (I expect that propensity_functions is up-to-date)
            if not sysmgr.propensity_functions:
                print("\n\nERROR: sysmgr.propensity_functions is:",
                      sysmgr.propensity_functions, " - ABORTING SIMULATION.\n\n")
            reaction_specs, propensity_functions = zip(*sysmgr.propensity_functions.items())
            assert set(reaction_specs) == set(sysmgr.possible_hybridization_reactions.keys())
            # reaction_specs: list of possible reaction_specs
            a = propensity_functions # propensity_functions[j]: propensity for reaction[j]
            a0_sum = sum(propensity_functions)

            # Step 2: Generate values for τ and j:  - easy.
            r1, r2 = random.random(), random.random()
            dt = ln(1/r1)/a0_sum
            if systime_max and self.sim_system_time + dt > systime_max:
                self.sim_system_time = systime_max
                print("\nSimulation system time reached systime_max = %s s !\n\n" % systime_max)
                return n_done
            # find j:
            breaking_point = r2*a0_sum
            j = 0 # python 0-based index: a[0] = j_1
            sum_j = a[j]
            while sum_j < breaking_point:
                j += 1
                sum_j += a[j]
            # [0.25, 0.25, 0.25, 0.25]
            # [0.25, 0.50, 0.75, 1.00]
            # bp = 0.9
            # j =
            # Propensity functions are only used to select which reaction path to fire.
            # Now that we have that, we just have to select an domain pair with doms_specs
            # is a doms_specs = {dom_spec1, dom_spec2} = {F₁, F₂}
            # Edit: If you want to group reactions by domain species, you need:
            #   reaction spec = ({domspec1, domspec2}, is_hybridizing, is_intracomplex)
            reaction_spec = reaction_specs[j]

            # 3. Effect the next reaction by updating time and counts (replacing t ← t + τ and x̄ ← x̄ + νj).
            # 3a: Update system time.
            # 3b: Hybridize/dehybridize domains and Update graphs/stats/etc

            # 3a: Update system time:
            self.sim_system_time += dt

            # 3b: Hybridize/dehybridize:
            c_j, is_hybridizing = sysmgr.possible_hybridization_reactions[reaction_spec]
            d1, d2 = sysmgr.react_and_process(reaction_spec, is_hybridizing)

            ## Post-hybridization assertions:
            try:
                # assert set(reaction_specs) == set(sysmgr.possible_hybridization_reactions.keys())
                # After processing, we do not expect the old reaction_specs to be the same as the new...
                assert set(sysmgr.propensity_functions.keys()) == set(sysmgr.possible_hybridization_reactions.keys())
            except AssertionError as e:
                print("\n\n", repr(e), sep="")
                print("set(reaction_specs):")
                pprint(set(reaction_specs))
                print("set(sysmgr.possible_hybridization_reactions.keys()):")
                pprint(set(sysmgr.possible_hybridization_reactions.keys()))
                print("set(sysmgr.propensity_functions.keys()):")
                pprint(set(sysmgr.propensity_functions.keys()))
                raise e

            # 3c: Dispatch the state change
            if self.dispatcher:
                # The state_change received by dispatch() is either a dict with keys:
                #     - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
                #     - forming: 1=forming, 0=eliminating
                #     - interaction: 1=backbone, 2=hybridization, 3=stacking (only for edge changes)
                #     - time: Current system time (might be subjective to rounding errors).
                #     - tau:  Change in system time between this directive and the former.
                #     - T:    Current system temperature.
                #     - multi: If True, interpret "nodes" as a list of pairs (u₁, u₂, ...) for NODE-altering directives,
                #                 or (u₁, v₁, u₂, v₂, ...) for EDGE-altering directives.
                #     - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
                directive = self.state_change_hybridization_template.copy()
                directive['forming'] = int(is_hybridizing)
                directive['time'] = self.sim_system_time
                directive['T'] = sysmgr.temperature
                directive['tau'] = dt
                directive['nodes'] = (d1, d2)  # Not sure if we have to give str representation here
                self.dispatcher.dispatch(directive, directive_is_list=False)


            answer = 'g' or input(("\nReaction complete. Type 'd' to enter debugger; 'g' to save plot to file; "
                                   "any other key to continue... "))
            if 'g' in answer:
                sysmgr.draw_and_save_graphs(n="%s_%0.03fs" % (n_done, self.sim_system_time))
            if 'd' in answer:
                pdb.set_trace()

            ## FINISH OF AND LOOP TO STEP (1)
            n_done += 1
            self.N_steps += 1
            self.system_stats['tau'] = dt
            self.system_stats['tau_deque'].append(dt)
            self.system_stats['tau_mean'] = sum(self.system_stats['tau_deque'])/len(self.system_stats['tau_deque'])
            self.timings['step_start_time'], self.timings['step_end_time'] = self.timings['step_end_time'], timer()
            self.timings['step_time'] = self.timings['step_start_time'] - self.timings['step_end_time']

            print(self.print_post_step_fmt.format(
                self=self, stats=self.system_stats, timings=self.timings, sysmgr=self.systemmgr),
                  end="")

            if n_done % 10000 == 0:
                print(("Simulated %s of %s steps at T=%s K (%0.0f C). "+
                       "%s state changes with %s selections in %s total steps.") %
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))

        # end while loop
        return n_done


    def anneal(self, T_start, T_finish, delta_T=-1, n_steps_per_T=None, time_per_T=None):
        """
        Simulate annealing repeatedly from T_start to T_finish,
        decreasing temperature by delta_T for every round,
        doing either:
            (a) at most n_steps number of steps at each temperature, or
            (b) at most time_per_T simulation time at each temperature.

        """

        # Range only useful for integers...
        n_steps_per_T = n_steps_per_T or self.N_steps_per_T
        time_per_T = time_per_T or self.time_per_T
        T = T_start
        assert delta_T != 0
        assert T_start > T_finish if delta_T < 0 else T_finish > T_start
        print(("\nStarting annealing - ramp is %s K to %s K in %s K increments, "
               "using %s steps or %s seconds simulation time at each temperature...") %
              (T_start, T_finish, delta_T, n_steps_per_T, time_per_T))
        while T >= T_finish if delta_T < 0 else T <= T_finish:
            print("\nSimulating at %s K for %s steps or %s seconds simulation time..." % \
                  (T, n_steps_per_T, time_per_T))
            self.simulate(T, n_steps_max=n_steps_per_T, simulation_time=time_per_T)
            T += delta_T
            self.save_stats_cache() # Save cache once per temperature
        print("Annealing complete! (%s)" % datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.print_setup()
