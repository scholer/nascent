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

# pylint: disable=C0103,W0212,W0142,W0141

"""

Module for

"""

import sys
import random
#from collections import defaultdict
from collections import deque
import copy
#import math
from math import log as ln # exp
import time
from datetime import datetime
from pprint import pprint
import pdb
import time
if sys.platform == "win32":
    timer = time.clock
else:
    timer = time.time

#import networkx as nx
import numpy as np


#from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
#from .thermodynamic_utils import thermodynamic_meltingcurve
from .simulator import Simulator
from .reactionmgr_grouped import ReactionMgrGrouped
# N_AVOGADRO in /mol,
# Universal Gas constant in cal/mol/K
from .constants import (N_AVOGADRO, R,
                        HYBRIDIZATION_INTERACTION,
                        PHOSPHATEBACKBONE_INTERACTION,
                        STACKING_INTERACTION)
from .debug import printd, pprintd

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
        domspec_pair = frozenset({domspec1, domspec2})

    Primary data structures:
        - possible_hybridization_reactions[domspec_pair] = propensity_constant cj, is_hybridizing
            Contains up-to-date hybridization and dehybridization reactions
            dict, keyed by domspec_pair = {domain1-state, domain2-state}
            This is "R" and "R_j" in Gillespie's formulas
        - propensity_functions[domspec_pair] = actual propensity for reaction between domain1 and domain2
            with domain-states domspec_pair = {domain1-state, domain2-state}
            This is "a" or "a(x)" in Gillespie's formulas.


    Overview of the different caches:
    To avoid doing a lot of re-calculation, I generally save all results.
    If the result involves only one domain, I store it as domspec for that domain.
    If the result involves two domains, I store it under the domspec_pair = {domain1-domspec, domain2-domspec}


    1-domain caches, state-dependent (domspec):
        - domain activity for inter-complex reactions.
            If the domains are not in part of the same complex, I can usually just take a "quick look around"
            and estimate steric interactions.

    2-domain caches, state-dependent (domspec_pair):
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
        ### Set up a "grouped" system manager for testing: ###
        # Make sure the strands, etc are all copied before you make a new system:
        strands = kwargs.pop('strands')
        if strands is not None:
            strands = copy.deepcopy(strands)
        self.reactionmgr_grouped = ReactionMgrGrouped(params=params, strands=strands, **kwargs)
        ### -- end making grouped reaction mgr --- ###
        self.time_per_T = params.get('time_per_T', 10) # seconds
        self.timings = {}  # Performance profiling of the simulator (real-world or cpu time)
        self.system_stats['tau_deque'] = deque(maxlen=10)
        self.step_sleep_factor = self.params.get('simulator_step_sleep_factor', 0)
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

        sysmgr = self.reactionmgr
        sysmgr_grouped = self.reactionmgr_grouped
        if T is None:
            T = sysmgr.temperature
        else:
            sysmgr.temperature = T
        sysmgr_grouped.temperature = T

        n_done = 0

        print("sysmgr.domain_pairs:")
        pprintd(sysmgr.domain_pairs)

        while n_done < n_steps_max:

            print("\n\n------ Step %03s -----------\n" % n_done)
            sysmgr.check_system()
            print(" - precheck OK.")

            # sysmgr.possible_hybridization_reactions is dict with  {domspec_pair: (c_j, is_hybridizing)}
            # sysmgr.propensity_functions is dict with              {domspec_pair: a_j}
            printd(" ---- sysmgr.print_reaction_stats(): ----")
            sysmgr.print_reaction_stats()

            # Step 1. (I expect that propensity_functions is up-to-date)
            if not sysmgr.propensity_functions:
                print("\n\nERROR: sysmgr.propensity_functions is:",
                      sysmgr.propensity_functions, " - ABORTING SIMULATION.\n\n")
                return
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
            # Now that we have that, we just have to select an domain pair with domspec_pair
            # is a domspec_pair = {dom_spec1, dom_spec2} = {F₁, F₂}
            # Edit: If you want to group reactions by domain species, you need:
            #   reaction spec = ({domspec1, domspec2}, is_hybridizing, is_intracomplex)
            reaction_spec = reaction_specs[j]

            # 3. Effect the next reaction by updating time and counts (replacing t ← t + τ and x̄ ← x̄ + νj).
            # 3a: Update system time.
            # 3b: Hybridize/dehybridize domains and Update graphs/stats/etc

            # 3a: Update system time:
            self.sim_system_time += dt

            # 3b: Hybridize/dehybridize:
            c_j = sysmgr.possible_hybridization_reactions[reaction_spec]
            # Determine reaction attributes:
            # We can currently distinguish grouped vs ungrouped based on whether reaction_spec is a tuple (vs frozenset)
            if isinstance(reaction_spec, tuple):
                is_hybridizing = reaction_spec[1]
                domspec_pair = tuple(reaction_spec[0])
            else:
                # Non-grouped, reaction_spec is just a domain_pair = frozenset((domain1, domain2))
                domain_pair = reaction_spec
                d1, d2 = tuple(domain_pair)
                if d1.partner == d2:
                    is_hybridizing = False  # already hybridized; must can only undergo de-hybridization reaction
                    is_intra = True
                    assert d2.partner == d1
                else:
                    is_hybridizing = True
                    is_intra = (d1.strand == d2.strand) or \
                               (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex)
                    assert d1.partner is None and d2.partner is None

            reaction_attr = sysmgr.reaction_attrs[domain_pair]
            if reaction_attr.is_hybridizing != is_hybridizing or reaction_attr.is_intra != is_intra:
                print("Expected reaction attr", reaction_attr,
                      "does not match actual reaction attributes is_hybridizing=%s, is_intra=%s" %
                      (is_hybridizing, is_intra))

            domspec_pair = frozenset((d1.domain_state_fingerprint(), d2.domain_state_fingerprint()))
            grouped_reaction_spec = (domspec_pair, is_hybridizing, is_intra)
            # reverted to reaction_spec = domain_pair for ungrouped reactions.
            # A domain instance should only have ONE way (is_hyb?, is_intra?) to react with another domain instance.
            dom1, dom2 = sysmgr.react_and_process(domain_pair, is_hybridizing, is_intra)
            sysmgr.check_system()
            print(" - post reaction_and_process check OK.")
            assert (d1, d2) == (dom1, dom2)

            # Repeat for grouped system mgr:
            mirror_grouped_system = False
            if mirror_grouped_system:
                print("\nRepeating reaction for grouped sysmgr; reaction_spec is:")
                pprint(reaction_spec)
                print("grouped_reaction_spec is:")
                pprint(grouped_reaction_spec)
                print("---- sysmgr reaction pathways: ----")
                sysmgr.print_reaction_stats()
                print("---- sysmgr_grouped reaction pathways: ----")
                sysmgr_grouped.print_reaction_stats()
                c_j_g = sysmgr_grouped.possible_hybridization_reactions[grouped_reaction_spec]
                print("c_j, is_hybridizing:", c_j, is_hybridizing)
                print("c_j_g:", c_j_g)
                assert np.isclose(c_j_g, c_j)
                print("React and process using grouped_reaction_spec %s:" % (grouped_reaction_spec,))
                d1_g, d2_g = sysmgr_grouped.react_and_process(*grouped_reaction_spec)
                # the grouped sysmgr might select different domains for hybridizing, but they should
                # have the same state, and the domspec distribution should be equivalent.
                try:
                    assert (d1_g.domain_state_fingerprint() == d1.domain_state_fingerprint() and \
                            d2_g.domain_state_fingerprint() == d2.domain_state_fingerprint()) or \
                           (d1_g.domain_state_fingerprint() == d2.domain_state_fingerprint() and \
                            d2_g.domain_state_fingerprint() == d1.domain_state_fingerprint())
                except AssertionError as e:
                    print("\nReacted domain species for grouped system does not match that for ungrouped:")
                    print("d1.domain_state_fingerprint():", d1.domain_state_fingerprint())
                    print("d1_g.domain_state_fingerprint():", d1_g.domain_state_fingerprint())
                    print("d2.domain_state_fingerprint():", d2.domain_state_fingerprint())
                    print("d2_g.domain_state_fingerprint():", d2_g.domain_state_fingerprint())
                    print("c_j, is_hybridizing:", c_j, is_hybridizing)
                    print("c_j_g:", c_j_g)
                    raise e


            ## Post-hybridization (and processing) assertions:
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
                # Dispatcher is responsible for translating domain objects to model-agnostic (str repr) graph events.
                directive['nodes'] = (d1, d2)
                self.dispatcher.dispatch(directive, directive_is_list=False)


            ## FINISH OF AND LOOP TO STEP (1)
            n_done += 1
            self.N_steps += 1
            self.system_stats['tau'] = dt
            self.system_stats['tau_deque'].append(dt)
            self.system_stats['tau_mean'] = sum(self.system_stats['tau_deque'])/len(self.system_stats['tau_deque'])
            # step_start_time: Real-world time when step was started.
            if self.step_sleep_factor:
                self.timings['step_end_time'] = timer()
            else:
                # Only call timer() once.
                self.timings['step_start_time'], self.timings['step_end_time'] = self.timings['step_end_time'], timer()
            self.timings['step_time'] = self.timings['step_start_time'] - self.timings['step_end_time']

            print(self.print_post_step_fmt.format(
                self=self, stats=self.system_stats, timings=self.timings, sysmgr=self.reactionmgr),
                  end="")

            if n_done % 10000 == 0:
                print(("Simulated %s of %s steps at T=%s K (%0.0f C). "+
                       "%s state changes with %s selections in %s total steps.") %
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))


            if self.step_sleep_factor:
                # Compensate for simulation calculation time:
                sleep_time = (dt-self.timings['step_time'])*self.step_sleep_factor
                if sleep_time > 0:
                    print("Sleeping for %s seconds before next step..." % sleep_time)
                    time.sleep(sleep_time)
                else:
                    print("Simulation step calculation time %s is larger than simulation reaction tau = %s" %
                          (self.timings['step_time'], dt))
                self.timings['step_start_time'] = timer()
            answer = 'c' or input(("\nReaction complete."
                            "Type 'd' to enter debugger;"
                            "'g' to save plot to file;"
                            "'q' to quit;"
                            "any other key to continue... "))
            if 'g' in answer:
                sysmgr.draw_and_save_graphs(n="%s_%0.03fs" % (n_done, self.sim_system_time))
                if mirror_grouped_system:
                    sysmgr_grouped.draw_and_save_graphs(prefix="grouped",
                                                        n="%s_%0.03fs" % (n_done, self.sim_system_time))
            if 'q' in answer:
                return
            if 'd' in answer:
                pdb.set_trace()


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
