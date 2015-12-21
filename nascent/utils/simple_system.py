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

Module with simplified systems for testing various simulation aspects.

"""

from math import log, exp
ln = log
import random
#import numpy as np



def random_choice(choices, weights=None):
    """ Select a random element in :choices: with probabilities from :weights:. """
    if weights is None:
        #weights = [1 for _ in range(len(choices))]
        return random.choice(choices)
    r2 = random.random() # default is the Mersenne Twister PRNG.
    breaking_point = r2*sum(weights)
    j = 0
    sum_j = weights[j]
    while sum_j < breaking_point:
        j += 1
        sum_j += weights[j]
    return choices[j]



def minimal_ssa(rate_constants, N_steps_max=100000, time_max=None, start_state=None):
    """ Minimal stochastic simulation algorithm. rate_constants: dict of dicts, [start][end] = k(start->end) """
    time, N_steps, stats = 0, 0, []
    state = start_state if start_state is not None else random.choice(list(rate_constants.keys()))
    while N_steps < N_steps_max and (time_max is None or time < time_max):
        end_states, a = zip(*rate_constants[state].items())
        tau = ln(1/random.random())/sum(a)  # or np.random.exponential(1/sum(a))
        stats.append({'state': state, 'tau': tau, 'time': time})
        state = random_choice(end_states, weights=a)
        time += tau
        N_steps += 1
    return stats


def throttled_ssa(rate_constants, throttle_func=None, N_steps_max=100000, time_max=None, start_state=None):
    """ Minimal stochastic simulation algorithm with throttle. """
    time, N_steps, stats = 0, 0, []
    state_count = {k: 0 for k in rate_constants.keys()}
    reaction_count = {(source, target): 0
                      for source, targets in rate_constants.items() for target in targets.keys()}
    cycle_count = {k: 0 for k in rate_constants.keys()}

    if throttle_func is None:
        throttle_func = lambda state, end_state, state_counts, cycle_counts, reaction_count: 1
    state = start_state if start_state is not None else random.choice(list(rate_constants.keys()))
    while N_steps < N_steps_max and (time_max is None or time < time_max):
        #end_states, a = zip(*rate_constants[state].items())
        end_states, a = zip(*[(end_state, k*throttle_func(state, end_state, state_count, cycle_count, reaction_count))
                              for end_state, k in rate_constants[state].items()])
        tau = ln(1/random.random())/sum(a)  # or np.random.exponential(1/sum(a))
        end_state = random_choice(end_states, weights=a)
        stats.append({'state': state, 'tau': tau, 'time': time})
        state_count[state] += 1
        reaction_count[(state, end_state)] += 1
        state = end_state
        time += tau
        N_steps += 1
    state_times = {state: sum(step['tau'] for step in stats if step['state'] == state)
                   for state in rate_constants.keys()}
    result = dict(zip("stats, cycle_count, state_count, time, state, state_times".split(", "),
                      [stats, cycle_count, state_count, time, state, state_times]))
    return result


def simulate_n_times(N_sims, rates, throttle_func=None, N_steps_max=20000, start_state="", states=None):
    """
    Invoke throttled_ssa(rates, throttled_ssa, N_steps_max, start_state) :N_sims: times.
    """
    if states is None:
        states = sorted(rates.keys())
    Ks = []
    all_stats = []
    for i in range(N_sims):
        stats = throttled_ssa(rates, throttle_func=throttle_func, N_steps_max=N_steps_max, start_state=start_state)
        all_stats.append(stats)
        state_partitions = {state: sum(step['tau'] for step in stats if step['state'] == state)
                            for state in states}
        K = sum(state_partitions[state] for state in states if state != "")/state_partitions[""]
        Ks.append(K)
        print(".", end="")
    print("")
    return all_stats, Ks


def cycle_test_1(k_AB, k1, k2, throttle_func=None, N_steps_max=1000):
    r"""
    Test throttling of the following cyclic reaction graph:

           ,----------.
     ,--- C1 --.       \
    /           \       \
  u              A ----  B
    \           /       /
     `--- C2 --´       /
           `----------´

    """
    stats = [] # list of dicts
    if throttle_func is None:
        throttle_func = lambda state, end_state, state_counts, cycle_counts, reaction_count: 1
    systime = 0
    cycle_count = {'N_AB': 0, 'N1': 0, 'N2': 0}
    rate_constants = {'A': {'B': k_AB},
                      'B': {'C1': k1, 'C2': k2},
                      'C1': {'A': 1},
                      'C2': {'A': 3},
                     }
    state_count = {k: 0 for k in rate_constants.keys()}
    reaction_count = {(source, target): 0
                      for source in rate_constants.keys()
                      for target in rate_constants[source].keys()}
    # use as k(A->B) = k['A']['B']

    state = "A"

    for _ in range(N_steps_max):
        r1, r2 = random.random(), random.random()
        a, end_states = zip(*[(k*throttle_func(state, end_state, state_count, cycle_count, reaction_count), end_state)
                              for end_state, k in rate_constants[state].items()])
        a0_sum = sum(a)
        # determine tau:
        tau = ln(1/r1)/a0_sum
        # determine reaction:
        breaking_point = r2*a0_sum
        j = 0
        sum_j = a[j]
        while sum_j < breaking_point:
            j += 1
            sum_j += a[j]

        state_count[state] += 1
        reaction_count[(state, end_states[j])] += 1
        if state == 'C1':
            cycle_count['N1'] += 1
        if state == 'C2':
            cycle_count['N2'] += 1
        if state == 'B':
            cycle_count['N_AB'] += 1

        step = {'state': state,
                'tau': tau,
                'systime': systime,
                'N1': cycle_count['N1'],
                'N2': cycle_count['N2']}

        systime += tau
        state = end_states[j]

        stats.append(step)

    state_times = {state: sum(step['tau'] for step in stats if step['state'] == state)
                   for state in rate_constants.keys()}
    result = dict(zip("stats, cycle_count, state_count, systime, state".split(", "),
                      [stats, cycle_count, state_count, systime, state]))
    result['state_times'] = state_times
    print("Stochastic equilibrium constant K = [B]/[A] = %0.05f" % (state_times['B']/state_times['A']))
    print(" -expected equilibrium constant K=k/(k1+k2) = %0.05f" % (k_AB/(k1+k2)))
    print("State counts:")
    print(state_count)
    print("Cycle counts:")
    print(cycle_count)

    return result


def state_partition_times(stats, states):
    """ Return the total time spent in each state. """
    return {state: sum(step['tau'] for step in stats if step['state'] == state) for state in states}

def plot_time_partitions(stats, states):
    """
    Plot state partition times.
    """
    from matplotlib import pyplot
    # state_partitions = state_times(stats, states)
    if states is None:
        states = sorted(set(step['state'] for step in stats))
    state_partitions = {state: sum(step['tau'] for step in stats if step['state'] == state) for state in states}
    times = [state_partitions[state] for state in states]
    indices = range(len(states))
    barwidth = 0.5
    ax = pyplot.bar(indices, times, width=barwidth, alpha=0.5)
    ticks = pyplot.xticks([i+barwidth/2 for i in indices], states)
    xlims = pyplot.xlim(xmin=-0.1)
    title = pyplot.title('Cumulative state time')
