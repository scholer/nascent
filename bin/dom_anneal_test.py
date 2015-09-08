#!/usr/bin/env python
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

# pylint: disable=W0142,C0103,C0301,W0141


"""

Test dom_anneal_sim module.

"""

import os
import sys
import yaml
from datetime import datetime

LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nascent_sim
from nascent.nascent_sim.dom_anneal_sim import Simulator
from nascent.nascent_sim.dom_utils import parse_strand_domains_file, check_strands

N_AVOGADRO = 6.022e23   # /mol


def get_max_oversampling_and_min_n_steps(volume, n_domains,
                                         target_coverage=1,
                                         max_selection_prob=1e-3):
    domain_conc = 1/(N_AVOGADRO*volume)
    max_oversampling = max_selection_prob*N_AVOGADRO*volume # = max_selection_prob/c
    # coverage = self.N_steps_per_T * oversampling_factor * c / self.N_domains
    min_n_steps = target_coverage*n_domains*N_AVOGADRO*volume/max_oversampling
    return max_oversampling, min_n_steps


if __name__ == '__main__':


    #strand_defs_file = os.path.join(os.path.dirname(__file__), "testfiles", "strand_defs01.txt")
    # LIBPATH,
    strand_defs_folder = os.path.join(os.path.dirname(nascent.nascent_sim.__file__), "testfiles")
    structure = "duplex1"
    #structure = "duplex1_2"
    n_strand_copies_default = 100

    # Load strand def and check the strands:
    strand_defs_file = os.path.join(strand_defs_folder, "strand_defs_{}.txt".format(structure))
    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file, n_clones_default=n_strand_copies_default)
    check_strands(input_oligos)

    statsfolder = os.path.join(strand_defs_folder, "simdata", structure)
    if not os.path.exists(statsfolder):
        os.mkdir(statsfolder)
    elif not os.path.isdir(statsfolder):
        raise OSError("Statsfolder already exists and is not a folder:", statsfolder)

    # Ensure that we have a unique filename:
    outputstatsfnfmt = os.path.join(statsfolder, "dom_anneal_stats_{uid}.txt")
    #outputstatsfile = os.path.expanduser("~/nascent_dom_anneal_test.txt")
    uid = 0
    outputstatsfile = os.path.expanduser(outputstatsfnfmt.format(uid=uid))
    while os.path.exists(outputstatsfile):
        uid += 1
        outputstatsfile = os.path.expanduser(outputstatsfnfmt.format(uid=uid))
    print("Using outputstatsfile:", outputstatsfile)
    # Make sure to touch file (in case you are running multiple simulations)
    with open(outputstatsfile, 'a'):
        os.utime(outputstatsfile)


    # Some calculations:
    # * Consider a volume of 1 um^3 = (1e-5 dm)^3 = 1e-15 L = 1 femto-liter.
    # * For c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
    #     N_Avogagro = 6.02e23 /mol

    n_domains = sum(len(s.Domains) for s in input_oligos)
    nl = 1e-15  # If volume = 1 nl, then each domain is approx 1 nM.
    al = 1e-18  # If volume = 1 al, then each domain is approx 1 uM.
    volume = al
    volume = volume*n_strand_copies_default
    n_steps_per_T = 20000
    #n_steps_per_T = 100000
    n_steps_per_T = 500000
    oversampling_factor = 100*n_strand_copies_default
    # Increase chance of selection or dehybridization by 1000

    oversampling_max, n_steps_min = get_max_oversampling_and_min_n_steps(
        volume, n_domains, target_coverage=0.1, max_selection_prob=1e-3)

    # Use maximum oversampling and minimal n_steps_per_T:
    oversampling_factor, n_steps_per_T = int(oversampling_max), int(n_steps_min)

    print("oversampling_factor:", oversampling_factor, "(max is %s)" % int(oversampling_max))
    if oversampling_factor > oversampling_max:
        print(" - oversampling_factor IS TOO HIGH FOR RELIABLE SIMULATION!")
    print("n_steps_per_T:", n_steps_per_T, "(min is %s)\n" % int(n_steps_min))
    if n_steps_per_T < n_steps_min:
        print(" - n_steps_per_T IS TOO LOW FOR RELIABLE SIMULATION!")

    adhoc_params = {'probablity_oversampling_factor': oversampling_factor,
                    'print_statsline_when_saving': True,
                    'n_steps_per_T': n_steps_per_T
                   }
    simulator = Simulator(volume=volume, strands=input_oligos, params=adhoc_params,
                          outputstatsfiles=outputstatsfile, verbose=0)

    #simulator.anneal(T_start=273+90, T_finish=273+20, n_steps_per_T=100000)
    # TODO:
    # - Melting (vs annealing) curve
    # -
    try:
        offset = 0 #0.3
        start = 40  # deg C
        #start = 60  # deg C
        stop = 90   # deg C
        #start, stop = stop, start   # invert
        step = -1 if start > stop else +1   # deg C
        T_start, T_finish = [273+val+offset for val in (start, stop)]
        with open(outputstatsfile, 'w') as fp:
            yaml.dump({'T_start': T_start, 'T_finish': T_finish, 'delta_T': step,
                       'volume': volume,
                       'params': adhoc_params,
                       'strand_defs_file': strand_defs_file,
                       'uid': uid,
                       'n_strand_copies_default': n_strand_copies_default,
                       'n_steps_per_T': n_steps_per_T,
                       'datetime': datetime.now().isoformat()
                      }, fp)
        simulator.anneal(T_start=273+start+offset, T_finish=273+stop+offset, delta_T=step)
    except KeyboardInterrupt:
        print("\n\nABORT: KeyboardInterrupt.\n\n")
        simulator.save_stats_cache()
    Tm_filename = os.path.splitext(outputstatsfile)[0] + ".energies.yaml"
    with open(Tm_filename, 'w') as fp:
        yaml.dump(simulator.Domain_dHdS, fp)

    reportfilename = os.path.splitext(outputstatsfile)[0] + ".domainstats.txt"
    simulator.save_domain_stats(reportfilename)

    reportfilename = os.path.splitext(outputstatsfile)[0] + ".report.txt"
    simulator.save_report(reportfilename)

    print("Simulation complete; stats are available in file", outputstatsfile)

"""

Debug errors:

 * Forgot to set domain.Partner = None when dehybridizing (nascent_dom_anneal_test.txt)

 * I was pairing the same domains with each other, e.g. pairing H3A with H3A instead of h3a.

"""
