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

LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nacent_sim

from nascent.nacent_sim.dom_anneal_sim import Simulator
from nascent.nacent_sim.dom_utils import parse_strand_domains_file, check_strands




if __name__ == '__main__':


    #strand_defs_file = os.path.join(os.path.dirname(__file__), "testfiles", "strand_defs01.txt")
    # LIBPATH,
    strand_defs_file = os.path.join(os.path.dirname(nascent.nacent_sim.__file__), "testfiles", "strand_defs01.txt")
    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file)
    check_strands(input_oligos)

    outputstatsfnfmt = "~/nascent_dom_anneal_test_{uid}.txt"
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

    #
    nl = 1e-15  # If volume = 1 nl, then each domain is approx 1 nM.
    al = 1e-18  # If volume = 1 al, then each domain is approx 1 uM.
    adhoc_params = {'probablity_oversampling_factor': 100, # Increase chance of selection or dehybridization by 1000
                    'print_statsline_when_saving': True,
                   }
    simulator = Simulator(volume=al, strands=input_oligos, params=adhoc_params,
                          outputstatsfiles=outputstatsfile, verbose=0)

    #simulator.anneal(T_start=273+90, T_finish=273+20, n_steps_per_T=100000)
    try:
        offset = 0.3
        start = 80  # deg C
        stop = 40   # deg C
        start, stop = stop, start   # invert
        step = -1 if start > stop else +1   # deg C
        simulator.anneal(T_start=273+start+offset, T_finish=273+stop+offset, delta_T=step, n_steps_per_T=50000)
    except KeyboardInterrupt:
        simulator.save_stats_cache()


"""

Debug errors:

 * Forgot to set domain.Partner = None when dehybridizing (nascent_dom_anneal_test.txt)

 * I was pairing the same domains with each other, e.g. pairing H3A with H3A instead of h3a.

"""
