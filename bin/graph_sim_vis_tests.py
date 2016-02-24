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
import glob
from datetime import datetime
from code import interact
import traceback

#try:
#    import rlcompleter
#    import readline
#except ImportError:
#    pass

LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nascent_sim
from nascent.graph_sim_nx.simulator_dm import DM_Simulator
from nascent.nascent_sim.dom_utils import parse_strand_domains_file, check_strands

N_AVOGADRO = 6.022e23   # /mol



def touch(filepath):
    with open(filepath, 'a'):
        os.utime(filepath)


def test(simulator):
    """
    Perform ad-hoc testing of simulator system.
    """

    #local = {'simulator': simulator}
    #ns = locals()
    readline.set_completer(rlcompleter.Completer(locals()).complete)
    readline.parse_and_bind("tab: complete")
    from importlib import reload
    from nascent.nascent_sim.moves import hybridize, dehybridize
    from nascent.nascent_sim import moves
    domA = next(d for d in simulator.Domains if d.Name == 'A')
    domB = next(d for d in simulator.Domains if d.Name == 'B')
    domC = next(d for d in simulator.Domains if d.Name == 'C')
    domb = next(d for d in simulator.Domains if d.Name == 'b')
    domD = next(d for d in simulator.Domains if d.Name == 'D')
    domE = next(d for d in simulator.Domains if d.Name == 'E')
    domF = next(d for d in simulator.Domains if d.Name == 'F')
    domc = next(d for d in simulator.Domains if d.Name == 'c')
    domG = next(d for d in simulator.Domains if d.Name == 'G')
    domH = next(d for d in simulator.Domains if d.Name == 'H')
    domg = next(d for d in simulator.Domains if d.Name == 'g')
    domK = next(d for d in simulator.Domains if d.Name == 'K')
    domi = next(d for d in simulator.Domains if d.Name == 'i')
    domL = next(d for d in simulator.Domains if d.Name == 'L')
    newc, oldc = moves.hybridize(domb, domB)
    interact(local=locals())



def cleanup(outputstatsfile):
    """ Remove all outputfiles related to outputstatsfile. """
    os.remove(outputstatsfile)  # clean up, if we are not saving.
    # TODO: You should remove all files explicitly. Some datafiles are appended, not overwritten.
    for fn in glob.glob(os.path.splitext(outputstatsfile)[0]+"*"):
        print("- Removing", fn)
        os.remove(fn)



def main():
    #strand_defs_file = os.path.join(os.path.dirname(__file__), "testfiles", "strand_defs01.txt")
    # LIBPATH,
    do_testing = False
    #do_testing = True
    strand_defs_folder = os.path.join(os.path.dirname(nascent.nascent_sim.__file__), "testfiles")
    #structure = "duplex1"
    #structure = "duplex1_2"
    #structure = "circ1"
    #structure = "circfb1a"
    #structure = "lin3s1"
    #structure = "polymer1h2s1"
    # TODO: The above structure should polymerize, but does not because the model lacks awareness of structural rigidity.
    structure = "catenane"

    #n_strand_copies_default = 400
    #n_strand_copies_default = 100
    #n_strand_copies_default = 40
    n_strand_copies_default = 10

    if do_testing:
        structure = "DistTest1"
        n_strand_copies_default = 1

    # Load strand def and check the strands:
    strand_defs_file = os.path.join(strand_defs_folder, "strand_defs_{}.txt".format(structure))
    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file, n_clones_default=n_strand_copies_default)
    check_strands(input_oligos)
    n_domains = sum(len(s.Domains) for s in input_oligos)

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


    # Some calculations:
    # * Consider a volume of 1 um^3 = (1e-5 dm)^3 = 1e-15 L = 1 femto-liter.
    # * For c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
    #     N_Avogagro = 6.02e23 /mol


    ### ADD SIMULATION NOTE/COMMENT/MESSAGE ###
    note = ""
    #note = "Explicit water during selection of domain2 - lower oversampling"
    #note = "Back to standard 1 al volume = 1.67 uM of each individual strand."
    #note = "First simulation with new probability model (repeat, now with r_hyb). Still selection based. Setting oversampling to 1."
    #note = "Second simulation with new probability model. Removing random.random() > math.exp(deltaG/(R*T)) before hybridizing selected strands."
    #note = "Third simulation with new probability model. Adding *0.5 to r_hyb to account for the fact that for every duplex we have two domains."
    #note = "Decreasing oversampling and n_strands."
    #note = "Adding an arbitrary factor 10 to p_hyb: p_hyb_old = 1 / (1 + math.exp(dG_std/(R*T))*self.Oversampling*10/compl_activity)."
    #note = "Still that arbitrary factor 10 to p_hyb, but ramping up oversampling to 2000."
    # That factor 10 helped... but why? Maybe I calculate the wrong dG_std?
    note = ("Re-introducing the intra-strand entropy penalty (only +2 cal/mol/K).")
    nl = 1e-15  # If volume = 1 nl, then each domain is approx 1.67 nM.
    al = 1e-18  # If volume = 1 al, then each domain is approx 1.67 uM.
    volume = al # / 1000
    volume = volume*n_strand_copies_default     # If you want to keep the strand concentration constant

    auto_oversample_and_steps = False

    # Increase chance of selection or dehybridization by 1000
    oversampling_max, n_steps_min = get_max_oversampling_and_min_n_steps(
        volume, n_domains, target_coverage=0.1, max_selection_prob=1e-4)


    if auto_oversample_and_steps:
        # Use maximum oversampling and minimal n_steps_per_T:
        oversampling_factor, n_steps_per_T = int(oversampling_max), int(n_steps_min)
    else:
        #n_steps_per_T = 20000
        n_steps_per_T = 100000
        #n_steps_per_T = 200000
        #n_steps_per_T = 400000
        #n_steps_per_T = 1000000
        #n_steps_per_T = 2000000
        #n_steps_per_T = 500000
        # oversampling_factor = 100*n_strand_copies_default
        oversampling_factor = 5000


    print("oversampling_factor:", oversampling_factor, "(max is %s)" % int(oversampling_max))
    if oversampling_factor > oversampling_max:
        print(" - oversampling_factor IS TOO HIGH FOR RELIABLE SIMULATION!")
    print("n_steps_per_T:", n_steps_per_T, "(min is %s)\n" % int(n_steps_min))
    if n_steps_per_T < n_steps_min:
        print(" - n_steps_per_T IS TOO LOW FOR RELIABLE SIMULATION!")


    ## SET UP SIMULATOR
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
    # Make and save thermodynamic melting curve:

    # Perform calculations and start simulation
    if do_testing:
        test(simulator)
    try:
        offset = 0 #0.3
        start = 40  # deg C
        #start = 60  # deg C
        #stop = 80   # deg C
        stop = 70   # deg C
        start, stop = stop, start   # invert => start hot
        step = -1 if start > stop else +1   # deg C
        T_start, T_finish = [273+val+offset for val in (start, stop)]

        cum_stats, domain_stats = simulator.thermodynamic_meltingcurve(T_start, T_finish)
        print("cum_stats:", len(cum_stats))
        print("domain_stats:", len(domain_stats))
        meltingfile = os.path.splitext(outputstatsfile)[0] + ".thermo_melting.yaml"
        with open(meltingfile, 'w') as fp:
            yaml.dump(cum_stats, fp)
            print("Thermodynamic melting curve save to file:", meltingfile)
        meltingfile = os.path.splitext(outputstatsfile)[0] + ".thermo_domainstats.yaml"
        with open(meltingfile, 'w') as fp:
            yaml.dump(domain_stats, fp)
            print("Thermodynamic melting curve domain stats save to file:", meltingfile)

        #raise KeyboardInterrupt    # Can be used to by-pass annealing
        #touch(outputstatsfile)
        with open(outputstatsfile, 'w') as fp:
            yaml.dump({'T_start': T_start, 'T_finish': T_finish, 'delta_T': step,
                       'volume': volume,
                       'params': adhoc_params,
                       'strand_defs_file': strand_defs_file,
                       'uid': uid,
                       'n_strand_copies_default': n_strand_copies_default,
                       'n_steps_per_T': n_steps_per_T,
                       'datetime': datetime.now().isoformat(),
                       'note': note
                      }, fp)

        ## INITIATE ANNEALING
        simulator.anneal(T_start=T_start, T_finish=T_finish, delta_T=step)

    except KeyboardInterrupt:
        print("\n\nABORT: KeyboardInterrupt.\n\n")
        print(traceback.format_exc()) # or print(sys.exc_info()[0])
        answer = input("Do you want to save the simulation data for this aborted run? [yes]/no  ")
        if answer and answer[0].lower() == 'n':
            cleanup(outputstatsfile)
            return
    except Exception as e:
        print("\n\nException during siumlation:", e)
        print(traceback.format_exc()) # or print(sys.exc_info()[0])
        answer = input("\nDo you want to enter interactive mode?  ")
        if answer and answer[0].lower() == 'y':
            interact(local=locals())
        answer = input("Do you want to save the simulation data for this aborted run? [yes]/no  ")
        if answer and answer[0].lower() == 'n':
            cleanup(outputstatsfile)
            return

    simulator.save_stats_cache()

    Tm_filename = os.path.splitext(outputstatsfile)[0] + ".energies.yaml"
    with open(Tm_filename, 'w') as fp:
        yaml.dump(simulator.Domain_dHdS, fp)

    reportfilename = os.path.splitext(outputstatsfile)[0] + ".domainstats.txt"
    simulator.save_domain_stats(reportfilename)

    reportfilename = os.path.splitext(outputstatsfile)[0] + ".report.txt"
    simulator.save_report(reportfilename)

    print("Simulation complete; stats are available in file", outputstatsfile)


if __name__ == '__main__':
    main()


"""

Debug errors:

 * Forgot to set domain.Partner = None when dehybridizing (nascent_dom_anneal_test.txt)

 * I was pairing the same domains with each other, e.g. pairing H3A with H3A instead of h3a.

"""
