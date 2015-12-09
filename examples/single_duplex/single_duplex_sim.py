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

Run a single simulation of a single duplex.

"""

import os
import sys
import yaml
import glob
from datetime import datetime
from code import interact
import traceback
from pprint import pprint as org_pprint
from collections import defaultdict, namedtuple
import pdb

import pprint
def info_pprint(*args, **kwargs):
    from inspect import currentframe, getframeinfo
    frameinfo = getframeinfo(currentframe().f_back)
    print(frameinfo.filename, frameinfo.lineno)
    org_pprint(*args, **kwargs)

pprint.pprint = info_pprint


# LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
# sys.path.insert(0, LIBPATH)
sys.path.insert(0, ".")
print("os.path.abspath('.'):", os.path.abspath('.'))

import nascent.nascent_sim
from nascent.graph_sim_nx.simulator_dm import DM_Simulator as Simulator
from nascent.graph_sim_nx.fileio import parse_strand_domains_file, check_strands
from nascent.graph_sim_nx.dispatcher import StateChangeDispatcher
from nascent.graph_sim_nx.constants import N_AVOGADRO #, R
from nascent.graph_sim_nx import debug

# Toggle debug printing:
debug.do_print = False


def touch(filepath):
    with open(filepath, 'a'):
        os.utime(filepath)


def main():
    """ Main test function. """
    enable_gephi_stream = False

    strand_defs_folder = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(nascent.__file__))), "testfiles")
    scriptdir = os.path.dirname(os.path.abspath(__file__))
    structure = "duplex2"

    #n_strand_copies_default = 400
    #n_strand_copies_default = 100
    # n_strand_copies_default = 40
    n_strand_copies_default = 20
    # n_strand_copies_default = 10
    # n_strand_copies_default = 2

    # Load strand def and check the strands:
    #strand_defs_file = os.path.join(strand_defs_folder, "strand_defs_{}.txt".format(structure))
    strand_defs_file = os.path.join(scriptdir, "{}.txt".format(structure))

    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file, n_clones_default=n_strand_copies_default)
    check_strands(input_oligos)
    n_domains = sum(len(s.domains) for s in input_oligos)

    run_directory = os.path.join(scriptdir,
                                 "simdata",
                                 structure,
                                 datetime.now().strftime("%Y-%m-%d %H%M%S"))
    if os.path.exists(run_directory) and not os.path.isdir(run_directory):
        print("Directory '" + run_directory + "' already exists BUT IS NOT A DIRECTORY. ABORTING!")
        raise OSError("Statsfolder/run_directory already exists and is not a folder:", run_directory)
    else:
        os.makedirs(run_directory)

    statsfolder = run_directory # os.path.join(strand_defs_folder, "simdata", structure)

    # Ensure that we have a unique filename:
    outputstatsfnfmt = os.path.join(statsfolder, "{statstype}.{ext}")
    def outputfn(**kwargs):
        return os.path.expanduser(outputstatsfnfmt.format(**kwargs))
    print("Using outputstatsfnfmt:", outputstatsfnfmt)


    # Some calculations:
    # * Consider a volume of 1 um^3 = (1e-5 dm)^3 = 1e-15 L = 1 femto-liter.
    # * For c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
    #     N_Avogagro = 6.02e23 /mol
    nl = 1e-15  # If volume = 1 nl, then each domain is approx 1.67 nM.
    al = 1e-18  # If volume = 1 al, then each domain is approx 1.67 uM.
    volume = al # / 1000
    volume = volume*n_strand_copies_default     # If you want to keep the strand concentration constant

    n_steps_per_T = 10000
    time_per_T = 100         # seconds

    print("n_steps_per_T:", n_steps_per_T)
    print("time_per_T:", time_per_T, "s")

    ## Dispatcher config:
    dispatcher_config = {'dispatcher_state_changes_fn': outputfn(statstype="graph_changes", ext="txt"),
                         'dispatcher_keep_file_open': True, # False,
                         'dispatcher_state_changes_file_add_header': True,
                         'dispatcher_cache_size': 100, # If output file is not kept open
                         # File output formatting:
                         # 'dispatcher_state_changes_line_fields': ['T', 'time', 'tau', 'change_type', 'forming', 'interaction']
                         #'dispatcher_state_changes_unpack_nodes': True,
                         #'dispatcher_state_changes_line_fmt': None,
                         # multi-directives are state change dicts with multiple nodes or multiple pairs, e.g.
                         # change['nodes'] = (node1, node2, ...) or (source1, target1, source2, target2, ...)
                         'dispatcher_multi_directive_support': True, # True, False
                         'dispatcher_graph_translation': 'domain-to-5p3p', # None, 'domain-to-5p3p', 'domain-to-strand'
                         #'livestreamer_graph_representation': '5p3p' # Not used.
                         'livestreamer_auto_apply_layout': 0,
                         'livestreamer_graph_layout_method': 'force-directed',
                         #None, 'gephi_mock', 'gephi_regular', 'gephi_multi'/'gephi', 'cytoscape'
                         'dispatcher_livestreamer': None, #'gephi_regular' if 0 else 'gephi_mock',
                         'visualization_use_websocket': False,  # For gephi or other websocket-enabled streamers.
                         'visualization_workspace': 'workspace0',
                         'dispatcher_debug_print': False, # print debug statements in dispatcher
                        }

    ## SET UP SIMULATOR
    params = {'time_per_T': time_per_T,
              'n_steps_per_T': n_steps_per_T,
              "working_directory": run_directory,
              # Sleep factor*tau after each DM simulation step: 0 = Do not wait, Higher = Wait longer.
              "simulator_step_sleep_factor": 0, # 0, 1, 0.5, 5...
              # A filename str or True to save to working directory.
              "save_invoked_reactions_to_file": True, # True or False
              # A filename str or True to save to working directory.
              "save_hybdehyb_to_file": True, # True or False
              # Enable or disable bending of single helices:
              "enable_helix_bending": False,
              # Merge complexes upon stacking and break if unstacked (if no hybridization interactions):
              "stacking_joins_complexes": True,
              # Allow ends in one complex to stack against ends in another complex:
              "enable_intercomplex_stacking": True,
              # Re-calculate changed domain reactions against all other (complementary) domains,
              # or only re-calculate for pairs against other changed domains?
              "reaction_update_pair_against_all": True,
              # Nric = normalized_reaction_invocation_count = sysmgr.reaction_invocation_count[reaction_spec]/len(sysmgr.domains_by_name[d.name])
              "reaction_throttle": False, # True: default to c_j_throttle_factor = exp(-Nric/10)
              "reaction_throttle_offset": 10,
              "reaction_throttle_reset_on_temperature_change": True,
              "dispatcher_enabled": False,  # True to use a dispatcher to save and visualize graph changes.
              "dispatcher_config": dispatcher_config,
              "stats_total_file": outputfn(statstype="time_totstats", ext="txt"),
              "stats_per_domain_file": outputfn(statstype="time_domain_stats", ext="txt"),
              "stats_per_strand_file": outputfn(statstype="time_strand_stats", ext="txt"),
             }
    with open(outputfn(statstype="config", ext="yaml"), 'w') as fp:
        yaml.dump(params, fp)
    simulator = Simulator(volume=volume, strands=input_oligos, params=params)
    # sysmgr = simulator.reactionmgr


    # Perform calculations and start simulation
    try:
        start = 40  # deg C
        #start = 60  # deg C
        #stop = 80   # deg C
        stop = 70   # deg C
        start, stop = stop, start   # invert => start hot
        step = -1 if start > stop else +1   # deg C
        T_start, T_finish = [273+val for val in (start, stop)]

        melt_cum_stats, domain_stats = simulator.thermodynamic_meltingcurve(T_start, T_finish)
        print("melt_cum_stats:", len(melt_cum_stats))
        print("domain_stats:", len(domain_stats))
        meltingfile = outputfn(statstype="thermo_melting", ext="yaml")
        with open(meltingfile, 'w') as fp:
            yaml.dump(melt_cum_stats, fp)
            print("Thermodynamic melting curve will be saved to file:", meltingfile)
        meltingfile = outputfn(statstype="thermo_domainstats", ext="yaml")
        with open(meltingfile, 'w') as fp:
            yaml.dump(domain_stats, fp)
            print("Thermodynamic melting curve domain stats save to file:", meltingfile)

        #raise KeyboardInterrupt    # Can be used to by-pass annealing
        #touch(outputstatsfile)
        infofn = outputfn(statstype="params", ext="yaml")
        with open(infofn, 'w') as fp:
            yaml.dump({'T_start': T_start, 'T_finish': T_finish, 'delta_T': step,
                       'volume': volume,
                       'params': params,
                       'strand_defs_file': strand_defs_file,
                       'n_strand_copies_default': n_strand_copies_default,
                       'n_steps_per_T': n_steps_per_T,
                       'datetime': datetime.now().isoformat()
                      }, fp)
        print("Basic simulation parameters saved to file:", infofn)
        ## INITIATE ANNEALING

        if 0:
            simulator.anneal(T_start=T_start, T_finish=T_finish, delta_T=step)
        else:
            ## Just simulate at a sigle temperature:
            do_profile = 0 # True, False
            if do_profile:
                import cProfile
                # run in context:
                profilefn = os.path.join(statsfolder, 'simulator.simulate.profile')
                print("Saving profiling stats to file:", profilefn)
                cProfile.runctx('simulator.simulate(T=330, n_steps_max=400, systime_max=1000)',
                                globals=None, locals=locals(), filename=profilefn)
                import pstats
                s = pstats.Stats(profilefn)
                print("Internal Time Top 10:")
                s.sort_stats('cumulative').print_stats(20)
                print("\nTotal Time Top 10:")
                s.sort_stats('time').print_stats(20)
            else:
                simulator.simulate(T=330, n_steps_max=10000, systime_max=100)

    except KeyboardInterrupt:
        print("\n\nABORT: KeyboardInterrupt.\n\n")
        print(traceback.format_exc()) # or print(sys.exc_info()[0])
        # answer = input("Do you want to enter debug mode? [yes]/no  ")
        # answer = input("Do you want to save the simulation data for this aborted run? [yes]/no  ")
        # if answer and answer[0].lower() == 'n':
        #     return
    # except Exception as e:
    #     print("\n\nException during siumlation:", e)
    #     print(traceback.format_exc()) # or print(sys.exc_info()[0])
    #     answer = input("\nDo you want to enter debug mode?  ")
    #     if answer and answer[0].lower() == 'y':
    #         #interact(local=locals())
    #         pdb.set_trace()
    #     answer = input("Do you want to save the simulation data for this aborted run? [yes]/no  ")
    #     # if answer and answer[0].lower() == 'n':
    #     #     cleanup(outputstatsfile)
    #     #     return

    # simulator.save_stats_cache()

    Tm_filename = outputfn(statstype="energies", ext="yaml")
    with open(Tm_filename, 'w') as fp:
        yaml.dump(simulator.reactionmgr.domain_dHdS, fp)

    reportfilename = outputfn(statstype="domainstats", ext="txt")
    simulator.save_domain_stats(reportfilename)

    reportfilename = outputfn(statstype="report", ext="txt")
    simulator.save_report(reportfilename)

    print("\n\nSimulation complete; stats are available in directory:", statsfolder)


if __name__ == '__main__':
    main()
