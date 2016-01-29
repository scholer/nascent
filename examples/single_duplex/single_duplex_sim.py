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
import logging
logger = logging.getLogger(__name__)
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
from nascent.utils.logging import init_logging

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
    structure = "duplex_16bp-d2" # "duplex2"

    #n_strand_copies_default = 400
    #n_strand_copies_default = 100
    # n_strand_copies_default = 40
    #n_strand_copies_default = 20
    # n_strand_copies_default = 10
    # n_strand_copies_default = 2
    n_strand_copies_default = 1

    # Load strand def and check the strands:
    #strand_defs_file = os.path.join(strand_defs_folder, "strand_defs_{}.txt".format(structure))
    strand_defs_file = os.path.join(scriptdir, "{}.txt".format(structure))

    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file, n_clones_default=n_strand_copies_default)
    check_strands(input_oligos)
    # Color the strands: (http://www.graphviz.org/doc/info/attrs.html#k:color)
    # Color names: https://en.wikipedia.org/wiki/X11_color_names
    colors = {"sienna": "#a0522d",
              "blue": "#0000FF",
              "dark magenta": "#8B008B",
              "crimson": "#DC143C",
              "turquoise": "#40e0d0",
             }
    # Or use colormgr?
    for strand, color in zip(input_oligos, colors.values()):
        for _, _, eattr in strand.edges(data=True):
            eattr['color'] = color
        for _, _, eattr in strand.ends5p3p_graph.edges(data=True):
            eattr['color'] = color

    #n_domains = sum(len(s.domains) for s in input_oligos)

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
    log_args = {'loglevel': logging.DEBUG,
                #'basic_logging': True,  # False = custom logging setup
                #'rotating': False,
               }
    init_logging(args=log_args, logdir=run_directory)

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
              "save_invoked_reactions_to_file": False, # True or False
              # A filename str or True to save to working directory.
              "save_hybdehyb_to_file": True, # True or False
              # Enable or disable bending of single helices:
              "enable_helix_bending": False,
              # Merge complexes upon stacking and break if unstacked (if no hybridization interactions):
              "stacking_joins_complexes": True,
              # Allow ends in one complex to stack against ends in another complex:
              "enable_intercomplex_stacking": False,
              # Re-calculate changed domain reactions against all other (complementary) domains,
              # or only re-calculate for pairs against other changed domains?
              "reaction_update_pair_against_all": True,
              # Nric = normalized_reaction_invocation_count = sysmgr.reaction_invocation_count[reaction_spec]/len(sysmgr.domains_by_name[d.name])
              "reaction_throttle": True, # True: default to c_j_throttle_factor = exp(-Nric/10)
              # Use a reaction throttle cache, decrementing the throttle when reaction is triggered, rather than calculated in calculate_c_j from Nric
              "reaction_throttle_use_cache": True,
              "reaction_throttle_per_complex": True,
              "reaction_throttle_offset": 0,
              "reaction_throttle_reset_on_temperature_change": True,
              "dispatcher_enabled": False,  # True to use a dispatcher to save and visualize graph changes.
              "dispatcher_config": dispatcher_config,
              "stats_total_file": outputfn(statstype="time_totstats", ext="txt"),
              "stats_per_domain_file": outputfn(statstype="time_domain_stats", ext="txt"),
              "stats_per_strand_file": outputfn(statstype="time_strand_stats", ext="txt"), #
              "stats_post_simulation_file": outputfn(statstype="post_simulation_stats", ext="yaml"), #
              # complexes - these are saved by reactionmgr as new complex assemblies (states) are encountered.
              "reaction_graph_complexes_directory": os.path.join(statsfolder, "complexes"), #
              # reaction_graph - are saved by statsmgr after every simulation.
              "reaction_graph_output_directory": os.path.join(statsfolder, "reaction_graph"), #
              "reaction_graph_output_fnfmt": "reaction_graph_{systime}.{ext}", #
              "reaction_graph_output_formats": "png",
              # How far back to look when looking for reaction microcycles:
              "reaction_microcycles_slice_size": 5,
             }
    with open(outputfn(statstype="config", ext="yaml"), 'w') as fp:
        yaml.dump(params, fp)
    simulator = Simulator(volume=volume, strands=input_oligos, params=params)
    # sysmgr = simulator.reactionmgr

    simulator.stats_writer.write_post_simulation_stats()

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
            from nascent.energymodels.biopython import hybridization_dH_dS, R, binary_state_probability
            from nascent.graph_sim_nx.constants import N_AVOGADRO, AVOGADRO_VOLUME_NM3
            from math import exp, log as ln
            T = 330
            # Either add volume energy here or pass a Q=<conc> parameter to binary_state_probability
            # conc = 1/(volume*N_AVOGADRO)
            volume_energy = -R*T*ln(volume*N_AVOGADRO/n_strand_copies_default)
            full_duplex_seq = "".join(domain.sequence for domain in input_oligos[0].domains)
            duplex_dH, duplex_dS = hybridization_dH_dS(full_duplex_seq) # does not account for concentration
            # kcal/mol for dH and cal/mol/K for dS.
            duplex_dG = duplex_dH*1000 - T * duplex_dS  # Units of cal/mol
            # Volume energy is for the un-bound state; subtract for duplex if unbound is defined as dG = 0.
            dG = duplex_dG - volume_energy
            print("Duplex sequence:", full_duplex_seq)
            print("Duplex energy:      %0.01f cal/mol" % duplex_dG)
            print(" - duplex_dH    : %0.01f" % (duplex_dH*1000))
            print(" - duplex_dS*T  : %0.01f" % (T * duplex_dS))
            print(" - volume_energy:   %0.01f" % volume_energy)
            print("dG incl volume:     %0.01f cal/mol" % dG)
            print("Duplex probability: %0.02f" % binary_state_probability(dG, T))
            # binary state probability just calculates:
            # print("dG/RT            = %0.01f" % (dG/(R*T)))
            # print("exp(-dG/RT)      = %0.04f" % (exp(-dG/(R*T))))
            # print("1/(1+exp(dG/RT)) = %0.04f" % (1/(1+exp(dG/(R*T)))))
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
                # pdb.set_trace()
                # simulator.simulate(T=330, n_steps_max=10000, systime_max=1000)
                #simulator.simulate(T=330, n_steps_max=1000, systime_max=10)
                simulator.simulate(T=330, n_steps_max=400, systime_max=1000)
                # simulator.simulate(T=330, n_steps_max=n_steps_per_T, systime_max=time_per_T)

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

    simulator.stats_writer.write_post_simulation_stats()
    simulator.stats_writer.close_all() # Close all open files...

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
