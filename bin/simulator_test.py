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
from pprint import pprint
from collections import defaultdict, namedtuple
import pdb

# try:
#     import rlcompleter
#     import readline
# except ImportError:
#     pass

LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nascent_sim
from nascent.graph_sim_nx.simulator_dm import DM_Simulator as Simulator
from nascent.graph_sim_nx.fileio import parse_strand_domains_file, check_strands
from nascent.graph_sim_nx.dispatcher import StateChangeDispatcher
from nascent.graph_sim_nx.constants import N_AVOGADRO #, R
from nascent.graph_sim_nx import debug

# Toggle debug printing:
debug.do_print = False

ReactionAttrs = namedtuple('ReactionAttrs', ['reaction_type', 'is_forming', 'is_intra'])

#N_AVOGADRO = 6.022e23   # /mol


def touch(filepath):
    with open(filepath, 'a'):
        os.utime(filepath)


def test(simulator, usepdb=False, testcase=None, execfile=None):
    """
    Perform ad-hoc testing of simulator system.
    """

    #local = {'simulator': simulator}
    #ns = locals()
    n_done = 0
    sysmgr = simulator.reactionmgr
    sysmgr_grouped = simulator.reactionmgr_grouped

    domains_by_duid = {}  # domain unique id => domain
    for domain in sysmgr.domains:
        assert domain.duid not in domains_by_duid
        domains_by_duid[domain.duid] = domain

    domA = next(d for d in sysmgr.domains if d.name == 'A')
    domB = next(d for d in sysmgr.domains if d.name == 'B')
    domC = next(d for d in sysmgr.domains if d.name == 'C')
    domb = next(d for d in sysmgr.domains if d.name == 'b')
    domD = next(d for d in sysmgr.domains if d.name == 'D')
    dome = next(d for d in sysmgr.domains if d.name == 'e')
    domE = next(d for d in sysmgr.domains if d.name == 'E')
    domF = next(d for d in sysmgr.domains if d.name == 'F')
    domc = next(d for d in sysmgr.domains if d.name == 'c')
    domG = next(d for d in sysmgr.domains if d.name == 'G')
    domH = next(d for d in sysmgr.domains if d.name == 'H')
    domg = next(d for d in sysmgr.domains if d.name == 'g')
    domK = next(d for d in sysmgr.domains if d.name == 'K')
    domi = next(d for d in sysmgr.domains if d.name == 'i')
    domI = next(d for d in sysmgr.domains if d.name == 'I')
    domL = next(d for d in sysmgr.domains if d.name == 'L')

    #execfile = r"C:\Users\scholer\Nascent\simdata\DistTest1\2015-11-10 120035\invoked_reactions.py"
    # testcase = None # None, 1, ...

    if execfile:
        # Execute the contents of file <execfile>, e.g. hyb_dehyb.py or invoked_reactions.py
        with open(execfile) as fp:
            code = compile(fp.read(), execfile, "exec") # exec or eval
        exec(code, globals(), {'sysmgr': sysmgr, 'domains_by_duid': domains_by_duid})

    elif testcase == 1:
        # sysmgr.hybridize(domB, domb)
        sysmgr.react_and_process(frozenset((domB, domb)), is_hybridizing=True, is_intra=False)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)

        # sysmgr.hybridize(domc, domC)
        sysmgr.react_and_process(frozenset((domc, domC)), is_hybridizing=True, is_intra=None)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)

        # sysmgr.hybridize(domi, domI)
        sysmgr.react_and_process(frozenset((domi, domI)), is_hybridizing=True, is_intra=None)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)

        # sysmgr.dehybridize(domi, domI)
        sysmgr.react_and_process(frozenset((domi, domI)), is_hybridizing=False, is_intra=None)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)
        # sysmgr.hybridize(domE, dome)

        sysmgr.react_and_process(frozenset((domE, dome)), is_hybridizing=True, is_intra=None)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)

        # sysmgr.dehybridize(domB, domb)
        sysmgr.react_and_process(frozenset((domB, domb)), is_hybridizing=False, is_intra=None)
        print("\n  - sysmgr.unhybridized_domains_by_name -")
        pprint(sysmgr.unhybridized_domains_by_name)
        sysmgr.hybridize(domG, domg)
        # --> KeyError in File "C:\Users\scholer\Dev\src-repos\na_strand_model\nascent\graph_sim_nx\systemmgr.py", line 380, in update_possible_reactions:
        #    self.unhybridized_domains_by_name[d1.name].remove(d1)

    else:
    # elif usepdb:
        # Use the python debugger (includes readline since 3.3):
        # import pdb
        #pdb.run('simulator.simulate(T=350, n_steps_max=1, systime_max=1)', locals=locals())
        # Or start main script from console using:
        # python -m pdb simulator_test.py
        # If you want to debug at a certain place in the code, use pdb.set_trace() to break out to the pdb.
        # To enter debug after experiencing an exception, use pdb.pm(). Debug starts from sys.last_traceback
        # To enter debug while you have a current exception or traceback, use pdb.post_morten
        # pdb.set_trace()
        # simulator.systemmgr.draw_and_save_graphs(n=1)
        # Or just run with python -m pdb <script>  (and press 'c' to continue until next error or breakpoint).
        # (Does not work if the exception is catched, of course...!)
        n_done = simulator.simulate(T=330, n_steps_max=100, systime_max=1000) #
        # try:
        #     # n_done = simulator.simulate(T=330, n_steps_max=10, systime_max=200)
        #     n_done = simulator.simulate(T=328, n_steps_max=20, systime_max=200) #
        #     # n_done = simulator.simulate(T=325, n_steps_max=20, systime_max=200) # Stable complex
        #     # n_done = simulator.simulate(T=320, n_steps_max=20, systime_max=200) # Stable complex
        #     # simulator.simulate(T=310, n_steps_max=10, systime_max=100) # Staple, irreversible.
        # except Exception as e:
        #     n_done = "failed"
        #     print("Exception:", type(e), e)
        #     #traceback.print_last()
        #     traceback.print_exc()
        #     enter_pdb = prompt_yes_no("Enter debugger?")
        #     if enter_pdb in ('y', True):
        #         # pdb.set_trace()
        #         exc_type, exc_value, tb = sys.exc_info() # get exception info
        #         traceback.print_exc()
        #         pdb.post_mortem(tb)
        # try:
        #     simulator.simulate(T=350, n_steps_max=10, systime_max=100)
        # except TypeError:
        #     pdb.post_mortem()  # use while handling an exception.
        #pdb.pm()  # only use after execution, when sys.last_traceback has been set.
    # else:
    #     readline.set_completer(rlcompleter.Completer(locals()).complete)
    #     readline.parse_and_bind("tab: complete")
    #     from importlib import reload
    #     #newc, oldc = moves.hybridize(domb, domB)
    #     # Simulate a single step:
    #     n_done = simulator.simulate(T=350, n_steps_max=10, systime_max=100)
    #     interact(local=locals())

    print("\n\nTEST RUN COMPLETE. n_done =", n_done)
    answer = input(("Type 'd' to enter debugger; "
                    "'g' to save plot to file; any other key to continue..."))
    if 'g' in answer:
        sysmgr.draw_and_save_graphs(n="end")
    if 'd' in answer:
        pdb.set_trace()


def run_repeatedly(simulator):
    sysmgr = simulator.reactionmgr
    answer = 'r'
    systime_max = 0
    T = 330
    n_steps_max = 30
    print("Starting simulation...")
    while 'q' not in answer:
        systime_max += 1000
        answer = input(("Type 'd' to enter debugger; " +
                        "'g' to save plot to file; " +
                        "'q' to quit; " +
                        "Type <temperature> [<max steps>] to change simulation parameters " +
                        ("(currently %s %s). " % (T, n_steps_max)) +
                        "Press any other key to run simulation again..."))
        if answer:
            if 'q' in answer:
                break
            if 'd' in answer:
                pdb.set_trace()
            elif 'g' in answer:
                sysmgr.draw_and_save_graphs(n="end")
            else:
                old_T = T
                try:
                    vals = answer.strip().split()
                    T = int(vals[0])
                    n_steps_max = int(vals[1])
                except TypeError:
                    pass
                except IndexError:
                    pass
                if T != old_T:
                    simulator.reactionmgr.reset_temperature(T)
        n_done = simulator.simulate(T=T, n_steps_max=n_steps_max, systime_max=systime_max)
        print("\nCompleted %s steps..." % n_done)




def cleanup(outputstatsfile):
    """ Remove all outputfiles related to outputstatsfile. """
    if os.path.exists(outputstatsfile):
        try:
            os.remove(outputstatsfile)  # clean up, if we are not saving.
        except IOError as e:
            print("Could not remove %s: %s" % (outputstatsfile, e))
    # TODO: You should remove all files explicitly. Some datafiles are appended, not overwritten.
    for fn in glob.glob(os.path.splitext(outputstatsfile)[0]+"*"):
        print("- Removing", fn)
        try:
            os.remove(fn)
            print("- Removed", fn)
        except IOError as e:
            print("- Unable to remove", fn, ":", e)


def prompt_yes_no(prompt, default=None):
    """ Prompt user for a yes/no answer, defaulting to :default:.
    If default is None, this function will continue to ask until a clear yes or no has been given."""
    answer = ""
    while answer not in ('y', 'n'):
        answer = input(prompt).strip().lower()
        if not answer:
            if default:
                answer = default
        else:
            answer = answer[0].lower()
        if answer not in ('y', 'n'):
            print("Answer must start with Y/y or N/n. Please try again.")
    return answer


def main():
    """ Main test function. """
    #strand_defs_file = os.path.join(os.path.dirname(__file__), "testfiles", "strand_defs01.txt")
    # LIBPATH,
    # adhoc_testing = False
    adhoc_testing = True
    usepdb = True
    do_run_repeatedly = True # False, True
    enable_gephi_stream = False # False, True
    strand_defs_folder = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(nascent.__file__))), "testfiles")
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
    #n_strand_copies_default = 10
    n_strand_copies_default = 2

    if adhoc_testing:
        structure = "DistTest1"
        n_strand_copies_default = 4

    # Load strand def and check the strands:
    strand_defs_file = os.path.join(strand_defs_folder, "strand_defs_{}.txt".format(structure))
    print("Loading strand/domain defs from file:", strand_defs_file)
    input_oligos = parse_strand_domains_file(strand_defs_file, n_clones_default=n_strand_copies_default)
    check_strands(input_oligos)
    n_domains = sum(len(s.domains) for s in input_oligos)

    run_directory = os.path.join(os.path.expanduser("~"),
                                 "Nascent",
                                 "simdata",
                                 structure,
                                 datetime.now().strftime("%Y-%m-%d %H%M%S"))
    if os.path.isdir(run_directory):
        clean = prompt_yes_no("Directory '" + run_directory + "' already exists. Remove this before proceeding? ")
        if clean == 'y':
            os.remove(run_directory)
    elif os.path.exists(run_directory):
        print("Directory '" + run_directory + "' already exists BUT IS NOT A DIRECTORY. ABORTING!")
        #raise OSError("Statsfolder/run_directory already exists and is not a folder:", run_directory)
        return
    else:
        os.makedirs(run_directory)

    statsfolder = run_directory # os.path.join(strand_defs_folder, "simdata", structure)

    # if not os.path.exists(statsfolder):
    #     os.mkdir(statsfolder)
    # elif not os.path.isdir(statsfolder):
    #     raise OSError("Statsfolder already exists and is not a folder:", statsfolder)

    # Ensure that we have a unique filename:
    outputstatsfnfmt = os.path.join(statsfolder, "graph_sim_nx.{uid}.{statstype}.{ext}")
    def outputfn(**kwargs):
        return os.path.expanduser(outputstatsfnfmt.format(**kwargs))
    #outputstatsfile = os.path.expanduser("~/nascent_dom_anneal_test.txt")
    uid = 0
    outputstatsfile = os.path.expanduser(outputstatsfnfmt.format(uid=uid, statstype="stats", ext="txt"))
    assert outputfn(uid=uid, statstype="stats", ext="txt") == outputstatsfile
    while os.path.exists(outputstatsfile):
        uid += 1
        outputstatsfile = outputfn(uid=uid, statstype="stats", ext="txt")
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

    n_steps_per_T = 100000
    time_per_T = 10         # seconds

    print("n_steps_per_T:", n_steps_per_T)
    print("time_per_T:", time_per_T, "s")


    ## SET UP SIMULATOR
    adhoc_params = {'time_per_T': time_per_T,
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
                    "reaction_throttle": True, # True: default to c_j_throttle_factor = exp(-Nric/10)
                    "reaction_throttle_offset": 10,
                    "reaction_throttle_reset_on_temperature_change": True,
                    "dispatcher_enabled": False, #enable_gephi_stream,
                    "stats_writer_endabled": True,
                   }
    simulator = Simulator(volume=volume, strands=input_oligos, params=adhoc_params,
                          outputstatsfiles=outputstatsfile)
    sysmgr = simulator.reactionmgr

    ## Hook up dispatcher
    dispatcher_config = {'dispatcher_state_changes_fn': outputfn(uid=uid, statstype="changes", ext="txt"),
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
    if enable_gephi_stream:
        dispatcher = StateChangeDispatcher(dispatcher_config)
        simulator.dispatcher = dispatcher
        dispatcher.init_graph(sysmgr.domain_graph)
        answer = input("Gephi initialized, press any key to continue (or 'q' to quit)...")
        if answer and answer[0] == 'q':
            return

    # Perform calculations and start simulation
    if do_run_repeatedly:
        run_repeatedly(simulator)
        return
    if adhoc_testing:
        test(simulator, usepdb=usepdb)
        return
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
        print("Basic simulation parameters saved to file:", outputstatsfile)
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
        yaml.dump(simulator.reactionmgr.domain_dHdS, fp)

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
