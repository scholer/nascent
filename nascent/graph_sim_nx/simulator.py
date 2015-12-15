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

# pylint: disable=C0103,W0142

"""

Module for

"""

import os
#import random
#from collections import defaultdict
import math
#from math import exp, log as ln
from datetime import datetime

#import networkx as nx
#import numpy as np


from nascent.energymodels.biopython import Tm_NN
# Consider using NN table that has units of R, R/K:
#from nascent.energymodels.biopython import energy_tables_in_units_of_R
#DNA_NN4 = energy_tables_in_units_of_R["DNA_NN4"]
from .thermodynamic_utils import thermodynamic_meltingcurve
from .reactionmgr import ReactionMgr
from .debug import printd, pprintd

# Module-level constants and variables
N_AVOGADRO = 6.022e23   # /mol
R = 1.987  # universal gas constant in cal/mol/K
VERBOSE = 0



def complex_sizes_hist(complexes):
    """
    Calculate complex size histogram:
        hist[N] => Count
    Where N is a complex size (int) and count is the number of complexes of this size.
    """
    hist = {}
    for c in complexes:
        N = len(c.strands)
        if N not in hist:
            hist[N] = 0
        hist[N] += 1
    return hist




class Simulator():
    """
    Simulator class to hold everything required for a single simulation.
    """

    def __init__(self, volume, strands, params, domain_pairs=None,
                 outputstatsfiles=None):
        """
        outputstatsfiles : A dict of <stat type>: <outputfilename>
        """
        self.params = params
        self.reactionmgr = ReactionMgr(volume=volume, strands=strands, params=params, domain_pairs=domain_pairs)
        self.sim_system_time = 0  # System time within the simulator, typically measured in micro-seconds.

        ## Simulation stats, counts:
        self.N_steps = 0    # Total number of steps
        self.N_changes = 0  # Number of state changes (hybridizations or de-hybridizations)
        self.N_selections = 0 # Total number of succeessfull selection of domain1 and domain2.
        self.N_steps_per_T = params.get('n_steps_per_T', 100000)
        self.system_stats = {}


        ## Data recording ##
        self.Default_statstypes = ("changesampling", )
        if isinstance(outputstatsfiles, str):
            base, _ = os.path.splitext(outputstatsfiles)
            outputstatsfiles = {k: base+"_"+k+'.csv' for k in self.Default_statstypes}
        self.Outputstatsfiles = outputstatsfiles    # dict with <stats type>: <outputfilepath> entries
        self.Record_stats = params.get("record_stats", True)    # Enable or disable stats recording
        self.stepsampling_frequency = self.params.get('stepsampling_frequency', 1)

        # Save stats to a cache and only occationally append them to file on disk.
        # Stats_cache: dict <stats type>: <stats>, where stats is a list of tuples:
        # [(Temperature, N_dom_hybridized, %_dom_hybridized, N_oligos_hybridized, %_oligos_hybridized), ...]
        self.Stats_cache = {k: [] for k in self.Default_statstypes}
        self.complex_size_stats = {k: [] for k in self.Default_statstypes}
        self.print_setup()
        self.print_statsline_when_saving = params.get("print_statsline_when_saving", True)

        ## Hooks ##
        self.Visualization_hook = None
        self.print_post_step_fmt = "{self.N_steps: 5} {self.N_selections: 5}"
        self.dispatcher = None


    def print_setup(self, fp=None):
        """ Print the simulator setup. """
        sysmgr = self.reactionmgr
        # concentration:
        c = 1/N_AVOGADRO/sysmgr.volume
        print("Simulator setup / parameters:", file=fp)
        print("Total concentration of each strand:", file=fp)
        print("\n".join("- {:10}: {:0.3g} uM (N={})".format(name, len(entries)*c*1e6, len(entries))
                        for name, entries in sysmgr.strands_by_name.items()), file=fp)
        print("Concentration of each domain is: %0.3g uM" % (c*1e6), file=fp)
        print("Domain pairing map:",
              #", ".join("->".join(kv) for kv in sysmgr.domain_pairs.items()), file=fp)
              sysmgr.domain_pairs, file=fp)
        print("Total number of strands:", sysmgr.N_strands, file=fp)
        print(" - number of different strands:", len(sysmgr.strands_by_name), file=fp)
        print("Total number of domains:", sysmgr.N_domains, file=fp)
        print(" - number of different domains:", len(sysmgr.domains_by_name), file=fp)
        print("Steps per T:", self.N_steps_per_T, file=fp)

    def save_report(self, reportfile):
        """ Save the output of print_setup to file. """
        with (open(reportfile, 'w') if isinstance(reportfile, str) else reportfile) as fp:
            self.print_setup(fp=fp)

    def print_domain_stats(self, fp=None, Tm_params=None):
        """ Print domain statistics (for each domain in each strand). """
        sysmgr = self.reactionmgr
        if Tm_params is None:
            Tm_params = {}
        #statdict = {}
        c_each = 1/(N_AVOGADRO*sysmgr.volume)
        for sname, strands in sysmgr.strands_by_name.items():
            c_strand_total = len(strands)*c_each
            # Perhaps also calculate the concentration of dnac2 explicitly
            # (insted of assuming it to be the same as c_strand_total)
            print("\n[ %s strand ]" % sname, file=fp)
            print("c_each: {:0.04g} uM".format(c_each*1e6), file=fp)
            print("n_species (copies):", len(strands), file=fp)
            print("c_strand_total: {:0.04g} uM".format(c_strand_total*1e6), file=fp)
            for domain in strands[0].domains:
                print("Domain:", domain.name, file=fp)
                print(" - total copy count of this domain:", len(sysmgr.domains_by_name[domain.name]), file=fp)
                print("\n".join(" - %s: %s" % (att, getattr(domain, att))
                                for att in ('sequence', )), file=fp)
                try:
                    deltaH, deltaS = sysmgr.domain_dHdS[domain.name]
                    print(" - deltaH, deltaS: {:.04g} kcal/mol, {:.04g} cal/mol/K".format(
                        *sysmgr.domain_dHdS[domain.name]), file=fp)
                except KeyError:
                    pass
                else:
                    # Tm_NN concentration unit is nM:
                    Q_melt = (c_strand_total - c_strand_total/2)
                    Tm = (1000 * deltaH) / (deltaS + (R * (math.log(Q_melt)))) - 273.15
                    print(" - Tm: {:.04g} (from deltaH, deltaS)".format(Tm), file=fp)
                Tm = Tm_NN(domain.sequence, dnac1=c_strand_total*1e9, dnac2=c_strand_total*1e9, **Tm_params)
                print(" - Tm: {:.04g} (from Tm_NN)".format(Tm), file=fp)
                Tm = Tm_NN(domain.sequence, dnac1=c_strand_total*1e9, dnac2=c_strand_total*1e9,
                           nn_table="DNA_NN4", **Tm_params)
                print(" - Tm: {:.04g} (from Tm_NN using DNA_NN4)".format(Tm), file=fp)


    def save_domain_stats(self, reportfile, Tm_params=None):
        """ Save domain stats to file. """
        with (open(reportfile, 'w') if isinstance(reportfile, str) else reportfile) as fp:
            self.print_domain_stats(fp=fp, Tm_params=Tm_params)


    def step(self, T):
        """ Overwrite in subclass (or overwrite simulate method entirely). """
        pass


    def simulate(self, T, n_steps_max=100000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        sysmgr = self.reactionmgr
        assert sysmgr.n_hybridized_domains() == sysmgr.N_domains_hybridized
        n_done = 0
        sysmgr.temperature = T
        while n_done < n_steps_max:
            try:
                self.step(T)
            except AssertionError as e:
                print("AssertionError:", e)
                print("self.n_hybridized_domains(), self.N_domains_hybridized = %s, %s" % \
                      (sysmgr.n_hybridized_domains(), sysmgr.N_domains_hybridized))
                raise e
            n_done += 1
            self.N_steps += 1
            if n_done % 10000 == 0:
                print(("Simulated %s of %s steps at T=%s K (%0.0f C). " +
                       "%s state changes with %s selections in %s total steps.") %
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))
            if self.Record_stats and (self.N_steps % self.stepsampling_frequency == 0):
                self.record_stats_snapshot(T, statstype="timesampling")
        assert sysmgr.n_hybridized_domains() == sysmgr.N_domains_hybridized


    def anneal(self, T_start, T_finish, delta_T=-1, n_steps_per_T=None):
        """
        Simulate annealing repeatedly from T_start to T_finish,
        decreasing temperature by delta_T for every round,
        doing at most n_steps number of steps at each temperature.
        # TODO: I am currently only capturing stats whenever the state changs (hybridization/melting).
        # This means that if we only have two states, those two states will be equally represented in the data,
        # even though one state might be favored and system spends the majority of the total time in this state.
        # It would probably be nice to also capture "time-interspersed" snapshots, so that I can distinguish this.
        """

        # Range only useful for integers...
        n_steps_per_T = self.N_steps_per_T or 100000
        T = T_start
        assert delta_T != 0
        assert T_start > T_finish if delta_T < 0 else T_finish > T_start
        while T >= T_finish if delta_T < 0 else T <= T_finish:
            print("\nSimulating at %s K for %s steps (ramp is %s K to %s K in %s K increments)" % \
                  (T, n_steps_per_T, T_start, T_finish, delta_T))
            self.simulate(T, n_steps_per_T)
            T += delta_T
            self.save_stats_cache() # Save cache once per temperature
        print("Annealing complete! (%s)" % datetime.now().strftime("%Y-%m-%d %H:%M"))
        self.print_setup()


    def save_stats_cache(self, outputfilenames=None):
        """ Save stats cache to outputfn. """
        sysmgr = self.reactionmgr
        if self.print_statsline_when_saving:
            print("| Total domain hybridization percentage: {:.0%} ({} of {})".format(
                sysmgr.N_domains_hybridized/sysmgr.N_domains,
                sysmgr.N_domains_hybridized, sysmgr.N_domains))
            #print(", ".join(str(i) for i in self.Stats_cache[-1]))
        if outputfilenames is None:
            outputfilenames = self.Outputstatsfiles
        if not outputfilenames:
            print("Unable to save stats cache: outputstatsfile is:", outputfilenames)
            return
        for statstype, outputfn in outputfilenames.items():
            print(" - Saving stats cache to file:", outputfn)
            with open(outputfn, 'a') as fp:
                # self.Stats_cache = [(Temperature, N_doms_hybridized, %_doms_hybridized), ...]
                fp.write("\n".join(", ".join(str(i) for i in line) for line in self.Stats_cache[statstype])+"\n")
            self.Stats_cache[statstype] = []    # Reset the cache
            c_size_fn = os.path.splitext(outputfn)[0] + '.complex_sizes.txt'
            if self.complex_size_stats[statstype]:
                with open(c_size_fn, 'a') as fp:
                    # dumping a list of dicts should be append'able.
                    # yaml ends dump with \n
                    #yaml.dump(self.complex_size_stats[statstype], fp, default_flow_style=False)
                    # yaml is not reliable for large data files, plus it is SUPER SLOW.
                    # list of dicts:
                    # [{T: {N: count of complexes with size N}}]
                    # T: size,count  size,count  ...
                    fp.write("\n".join("%s: " % T + "\t".join("%s,%s" % tup for tup in hist.items())
                                       for entry in self.complex_size_stats[statstype]
                                       for T, hist in entry.items()) + "\n")
                self.complex_size_stats[statstype] = [] # reset cache


    def record_stats_snapshot(self, T, statstype="changesampling"):
        """ Save stats snapshot to stats cache. """
        sysmgr = self.reactionmgr
        self.Stats_cache[statstype].append((T,
                                            sysmgr.N_domains_hybridized,
                                            sysmgr.N_domains_hybridized/sysmgr.N_domains,
                                            sysmgr.N_strands_hybridized,
                                            sysmgr.N_strands_hybridized/sysmgr.N_strands
                                           ))
        # TODO: Consider option to also count all "non-complexed strands (as complexes of size 1)."
        self.complex_size_stats[statstype].append({T: complex_sizes_hist(sysmgr.complexes)})


    def thermodynamic_meltingcurve(self, T_start, T_finish, delta_T=None, volume=None):
        """
        Calculate thermodynamic melting curve.
        This uses a lot of parameters stored in self:
        - self.volume
        - self.domains_by_name  # including the number of domains for concentration
        - self.domain_pairs
        - self.domain_dHdS
        Returns a list/ordered dict of
            temperature : [c_hybridized, c_non_hybridized, total_conc, domains_total]
            where c_hybridized, c_non_hybridized, and total_conc are in units of M = mol/L.
        """
        sysmgr = self.reactionmgr
        if delta_T is None:
            delta_T = -0.5 if T_start > T_finish else +0.5
        if volume is None:
            volume = sysmgr.volume
        return thermodynamic_meltingcurve(T_start, T_finish, delta_T, volume,
                                          sysmgr.domains_by_name, sysmgr.domain_pairs,
                                          sysmgr.domain_dHdS)
