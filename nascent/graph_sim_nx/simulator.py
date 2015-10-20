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

# pylint: disable=C0103

"""

Module for

"""

import os
import random
from collections import defaultdict
import math
from math import exp, log as ln
from datetime import datetime

import networkx as nx
import numpy as np


from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
from .thermodynamic_utils import thermodynamic_meltingcurve



# Module-level constants and variables
N_AVOGADRO = 6.022e23   # /mol
R = 1.987  # universal gas constant in cal/mol/K
VERBOSE = 0



def complex_sizes_hist(complexes):
    hist = {}
    for c in complexes:
        N = len(c.Strands)
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
        self.simulation_time = 0
        self.temperature = params.get('temperature', 300)
        self.params = params
        self.volume = volume
        self.complexes = []
        self.removed_complexes = []
        self.strands = strands
        self.strands_by_name = defaultdict(list)
        for strand in strands:
            self.strands_by_name[strand.Name].append(strand)
        print("Strands in self.strands_by_name:")
        print("\n".join("- %10s: %s species" % (sname, len(strands))
                        for sname, strands in self.strands_by_name.items()))
        self.domains_list = [domain for strand in strands for domain in strand.Domains]
        self.domains = set(self.domains_list)

        # Stats - counts
        self.N_domains = len(self.domains)
        self.N_strands = len(self.strands)
        self.N_domains_hybridized = sum(1 for domain in self.domains_list if domain.Partner)
        self.N_strands_hybridized = sum(1 for oligo in self.strands if oligo.is_hybridized())
        self.N_steps = 0    # Total number of steps
        self.N_changes = 0  # Number of state changes (hybridizations or de-hybridizations)
        self.N_selections = 0 # Total number of succeessfull selection of domain1 and domain2.
        self.N_steps_per_T = params.get('n_steps_per_T', 100000)

        self.domains_by_name = defaultdict(list)
        for d in self.domains:
            self.domains_by_name[d.Name].append(d)
        print("Domains in self.domains_by_name:")
        print("\n".join("- %10s: %s species" % (dname, len(domains))
                        for dname, domains in self.domains_by_name.items()))
        if domain_pairs is None:
            # mapping: dom_a -> dom_A, dom_A -> dom_a
            # TODO: This could perhaps be a list, if you want to have different types of domains interacting,
            # E.g. dom_a could be perfect match for dom_A, while dom_ax has 1 mismatch:
            # domain_pairs[dom_A.Name] = [dom_a.Name, dom_ax.Name]
            # Or it could be a set of sets: {{da, dA}, {dA, dax}} and then generate partners by:
            # partners_species = set(chain(pair for pair in domain_pairs if dA in pair)) - {dA}
            # However, might as well only do this once and save the list!
            # Also, if you change domain_pairs mapping, remember to adjust Domain_dHdS cache as well.
            domain_pairs = {d.Name: d.name.lower() if d.name == d.name.upper() else d.name.upper()
                            for d in self.domains_list}
            # remove pairs without partner:
            domain_pairs = {d1name: d2name for d1name, d2name in domain_pairs.items()
                            if d2name in self.domains_by_name}
        # allow self-complementarity?
        assert not any(k == v for k, v in domain_pairs.items())
        self.domain_pairs = domain_pairs
        #self.upper_case_domains = [d for d in self.domains if d == d.upper()]
        #self.lower_case_domains = [d for d in self.domains if d == d.lower()]


        ### System graphs ###
        # Not sure whether we should have system-level graphs or a graph for each complex.


        ### Caches: ###
        # Standard enthalpy and entropy of hybridization,
        # indexed as [frosenset({d1.Name, d2.Name})][0 for enthalpy, 1 for entropy]
        # - Note: index with frozenset((a,b)) or just cache[a, b] = cache[b, a] = value? # d[1,2] same as d[(1,2)]
        # --> Creating a set() or frozenset() takes about 10x longer than to make tuple.
        # --> dict assignment with frozenset is 0.4/0.5 us vs 0.17/0.05 for the "double tuple entry" (python/pypy),
        #       so if you have the memory for it, go for the faster tuples which takes 2x memory.
        # Note: Whenever I use "name", I'm refering to the name of the specie - domain or strand specie.
        self.domain_dHdS = {}

        ## Data recording ##
        self.Default_statstypes = ("changesampling", )
        if isinstance(outputstatsfiles, str):
            base, ext = os.path.splitext(outputstatsfiles)
            outputstatsfiles = {k: base+"_"+k+'.csv' for k in self.Default_statstypes}
        self.Outputstatsfiles = outputstatsfiles    # dict with <stats type>: <outputfilepath> entries
        self.Record_stats = params.get("record_stats", True)    # Enable or disable stats recording

        # Save stats to a cache and only occationally append them to file on disk.
        # Stats_cache: dict <stats type>: <stats>, where stats is a list of tuples:
        # [(Temperature, N_dom_hybridized, %_dom_hybridized, N_oligos_hybridized, %_oligos_hybridized), ...]
        self.Stats_cache = {k: [] for k in self.Default_statstypes}
        self.Complex_size_stats = {k: [] for k in self.Default_statstypes}
        print("Simulator initiated at V=%s with %s strands spanning %s domains." \
              % (self.volume, len(self.strands), len(self.domains)))
        self.print_setup()


        ## Hooks ##
        self.Visualization_hook = None


    def print_setup(self, fp=None):
        """ Print the simulator setup. """
        c = concentration = 1/N_AVOGADRO/self.volume
        print("Simulator setup / parameters:", file=fp)
        print("Total concentration of each strand:", file=fp)
        print("\n".join("- {:10}: {:0.3g} uM (N={})".format(name, len(entries)*c*1e6, len(entries))
                        for name, entries in self.strands_by_name.items()), file=fp)
        print("Concentration of each domain is: %0.3g uM" % (c*1e6), file=fp)
        print("Domain pairing map:", ", ".join("->".join(kv) for kv in self.domain_pairs.items()), file=fp)
        print("Total number of strands:", self.N_strands, file=fp)
        print(" - number of different strands:", len(self.strands_by_name), file=fp)
        print("Total number of domains:", self.N_domains, file=fp)
        print(" - number of different domains:", len(self.domains_by_name), file=fp)
        print("Steps per T:", self.N_steps_per_T, file=fp)
        oversampling_factor = self.params['probablity_oversampling_factor']
        print("Probability oversampling factor:", oversampling_factor, file=fp)
        print("Domain concentration * oversampling factor: {:.02g}\n - less than 1e-3: {}".format(
            c*oversampling_factor, c*oversampling_factor < 1e-3), file=fp)
        print("n_steps_per_T * domain concentration * oversampling factor: {:.02g}".format(
            self.N_steps_per_T * oversampling_factor * c), file=fp)
        coverage = self.N_steps_per_T * oversampling_factor * c / self.N_domains
        print("Coverage: {:.02g}\n - larger than 1: {}".format(
            coverage, coverage > 1), file=fp)

    def save_report(self, reportfile):
        with (open(reportfile, 'w') if isinstance(reportfile, str) else reportfile) as fp:
            self.print_setup(fp=fp)

    def print_domain_stats(self, fp=None, Tm_params=None):
        if Tm_params is None:
            Tm_params = {}
        #statdict = {}
        c_each = 1/(N_AVOGADRO*self.volume)
        for sname, strands in self.strands_by_name.items():
            c_strand_total = len(strands)*c_each
            # Perhaps also calculate the concentration of dnac2 explicitly
            # (insted of assuming it to be the same as c_strand_total)
            print("\n[ %s strand ]" % sname, file=fp)
            print("c_each: {:0.04g} uM".format(c_each*1e6), file=fp)
            print("n_species (copies):", len(strands), file=fp)
            print("c_strand_total: {:0.04g} uM".format(c_strand_total*1e6), file=fp)
            for domain in strands[0].Domains:
                print("Domain:", domain.Name, file=fp)
                print(" - total copy count of this domain:", len(self.domains_by_name[domain.Name]), file=fp)
                print("\n".join(" - %s: %s" % (att, getattr(domain, att))
                                for att in ('Sequence', )), file=fp)
                try:
                    deltaH, deltaS = self.domain_dHdS[domain.Name]
                    print(" - deltaH, deltaS: {:.04g} kcal/mol, {:.04g} cal/mol/K".format(
                        *self.domain_dHdS[domain.Name]), file=fp)
                except KeyError:
                    pass
                else:
                    # Tm_NN concentration unit is nM:
                    Q_melt = (c_strand_total - c_strand_total/2)
                    Tm = (1000 * deltaH) / (deltaS + (R * (math.log(Q_melt)))) - 273.15
                    print(" - Tm: {:.04g} (from deltaH, deltaS)".format(Tm), file=fp)
                Tm = Tm_NN(domain.Sequence, dnac1=c_strand_total*1e9, dnac2=c_strand_total*1e9, **Tm_params)
                print(" - Tm: {:.04g} (from Tm_NN)".format(Tm), file=fp)
                Tm = Tm_NN(domain.Sequence, dnac1=c_strand_total*1e9, dnac2=c_strand_total*1e9, nn_table=DNA_NN4, **Tm_params)
                print(" - Tm: {:.04g} (from Tm_NN using DNA_NN4)".format(Tm), file=fp)


    def save_domain_stats(self, reportfile, Tm_params=None):
        with (open(reportfile, 'w') if isinstance(reportfile, str) else reportfile) as fp:
            self.print_domain_stats(fp=fp, Tm_params=Tm_params)



    def n_hybridized_domains(self):
        """ Count the number of hybridized domains. """
        count = sum(1 for domain in self.domains_list if domain.Partner)
        if not count % 2 == 0:
            print("Weird - n_hybridized_domains counts to %s (should be an even number)" % count)
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.domains_list if domain.Partner))
        return count

    def n_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for oligo in self.strands if oligo.is_hybridized())



    def simulate(self, T, n_steps_max=100000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        assert self.n_hybridized_domains() == self.N_domains_hybridized
        n_done = 0
        self.temperature = T
        while n_done < n_steps_max:
            try:
                self.step(T)
            except AssertionError as e:
                print("AssertionError:", e)
                print("self.n_hybridized_domains(), self.N_domains_hybridized = %s, %s" % \
                      (self.n_hybridized_domains(), self.N_domains_hybridized))
                raise(e)
            n_done += 1
            self.N_steps += 1
            if n_done % 10000 == 0:
                print("Simulated %s of %s steps at T=%s K (%0.0f C). %s state changes with %s selections in %s total steps." % \
                      (n_done, n_steps_max, T, T-273.15, self.N_changes, self.N_selections, self.N_steps))
            if self.Record_stats and (self.N_steps % self.Timesampling_frequency == 0):
                self.record_stats_snapshot(T, statstype="timesampling")
        assert self.n_hybridized_domains() == self.N_domains_hybridized


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
        if self.Print_statsline_when_saving:
            print("| Total domain hybridization percentage: {:.0%} ({} of {})".format(
                self.N_domains_hybridized/self.N_domains,
                self.N_domains_hybridized, self.N_domains))
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
            if self.Complex_size_stats[statstype]:
                with open(c_size_fn, 'a') as fp:
                    # dumping a list of dicts should be append'able.
                    # yaml ends dump with \n
                    #yaml.dump(self.Complex_size_stats[statstype], fp, default_flow_style=False)
                    # yaml is not reliable for large data files, plus it is SUPER SLOW.
                    # list of dicts:
                    # [{T: {N: count of complexes with size N}}]
                    # T: size,count  size,count  ...
                    fp.write("\n".join("%s: " % T + "\t".join("%s,%s" % tup for tup in hist.items())
                                       for entry in self.Complex_size_stats[statstype]
                                       for T, hist in entry.items()) + "\n")
                self.Complex_size_stats[statstype] = [] # reset cache


    def record_stats_snapshot(self, T, statstype="changesampling"):
        """ Save stats snapshot to stats cache. """
        self.Stats_cache[statstype].append((T,
                                            self.N_domains_hybridized,
                                            self.N_domains_hybridized/self.N_domains,
                                            self.N_strands_hybridized,
                                            self.N_strands_hybridized/self.N_strands
                                           ))
        # TODO: Consider option to also count all "non-complexed strands (as complexes of size 1)."
        self.Complex_size_stats[statstype].append({T: complex_sizes_hist(self.complexes)})


    def thermodynamic_meltingcurve(self, T_start, T_finish, delta_T=None, volume=None):
        """
        Calculate thermodynamic melting curve.
        This uses a lot of parameters stored in self:
        - self.volume
        - self.domains_by_name  # including the number of domains for concentration
        - self.domain_pairs
        - self.domain_dHdS
        """
        if delta_T is None:
            delta_T = -0.5 if T_start > T_finish else +0.5
        if volume is None:
            volume = self.volume
        return thermodynamic_meltingcurve(T_start, T_finish, delta_T, volume,
                                          self.domains_by_name, self.domain_pairs,
                                          self.domain_dHdS)
