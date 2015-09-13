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

Concept:
* Step-wise statistical thermodynamic annealing simulation.
* Simular to multi-strand, but using *domains* rather than bases.
* The hybridization energy (dG) of a domain is determined by a finite number of state permutations:
** For each of the two strands:
** Whther the domain has a neighboring domain (danling penalty).
** Whether the neighboring domain is hybridized (stacking gain).
    Note: The stacking gain can only be used once per side: Even if both strands has both neighboring domains
    hybridized, you only get one stacking gain for each side.
** The change in system entropy caused by the hybridization:
    Combining two complexes to one will decrease the combined entropy of the two complexes.
    Joining the ends of a complex might also decrease rotational freedom of the complex -> entropy reduction.
* Molar concentrations are not an explicit part of the thermodynamic calculations, but are accounted for
    through the interaction probability: If molecule 1 is present in 5 copies and molecule 2 only present in 1 copy,
    and both has the same interaction energy, then molecule 1 has 5 times higher probability than molecule 2
    of being selected as the next "winner".

Questions:
* How to account for intra-molecular reactions? - Hybridizations within a complex?
** Effect on hybridization energy: Depends on whether hybridization will restrict rotational freedom, etc.
** Effect on interaction probability? Increase "effective" concentration vs increase in -dG?
    I think increase in -dG is the better option as this allows us to treat all interactions "equally".
    (I.e. not having to say: If this domain interacts with that domain, then the domain is in 1000x concentration,
    but otherwise it is in 1x concentration.)


Calculating domain hybridization energy for a specific strand:
* Find the "global" dG lookup table for the domain, indexed by the 8-bit-tuple:
    (strand1: (<5p-domain>, <5p-domain-is-hybridized>, <3p-domain>, <3p-domain-is-hybridized),
     strand1: (<5p-domain>, <5p-domain-is-hybridized>, <3p-domain>, <3p-domain-is-hybridized))
*


How to calculate interaction probability?
* Statistical mechanics/physics
* Partition function: Z = sum(exp(-E(q)/(k*T)) for all q)
* Probability of state q: P(q) = exp(-E(q)/(k*T) / Z

How to select which domain will hybridize next?
a) First select a domain based on interaction probability, then select a strand with that domain?
b) Enter all domains of all strands into the lottery, weighted by -dG.
    Perhaps consider all "events" that could possibly happen.
    Uh, if you have 10 strands with domain A and 5 strands with the complementary domain a, that's 50 possible events..
    It might scale considerably better if you just say "all identical domains are treated as one".

Edit: New random model
1) Select a random domain.
2) If the domain is not hybridized: select a rcomplement domain at random.
    The selection has to be at least partially biased against the number/concentration of rcomplement domains.
    Intra-complex reactions should be biased higher.
    In particular, there should be a chance that no rcomplement domain is selected (and nothing will happen).
        Hybridization rate = k * [A] * [a]
21) An alternative way is to select two domains at random and check whether they are complementary.
    However, that is expensive (lots of failed picks).
    Actually, since selection of the rcomplement domain above must be scaled, that might be the same...
        - Except that intra-complex reactions will have higher chance of success.
            And taking that into account is easier if we have already selected the first one..
    Hybridized domains could still be checked whenever they are picked.
3) Calculate the energy difference and partition function for hybridized vs melted state
    of the domain and it's partner.
        E = (look up energy in domain dG table or calculate anew)
    Q: How to properly account for entropy?
        We have already done the random selection, which affects the k_on rate.
    Q: For intra-complex reactions, how much entropy is accounted for by
        "increased effective concentration" and how much must be adjusted in the energy function?
        - The k_off rate must be correct. k_off does not depend on changes to effective concentration.
4) Roll a dice, and change state of the domain and its partner if the dice says so.




Are you considering domain A and it complement, a, seperately?
* YES: Because every strand with domain A is considered separately, thus A and a must also be considered separately.
* NO : A domain hybridization is one event, A hybridizing to a is the same as a hybridizing to A.
    You could say that you are only considering all capital A and their interactions with all lower-case a,
    and not considering all lowercase a and their interaction with capital A.


How to account for events that would occour on different time scales?
* Multistrand selects an event, then forwards the clock by the characteristic time of that event..
    Not sure this would be realistic when we have multiple species in the reaction tube.
* Maybe just consider "sufficiently large" time steps?
    This might dis-favour formations that rely on a lot of fast steps.
* Perhaps not worry too much about this just yet...


Reversability:
* Should I allow domains to melt and dissociate?
* Perhaps add this through a modification to selection method (b) above: consider domain
    melting as an "event", similar to hybridization.


Some calculations:
* Consider a volume of 1 um^3 = (1e-5 dm)^3 = 1e-15 L = 1 femto-liter.
* For c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
    N_Avogagro = 6.02e23 /mol



TODOS:
* TODO: Implement complex manager to keep track of complexes.
* TODO: Implement complex graph drawing output.
* TODO: Determine the best way to draw the graph network.
        - Strands are nodes
        - Domains are nodes
        - Domain ends (5p, 3p) are nodes.
* TODO: Equilibrium simulation and plotting
        - Goal: To see how long it takes to reach equilibrium depending on:
            - Number of strands/domains
            - Volume
            - Oversampling factor
        - Plot: hyb% vs N_steps
* TODO: Evaluate the effect of the number of strands/clones.
        - Make sure you have the same concentrations when comparing!
* TODO: Compare simulated Tm vs NN calculated ones -- they should match for simple duplexes.
* TODO: Compare simulation against NuPack data.
* TODO: Go back to "finite number" statistical mechanics and check that
        1) the formulas are still valid, (e.g. Stirlin's approximation)
        2) we produce the correct result.
* TODO: (Advanced) Detech DPX/DAX and other motifs that produce rigitidy
        in the structure (preventing ends from joining/circularizing)
* TODO: Make "cadnano json -> strand+domains defs" converter script.
* TODO: Consider making use of SimPy or other python simulation library
        http://simpy.readthedocs.org/, http://www.grotto-networking.com/DiscreteEventPython.html
* TODO: Determine the relationship between
        - n_steps_per_T
        - n_domains_total: total number of domains
        - concentration / volume
        - oversampling_factor
        For equilibrium, my gut feeling is they they are related as:
        - n_steps_per_T -> is the thing that needs to be large enough.
        - concentration -> lower c means that p_on is lower -> more steps are required (for the same oversampling factor)
        - total number of domains -> Even for a constant concentration, more strands means that more steps are required to reach equilibrium.
        - oversampling_factor -> Higher oversampling means that less steps are required.

    n_steps_per_T  * oversampling * domain_concentration / n_domains_total >= coverage
    n_steps_per_T * oversampling / (n_domains_total * volume * N_Avogadro ) >= coverage
    Example: n_steps_per_T = 500000, volume = 1e-18, oversampling=100, n_domains_total=200
        1e5 * 1e3 / (20 * 1e-18 * 6e23) = 1e8 * 2e1 * 1e18 * 1.6e-24
        = 2*1.6*1e(8+1+18-24) = 3.2*1e3

    Note that typically I want to keep the oligo concentration constant, so:
        volume = volume_default * n_strand_copies
    But also:
        n_domains_total = n_domains_per_strand * n_strand_copies
    So at constant oligo concentration, the coverage increases by by the SQUARE of n_strands_copies!

    Also, we want the number of strands to be sufficiently high:
        n_strand_copies > 10 ?
    We also can't bump up the oversampling too much, because that would
    make the random domain selection too high:
        domain_concentration * oversampling_factor < 1e-3 ?
        oversampling_factor_max = 1000 / domain_concentration
                                = 1000 * N_Avogadro * volume

* TODO: Investigate the current issue (Sep 7th) where a simple duplex (duplex1) simulation
        yields an annealing curve with a prominent shoulder between T=80C and T=64 C (Tm=62C).
        - also present for melting curve (duplex1/#15).
        - more n_steps_per_T makes it worse.
        - high oversampling_factor makes it worse.

        It definitely seems to be an issue with oversampling.

* TODO: Review domain selection vs thermodynamic probability.
        Currently, probability selection of a complementary domain is:
            p_hyb = k * [compl_domain]
        Which is similar to
            Q = [duplex]/([domain]*[complement])
        and:
            p_i = 1 / (1 + math.exp(deltaE/(R*T))*Q)
        The idea with selection was that
            p_hyb = p_selection * (p_hyb|selection)
        But that doesn't sound like it would give the same as above:




* TODO: Implement option to do hybridization selection through thermodynamic calculation
        rather than during selection.


* TODO: Go back and review simulation strategy and model with respect
        to v_hyb, v_melt and K = k_hyb/k_melt

* TODO: Implement graph-based visualization

* TODO: Implement complex manager/tracker.

* TODO: Add complex size histogram data to output stats.


* Done: Calculate the thermodynamicly correct melting curve for reference.
        This is just calculating the partition functions for the two states for every T in the range.
        Probably even just use
            K = 1 / ([strandB]_init - [strandA]_init/2)    # Equilibrium constant, because we are at equilibrium.
            ΔG° = ΔH° - T*ΔS°
            binary_state_probability_cal_per_mol(ΔG°, T, Q=K)

Done:
* Done: Plotting of multiple simulation datasets.
* Done: Re-check probability functions
* Done: Re-check selection function.
* Done: Check the quality of random and numpy.random - both seem good.
* Done: thermodynamic_meltingcurve calculation -- based on K from deltaG from deltaH and deltaS.
* Done: plotting: plot_thermodynamic_meltingcurve
* Done: plot_nupack_f_complexed_strands_vs_T - using downloaded job data from NuPack.
* Done: Plotting of f_hyb_vs_N_steps, with temperature range selection.


==== Comparison with NuPack ====

NuPack is generally oriented towards calculating the exact mathematical solution
rather than using simulation.

"Thermodynamic Analysis of Interacting Nucleic Acid Strands"
* Considers a test tube with complexes of strands. (section 3)
* Introduces a "partition function of the box" (section 3.2)
* "equilibrium concentration for each species of complex in the thermodynamic
    limit of large populations" (section 3.3)


==== Comparison with Multistrand ====

Multistrand can to some degree be considered as a "simulation extension" to NuPack.
- It has a strandcomplex list that is capable of making random choices.
- Complexes can be joined (merged), similar to my Complex.merge.
- Choices gets a random choice input value ('rchoice' = "random choice")
- The controlling simulations system (ssystem.cc) makes heavy use of random number generators (RNGs).
- Selection of complexes and moves is based on "rates" or "flux".
    e.g. move->getRate(), move->totalrate,
    loop->returnFlux() - returns move->getRate()
- Loop seems be be the basic structural unit. Sort of similar ot my domains.
- There are several types of moves (join; delete; )
- OMG: There's over 6k LOC in loop.cc !
- Python interface (in the 'interface' directory) -
- Uses graphviz for drawing ()



===== Structures to simulate ======

NuPack-compatible (un-pseudoknotted) structures:
* Simple duplex
** Duplex separated into two consecutive domains
* Duplex with small, symmetric bulge
* Duplex with large, asymetric bulge
* Threeway junction
* Fourway junction
* Fiveway junction
* Un-pseudoknotted tile with single crossovers

2D tile structures to simulate:
* DAE
* DPE

Other non-scaffolded structures:
* Tetrahedrons

Small scaffolded structures:
* WC-long catenane
* Small tetrahedron

Large scaffolded structures (origamis)
* Rothemund's Tall Rectangle
* Rothemund's Smiley-face
* 6HB
* 24HB

"""

#import sys
import os
import random
#from random import choice
from collections import defaultdict
import math
from datetime import datetime
#from math import log, exp

try:
    # Only used for the weighted choice
    import numpy as np
except ImportError:
    print("Could not import numpy. Will use slower native alternative.")
    np = None

# import yaml
import glob

from .dom_anneal_models import Tube #, Strand, Domain, Complex
from .energymodels.biopython import (binary_state_probability_cal_per_mol,
                                     hybridization_dH_dS, Tm_NN, DNA_NN3, DNA_NN4)


N_AVOGADRO = 6.022e23   # /mol
R = 1.987  # universal gas constant in cal/mol/K


def get_NN_stacks(domain1, domain2):
    """
    Returns two NN stack keys for two hybridized domains.
    The first NN stack is at the 5p end of domain1, (equivalent to the 3p end of domain2)
    the second at the 3p end of domain1 (equivalent to the 5p end of domain2).
    """
    ## "GA/CT" NN stack is  5'-G:A-3'    i.e. domain1.domain5p().Sequence[-1] : domain1.Sequence[ 0]
    ##                      3'-C:T-5'    i.e. domain2.domain3p().Sequence[ 0] : domain2.Sequence[-1]
    ## Assuming domain1.domain5p() == domain2.domain3p().Partner
    ## And the same for the other side at domain1.domain3p
    neighbors = [domain1.domain5p(), domain1,
                 domain2.domain3p(), domain2]
    seqidxs = [-1, 0, 0, -1]
    bases5p = "".join([d.Sequence[seqidxs[i]] if d else "." for i, d in enumerate(neighbors)])

    neighbors = [domain1, domain1.domain3p(),
                 domain2, domain2.domain5p()]
    seqidxs = [-1, 0, 0, -1]
    bases3p = "".join([d.Sequence[seqidxs[i]] if d else "." for i, d in enumerate(neighbors)])

    NN_keys = [bases[0:2]+'/'+bases[2:] for bases in (bases5p, bases3p)]

    return NN_keys


# To use the system's random number generator rather than Python's Mersenne Twister algorithm
# Note: The Mersenne Twister algorithm is very widely used (e.g. GLib, Matlab, C++, etc),
# so it is probably "good enough". https://en.wikipedia.org/wiki/Mersenne_Twister
# NumPy also uses the Mersenne Twister, and should be just as random as python's random.
# sys_random = random.SystemRandom()


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

    def __init__(self, volume, strands, params, verbose=0, domain_pairs=None,
                 outputstatsfiles=None, nn_table=None):
        """
        outputstatsfiles : A dict of <stat type>: <outputfilename>
        """
        ## TODO: Implement a ComplexManager to keep track of all the complexes.
        ## Note: At the moment we rely on complexes being garbage collected when it is not bound by any strands
        ## via strand.Complex attribute.
        self.NN_Table = nn_table if nn_table else DNA_NN4
        self.Params = params
        self.Tube = Tube(volume, strands)
        self.Volume = volume
        self.Complexes = []
        self.Removed_complexes = []
        self.Strands = strands
        self.Strands_by_name = defaultdict(list)
        for strand in strands:
            self.Strands_by_name[strand.Name].append(strand)
        print("Strands in self.Strands_by_name:")
        print("\n".join("- %10s: %s species" % (sname, len(strands))
                        for sname, strands in self.Strands_by_name.items()))
        self.VERBOSE = verbose
        self.Print_statsline_when_saving = params.get('print_statsline_when_saving', False)
        # random.choice requires a list, not a set.
        self.Domains_list = [domain for strand in strands for domain in strand.Domains]
        self.Domains = set(self.Domains_list)
        self.N_domains = len(self.Domains)
        self.N_strands = len(self.Strands)
        self.N_domains_hybridized = sum(1 for domain in self.Domains_list if domain.Partner)
        self.N_strands_hybridized = sum(1 for oligo in self.Strands if oligo.is_hybridized())
        self.Domains_by_name = defaultdict(list)
        for d in self.Domains:
            self.Domains_by_name[d.Name].append(d)
        print("Domains in self.Domains_by_name:")
        print("\n".join("- %10s: %s species" % (dname, len(domains))
                        for dname, domains in self.Domains_by_name.items()))
        if domain_pairs is None:
            # mapping: dom_a -> dom_A, dom_A -> dom_a
            # TODO: This could perhaps be a list, if you want to have different types of domains interacting,
            # E.g. dom_a could be perfect match for dom_A, while dom_ax has 1 mismatch.
            # But then you would also have to adjust Domain_dHdS to account for different types of interacting domains.
            domain_pairs = {d.Name: d.lower() if d.Name == d.upper() else d.upper() for d in self.Domains_list}
        assert not any(k == v for k, v in domain_pairs.items())
        self.Domain_pairs = domain_pairs
        self.Uppers = [d for d in self.Domains if d == d.upper()]
        self.Lowers = [d for d in self.Domains if d == d.lower()]
        # Standard enthalpy and entropy of hybridization,
        # indexed as [<domain-name-uppercase>][0 or 1]
        self.Domain_dHdS = {}
        #self.Visualization_hook = self.print_domain_hybridization_percentage
        #self.Visualization_hook = self.randomly_print_stats
        self.Compl_dom_selection = "conc"
        self.Default_statstypes = ("timesampling", "changesampling")
        # Visualization_hook - can easily become rather obtrusive and verbose.
        self.Visualization_hook = None
        # viz-hook options: self.save_stats_if_large
        #self.Visualization_hook = self.print_complexes  # print the complexes at every change.

        if isinstance(outputstatsfiles, str):
            base, ext = os.path.splitext(outputstatsfiles)
            outputstatsfiles = {k: base+"_"+k+'.csv' for k in self.Default_statstypes}
        self.Outputstatsfiles = outputstatsfiles    # dict with <stats type>: <outputfilepath> entries
        self.Record_stats = params.get("record_stats", True)    # Enable or disable stats recording
        # Record stats every N number of steps.
        self.Timesampling_frequency = self.Params.get('timesampling_frequency', 10)
        self.Oversampling = self.Params['probablity_oversampling_factor']
        self.N_steps_per_T = params.get('n_steps_per_T', 100000)
        self.N_steps = 0    # Total number of steps
        self.N_changes = 0  # Number of state changes (hybridizations or de-hybridizations)
        self.N_selections = 0 # Total number of succeessfull selection of domain1 and domain2.
        # Save stats to a cache and only occationally append them to file on disk.
        # Stats_cache: dict <stats type>: <stats>, where stats is a list of tuples:
        # [(Temperature, N_dom_hybridized, %_dom_hybridized, N_oligos_hybridized, %_oligos_hybridized), ...]
        self.Stats_cache = {k: [] for k in self.Default_statstypes}
        self.Complex_size_stats = {k: [] for k in self.Default_statstypes}
        print("Simulator initiated at V=%s with %s strands spanning %s domains." \
              % (self.Volume, len(self.Strands), len(self.Domains)))
        self.print_setup()


    def print_setup(self, fp=None):
        """ Print the simulator setup. """
        c = concentration = 1/N_AVOGADRO/self.Volume
        print("Simulator setup / parameters:", file=fp)
        print("Total concentration of each strand:", file=fp)
        print("\n".join("- {:10}: {:0.3g} uM (N={})".format(name, len(entries)*c*1e6, len(entries))
                        for name, entries in self.Strands_by_name.items()), file=fp)
        print("Concentration of each domain is: %0.3g uM" % (c*1e6), file=fp)
        print("Domain pairing map:", ", ".join("->".join(kv) for kv in self.Domain_pairs.items()), file=fp)
        print("Total number of strands:", self.N_strands, file=fp)
        print(" - number of different strands:", len(self.Strands_by_name), file=fp)
        print("Total number of domains:", self.N_domains, file=fp)
        print(" - number of different domains:", len(self.Domains_by_name), file=fp)
        print("Steps per T:", self.N_steps_per_T, file=fp)
        oversampling_factor = self.Params['probablity_oversampling_factor']
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
        c_each = 1/(N_AVOGADRO*self.Volume)
        for sname, strands in self.Strands_by_name.items():
            c_strand_total = len(strands)*c_each
            # Perhaps also calculate the concentration of dnac2 explicitly
            # (insted of assuming it to be the same as c_strand_total)
            print("\n[ %s strand ]" % sname, file=fp)
            print("c_each: {:0.04g} uM".format(c_each*1e6), file=fp)
            print("n_species (copies):", len(strands), file=fp)
            print("c_strand_total: {:0.04g} uM".format(c_strand_total*1e6), file=fp)
            for domain in strands[0].Domains:
                print("Domain:", domain.Name, file=fp)
                print(" - total copy count of this domain:", len(self.Domains_by_name[domain.Name]), file=fp)
                print("\n".join(" - %s: %s" % (att, getattr(domain, att))
                                for att in ('Sequence', )), file=fp)
                try:
                    deltaH, deltaS = self.Domain_dHdS[domain.Name]
                    print(" - deltaH, deltaS: {:.04g} kcal/mol, {:.04g} cal/mol/K".format(
                        *self.Domain_dHdS[domain.Name]), file=fp)
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
        count = sum(1 for domain in self.Domains_list if domain.Partner)
        if not count % 2 == 0:
            print("Weird - n_hybridized_domains counts to %s (should be an even number)" % count)
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.Domains_list if domain.Partner))
        return count

    def n_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for oligo in self.Strands if oligo.is_hybridized())


    def get_complexes_and_strands(self):
        """ Return all complexes and strands that are in a complex. """
        complexed_strands = [strand for strand in self.Strands if strand.Complex]
        complexes = {strand.Complex for strand in complexed_strands}
        return complexed_strands, complexes

    def get_potential_partners(self, dom1):
        """
        Return a list of possible binding partners for domain <dom1>.
        The activity or biding energy of the individual domains is not necessarily the same.
        """
        candidates = [d for d in self.Domains_by_name[self.Domain_pairs[dom1.Name]] if not d.Partner]
        return candidates

    def select_event_domains(self):
        """
        Returns a tuple with (domain1, domain2, is_hybridized, compl_activity).
        If domain2 is None, then no domain was found for domain1 and nothing should happen.
        Notes:
            For strandA + strandB --> duplex reaction,
            at T=Tm (for strandA), K = [duplex]/([strandA]*[strandB] = 1/[strandB]
            since [duplex] = [strandA] = [strandA]_init/2
            (provided the initial amount of strandB is equal to or greater than strandA)

        """
        is_hybridized = None
        dom1 = random.choice(self.Domains_list)
        if dom1.Partner:
            # TODO: Consider flipping a coin here, since [duplex] = 0.5 * [duplexed domain or complement]
            # Note: This will decrease p_melt
            # dom2 = dom1.Partner if random.random() < 0.5 else None
            return (dom1, dom1.Partner, True, 1)

        ## Note: The steps speed up considerably when reaching low temperatures where most domains are hybridized.
        ## That indicates that the steps below takes up considerable computation time.
        ## I imagine calling d.effective_activity when calculating domain_weights takes considerable time.
        ## However, this should only be slow when the two domains are in the same complex and we have to calculate
        ## the distance between the two domains.
        ## Note sure about np.random.choice(...), but I imagine that is fairly fast.
        ## OK, think I optimized it a bit. Now simulations slows down when there are lots of state *changes*.

        # Only consider complementary domains that are not already paired up:
        candidates = [d for d in self.Domains_by_name[self.Domain_pairs[dom1.Name]] if not d.Partner]
        # Whether we find a candidate depends on (1) the number of candidates and (2) the volume.
        # empty_volume_number should be scaled with the tube volume.
        # empty_volume_number = self.Params.get("empty_volume_number", 100)
        # Edit: Using volume and numbers to get an actual concentration:
        # c2 = 1/N_AVOGADRO/self.Volume     # Edit2: concentration is determined for each domain individually.

        # Need to scale intra-complex domain higher: (higher effective concentration)
        # Doing this by making a list with repeated entries:
        # (This is only able to scale by integers, but should suffice for now...)
        #if dom1.Complex and any(dom1.Complex == d.Complex for d in candidates):
        #    candidates = sum(([d]*d.effective_conc(dom1) if d.Complex == dom1.Complex else [d]
        #                      for d in candidates), []) # Need to give a list as starting value to sum

        # Simple way: just add the two together (then we won't have to worry about len(candidates) >= empty_volume_number)
        # If empty_volume_number = 0, there should be 100% chance of selecting a candidate, so subtract 1:
        #n_cands = len(candidates)
        #event_number = random.randint(0, empty_volume_number + len(candidates) - 1)
        # More complex: empty_volume_number is the "total volume" of which the domains can occupy a certain portion.
        #if len(candidates) < empty_volume_number:
        #    event_number = random.randint(0, empty_volume_number)
        #    # event_number = random.randint(0, empty_volume_number - len(candidates)) # No...?
        #else:
        #    # Will most likely never happen
        #    event_number = random.randint(0, len(candidates) - 1)
        #if event_number >= len(candidates):
        #    # Nothing should happen, dom1 did not find a partner in the "allotted time".
        #    return (dom1, None)
        #dom2 = candidates[event_number]


        ### New concentration-based weighting and selection:
        # Note: oversampling is only used if the domains are not in the same complex:
        # return oversampling*c, c = 1/(N_AVOGADRO*volume)
        # Are we sure oversampling can be applied linearly?
        domain_weights = [d.effective_activity(dom1, volume=self.Volume, oversampling=self.Oversampling)
                          for d in candidates]

        # Adding constant water solvent option:
        #candidates.append(None)
        #domain_weights.append(1)    # 1 or maybe even 56 = the molarity of water?
        domain_weights_sum = sum(domain_weights)

        # Do you add a specific "water/None" option?  (yes, if domain_weights_sum < 1)
        # Or do you first calculate "is a domain gonna be selected?" and then select which domain is selected?
        # Adding a specific "None" option as "water":
        ## Consider maybe always adding the "water/solvent" option...?
        if domain_weights_sum < 1:
            # Maybe insert at pos 0 instead of the end, depending on the implementation of weighted choice?
            # Result: Doesn't make much difference.
            #domain_weights.append(1-domain_weights_sum)
            #candidates.append(None)
            domain_weights.insert(0, 1-domain_weights_sum)
            candidates.insert(0, None)
        if domain_weights_sum > 1:
            # Normalize...:
            domain_weights = [w/domain_weights_sum for w in domain_weights]
        if np:
            # http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.choice.html
            dom2 = np.random.choice(candidates, p=domain_weights)
        else:
            # If you want to do the same without using numpy, you could do something like:
            rvalue = random.random()
            dom2 = None
            # Consider sorting the candidates by weight to increase # Nope, sorting takes considerable time.
            #domain_weights, candidates = zip(*reversed(sorted(zip(domain_weights, candidates)))) # Only if domains are orderable
            #candidates.sort(key=dict(zip(candidates, domain_weights)).get)  # First sort candidates by their weight,
            #domain_weights.sort()                                           # Then sort the weights.
            for i, w in enumerate(domain_weights):    # are normalized
                if rvalue <= w:
                    dom2 = candidates[i]
                    break
                rvalue -= w
            else:
                # Should not happen since sum(w) == 1.
                raise ValueError("No domain found, sum(w) does not add up to 1.")
        #print("Candidates and weights: (dom1=%s)" % dom1)
        #print("\n".join("{:<30}: {:.9g}".format(str(d), w) for d, w in zip(candidates, domain_weights)))
        return dom1, dom2, is_hybridized, domain_weights_sum


    def hybridization_energy(self, domain1, domain2, T):
        """
        Calculate standard Gibbs free energy of hybridization between domain1 and domain2,
        assuming that domain1 and domain2 does not form any other structures in their single-stranded state.
        Note that this assumption does often not hold true.

        Is there a linear or just simple dependence for dG on temperature?
        Note: The values of DNA_NN tables are in
            kcal/mol for enthalpy  (the first entry)
            cal/mol/K for entropy  (the second entry)
        The Nearest Neighbor model calculates ΔG° as:
            ΔG°(total) = ΔG°(init) + ΔG°(symmetry) + sum(ΔG°(stacking)) + ΔG°(AT terminal)
        The stacking parameters are given as DNA_NN3:
            e.g. AA/TT is for the stacking between AA
                                                   TT
            The directionality does not seem to matter, e.g. 5'-AC is same as 5'-CA.
        Where the Gibbs free energy is (assuming ΔCp° = 0, so ΔH° and ΔS° are temperature independent):
            ΔG° = ΔH° - T*ΔS°
        This means you can contract all entropy terms:
            ΔG°(total) = ΔH°(init) + ΔH°(symmetry) + sum(ΔH°(stacking)) + ΔH°(AT terminal)
                        - T * (ΔS°(init) + ΔS°(symmetry) + sum(ΔS°(stacking)) + ΔS°(AT terminal))
                       = ΔH°(total) - T*ΔS°(total)
        If you have saved ΔH°(total) and ΔS°(total), then calculating ΔG°(total) is fast.
        It would, perhaps, be nice to have the updated values from David Zhang..
        """
        # returned is (5p domain-exists, 5p domain-hybridized, 3p domain-exists, 3p domain-hybridized)
        # d1_params, d2_params = domain1.neighbor_state_tuple(), domain2.neighbor_state_tuple()
        # param_tuple = d1_params + d2_params
        # indexed as [<domain-name-uppercase>][<8-bit-tuple>][0 or 1]
        # Edit: Currently just saving the dH_dS result as [<domain-name-uppercase>]
        # and performing neighbor-induced adjustments on the fly...
        # Do not use dict.setdefault(key, energy-if-no-key()) for this!
        # This will call energy-if-no-key() before checking whether keys is in the dict.
        if domain1.Name not in self.Domain_dHdS:
            print(domain1.Name, "domain not in self.Domain_dHdS - calculating...")
            self.Domain_dHdS[domain1.Name] = deltaH, deltaS = \
                hybridization_dH_dS(domain1.Sequence, domain2.Sequence, nn_table=self.NN_Table)
            print(domain1.Name, "deltaH, deltaS = {:.04g}, {:.04g}".format(deltaH, deltaS),
                  "({} entries in self.Domain_dHdS".format(len(self.Domain_dHdS)))
        else:
            deltaH, deltaS = self.Domain_dHdS[domain1.Name]

        deltaH_NN_domain, deltaS_NN_domain = None, None
        NN_domains = [False, False]
        deltaH_corr, deltaS_corr = (0, 0)
        #return deltaH*1000 - T * deltaS, deltaH, deltaS, deltaH_corr, deltaS_corr

        ## TODO: If we split a large domain into two half-sized domains, we should get the same result.
        if domain1.domain5p() and domain2.domain3p() \
            and domain1.domain5p().Partner and domain2.domain3p().Partner \
            and domain1.domain5p().Partner == domain2.domain3p():
            # Need to add another nearest-neighbour contribution.
            # (but then somehow not do the dangling ends contributions below)
            # 'GA/CT' is equivalent to the stacking between 5'-G:A-(domain1)-3'
            # and                                           3'-C:T-(domain2)-5'
            # In this if we are checking on the 5p end of domain1:
            NN_stack = "{}{}/{}{}".format(domain1.domain5p().Sequence[-1],
                                          domain1.Sequence[0],
                                          domain2.domain3p().Sequence[0],
                                          domain2.Sequence[-1])
            deltaH_NN_domain, deltaS_NN_domain = self.NN_Table[NN_stack]
            NN_domains[0] = True
            deltaH_corr += deltaH_NN_domain
            deltaS_corr += deltaS_NN_domain
        if domain2.domain5p() and domain1.domain3p() \
            and domain2.domain5p().Partner and domain1.domain3p().Partner \
            and domain2.domain5p().Partner == domain1.domain3p():
            NN_stack = "{}{}/{}{}".format(domain2.domain5p().Sequence[-1],
                                          domain2.Sequence[0],
                                          domain1.domain3p().Sequence[0],
                                          domain1.Sequence[-1])
            deltaH_NN_domain, deltaS_NN_domain = self.NN_Table[NN_stack]
            NN_domains[1] = True
            deltaH_corr += deltaH_NN_domain
            deltaS_corr += deltaS_NN_domain

        # Unless all nearest neighbor domain are forming regular bp stack,
        # we have to calculate corrections to deltaH and deltaS:
        if not all(NN_domains):
            ### Dangling adjustment from electrostatic repulsion ###
            # Seems like this should be sequence dependent, based on DNA_DE1 lookup table.
            # There is even a difference whether a nucleotide is on the 5p end or 3p... sigh...
            # Use dG = dH - dS*T with T = 330 K (57 C) to get a feeling for the contributions.
            # According to the DNA_DE1 table, dangling ends actually repel primarily via entropic effects?
            # - both dH and dS are negative.
            # Actually, it seems most dangling ends contribute a negative ΔΔG to hybridization energy:
            # DE_dG = {k: (v[0]-0.300*v[1], v[0]-0.360*v[1]) for k, v in DNA_DE1.items()}   # from about 30 C to 90 C:
            # [pyplot.plot([300, 360], v) for v in DE_dG.values()]  # followed by display(*getfigs()), with %pylab enabled
            # In all cases the ΔΔG is between -1 and +1 kcal/mol.
            # FWIW, -1 kcal/mol corresponds to a change in K of exp(-ΔG/RT) = exp(ΔS/R-ΔH/RT) = 4.6
            # At 60 C: RT = 1.987 cal/mol/K * 330 K = 0.65 kcal/mol.
            # Hmm... It also seems to me that the effect of electrostatic repulsion should be salt dependent.
            N_neighbors = sum(bool(v) for v in (
                (domain1.domain5p() and not NN_domains[0],
                 domain2.domain3p() and not NN_domains[0],
                 domain2.domain5p() and not NN_domains[1],
                 domain1.domain3p() and not NN_domains[1])))
            deltaH_corr = -3.0 * N_neighbors    #  -3 cal/mol   electrostatic repulsion for each neighboring domain
            deltaS_corr = -10.0 * N_neighbors   # -10 cal/mol/K electrostatic repulsion for each neighboring domain
            # Domain repulsion: ΔΔG of 0 kcal/mol at 300 K (30 C) and +0.4 kcal/mol at 340 K (70 C) - seems reasonable.

            ### Stacking interactions when the neighboring domain is hybridized:
            # Again, this should probably be sequence dependent:
            # TODO: Consider whether this can simply be the same as the domain NN stacking contributions above:
            N_stacking = (domain1.domain5p_is_hybridized() or domain2.domain3p_is_hybridized() +
                          domain2.domain5p_is_hybridized() or domain1.domain3p_is_hybridized())
            # A CG/GC NN has dH, dS: 'CG/GC': (-10.6, -27.2), while AT/TA has: 'AT/TA': (-7.2, -20.4)
            deltaH_corr += -7.0 * N_stacking    # -10 cal/mol/K for each stacking interaction
            deltaS_corr += -20.0 * N_stacking   # -20 cal/mol   for each stacking interaction
            # Stacking: ΔΔG of -1 kcal/mol at 300 K (30 C) and -0.2 kcal/mol at 340 K (70 C) - seems reasonable.

            # It is probably a lot better to compare the two states directly using NuPack...!
            # (or at least see what they do and do the same...)

            # ΔG° = ΔH° - T*ΔS°             # Standard Gibbs free energy
            # ΔG = ΔG° + RT ln(Q)           # Gibbs free energy at non-equilibrium
            # ΔG° = ΔH° - T*ΔS° = -RT ln(K) # At equilibrium ΔG=0. Notice the minus in front of RT!

            # You may want to correct the deltaS depending on how the hybridization connects the strands.
            # I.e. how domain hybridization will affect the overall entropy of the complex.
            #  - strand hybridization will reduce the conformational freedom of the complex.
            # This is a dynamic function of the current complex structure and cannot be stored for later use.
            # Obviously, this is only applicable if the two domains are already in the same complex.
            if domain1.Complex and domain1.Complex == domain2.Complex:
                # TODO: Calculate a better entropy correction for intra-complex hybridizations
                #pass
                deltaS_corr += 2 # exp(4/R) = 7.4 !

        ### Add corrigating factors
        deltaH += deltaH_corr
        deltaS += deltaS_corr

        # The deltaH and deltaS calculated by melting_dH_dS is based on SantaLucia1997 - DNA_NN3.
        # The initial entropy of formation are rather negligible.
        deltaG = deltaH*1000 - T * deltaS    # Values in cal/mol - rather than kcal/mol

        return deltaG, deltaH, deltaS, deltaH_corr, deltaS_corr






    def step(self, T):
        """
        Perform a single step in the simulation at temperature T.
        """
        domain1, domain2, is_hybridized, compl_activity = \
            self.select_event_domains()
        if not domain2:
            # If we are selecting for an activity of 1 (or equal to the duplex activity, whatever that is),
            # then there must be a very high probability of NOT finding a domain2. Like, 1 in a million..
            # And, similarly, if we have a duplex, there is a 100% chance of selecting that duplex, but
            # then, even at T=Tm, there is only a 1 in a million chance that the duplex will melt.
            # In other words, we have to perform a lot of "selections" before two strands will hybridize,
            # and we have to perform a lot of p_hyb < random.random() tests before a duplex will melt.
            # This may or may not be suitable for simulation.
            # Also, as you can see, it becomes important that the selection probability is correct
            # (given that p_hyb is close to 1 once two strands have been selected).
            # Instead of trying to get the selection probability right, it might be more reliable
            # to delegate this part to hybridization_probability() using Q to correct for non-unity activity.
            # This will still get a lot of failed attempts, failing at p < random.random() instead of during
            # duplex2 selection, but it is probably a lot easier to get right for a novice like me :)
            # https://en.wikipedia.org/wiki/Thermodynamic_activity
            #print("Failed to select a 2nd domain for domain 1.")
            #sys.stdout.write(".")  # No, even this is too much.
            return

        # Assertion check (at least while we are still debugging and ensuring compliancy).
        # Eventually we might allow self-complementary domains... (I.e. two different domain objects but with same name)
        assert domain1 != domain2 and domain1.Name != domain2.Name

        # Selection takes care of the concentration-dependent probability of bringing non-complexed
        # strands together to Q=1:
        #p_hyb = self.hybridization_probability(domain1, domain2, T)

        self.N_selections += 1
        state_change = False
        # r_hyb = 0.1         # If two domains are selected and not hybridized, what is the chance they will hybridize?
        # deltaG, deltaH, deltaS, deltaH_corr, deltaS_corr
        # dG_std, dH_std, dS_std, dH_corr, dS_corr
        #dG_std, *rest = \
        dG_std, dH_std, dS_std, dH_corr, dS_corr = \
            self.hybridization_energy(domain1, domain2, T)    # return value in cal/mol
        # The returned dH_std and dS_std are *after* adding dH_corr or dS_corr respectively.
        #if dH_corr or dS_corr or dG_std != dH_std - T*dS_std:
        #    print("Unexpected values for dG_std, dH_std, dS_std, dH_corr, dS_corr:")
        #    print(", ".join("%0.5g" % val for val in (dG_std, dH_std, dS_std, dH_corr, dS_corr)))
        #    print("dH_std - T*dS_std = %0.6g" % (dH_std*1000 - T*dS_std))
        #    raise AssertionError("Unexpected event")
        # Old p_hyb model:
        # p_hyb_old = 1 / (1 + math.exp(dG_std/(R*T))*Q)    # always used Q=1 after selection.
        # compl_activity ~0...1 so dividing by compl_activity makes the denominator larger = smaller p_hyb.
        # oversampling also makes the denominator larger -> smaller p_hyb
        #p_hyb_old = 1 / (1 + math.exp(dG_std/(R*T))/compl_activity*self.Oversampling) # oversampling might be included in compl_activity.
        # exp(dG/(R*T) = exp(dH/RT-dS/R) = exp(dH/RT)/exp(dS/R)
        #p_hyb_old = 1 / (1 + math.exp(dG_std/(R*T))*self.Oversampling*10/compl_activity)
        p_hyb_old = 1 / (1 + math.exp(dH_std*1000/(R*T)-dS_std/R)*self.Oversampling/compl_activity)
        # We have shown that p_mel = 1 - p_hyb and p_mel = exp(+ΔG°/RT) * p_hyb
        # is in agreement with p_hyb = x * 1/(1+exp(+ΔG°/RT)
        if is_hybridized:
            # Are we sure that a linear 0...1 choice is the best probability distribution?
            # Since math.exp(dG_std/(R*T)) can be > 1
            # r_mel = exp(+ΔG°/RT) * r_hyb = exp(+ΔG°/RT) * 0.1
            # Here we are looking at melting probability/rate.
            #if random.random() < math.exp(dG_std/(R*T))*r_hyb*0.5*self.Oversampling:
            if random.random() >= p_hyb_old:
                new_dist, new_complexes, obsolete_complexes = domain1.dehybridize(domain2)
                state_change = True
        else:
            #
            # possibly add random.random() > exp(+ΔG°/RT) dependency ?
            # Usually, I would say that if two domains have been selected, they should (almost) always be hybridized.
            # But, that introduces a "sampling bias", since the state might be very short lived.
            # Or maybe add check for compl_activity (Q) ?
            # You could also add oversampling to
            # K = exp(dG/RT) = exp(dG/(RT*oversampling)) = exp(dG/RT)^(1/oversampling)
            # This would make K larger for negative dG and smaller for positive dG.
            #if random.random() > math.exp(dG_std/(R*T))*r_hyb*self.Oversampling: # sim #41
            if random.random() < p_hyb_old:
                new_dist, new_complexes, obsolete_complexes = domain1.hybridize(domain2)
                state_change = True

        if state_change:
            self.N_changes += 1
            self.N_domains_hybridized += -2 if is_hybridized else 2

            if new_complexes:
                for c in new_complexes:
                    if c in self.Complexes:
                        print("WEIRD: new complex %s is already in self.Complexes:" % c)
                        print(self.Complexes)
                    self.Complexes.append(c)
            if obsolete_complexes:
                for c in obsolete_complexes:
                    if self.Complexes.count(c) > 1:
                        print("WEIRD: obsolete complex %s is present %s times in self.Complexes:" \
                              % (c, self.Complexes.count(c)))
                        print(self.Complexes)
                    try:
                        self.Complexes.remove(c)
                    except ValueError:
                        print("Error removing complex", str(c))
                        print(" - Is complex in self.Removed_complexes? -", c in self.Removed_complexes)
                        if self.VERBOSE > 1:
                            print(" - domain1:", domain1)
                            print(" - domain2:", domain2)
                            print(" - len(self.Complexes):", len(self.Complexes))
                            print(" - self.Complexes:", self.Complexes)
                            print(" - len(obsolete_complexes):", len(obsolete_complexes))
                            print(" - obsolete_complexes:", obsolete_complexes)
                            print(" - new_complexes:", new_complexes)
                            print(" - is_hybridized:", is_hybridized)
                            print(" - new_dist:", new_dist)
                            if c:
                                print("c.Strands:", c.Strands)
                                print("c.Connections:", c.Connections)
                                print("c.N_strand_changes:", c.N_strand_changes)
                                print("c.Strands_changes:", c.Strands_changes)
                                print("c.Strands_history:", c.Strands_history)
                            # if is_hybridized, then we've used domain1.dehybridize(domain2)
                    else:
                        self.Removed_complexes.append(c)
            if self.Record_stats:
                self.record_stats_snapshot(T)
            if self.Visualization_hook:
                self.Visualization_hook(updated_domains=(domain1, domain2))

        return


    def simulate(self, T, n_steps_max=100000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        assert self.n_hybridized_domains() == self.N_domains_hybridized
        n_done = 0
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

    def save_stats_if_large(self, **kwargs):
        """ Save stats cache when it is sufficiently large. """
        if any(len(cache) > 10000 for cache in self.Stats_cache.values()):
            self.save_stats_cache()

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
        self.Complex_size_stats[statstype].append({T: complex_sizes_hist(self.Complexes)})


    def randomly_print_stats(self, **kwargs):
        """ Print stats at random intervals. """
        if random.random() < 0.05:
            self.print_domain_hybridization_stats(**kwargs)


    def print_complexes(self, *args, **kwargs):
        """
        For use as visualization hook.
        Note: This one does become rather verbose, especially just below T=Tm.
        """
        p_update = self.Params.get('viz_hook_update_probability', 1)
        if p_update < 1 and random.random > p_update:
            return
        print("- viz-hook: Complexes = ", self.Complexes)

    def print_N_complexes(self, *args, **kwargs):
        """
        For use as visualization hook.
        Note: This one does become rather verbose, especially just below T=Tm.
        """
        p_update = self.Params.get('viz_hook_update_probability', 1)
        if p_update < 1 and random.random > p_update:
            return
        print("- viz-hook: N Complexes = ", len(self.Complexes))


    def print_domain_hybridization_stats(self, updated_domains=None):
        """ Print domain hybridization stats """
        N_domains = len(self.Domains)
        N_hybridized = sum(1 for domain in self.Domains if domain.Partner)
        if self.VERBOSE > 1:
            print("| - Updated domains %s and %s" % updated_domains)
        print("| Total domain hybridization percentage: {:.0%} ({} of {})".format(
            N_hybridized/N_domains, N_hybridized, N_domains))


    def thermodynamic_meltingcurve(self, T_start, T_finish, delta_T=None, volume=None):
        """ Calculate thermodynamic melting curve. """
        if delta_T is None:
            delta_T = -0.5 if T_start > T_finish else +0.5
        if volume is None:
            volume = self.Volume

        # from collections import OrderedDict
        cum_stats = {} #OrderedDict()       # Has T: <fraction of duplexed domains>
        domain_stats = {} #OrderedDict()    # Same as cum_stats, but for individual domains
        concentrations = {dname: len(domains)/(volume*N_AVOGADRO)
                          for dname, domains in self.Domains_by_name.items()}
        T = T_start
        while T >= T_finish if delta_T < 0 else T <= T_finish:
            hybridized = 0      # can be a fraction
            non_hybridized = 0
            total_conc = 0
            domains_total = 0
            domain_stats[T] = {}
            for dname, domains in self.Domains_by_name.items():
                if dname not in self.Domain_pairs or self.Domain_pairs[dname] not in concentrations:
                    # domain has no complementary domains:
                    continue
                if dname not in self.Domain_dHdS:
                    self.Domain_dHdS[dname] = hybridization_dH_dS(domains[0].Sequence)
                # standard-condition energies:
                deltaH, deltaS = self.Domain_dHdS[dname]
                # deltaH in kcal/mol, deltaS in cal/mol/K:
                deltaG = deltaH*1000 - T*deltaS
                R = 1.987  # universal gas constant in cal/mol/K
                try:
                    # If deltaG is really large, we might get overflow error
                    K = math.exp(-deltaG/(R*T))
                except OverflowError:
                    print("Warning: OverflowError while calculating K = math.exp(-deltaG/(R*T)) " +
                          "= math.exp(-{:0.03g}/({}*{})). Setting K=1e10.".format(deltaG, R, T))
                    K = 1e10
                # We assume all partnering domains are perfect complements:
                c_domain = concentrations[dname]
                c_complement = concentrations[self.Domain_pairs[dname]]
                # f_hyb = binary_state_probability_cal_per_mol(deltaG)
                if c_complement > c_domain:
                    # Ai must be the higher of the two concentrations:
                    D, A, B, K2 = solvetwocomponent(c_complement, c_domain, K)
                    x = D/B
                else:
                    #
                    D, A, B, K2 = solvetwocomponent(c_domain, c_complement, K)
                    x = D/A
                if abs(K2-K)/K > 0.1:
                    print(("- Note: solvetwocomponent output K={:0.04g} is different " +
                          "from input K={:0.04g} (T={}, dname={}").format(K, K2, T, dname))
                # x is the fraction of domain that is hybridized:
                # Actually, we just want the concentration of duplexed ?
                # TODO: We could scale by number of bases in each domain...
                hybridized += D         # is a concentration
                non_hybridized += c_domain - D
                total_conc += c_domain
                domains_total += len(domains)
                domain_stats[T]['dname'] = [D, c_domain - D, c_domain, len(domains), x]
            cum_stats[T] = (hybridized, non_hybridized, total_conc, domains_total)
            T += delta_T
        return cum_stats, domain_stats


def solvepol2(a, b, c):
    det = math.sqrt(b**2-4*a*c)
    x1, x2 = (-b + det) / (2*a), (-b - det) / (2*a)
    return x1, x2, det

def report_weird(x1, x2, Ai, Bi):
    print("(x1, x2):", (x1, x2))
    for x in (x1, x2):
        D = x
        A = Ai - D
        B = Bi - D
        print("\nx = D =", D)
        print("A =", A)
        print("B =", B)


def solvetwocomponent(Ai, Bi, K):
    """
     0  = D^2 + (-Bi-Ai-1/K)*D + Ai*Bi
    """
    a = 1
    b = (-Bi-Ai-1/K)
    c = Ai*Bi
    #print("a=%s, b=%s, c=%s" % (a, b, c))
    x1, x2, _ = solvepol2(a, b, c)
    D = None
    if not any(M-x1 < 0 for M in (Ai, Bi)):
        # x1 is a viable solution
        D = x1
    if not any(M-x2 < 0 for M in (Ai, Bi)):
        # x1 is a viable solution
        if D:
            print("Both x1 and x2 seems to be good solutions?!?!")
            report_weird(x1, x2, Ai, Bi)
        D = x2
    if D is None:
        print("Neither x1 nor x2 are viable solutions?")
        report_weird(x1, x2, Ai, Bi)
    A = Ai-D
    B = Bi-D
    K = D/(A*B)
    return D, A, B, K



if __name__ == '__main__':

    import os
    from .dom_utils import parse_strand_domains_file

    strand_defs_file = os.path.join(os.path.dirname(__file__), "testfiles", "strand_defs01.txt")
    input_oligos = parse_strand_domains_file(strand_defs_file)

    # Some calculations:
    # * Consider a volume of 1 um^3 = (1e-5 dm)^3 = 1e-15 L = 1 femto-liter.
    # * For c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
    #     N_Avogagro = 6.02e23 /mol

    #
    adhoc_params = {}
    simulator = Simulator(volume=1e-15, strands=input_oligos, params=adhoc_params)
