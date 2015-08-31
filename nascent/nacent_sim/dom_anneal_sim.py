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


"""


import random
#from random import choice
from collections import defaultdict
from math import log, exp

from .dom_anneal_models import Tube, Strand, Domain, Complex
from .energymodels.biopython import binary_state_probability, hybridization_dH_dS


N_Avogadro = 6.022e23


class Simulator():

    def __init__(self, volume, strands, params):
        self.Tube = Tube(volume, strands)
        self.Volume = volume
        self.Strands = strands
        self.Domains = [strand.Domain for strand in strands]
        self.Domains_by_name = defaultdict(list)
        for d in self.Domains:
            self.Domains_by_name[d.Name].append(d)
        self.Uppers = [d for d in self.Domains if d == d.upper()]
        self.Lowers = [d for d in self.Domains if d == d.lower()]
        self.Params = params
        #self.Domain_dG = {}  # currently indexed as [<domain-name-uppercase>][<temp>][<8-bit-tuple>]
        # Standard enthalpy and entropy, indexed as [<domain-name-uppercase>][<8-bit-tuple>][0 or 1]
        self.Domain_dHdS = {}

    def select_event_domains(self):
        """
        Returns a tuple with (domain1, domain2).
        If domain2 is None, then no
        """
        dom1 = random.choice(self.Domains)
        if dom1.Partner:
            return (dom1, dom1.Partner)
        # Only consider domains that are not already paired up:
        candidates = [d for d in self.Domains_by_name[dom1.Name] if not d.Partner]
        # Whether we find a candidate depends on (1) the number of candidates and (2) the volume.
        # empty_volume_number should be scaled with the tube volume.
        empty_volume_number = self.Params.get("empty_volume_number", 100)

        # Need to scale intra-complex domain higher: (higher effective concentration)
        # Doing this by making a list with repeated entries:
        # (This is only able to scale by integers, but should suffice for now...)
        if dom1.Complex and any(dom1.Complex == d.Complex for d in candidates):
            candidates = sum(([d]*d.effective_conc(dom1) if d.Complex == dom1.Complex else [d]
                              for d in candidates), []) # Need to give a list as starting value to sum

        # Simple way: just add the two together (then we won't have to worry about len(candidates) >= empty_volume_number)
        # If empty_volume_number = 0, there should be 100% chance of selecting a candidate, so subtract 1:
        n_cands = len(candidates)
        event_number = random.randint(0, empty_volume_number + len(candidates) - 1)
        # More complex: empty_volume_number is the "total volume" of which the domains can occupy a certain portion.
        if len(candidates) < empty_volume_number:
            event_number = random.randint(0, empty_volume_number)
            # event_number = random.randint(0, empty_volume_number - len(candidates)) # No...?
        else:
            # Will most likely never happen
            event_number = random.randint(0, len(candidates) - 1)
        if event_number >= len(candidates):
            # Nothing should happen, dom1 did not find a partner in the "allotted time".
            return (dom1, None)
        dom2 = candidates[event_number]
        return (dom1, dom2)


    def calculate_energy(self, domain, configuration):
        pass


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
        d1_params, d2_params = domain1.neighbor_state_tuple(), domain2.neighbor_state_tuple()
        param_tuple = d1_params + d2_params

        # indexed as [<domain-name-uppercase>][<8-bit-tuple>][0 or 1]
        # Edit: Currently just saving the dH_dS sequence result as [<domain-name-uppercase>]
        # Will perform neighbor adjustments on the fly...
        deltaH, deltaS = self.Domain_dHdS.setdefault(domain1.Name.upper(), #{})\
                                         #.setdefault(param_tuple,       # WAIT, this does not capture or save any dependence on the param_tuple!!
                                                     melting_dH_dS(domain1.Sequence,
                                                                    domain2.Sequence))
        # Perform neighbor adjustments on the fly...
        # Dangling adjustment... (electrostatic repulsion)
        # Seems like this should be sequence dependent, based on DNA_DE1 lookup table.
        # There is even a difference whether a nucleotide is on the 5p end or 3p... sigh...
        deltaH_corr =  -3.0 * sum(param_tuple[idx] for idx in range(0, 8, 2))
        deltaS_corr = -10.0 * sum(param_tuple[idx] for idx in range(0, 8, 2))
        # Stacking interactions when the neighboring domain is hybridized:
        deltaH_corr += -10.0 * sum(param_tuple[idx] for idx in range(1, 8, 2))
        deltaS_corr += -20.0 * sum(param_tuple[idx] for idx in range(1, 8, 2))

        # It is probably a lot better to compare the two states directly using NuPack...!


        # ΔG° = ΔH° - T*ΔS°             # Standard Gibbs free energy
        # ΔG = ΔG° + RT ln(Q)           # Gibbs free energy at non-equilibrium
        # ΔG° = ΔH° - T*ΔS° = -RT ln(K) # At equilibrium ΔG=0. Notice the minus in front of RT!

        # You may want to correct the deltaS depending on how the hybridization connects the strands.
        # This is a dynamic function of the current complex structure and cannot be stored for later use.
        deltaS_corr += 0

        # Add corrigating factors
        deltaH += deltaH_corr
        deltaS += deltaS_corr

        # The deltaH and deltaS calculated by melting_dH_dS is based on SantaLucia1997 - DNA_NN3.
        # The initial entropy of formation are rather negligible.
        deltaG = deltaH*1000 - T * deltaS    # Values in cal/mol - rather than kcal/mol

        return deltaG, deltaH, deltaS





    def hybridization_probability(self, domain1, domain2, T, Q=None):
        """
        Calculate the partition

        SantaLucia, Dirks (Ann Rev Biophys, 2004): "We note that the database presented is not appropriate
        for partition function computations (50; J. SantaLucia, unpublished results)." - WTF?
        Ahh... they are only good for determining Tm.. Not for partition functions far from Tm.
            - The most likely reason is that either there's some ensemble effects that are not properpy captured
            by the NN models, or the values are just not accurate enough for partition calculations.

        SantaLucia, Dirks (Ann Rev Biophys, 2004): "Note that many duplexes have competing single-strand structure,
        and this compromises the validity of the two-state approximation and results in systematically lower TMs than
        would be predicted by Equation 3 [the equation used to predict Tm]."

        Perhaps I should look at how NuPack calculates its partition functions...
        NuPack refs:
            http://www.nupack.org/home/model
            "The free energy of an unpseudoknotted secondary structure is calculated using nearest-neighbor empirical
            parameters for RNA in 1M Na+ (Serra and Turner, 1995; Mathews et al., 1999) or DNA in user-specified
            Na+ and Mg++ concentrations (SantaLucia, 1998; SantaLucia and Hicks, 2004; Koehler and Peyret, 2005)"
        """
        deltaG, deltaH, deltaS = self.hybridization_energy(domain1, domain2, T)    # return value in cal/mol
        #R = 1.987  # universal gas constant in cal/mol/K
        ## Note: R = k * NA , where k is the boltzmann constant and NA is Avogadro's constant.
        #K = exp(-deltaG/(R*T))
        ## Not really sure how to convert from equilibrium constant K to probability that they will hybridize
        ##   the deltaG used to calculate the equilibrium constant is at standard conditions (1 M) !
        ## Probably refer to NuPack or similar to make sure you get it right...
        ## Also remember that K = k_on/k_off
        ##                         k_on              [1/M/s]      # v_on = k_on [strandA] [strandB]
        ##     domainA + domainB --------> duplex
        ##                       <--------
        ##                         k_off             [1/s]       # v_off = k_off [duplex]
        #if domain1.Complex and domain2.Complex and domain1.Complex == domain2.Complex:
        #    # What if we are breaking the bond, and the bond is the only bond holding the domains in the same complex?
        #    # Well, that is a k_off rate, that shouldn't be affected by concentrations anyways.
        #    c = 1
        #else:
        #    # Different complexes
        #    # 1 nM is about 1 molecule per femto-liter: c = 1 nM: n = 1e-15 L * 1e-9 mol/L = 1e-24 mol = 0.6
        #    # N_Avogadro = 6.022e23/mol
        #    c = 1/N_Avogadro/self.Volume
        #
        ## Fuck that, I can't see how to include the concentration in the probability.
        ## As I see it, it should already be included
        #p = K/(1+K)     # If K = 7/2, then p = 7/(7+2) = K/(K+1) = 1/(1+1/K)
        #
        ## Edit: How to include concentration in probability:
        ## 1) Calculate non-standard ΔG: ΔG = ΔG° + RT ln(Q)
        ##       where Q = [duplex]/([strandA] [strandB])    # If concentrations are around 1 uM, then Q = 1e6.

        # Note that the concentration is obviously very important. If you have a duplex at T=Tm (so p_duplex = 50%),
        # at a concentration of 1 uM. If you then increase the concentration/activity to 1 M, then suddenly
        # the p_duplex is something like 99.9998 % !
        # Another way to see this is to compare Tm at c = 25 nM with Tm at c = 1 M.
        # -- typically increases about 30 degC !!
        # If you want to simulate this, then you need to incorporate this at the earlier "selection" step.
        # This essentially boils down to setting a very high empty_volume_number in select_event_domains() above.
        # Alternatively, you can use a lower empty_volume_number (and thus get more selections), but
        # then you have to compensate using a correcting Q.
        # For instance, instead of saying "we put the strands so close together that they have an activity of 1",
        # you could say, "we put them within a distance so their concentrations are 1 mM", and then give Q = 1e3.
        # Note: Requiring an activity of 1 before hybridization can occour might be way too high:
        # DNA strands has a certain extend, and can interact at great distances.
        # Of course, they have to eventually end up in a duplex state, where we can assume an activity of 1...
        # But that might still yield kinetics that are quite wrong.. Maybe they don't have an activity of 1 in the
        # duplex state?
        # Note: There is also a certain probability that the standard ΔH°, ΔS°, and ΔG° values are not suited for
        # these types of partition calculations (they've stated that them self).
        # Again: Try to look at what NuPack does.
        # Edit: Essentially, we just need to assure that at T=Tm, p_on == p_off.
        # p_on is defined at two stages: selection and hybridization (after it has been selected).
        #       p_on = p_selection * p_hyb
        # There are only two easy ways to guarantee that you get the correct p_on:
        # 1) set p_selection = 1 -- and use concentration to calculate Q which is then used when calculating p_hyb.
        # 2) Use concentration to calculate p_selection, and use the same p_hyb for both hybridization and melting.
        #
        # Another concern is that for low concentrations (1 uM - 1 nM), p_on and p_off will both be very, very
        # low -- like 1 in a millionth to 1 in a billionth. This requires a lot of fruitless "dice rolling",
        # before the system actually changes state. And that is at T=Tm...
        # One way to mitigate this is to do oversampling (or whatever it should be called):
        # Multiplying both p_on and p_off by a certain probablistic_oversampling_constant, to ensure that something
        # actually happens within a reasonable timeframe.
        # The thing here is that you need to take state into account:
        # If the strand/domain is currently hybridized, you need to increase the chance that it will melt.
        # If it is not currently hybridized, then increase the chance that it (1) find a partner,
        # (2) to which it will hybridize.
        if Q is None:
            Q = 1   # Do selection in select_event_domains()
            #Q = 1e3 # select_event_domains() selects to 1 mM.
            #Q = 1e6 # select_event_domains() selects to 1 uM.
            # Q = (concentration) # If p_selection = 1

        p_hyb = binary_state_probability(deltaG, T, Q=Q)

        return p_hyb


    def step(self, T):
        """
        Perform a single step in the simulation at temperature T.
        """
        domain1, domain2 = self.select_event_domains()
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
            print("Failed to select a 2nd domain for domain 1.")
            return
        print("Selected domain1 and 2:", domain1, domain2)
        # Probability that the two strands are in hybridized state:
        p = self.hybridization_probability(domain1, domain2, T)
        print("Hybridization probability:", p)
        if p < random.random():
            if domain1.Partner == domain2:
                print("De-hybridizing domain1 and domain2.")
                ## TODO: disassociate the two domains, splitting the complex if required.
                ## TODO: Update rendering/visualization
            else:
                print("Domain 1 did not hybridize to domain2 (p=%s)" % p)
        else:
            if domain1.Partner == domain2:
                print("Domain1 remains hybridized to domain2")
            else:
                # HYBRIDIZE:  (separate this logic out)
                domain1.Partner = domain2
                domain2.Partner = domain1
                # Make complex if required:
                # if domain1.Complex and domain2.Complex
                ## TODO: Update rendering/visualization



    def simulate(self, T, n_steps=1000):
        """
        Simulate at most n_steps number of rounds at temperature T.
        """
        n_done = 0
        while n_done < n_steps:
            self.step(T)


    def anneal(self, T_start, T_finish, delta_T, n_steps_per_T=1000):
        """
        Simulate annealing repeatedly from T_start to T_finish,
        decreasing temperature by delta_T for every round,
        doing at most n_steps number of steps at each temperature.
        """

        for T in range(T_start, T_finish, delta_T):
            self.simulate(T, n_steps_per_T)






if __name__ == '__main__':
    pass
