#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2004-2008 by Sebastian Bassi.
# Copyright 2013 by Markus Piotrowski.
# Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com

# Most of this file is copied from the Biopython distribution,
# and this file is distributed under the official biopython license.
# See the biopython license file for license information.
# Rasmus Sorensen further permits the use and distribution of his specific
# commits under the GNU Public License version 3.

# pylint: disable=C0111,C0103,R0913,R0914



import math


WC = dict(zip("ATGC", "TACG"))
R = 1.987 # cal/mol/K

def compl(seq):
    return "".join(WC[b] for b in seq)
def rcompl(seq):
    return "".join(reversed([WC[b] for b in seq]))

def gc_percent(seq):
    seq = seq.upper()
    return (seq.count("G")+seq.count("C"))/len(seq)


# Values are (deltaH, deltaS) pairs = (kcal/mol, cal/mol/K)

# Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594
DNA_NN3 = {
    'init': (0, 0), 'init_A/T': (2.3, 4.1), 'init_G/C': (0.1, -2.8),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT': (-7.9, -22.2), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -21.3),
    'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0),
    'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4),
    'GG/CC': (-8.0, -19.9)}

# SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
DNA_NN4 = {
    'init': (0.2, -5.7), 'init_A/T': (2.2, 6.9), 'init_G/C': (0, 0),
    'init_oneG/C': (0, 0), 'init_allA/T': (0, 0), 'init_5T/A': (0, 0),
    'sym': (0, -1.4),
    'AA/TT': (-7.6, -21.3), 'AT/TA': (-7.2, -20.4), 'TA/AT': (-7.2, -20.4),
    'CA/GT': (-8.5, -22.7), 'GT/CA': (-8.4, -22.4), 'CT/GA': (-7.8, -21.0),
    'GA/CT': (-8.2, -22.2), 'CG/GC': (-10.6, -27.2), 'GC/CG': (-9.8, -24.4),
    'GG/CC': (-8.0, -19.0)}


# Internal mismatch and inosine table (DNA)
# Allawi & SantaLucia (1997), Biochemistry 36: 10581-10594
# Allawi & SantaLucia (1998), Biochemistry 37: 9435-9444
# Allawi & SantaLucia (1998), Biochemistry 37: 2170-2179
# Allawi & SantaLucia (1998), Nucl Acids Res 26: 2694-2701
# Peyret et al. (1999), Biochemistry 38: 3468-3477
# Watkins & SantaLucia (2005), Nucl Acids Res 33: 6258-6267
DNA_IMM1 = {
    'AG/TT': (1.0, 0.9), 'AT/TG': (-2.5, -8.3), 'CG/GT': (-4.1, -11.7),
    'CT/GG': (-2.8, -8.0), 'GG/CT': (3.3, 10.4), 'GG/TT': (5.8, 16.3),
    'GT/CG': (-4.4, -12.3), 'GT/TG': (4.1, 9.5), 'TG/AT': (-0.1, -1.7),
    'TG/GT': (-1.4, -6.2), 'TT/AG': (-1.3, -5.3), 'AA/TG': (-0.6, -2.3),
    'AG/TA': (-0.7, -2.3), 'CA/GG': (-0.7, -2.3), 'CG/GA': (-4.0, -13.2),
    'GA/CG': (-0.6, -1.0), 'GG/CA': (0.5, 3.2), 'TA/AG': (0.7, 0.7),
    'TG/AA': (3.0, 7.4),
    'AC/TT': (0.7, 0.2), 'AT/TC': (-1.2, -6.2), 'CC/GT': (-0.8, -4.5),
    'CT/GC': (-1.5, -6.1), 'GC/CT': (2.3, 5.4), 'GT/CC': (5.2, 13.5),
    'TC/AT': (1.2, 0.7), 'TT/AC': (1.0, 0.7),
    'AA/TC': (2.3, 4.6), 'AC/TA': (5.3, 14.6), 'CA/GC': (1.9, 3.7),
    'CC/GA': (0.6, -0.6), 'GA/CC': (5.2, 14.2), 'GC/CA': (-0.7, -3.8),
    'TA/AC': (3.4, 8.0), 'TC/AA': (7.6, 20.2),
    'AA/TA': (1.2, 1.7), 'CA/GA': (-0.9, -4.2), 'GA/CA': (-2.9, -9.8),
    'TA/AA': (4.7, 12.9), 'AC/TC': (0.0, -4.4), 'CC/GC': (-1.5, -7.2),
    'GC/CC': (3.6, 8.9), 'TC/AC': (6.1, 16.4), 'AG/TG': (-3.1, -9.5),
    'CG/GG': (-4.9, -15.3), 'GG/CG': (-6.0, -15.8), 'TG/AG': (1.6, 3.6),
    'AT/TT': (-2.7, -10.8), 'CT/GT': (-5.0, -15.8), 'GT/CT': (-2.2, -8.4),
    'TT/AT': (0.2, -1.5),
    'AI/TC': (-8.9, -25.5), 'TI/AC': (-5.9, -17.4), 'AC/TI': (-8.8, -25.4),
    'TC/AI': (-4.9, -13.9), 'CI/GC': (-5.4, -13.7), 'GI/CC': (-6.8, -19.1),
    'CC/GI': (-8.3, -23.8), 'GC/CI': (-5.0, -12.6),
    'AI/TA': (-8.3, -25.0), 'TI/AA': (-3.4, -11.2), 'AA/TI': (-0.7, -2.6),
    'TA/AI': (-1.3, -4.6), 'CI/GA': (2.6, 8.9), 'GI/CA': (-7.8, -21.1),
    'CA/GI': (-7.0, -20.0), 'GA/CI': (-7.6, -20.2),
    'AI/TT': (0.49, -0.7), 'TI/AT': (-6.5, -22.0), 'AT/TI': (-5.6, -18.7),
    'TT/AI': (-0.8, -4.3), 'CI/GT': (-1.0, -2.4), 'GI/CT': (-3.5, -10.6),
    'CT/GI': (0.1, -1.0), 'GT/CI': (-4.3, -12.1),
    'AI/TG': (-4.9, -15.8), 'TI/AG': (-1.9, -8.5), 'AG/TI': (0.1, -1.8),
    'TG/AI': (1.0, 1.0), 'CI/GG': (7.1, 21.3), 'GI/CG': (-1.1, -3.2),
    'CG/GI': (5.8, 16.9), 'GG/CI': (-7.6, -22.0),
    'AI/TI': (-3.3, -11.9), 'TI/AI': (0.1, -2.3), 'CI/GI': (1.3, 3.0),
    'GI/CI': (-0.5, -1.3)}


# Terminal mismatch table (DNA)
# SantaLucia & Peyret (2001) Patent Application WO 01/94611
DNA_TMM1 = {
    'AA/TA': (-3.1, -7.8), 'TA/AA': (-2.5, -6.3), 'CA/GA': (-4.3, -10.7),
    'GA/CA': (-8.0, -22.5),
    'AC/TC': (-0.1, 0.5), 'TC/AC': (-0.7, -1.3), ' CC/GC': (-2.1, -5.1),
    'GC/CC': (-3.9, -10.6),
    'AG/TG': (-1.1, -2.1), 'TG/AG': (-1.1, -2.7), 'CG/GG': (-3.8, -9.5),
    'GG/CG': (-0.7, -19.2),
    'AT/TT': (-2.4, -6.5), 'TT/AT': (-3.2, -8.9), 'CT/GT': (-6.1, -16.9),
    'GT/CT': (-7.4, -21.2),
    'AA/TC': (-1.6, -4.0), 'AC/TA': (-1.8, -3.8), 'CA/GC': (-2.6, -5.9),
    'CC/GA': (-2.7, -6.0), 'GA/CC': (-5.0, -13.8), 'GC/CA': (-3.2, -7.1),
    'TA/AC': (-2.3, -5.9), 'TC/AA': (-2.7, -7.0),
    'AC/TT': (-0.9, -1.7), 'AT/TC': (-2.3, -6.3), 'CC/GT': (-3.2, -8.0),
    'CT/GC': (-3.9, -10.6), 'GC/CT': (-4.9, -13.5), 'GT/CC': (-3.0, -7.8),
    'TC/AT': (-2.5, -6.3), 'TT/AC': (-0.7, -1.2),
    'AA/TG': (-1.9, -4.4), 'AG/TA': (-2.5, -5.9), 'CA/GG': (-3.9, -9.6),
    'CG/GA': (-6.0, -15.5), 'GA/CG': (-4.3, -11.1), ' GG/CA': (-4.6, -11.4),
    'TA/AG': (-2.0, -4.7), 'TG/AA': (-2.4, -5.8),
    'AG/TT': (-3.2, -8.7), 'AT/TG': (-3.5, -9.4), 'CG/GT': (-3.8, -9.0),
    'CT/GG': (-6.6, -18.7), 'GG/CT': (-5.7, -15.9), 'GT/CG': (-5.9, -16.1),
    'TG/AT': (-3.9, -10.5), 'TT/AG': (-3.6, -9.8)}


# Dangling ends table (DNA)
# Bommarito et al. (2000), Nucl Acids Res 28: 1929-1934
DNA_DE1 = {
    'AA/.T': (0.2, 2.3), 'AC/.G': (-6.3, -17.1), 'AG/.C': (-3.7, -10.0),
    'AT/.A': (-2.9, -7.6), 'CA/.T': (0.6, 3.3), 'CC/.G': (-4.4, -12.6),
    'CG/.C': (-4.0, -11.9), 'CT/.A': (-4.1, -13.0), 'GA/.T': (-1.1, -1.6),
    'GC/.G': (-5.1, -14.0), 'GG/.C': (-3.9, -10.9), 'GT/.A': (-4.2, -15.0),
    'TA/.T': (-6.9, -20.0), 'TC/.G': (-4.0, -10.9), 'TG/.C': (-4.9, -13.8),
    'TT/.A': (-0.2, -0.5),
    '.A/AT': (-0.7, -0.8), '.C/AG': (-2.1, -3.9), '.G/AC': (-5.9, -16.5),
    '.T/AA': (-0.5, -1.1), '.A/CT': (4.4, 14.9), '.C/CG': (-0.2, -0.1),
    '.G/CC': (-2.6, -7.4), '.T/CA': (4.7, 14.2), '.A/GT': (-1.6, -3.6),
    '.C/GG': (-3.9, -11.2), '.G/GC': (-3.2, -10.4), '.T/GA': (-4.1, -13.1),
    '.A/TT': (2.9, 10.4), '.C/TG': (-4.4, -13.1), '.G/TC': (-5.2, -15.0),
    '.T/TA': (-3.8, -12.6)}

# Hmm... dangling ends observation: ΔH° is almost always equal to T*ΔS° for T = 325 K.
#  scipy.stats.linregres(dHs, [325*dS/1000 for dS in dSs]) - slope = 1.0016, intersect = 0.277.



energy_tables_in_units_of_R = {}
for tbl_name in "DNA_NN3, DNA_NN4, DNA_IMM1, DNA_TMM1, DNA_DE1".split(", "):
    tbl = globals()[tbl_name]
    energy_tables_in_units_of_R[tbl_name] = {
        # R is in cal/mol/K
        k: (H*1000/R, S/R) for k, (H, S) in tbl.items()
    }

# To use energies in unit of R:
# K = exp(-dG/T) = exp(dS-dH/T)  # Very simple to interpret and use, right?


def canonical_partitions(energies, T, unit='cal/mol'):
    """
    Return list [exp(-E_s / kT) for s in energies microstates]
    """
    k = {'J/K': 1.380e-23, 'cal/K': 0.329e-23}
    R = {'J/mol/K': 8.314, 'cal/mol/K': 1.987}
    if 'mol' in unit:
        beta = 1/(R['cal/mol/K' if 'cal' in unit else 'J/mol/K']*T)
    else:
        beta = 1/(k['cal/K' if 'cal' in unit else 'J/K']*T)
    if unit[0] == 'k':
        beta = beta / 1000
    try:
        return [math.exp(-beta*E_s) for E_s in energies]
    except TypeError:
        return math.exp(-beta*energies)

def probabilities(energies, T, unit='cal/mol'):
    """
    Return a list of probabilities (probability distribution) given a list of energies,
    classical discrete canonical partition function, assuming fixed temperature, volume, and number.
    """
    partitions = canonical_partitions(energies, T, unit)
    Z = sum(partition for partition in partitions)      # Partition function
    try:
        # p_s = exp(-E_s / kT) / Z # Probability of system being in microstate s
        return [partition/Z for partition in partitions]
    except TypeError:
        # Assume we only got one energy and we are calculating its probability in a binary state system,
        # i.e. a ΔE, ΔG or similar: Set E_other_state = 0, so part_other_state = 1:
        # prob_s = partition_s / (partiton_s + 1)
        return partitions / (partitions + 1)

def binary_state_probability(energy, T, unit='cal/mol', Q=1):
    """
    Consideration: Is Stirling's approximation applicable here?
        ln(N!) ≈ N ln(N) - N        # When N is very large, N! = factorial(N) = N*(N-1)*(N-2)*...*(N-(N-1))
    Actually, N ln(N) - N + ln(N)/2 + 1/ln(3) is a slight improvement of the approximation.

    Q is the reaction fraction in ΔG = ΔG° + RT ln(Q).
    This can be used to get the state probability when the transition does not have activity 1.
    E.g. for duplex formation:
        Q = [duplex] / ([strandA] [strandB])        # Q will be > 1 for typical concentrations.
    If you use a virtual concentration of say 1e-3 M for all parties, then:
        Q = [duplex] / ([strandA] [strandB]) = 1 / 1e-3 = 1e3.
    At T = 330 K, this will give a correction factor of
        RT ln(Q) = 1.987 cal/mol/K * 330 K * ln(1e3) = 4529 cal/mol.

    TODO: Notice - Pierce et al often use Q to denote partition functions rather than
        reaction fractions. Consider using "K" instead of Q.
    """
    if Q is None:
        Q = 1
    # TODO: Optimize this to make it more efficient.
    k = {'J/K': 1.380e-23, 'cal/K': 0.329e-23}  # Boltzman. For single particle states.
    R = {'J/mol/K': 8.314, 'cal/mol/K': 1.987}  # Gas constant. For molar quantities.
    if 'mol' in unit:
        beta = 1/(R['cal/mol/K' if 'cal' in unit else 'J/mol/K']*T)
    else:
        beta = 1/(k['cal/K' if 'cal' in unit else 'J/K']*T)
    if unit[0] == 'k':
        beta = beta / 1000
    if Q is not None and Q != 1:
        energy = energy + math.log(Q)/beta # ΔG = ΔG° + RT ln(Q), beta = 1/RT
    p_i = 1 / (1 + math.exp(beta*energy))   # Same as e^(ΔE/kT)...
    # The minus in part = e^(-ΔE/kT) is lost during algebraic transformation:
    # p_i = exp(-ΔE/kT) / (exp(-ΔE/kT) + 1) = 1/(1+exp(ΔE/kT))
    return p_i

def binary_state_probability_cal_per_mol(deltaE, T, Q=1):
    """
    Optimized calculation of binary state probability from deltaE, T and Q.
    Gives the probability of being in state f, considering the a system of two
    states, i and f and the state change: i -> f where deltaE is E_f - E_i
    If the initial and final states have multiple components, and the given
    deltaE is for another system, you can include a proper Q, which is the
    molecualar ratio: Q = [f_1]*[f_2]*(...) / [i_1]*[i_2]*(...)
    For a hybridization reaction, the reaction quotient is
        Q = [duplex] / ([strandA]*[strandB])
    """
    #if Q and Q != 1:
    #    deltaE = deltaE + math.log(Q)/beta
    # e^(ΔE/RT) = e^((ΔE° + RT ln(Q))/RT) =
    #           = e^(ΔE°/RT + ln(Q)) = e^(ΔE°/RT)*Q
    # Equivalently,
    #   e^(-ΔE/RT) = e^(-ΔE°/RT) * 1/Q
    R = 1.987 # cal/mol/K
    p_i = 1 / (1 + math.exp(deltaE/(R*T))*Q)    # Verified, is OK.
    return p_i


def hybridization_dH_dS(seq, check=True, c_seq=None, shift=0, nn_table=DNA_NN4,
                        tmm_table=DNA_TMM1, imm_table=DNA_IMM1, de_table=DNA_DE1,
                        selfcomp=False, Na=50, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5,
                        verbose=False):
    """
    Calculate standard enthalpy and entropy of hybridization, i.e. the reaction

            strandA + strandB -> duplex(strandA, strandB)

    *without* accounting for intra-strand base-pairs that might stabilize the
    ensemble state of the single strands in their non-hybridized state.
    (this might be very important!)

    Return values are in units of kcal/mol for dH and cal/mol/K for dS.

    Arguments:
     - seq: The primer/probe sequence as string or Biopython sequence object. For
       RNA/DNA hybridizations seq must be the RNA sequence.
     - c_seq: Complementary sequence. The sequence of the template/target in
       3'->5' direction. c_seq is necessary for mismatch correction and
       dangling-ends correction. Both corrections will automatically be
       applied if mismatches or dangling ends are present. Default=None.
     - shift: Shift of the primer/probe sequence on the template/target sequence,
       e.g.::
                           shift=0       shift=1        shift= -1
        Primer (seq):      5' ATGC...    5'  ATGC...    5' ATGC...
        Template (c_seq):  3' TACG...    3' CTACG...    3'  ACG...

       The shift parameter is necessary to align seq and c_seq if they have
       different lengths or if they should have dangling ends. Default=0


    Note: The energies calculated are indeed for *hybridization*:
        We get a negative enthalpy and a negative entropy (for complementary strands).

    To calculate ΔG°:
        ΔG° = ΔH° - T*ΔS°
    To calculate Gibbs free energy at non-equilibrium:
        ΔG = ΔG° + RT ln(Q)     # Where Q = [duplex]/([strandA] [strandB])
    Or, to calculate equilibrium constant (at equilibrium: ΔG = 0):
        ΔG° = ΔH° - T*ΔS° = -RT ln(K)
        K = exp(-ΔG°/RT) = exp(ΔS°/R - ΔH°/RT)

    Note that K = [strandA]_tot/2 at T = Tm, if [strandA] = [strandB].
    It might be nice to use a "extend of hybridization" fraction that goes from
    0.0 - 1.0 (0% to 100% hybridized):
        X_A = 1 - [StrandA] / [strandA]_tot = 1 - [StrandA] / (K*2) # Conversion of strandA (consumed)
        Y_A = [StrandA] / [strandA]_tot = [StrandA] / (K*2)         # Yield of duplex formation

    The Tm_NN function below inverts the K_melt (watch out for the minus in the equations)
        ln(K_hyb) = -ln(K_melt)
        Tm = ΔH°/(ΔS° + R ln(K))    # K when 50% of strand 1 is melted.
    Tm (in degC) is calculated as:
        Tm = (1000 * deltaH) / (deltaS + (R * (math.log(k)))) - 273.15

    How to get partition function?
    Two scenarios:
    a) If you have to account for concentration dependence:
            ??

    b) If you do not have to account for concentration, i.e. if you have already done so by other statistical means
        and you can assume an activity of 1.
        From statistical mechanics:
            Z   = sum(exp(-E_s / kT) for s in microstates)   # Partition function for classical canonical ensemble
            p_s = exp(-E_s / kT) / Z                         # Probability of system being in microstate s
        We usually only consider two states, hybridized and non-hybridized
            Z   = exp(-E_hyb / kT) + exp(-E_melt / kT)
        We can set the energy offset arbitrarily, so set E_melt = 0:
            Z   = exp(-E_hyb / kT) + exp(0 / kT) = exp(-E_hyb / kT) + 1
            p_h = exp(-E_hyb / kT) / (exp(-E_hyb / kT) + 1)
                = 1 / ((A+1)/A) = 1 / (1 + 1/A)     # Using A = exp(-E_hyb / kT) as abbreviation
                = 1 / (1+exp(ΔE_hyb/kT))            # ΔE_hyb = E_hyb if E_melt = 0
        Examples:
            ΔE_hyb = 0 kT: p_hyb = 1 /  (1+exp( 0kT/kT)) = 1 / (1+1) = 0.5          # Both states equally probable
            ΔE_hyb = 1 kT: p_hyb = 1 /  (1+exp( 1kT/kT)) = 1 / (1+2.72) = 0.269     # Melting is slightly favored
            ΔE_hyb = 2 kT: p_hyb = 1 /  (1+exp( 2kT/kT)) = 1 / (1+7.39) = 0.119     # Melting is significantly favored
            ΔE_hyb = 3 kT: p_hyb = 1 /  (1+exp( 3kT/kT)) = 1 / (1+20.1) = 0.047     # Melting is strongly favored
            ΔE_hyb = -1 kT: p_hyb = 1 / (1+exp(-1kT/kT)) = 1 / (1+0.37) = 0.731     # Hybridization is slightly favored
            ΔE_hyb = -2 kT: p_hyb = 1 / (1+exp(-2kT/kT)) = 1 / (1+0.13) = 0.881     # Hybridization is signifly favored
            ΔE_hyb = -3 kT: p_hyb = 1 / (1+exp(-3kT/kT)) = 1 / (1+0.05) = 0.953     # Hybridization is strongly favored

        Note: If ΔE_hyb is in values of kcal/mol or kJ/mol, then use the Gas constant R = k N_A instead of k:
            k = 1.380e−23 J/K = 0.329e-23 cal/K = 8.617e-5 eV/K
            R = 8.314 J/mol/K = 1.987 cal/mol/K

    Note:
        If you have calculated ΔG° = ΔH° - ΔS°*T  at a certain Tm,
        where the p_hyb is 50% with concentrations of say 1 uM. Then at a concentration/activity of 1 M, that would
        probability is suddenly increased to something like 99.9998 %. In other words, at micro-molar concentrations,
        k_on is low because of the very large volume that the strand needs to search before finding a partner.
        If our simulation is at T=Tm, we need k_on to equal k_off. If we assume an activity of 1 M when we decide
        whether an oligo forms a duplex or not, then the "oligo selection" probability has to be very low,
        compared to the selection chance when the oligos are hybridized.
        Alternatively, we could say that "selection" means bringing the two strands to within some pre-defined
        distance (=concentration), and then using
            ΔG = ΔG° + RT ln(Q) = ΔG° + RT ln([duplex] / ([strandA][strandB]))

        When Q is < 1: RT ln(Q) is negative... Making the reaction even more favorable??
        Yes, Q < 1  <=>  1/[strandA] < 1  <=> [strandA] > 1.

        If we are only considering those two strands, then :
            ΔG = ΔG° + RT ln(1/[strandA])

        It might be nice to calculate a ΔG vs ln(Q) plot ()
    """
    comb_table = nn_table.copy()
    comb_table.update(tmm_table)
    comb_table.update(imm_table)
    comb_table.update(de_table)
    # Also add the "reversed" neighbourgh pair:
    comb_table.update({k[::-1]: v for k, v in comb_table.items()})
    # Neither GA/TA nor AT/AG is present?

    seq = str(seq)
    if c_seq is None:
        # Wait... are we taking the rcompl or compl? - compl. c_seq is 3'-to-5'.
        c_seq = compl(seq)
    if check:
        baseset = ('A', 'C', 'G', 'T', 'I')
        seq = ''.join([base for base in seq if base in baseset])
    tmpseq = seq
    tmp_cseq = c_seq
    deltaH = 0
    deltaS = 0
    dH = 0  # Names for indexes
    dS = 1  # 0 and 1

    # Dangling ends?
    if shift or len(seq) != len(c_seq):
        # Align both sequences using the shift parameter
        if shift > 0:
            tmpseq = '.' * shift + seq
        if shift < 0:
            tmp_cseq = '.' * abs(shift) + c_seq
        if len(tmp_cseq) > len(tmpseq):
            tmpseq += (len(tmp_cseq) - len(tmpseq)) * '.'
        if len(tmp_cseq) < len(tmpseq):
            tmp_cseq += (len(tmpseq) - len(tmp_cseq)) * '.'
        if verbose:
            print("Aligned sequences:\n  %s\n  %s" % (tmpseq, tmp_cseq))
        # Remove 'over-dangling' ends
        while tmpseq.startswith('..') or tmp_cseq.startswith('..'):
            tmpseq = tmpseq[1:]
            tmp_cseq = tmp_cseq[1:]
        while tmpseq.endswith('..') or tmp_cseq.endswith('..'):
            tmpseq = tmpseq[:-1]
            tmp_cseq = tmp_cseq[:-1]
        # Now for the dangling ends
        if tmpseq.startswith('.') or tmp_cseq.startswith('.'):
            left_de = tmpseq[:2] + '/' + tmp_cseq[:2]
            deltaH += de_table[left_de][dH]
            deltaS += de_table[left_de][dS]
            if verbose:
                print("Adding left dangling end contribution for %s: dH+=%s, dS+=%s" %
                      (left_de, de_table[left_de][dH], de_table[left_de][dS]))
            tmpseq = tmpseq[1:]
            tmp_cseq = tmp_cseq[1:]
        if tmpseq.endswith('.') or tmp_cseq.endswith('.'):
            right_de = tmp_cseq[-2:][::-1] + '/' + tmpseq[-2:][::-1]
            deltaH += de_table[right_de][dH]
            deltaS += de_table[right_de][dS]
            if verbose:
                print("Adding right dangling end contribution for %s: dH+=%s, dS+=%s" %
                      (right_de, de_table[right_de][dH], de_table[right_de][dS]))
            tmpseq = tmpseq[:-1]
            tmp_cseq = tmp_cseq[:-1]

    # Now for terminal mismatches
    left_tmm = tmp_cseq[:2][::-1] + '/' + tmpseq[:2][::-1]
    if left_tmm in tmm_table:
        deltaH += tmm_table[left_tmm][dH]
        deltaS += tmm_table[left_tmm][dS]
        tmpseq = tmpseq[1:]
        tmp_cseq = tmp_cseq[1:]
        if verbose:
            print("Adding left terminal mismatch contribution for %s: dH+=%s, dS+=%s" %
                  (left_tmm, tmm_table[left_tmm][dH], tmm_table[left_tmm][dS]))
    right_tmm = tmpseq[-2:] + '/' + tmp_cseq[-2:]
    if right_tmm in tmm_table:
        deltaH += tmm_table[right_tmm][dH]
        deltaS += tmm_table[right_tmm][dS]
        tmpseq = tmpseq[:-1]
        tmp_cseq = tmp_cseq[:-1]
        if verbose:
            print("Adding right terminal mismatch contribution for %s: dH+=%s, dS+=%s" %
                  (right_tmm, tmm_table[right_tmm][dH], tmm_table[right_tmm][dS]))

    # Now everything 'unusual' at the ends is handled and removed and we can
    # look at the initiation.
    # One or several of the following initiation types may apply:

    # Type: General initiation value
    deltaH += nn_table['init'][dH]
    deltaS += nn_table['init'][dS]

    # Type: Duplex with no (allA/T) or at least one (oneG/C) GC pair
    if seq.count("C") + seq.count("C") == 0:
        deltaH += nn_table['init_allA/T'][dH]
        deltaS += nn_table['init_allA/T'][dS]
    else:
        deltaH += nn_table['init_oneG/C'][dH]
        deltaS += nn_table['init_oneG/C'][dS]

    # Type: Penalty if 5' end is T (for both strands)
    if seq.startswith('T'):
        deltaH += nn_table['init_5T/A'][dH]
        deltaS += nn_table['init_5T/A'][dS]
    if seq.endswith('A'):
        deltaH += nn_table['init_5T/A'][dH]
        deltaS += nn_table['init_5T/A'][dS]

    # Type: Different values for G/C or A/T terminal basepairs
    ends = seq[0] + seq[-1]
    AT = ends.count('A') + ends.count('T')
    GC = ends.count('G') + ends.count('C')
    deltaH += nn_table['init_A/T'][dH] * AT
    deltaS += nn_table['init_A/T'][dS] * AT
    deltaH += nn_table['init_G/C'][dH] * GC
    deltaS += nn_table['init_G/C'][dS] * GC

    # Finally, the 'zipping'
    for basenumber in range(len(tmpseq) - 1):
        neighbors = tmpseq[basenumber:basenumber + 2] + '/' + \
            tmp_cseq[basenumber:basenumber + 2]
        #try:
        #    deltaH += comb_table[neighbors][dH]
        #    deltaS += comb_table[neighbors][dS]
        #except KeyError as e:
        #    print("Key not found: ", e)
        #    print("comb_table keys: ", ", ".join(comb_table.keys()))
        #    raise e
        if neighbors in imm_table:
            deltaH += imm_table[neighbors][dH]
            deltaS += imm_table[neighbors][dS]
        elif neighbors[::-1] in imm_table:
            deltaH += imm_table[neighbors[::-1]][dH]
            deltaS += imm_table[neighbors[::-1]][dS]
        elif neighbors in nn_table:
            deltaH += nn_table[neighbors][dH]
            deltaS += nn_table[neighbors][dS]
        elif neighbors[::-1] in nn_table:
            deltaH += nn_table[neighbors[::-1]][dH]
            deltaS += nn_table[neighbors[::-1]][dS]
        else:
            raise KeyError("Could not find key %s when zipping" % neighbors)

    if selfcomp:
        # Symmetry addition:
        deltaH += nn_table['sym'][dH]
        deltaS += nn_table['sym'][dS]
    if saltcorr and saltcorr == 5:
        # Only salt correction method #5 considers energy, the others applies to Tm.
        corr = salt_correction(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                               method=saltcorr, seq=seq)
        deltaS += corr
    return deltaH, deltaS



def Tm_NN(seq, check=True, c_seq=None, shift=0, nn_table=DNA_NN3,
          tmm_table=DNA_TMM1, imm_table=DNA_IMM1, de_table=DNA_DE1,
          dnac1=25, dnac2=25, selfcomp=False, Na=50, K=0, Tris=0, Mg=0,
          dNTPs=0, saltcorr=5, molar_unit='nM'):
    """
    Return the Tm using nearest neighbor thermodynamics.

    Arguments:
     - seq: The primer/probe sequence as string or Biopython sequence object. For
       RNA/DNA hybridizations seq must be the RNA sequence.
     - c_seq: Complementary sequence. The sequence of the template/target in
       3'->5' direction. c_seq is necessary for mismatch correction and
       dangling-ends correction. Both corrections will automatically be
       applied if mismatches or dangling ends are present. Default=None.
     - shift: Shift of the primer/probe sequence on the template/target sequence,
       e.g.::
                           shift=0       shift=1        shift= -1
        Primer (seq):      5' ATGC...    5'  ATGC...    5' ATGC...
        Template (c_seq):  3' TACG...    3' CTACG...    3'  ACG...

       The shift parameter is necessary to align seq and c_seq if they have
       different lengths or if they should have dangling ends. Default=0
     - table: Thermodynamic NN values, eight tables are implemented:
       For DNA/DNA hybridizations:

        - DNA_NN1: values from Breslauer et al. (1986)
        - DNA_NN2: values from Sugimoto et al. (1996)
        - DNA_NN3: values from Allawi & SantaLucia (1997) (default)
        - DNA_NN4: values from SantaLucia & Hicks (2004)

       For RNA/RNA hybridizations:

        - RNA_NN1: values from Freier et al. (1986)
        - RNA_NN2: values from Xia et al. (1998)
        - RNA_NN3: valuse from Chen et al. (2012)

       For RNA/DNA hybridizations:

        - R_DNA_NN1: values from Sugimoto et al. (1995)

       Use the module's maketable method to make a new table or to update one
       one of the implemented tables.
     - tmm_table: Thermodynamic values for terminal mismatches. Default: DNA_TMM1
       (SantaLucia & Peyret, 2001)
     - imm_table: Thermodynamic values for internal mismatches, may include
       insosine mismatches. Default: DNA_IMM1 (Allawi & SantaLucia, 1997-1998;
       Peyret et al., 1999; Watkins & SantaLucia, 2005)
     - de_table: Thermodynamic values for dangling ends:

        - DNA_DE1: for DNA. Values from Bommarito et al. (2000). Default
        - RNA_DE1: for RNA. Values from Turner & Mathews (2010)

     - dnac1: Concentration of the higher concentrated strand [nM]. Typically this
       will be the primer (for PCR) or the probe. Default=25.
     - dnac2: Concentration of the lower concentrated strand [nM]. In PCR this is
       the template strand which concentration is typically very low and may
       be ignored (dnac2=0). In oligo/oligo hybridization experiments, dnac1
       equals dnac1. Default=25.
       MELTING and Primer3Plus use k = [Oligo(Total)]/4 by default. To mimic
       this behaviour, you have to divide [Oligo(Total)] by 2 and assign this
       concentration to dnac1 and dnac2. E.g., Total oligo concentration of
       50 nM in Primer3Plus means dnac1=25, dnac2=25.
     - selfcomp: Is the sequence self-complementary? Default=False. If 'True' the
       primer is thought binding to itself, thus dnac2 is not considered.
     - Na, K, Tris, Mg, dNTPs: See method 'Tm_GC' for details. Defaults: Na=50,
       K=0, Tris=0, Mg=0, dNTPs=0.
     - saltcorr: See method 'Tm_GC'. Default=5. 0 means no salt correction.

    """

    deltaH, deltaS = hybridization_dH_dS(seq, check, c_seq, shift,
                                         nn_table, tmm_table, imm_table, de_table,
                                         selfcomp, Na, K, Tris, Mg, dNTPs, saltcorr)

    units = {'M': 1, 'm': 1e-3, 'u': 1e-6, 'n': 1e-9, 'p': 1e-12, 'f': 1e-15, 'a': 1e-18}
    molar_unit = units[molar_unit[0]]

    if selfcomp:
        k = dnac1 * molar_unit    # equilibrium constant for bi-molecular homo-dimerization reaction
        # We are looking at *melting* temperature, and all thermodynamic data are for the melting reaction:
        #   duplex -> strandA + strandB
        # K = [strandA]*[strandA] / [hybridized] = [strandA]    # because [strandA] = [hybridized] at Tm.
        # K = exp(-dG/RT) = 1
    else:
        # "x equals 4 for nonself-complementary duplexes and equals 1 for self-complementary duplexes."
        # For bi-molecular hybridization reaction:
        # K = [strandA]*[strandB]/[duplex]
        # StrandB (in excess): [strandB] = [strandB]_init - [duplex]
        # At T=Tm (half of strandA is in duplex), [duplex] = [strandA] = [strandA]_init/2
        # and let dnac2 = [strandA]_init/2, dnac1 = [strandB]_init:
        # At T=Tm: K = [strandB]*[strandA]/[duplex] = [strandB] # since [duplex] = [strandA]
        #            = [strandB]_init - [duplex]                # since [strandB] = [strandB]_init - [duplex]
        #            = [strandB]_init - [strandA]_init/2        # since [duplex] = [strandA]_init/2
        #            = dnac1 - dnac2 / 2
        k = (dnac1 - (dnac2 / 2.0)) * molar_unit    # Note: This is not for hybridization but for MELTING!
        # I don't get it... If dnac2 > 2 dnac1 then this is negative...?
        # dnac2 MUST be the lower of the two concentrations.
        # if dnac2 = dnac1, then dnac1 - (dnac2 / 2) = (2*dnac1/2 - (dnac1 / 2)) = dnac1/2

    # K = exp(-dG/RT) = 1
    #
    # See e.g. SantaLucia & Hicks (2004) eq (3).
    R = 1.987  # universal gas constant in Cal/degrees C*Mol
    # deltaH is in kcal/mol, deltaS is in cal/mol/K
    Tm = (1000 * deltaH) / (deltaS + (R * (math.log(k)))) - 273.15
    # ΔG° = ΔH° - T*ΔS°
    # -ΔG°/RT = (ΔS°-ΔH°/T)/R
    # Non-equilibrium: ΔG = ΔG° + RT ln(Q)  ,  Q = reaction quotient, also sometimes denoted 'Y'.
    # At equilibrium, ΔG = 0, so ΔG° = -RT ln(Q) and K = Q(equilibrium) = exp(-ΔG°/RT)
    # At equilibrium:
    #   ΔG = ΔH - T*ΔS = 0
    #   T = ΔH/ΔS   # but we usually do not know ΔH or ΔS , we know ΔH° and ΔS° (at standard conditions).
    #   ΔG° = ΔH° - T*ΔS° = -RT ln(K) <=> ΔH° = (T*ΔS°-RT ln(K)) = T * (ΔS°-R ln(K))
    #       <=> T = ΔH° / (ΔS° - R ln(K))   # Note the change in -R ln(K)
    # The energies we are calculating are base-pairing energies: strandA + strandB -> strandA:strandB
    # However, the equilibrium constant is for the melting reaction.
    # Letting K_melt = 1/K_hyb changes the minus to a plus in the equation above:
    #       ln(K_hyb) = ln(1/K_melt) = -ln(K_melt)
    # We want K = Q = [products]/[reagents] = 1
    # Note: This is based on a two-state approximation, where the melted strands
    # only has a single un-basepaired structure and does not form any secondary structures.
    # if saltcorr == 5, then e have already done salt correction when calculating deltaS above:
    if saltcorr and saltcorr != 5:
        corr = salt_correction(Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
                               method=saltcorr, seq=seq)
    if saltcorr in (1, 2, 3, 4):
        Tm += corr
    if saltcorr in (6, 7):
        # Tm = 1/(1/Tm + corr)
        Tm = (1 / (1 / (Tm + 273.15) + corr) - 273.15)

    return Tm



def salt_correction(Na=0, K=0, Tris=0, Mg=0, dNTPs=0, method=1, seq=None):
    """Calculate a term to correct Tm for salt ions.

    Depending on the Tm calculation, the term will correct Tm or entropy. To
    calculate corrected Tm values, different operations need to be applied:

     - methods 1-4: Tm(new) = Tm(old) + corr
     - method 5: deltaH(new) = deltaH(old) + corr
     - methods 6+7: Tm(new) = 1/(1/Tm(old) + corr)

    Parameters:
     - Na, K, Tris, Mg, dNTPS: Millimolar concentration of respective ion. To have
       a simple 'salt correction', just pass Na. If any of K, Tris, Mg and
       dNTPS is non-zero, a 'sodium-equivalent' concentration is calculated
       according to von Ahsen et al. (2001, Clin Chem 47: 1956-1961):
       [Na_eq] = [Na+] + [K+] + [Tris]/2 + 120*([Mg2+] - [dNTPs])^0.5
       If [dNTPs] >= [Mg2+]: [Na_eq] = [Na+] + [K+] + [Tris]/2
     - method: Which method to be applied. Methods 1-4 correct Tm, method 5
       corrects deltaS, methods 6 and 7 correct 1/Tm. The methods are:

       1. 16.6 x log[Na+]
          (Schildkraut & Lifson (1965), Biopolymers 3: 195-208)
       2. 16.6 x log([Na+]/(1.0 + 0.7*[Na+]))
          (Wetmur (1991), Crit Rev Biochem Mol Biol 126: 227-259)
       3. 12.5 x log(Na+]
          (SantaLucia et al. (1996), Biochemistry 35: 3555-3562
       4. 11.7 x log[Na+]
          (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465
       5. Correction for deltaS: 0.368 x (N-1) x ln[Na+]
          (SantaLucia (1998), Proc Natl Acad Sci USA 95: 1460-1465)
       6. (4.29(%GC)-3.95)x1e-5 x ln[Na+] + 9.40e-6 x ln[Na+]^2
          (Owczarzy et al. (2004), Biochemistry 43: 3537-3554)
       7. Complex formula with decision tree and 7 empirical constants.
          Mg2+ is corrected for dNTPs binding (if present)
          (Owczarzy et al. (2008), Biochemistry 47: 5336-5353)

    Examples:

        >>> from Bio.SeqUtils import MeltingTemp as mt
        >>> print('%0.2f' % mt.salt_correction(Na=50, method=1))
        -21.60
        >>> print('%0.2f' % mt.salt_correction(Na=50, method=2))
        -21.85
        >>> print('%0.2f' % mt.salt_correction(Na=100, Tris=20, method=2))
        -16.45
        >>> print('%0.2f' % mt.salt_correction(Na=100, Tris=20, Mg=1.5,
        ...                                    method=2))
        -10.99

    """
    if method in (5, 6, 7) and not seq:
        raise ValueError('sequence is missing (is needed to calculate ' +
                         'GC content or sequence length).')
    if seq:
        seq = str(seq)
    corr = 0
    if not method:
        return corr
    Mon = Na + K + Tris / 2.0  # Note: all these values are millimolar
    mg = Mg * 1e-3             # Lowercase ions (mg, mon, dntps) are molar
    # Na equivalent according to von Ahsen et al. (2001):
    if sum((K, Mg, Tris, dNTPs)) > 0 and not method == 7:
        if dNTPs < Mg:
            # dNTPs bind Mg2+ strongly. If [dNTPs] is larger or equal than
            # [Mg2+], free Mg2+ is considered not to be relevant.
            Mon += 120 * math.sqrt(Mg - dNTPs)
    mon = Mon * 1e-3
    # Note: math.log = ln(), math.log10 = log()
    if method in range(1, 7) and not mon:
        raise ValueError('Total ion concentration of zero is not allowed in ' +
                         'this method.')
    if method == 1:
        corr = 16.6 * math.log10(mon)
    if method == 2:
        corr = 16.6 * math.log10((mon) / (1.0 + 0.7 * (mon)))
    if method == 3:
        corr = 12.5 * math.log10(mon)
    if method == 4:
        corr = 11.7 * math.log10(mon)
    if method == 5:
        corr = 0.368 * (len(seq) - 1) * math.log(mon)
    if method == 6:
        corr = (4.29 * gc_percent(seq) / 100 - 3.95) * 1e-5 * math.log(mon) + \
            9.40e-6 * math.log(mon) ** 2
    if method == 7:
        a, b, c, d = 3.92, -0.911, 6.26, 1.42
        e, f, g = -48.2, 52.5, 8.31
        if dNTPs > 0:
            dntps = dNTPs * 1e-3
            ka = 3e4  # Dissociation constant for Mg:dNTP
            # Free Mg2+ calculation:
            mg = (-(ka * dntps - ka * mg + 1.0) +
                  math.sqrt((ka * dntps - ka * mg + 1.0) ** 2 + 4.0 * ka * mg)) / (2.0 * ka)
        if Mon > 0:
            R = math.sqrt(mg) / mon
            if R < 0.22:
                corr = (4.29 * gc_percent(seq) / 100 - 3.95) * \
                    1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2
                return corr
            elif R < 6.0:
                a = 3.92 * (0.843 - 0.352 * math.sqrt(mon) * math.log(mon))
                d = 1.42 * (1.279 - 4.03e-3 * math.log(mon) -
                            8.03e-3 * math.log(mon) ** 2)
                g = 8.31 * (0.486 - 0.258 * math.log(mon) +
                            5.25e-3 * math.log(mon) ** 3)
        corr = (a + b * math.log(mg) + (gc_percent(seq) / 100) *
                (c + d * math.log(mg)) + (1 / (2.0 * (len(seq) - 1))) *
                (e + f * math.log(mg) + g * math.log(mg) ** 2)) * 1e-5
    if method > 7:
        raise ValueError('Allowed values for parameter \'method\' are 1-7.')
    return corr


def chem_correction(Tm, DMSO=0, fmd=0, DMSOfactor=0.75, fmdfactor=0.65,
                    fmdmethod=1, GC=None):
    """Correct a given Tm for DMSO and formamide.

    Please note that these corrections are +/- rough approximations.

    Arguments:
     - Tm: Melting temperature.
     - DMSO: Percent DMSO.
     - fmd: Formamide concentration in %(fmdmethod=1) or in molar (fmdmethod=2).
     - DMSOfactor: How much should Tm decreases per percent DMSO. Default=0.65
       (von Ahsen et al. 2001). Other published values are 0.5, 0.6 and 0.675.
     - fmdfactor: How much should Tm decrease per percent formamide. Default=0.65.
       Several papers report factors between 0.6 and 0.72.
     - fmdmethod:

         1. Tm = Tm - factor(%formamide) (Default)
         2. Tm = Tm + (0.453(f(GC)) - 2.88) x [formamide]

       Here f(GC) is fraction of GC.
       Note (again) that in fmdmethod=1 formamide concentration is given in %,
       while in fmdmethod=2 it is given in molar.
     - GC: GC content in percent.

    Examples:

        >>> from Bio.SeqUtils import MeltingTemp as mt
        >>> mt.chem_correction(70)
        70
        >>> print('%0.2f' % mt.chem_correction(70, DMSO=3))
        67.75
        >>> print('%0.2f' % mt.chem_correction(70, fmd=5))
        66.75
        >>> print('%0.2f' % mt.chem_correction(70, fmdmethod=2, fmd=1.25,
        ...                                    GC=50))
        66.68

    """
    if DMSO:
        Tm -= DMSOfactor * DMSO
    if fmd:
        # McConaughy et al. (1969), Biochemistry 8: 3289-3295
        if fmdmethod == 1:
            Tm -= fmdfactor * fmd                 # Note: fmd is percent
        # Blake & Delcourt (1996), Nucl Acids Res 11: 2095-2103
        if fmdmethod == 2:
            if GC is None or GC < 0:
                raise ValueError('\'GC\' is missing or negative')
            Tm += (0.453 * (GC / 100.0) - 2.88) * fmd   # Note: fmd is molar
        if fmdmethod not in (1, 2):
            raise ValueError('\'fmdmethod\' must be 1 or 2')
    return Tm
