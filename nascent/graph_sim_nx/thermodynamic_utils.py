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

# pylint:disable=C0103,R0914

"""

Module for

"""

import math

from .constants import N_AVOGADRO, R
from nascent.energymodels.biopython import hybridization_dH_dS




def thermodynamic_meltingcurve(T_start, T_finish, delta_T, volume,
                               domains_by_name, domain_pairs, domain_dHdS):
    """
    This function currently only supports one mapping in domain_pairs:
    Returns a list/ordered dict of
        temperature : [hybridized, non_hybridized, total_conc, domains_total]
    """
    # from collections import OrderedDict
    cum_stats = {} #OrderedDict()       # Has T: <fraction of duplexed domains>
    domain_stats = {} #OrderedDict()    # Same as cum_stats, but for individual domains
    concentrations = {dname: len(domains)/(volume*N_AVOGADRO)
                      for dname, domains in domains_by_name.items()}
    T = T_start

    while T >= T_finish if delta_T < 0 else T <= T_finish:
        hybridized = 0      # can be a fraction
        non_hybridized = 0
        total_conc = 0
        domains_total = 0
        domain_stats[T] = {}
        for dname, domains in domains_by_name.items():
            if len(domains) == 0:
                continue
            try:
                compl_name = domain_pairs[dname]
            except KeyError:
                # domain has no complementary domains:
                continue
            if not isinstance(compl_name, str):
                # domain_pairs[dname] is a list
                try:
                    compl_name = sorted([name for name in compl_name if name in concentrations])[0]
                except IndexError:
                    # compl_name list is empty after filtering
                    continue
            if dname not in domain_dHdS:
                domain_dHdS[dname] = hybridization_dH_dS(domains[0].sequence)
            # standard-condition energies:
            deltaH, deltaS = domain_dHdS[dname]
            # deltaH in kcal/mol, deltaS in cal/mol/K:
            deltaG = deltaH*1000 - T*deltaS
            try:
                # If deltaG is really large, we might get overflow error
                K = math.exp(-deltaG/(R*T))
            except OverflowError:
                print("Warning: OverflowError while calculating K = math.exp(-deltaG/(R*T)) " +
                      "= math.exp(-{:0.03g}/({}*{})). Setting K=1e10.".format(deltaG, R, T))
                K = 1e10
            # We assume all partnering domains are perfect complements:
            conc_domain = concentrations[dname]
            conc_compl = concentrations[compl_name]
            # f_hyb = binary_state_probability_cal_per_mol(deltaG)
            if conc_compl > conc_domain:
                # Ai must be the higher of the two concentrations:
                D, A, B, K2 = solvetwocomponent(conc_compl, conc_domain, K)
                x = D/B
            else:
                #
                D, A, B, K2 = solvetwocomponent(conc_domain, conc_compl, K)
                x = D/A
            if abs(K2-K)/K > 0.1:
                print(("- Note: solvetwocomponent output K={:0.04g} is different " +
                       "from input K={:0.04g} (T={}, dname={}").format(K, K2, T, dname))
            # x is the fraction of domain that is hybridized:
            # Actually, we just want the concentration of duplexed ?
            # TODO: We could scale by number of bases in each domain...
            hybridized += D         # is a concentration
            non_hybridized += conc_domain - D
            total_conc += conc_domain
            domains_total += len(domains)
            domain_stats[T]['dname'] = [D, conc_domain - D, conc_domain, len(domains), x]
        cum_stats[T] = [hybridized, non_hybridized, total_conc, domains_total]
        T += delta_T
    return cum_stats, domain_stats


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
