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

# pylint: disable=C0103,W0212

"""
Module for analysing reaction graph.

"""

import math
import yaml
import networkx
nx = networkx
from networkx.classes.digraph import DiGraph

from .constants import STACKING_INTERACTION, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION


def load_reaction_graph(fn):
    with open(fn) as fp:
        graph = yaml.load(fp)
    return graph


def graph_state_partitions(g, T, update=False, unit='R'):
    """
    Calculate dG and partitions for reaction graph nodes (states).
    params:
    :g:         Graph
    :T:         Temperature in K
    :update:    If True, will update graph nodes with result.
    :unit:      The units of dH and dS. Currently only units of R*K and R are supported.
    """
    dHdS = {node: nattr.get('dH_dS_first', (0.0, 0.0)) for node, nattr in g.node.items()}
    dG = {node: dH - T*dS for node, (dH, dS) in dHdS.items()}
    # Universe entropy, dS_uni = dS_surroundings + dS_system, in units of R:
    dSuni_R = {node: dS-dH/T for node, (dH, dS) in dHdS.items()}
    partitions = {node: math.exp(dS) for node, dS in dSuni_R.items()}
    Z = sum(partitions.values())
    partitions = {node: part/Z for node, part in partitions.items()}
    if update:
        for node, nattr in g.node.items():
            nattr['partition'] = partitions[node]
            nattr['dSuni'] = dSuni_R[node]
            nattr['dG'] = dG[node]
    return partitions


def reaction_attr_str(reaction_attr):
    """
    Return
        h+ for forming  hybridization reactions,
        h- for breaking hybridization reactions,
        s+ for forming  stacking reactions,
        s- for breaking stacking reactions,
    Appends an asterix (*) for inter-molecular reactions.
    """
    return "".join((reaction_attr.reaction_type,
                    "+" if reaction_attr.is_forming else "-",
                    " " if reaction_attr.is_intra else "*"))

def reaction_spec_pair_str(reaction_spec_pair):
    """
    Return
    Appends an asterix (*) for inter-molecular reactions.
    """
    return ", ".join(sorted(reaction_spec_pair))


def reaction_str(reaction_spec_pair, reaction_attr):
    """
    Return
        h+*: s1_A > < s2_a  (0, 0)  for inter-complex hybridization
        h+ : s1_A >_< s2_a          for intra-complex hybridization
        h- : s1_A-<_>-s2_a          for de-hybridization (always intra-complex)
    for hybridization reactions, and
        s+ : s1_A3p >_< s2_B5p / s3_a5p >_< s4_b3p   (fp, fp, fp, fp)   for intra-complex stacking
        s- : s1_A3p:<_>:s2_B5p / s3_a5p:<_>:s4_b3p   (fp, fp, fp, fp)   for intra-complex stacking
    for stacking reactions, and
        b+ : s1_A3p >_< s2_B5p  (fp, fp)   for intra-complex ligation
        b- : s1_A3p.<_>.s2_B5p  (fp, fp)   for intra-complex backbone hydrolysis
    for backbon reactions.
    Some of the information is intentionally redundant to verify that everything is as it should be:
    E.g. the reaction_attr_str before ':' can be extracted from the reaction_spec_pair after the colon.
    Also the current domain hybridization state (indicated by dash after domain name in -<_>-) and the
    current domain-end stacking state (indicated by colon after end in 3p:<_>:) is (or should be)
    redundant with the indication of whether the reaction is forming (indicated by "> <" for hybridization/stacking
    reactions) or breaking (indicated by "< >" for dehybridization/unstacking reactions).
    """
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        d1, d2 = {}, {}
        d1['fp'], d2['fp'] = sorted(tuple(reaction_spec_pair))
        for d in (d1, d2):
            dspecie, is_hyb, cstate, icid = d['fp']
            d['dspecie'] = "_".join(dspecie)     # dspecie = (strand.name, domain.name)
            d['hyb'] = "-" if is_hyb else " "   # whether domain is currently hybridized
            d['cstate'] = str(cstate % 100000)   # complex state (fingerprint)
            d['icid'] = "" if icid == 0 else str(icid)  # in-complex identifier (for symmetric complexes)
        cstates = []
        for d in (d1, d2):
            if d['cstate'] not in cstates:
                cstates.append(d['cstate'])
        states_str = ", ".join(cstates)
        form_str = (">%s<" if reaction_attr.is_forming else "<%s>") % ("_" if reaction_attr.is_intra else " ")
        fmt = "{ra}:  {d1[dspecie]}{d1[hyb]}{form_str}{d2[hyb]}{d2[dspecie]}\t({states_str})"
        return fmt.format(ra=reaction_attr_str(reaction_attr), form_str=form_str, d1=d1, d2=d2, states_str=states_str)
    elif reaction_attr.reaction_type is STACKING_INTERACTION:
        h1end3p, h2end5p, h2end3p, h1end5p = {}, {}, {}, {}
        for d, efp in zip((h1end3p, h2end5p, h2end3p, h1end5p),
                          (efp for pair in tuple(reaction_spec_pair) for efp in pair)):
            # DomainEnd fingerprint = (domain.state_fingerprint, self.end,self.is_stacked())
            dfp, d['end'], is_stacked = efp
            # domain attributes:
            dspecie, is_hyb, cstate, icid = dfp
            d['dspecie'] = "_".join(dspecie)     # dspecie = (strand.name, domain.name)
            d['stacked'] = ":" if is_stacked else " "   # whether domain is currently hybridized
            d['cstate'] = str(cstate % 100000)   # complex state (fingerprint)
            d['icid'] = "" if icid == 0 else str(icid)  # in-complex identifier (for symmetric complexes)
        cstates = []
        for d in (h1end3p, h2end5p, h2end3p, h1end5p):
            if d['cstate'] not in cstates:
                cstates.append(d['cstate'])
        states_str = ", ".join(cstates)
        form_str = (">%s<" if reaction_attr.is_forming else "<%s>") % ("_" if reaction_attr.is_intra else " ")
        fmt = ("{ra}:"
               "{h1end3p[dspecie]}{h1end3p[end]}{h1end3p[stacked]}{form_str}"
               "{h1end5p[stacked]}{h1end5p[dspecie]}{h1end5p[end]} / "
               "{h1end3p[dspecie]}{h1end3p[end]}{h1end3p[stacked]}{form_str}"
               "{h2end5p[stacked]}{h2end5p[dspecie]}{h2end5p[end]} \t"
               #"({h1end3p[cstate]}, {h1end5p[cstate]}, {h2end3p[cstate]}, {h2end5p[cstate]})"
               "({states_str})"
              )
        return fmt.format(ra=reaction_attr_str(reaction_attr), form_str=form_str, states_str=states_str,
                          h1end3p=h1end3p, h2end5p=h2end5p, h2end3p=h2end3p, h1end5p=h1end5p)


def reaction_eattrs(reaction_attr, activity, c_j, throttle_factor, dHdS):
    """
    """
    return


def reaction_edge_label(eattr):
    """
    Add a suitable label to reaction edge attribute dict :eattr:.
    eattr dict should already contain the following entries:
    * reaction_attr tuple values:
    * reaction_type, is_forming, is_intra
    * is_forming_str
    * is_intra_str
    * dS
    * dH
    * loop type enum (0=no loop, 1, 2)
    * activity
    * c_j
    * throttle_factor
    """
    return "{reaction_type}{is_forming_str}{is_intra_str} {c_j:0.1e}"

