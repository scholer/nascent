# -*- coding: utf-8 -*-
##    Copyright 2015-2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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

# pylint: disable=C0103,W0142,W0212,R0902

"""

Module with various utility functions for reactions.

"""
from nascent.graph_sim_nx.constants import HYBRIDIZATION_INTERACTION, STACKING_INTERACTION
from .constants import HYBRIDIZATION_INTERACTION, STACKING_INTERACTION

def get_reaction_spec_pair(elem1, elem2, reaction_attr):
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
    elif reaction_attr.reaction_type is STACKING_INTERACTION:
        reaction_spec_pair = frozenset(((elem1[0].state_fingerprint(), elem2[1].state_fingerprint()),
                                        (elem2[0].state_fingerprint(), elem1[1].state_fingerprint())))
    return reaction_spec_pair


def reaction_attr_to_str(reaction_attr):
    """
    Alternative name: reaction_attr_repr ?   ("to_str" vs "repr")
    Would be nice to have a generic, overloadable, Julia-style multiple-dispatch, global "repr" function...
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


def reaction_spec_pair_to_str(reaction_spec_pair):
    """
    Alternative name: reaction_spec_pair_repr ?
    Return

    Appends an asterix (*) for inter-molecular reactions.
    """
    return ", ".join(str(val) for val in sorted(reaction_spec_pair))


def reaction_to_str(reaction_spec_pair, reaction_attr, reaction_pair=None):
    """
    Alternative names: reaction_repr, reaction_as_str, reaction_str_repr ?
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
    if reaction_spec_pair is None:
        elem1, elem2 = tuple(reaction_pair)
        reaction_spec_pair = get_reaction_spec_pair(elem1, elem2, reaction_attr)
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        d1, d2 = {}, {}
        d1['fp'], d2['fp'] = sorted(tuple(reaction_spec_pair))
        for d in (d1, d2):
            dspecie, is_hyb, cstate, icid = d['fp']
            d['dspecie'] = "_".join(dspecie)     # dspecie = (strand.name, domain.name)
            d['hyb'] = "-" if is_hyb else " "   # whether domain is currently hybridized
            d['cstate'] = str(cstate % 100000) if isinstance(cstate, int) else cstate  # complex state (fingerprint)
            d['icid'] = "" if icid == 0 else str(icid)  # in-complex identifier (for symmetric complexes)
        cstates = []
        for d in (d1, d2):
            if d['cstate'] not in cstates:
                cstates.append(d['cstate'])
        states_str = ", ".join(cstates)
        form_str = (">%s<" if reaction_attr.is_forming else "<%s>") % ("_" if reaction_attr.is_intra else " ")
        fmt = "{ra}:  {d1[dspecie]}{d1[hyb]}{form_str}{d2[hyb]}{d2[dspecie]}    ({states_str})"
        return fmt.format(ra=reaction_attr_to_str(reaction_attr), form_str=form_str, d1=d1, d2=d2, states_str=states_str)
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
        #             h1end3p         h1end5p
        # Helix 1   ----------3' : 5'----------
        # Helix 2   ----------5' : 3'----------
        #             h2end5p         h2end3p
        # Note: Biopython's "stack_string" is h1end3p h1end5p/h2end5p h2end3p
        # But the format below is in format   h1end3p h1end5p/h2end3p h2end5p
        fmt = ("{ra}:"
               "{h1end3p[dspecie]}{h1end3p[end]}{h1end3p[stacked]}{form_str}"
               "{h1end5p[stacked]}{h1end5p[dspecie]}{h1end5p[end]} / "
               "{h2end3p[dspecie]}{h2end3p[end]}{h2end3p[stacked]}{form_str}"
               "{h2end5p[stacked]}{h2end5p[dspecie]}{h2end5p[end]}    "
               #"({h1end3p[cstate]}, {h1end5p[cstate]}, {h2end3p[cstate]}, {h2end5p[cstate]})"
               "({states_str})"
              )
        return fmt.format(ra=reaction_attr_to_str(reaction_attr), form_str=form_str, states_str=states_str,
                          h1end3p=h1end3p, h2end5p=h2end5p, h2end3p=h2end3p, h1end5p=h1end5p)


