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

# import pytest
import sys
# from random import choice
# from pprint import pprint
import pdb
from math import log
ln = log
try:
    from numpy import isclose
    import numpy as np
except ImportError:
    def isclose(a, b, rtol=1e-5, atol=1e-8):
        """ Default rtol, atol is actually pretty lenient, only compares to about 5 significant digits..."""
        return abs(a-b) <= (atol + rtol*abs(b))
sys.path.insert(0, ".")


from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
# from nascent.graph_sim_nx.graph_manager import GraphManager
# from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx.reactionmgr import ReactionMgr
# from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.reaction_utils import get_reaction_spec_pair, reaction_to_str
from nascent.graph_sim_nx.constants import STACKING_INTERACTION, HYBRIDIZATION_INTERACTION, RA_UNSTACK_INTRA
from nascent.graph_sim_nx.constants import ReactionAttrs, RA_HYB_INTRA, RA_HYB_INTER, RA_DEHYB_INTRA, RA_STACK_INTRA
from nascent.graph_sim_nx.constants import I_DH, I_DS, I_HYBRIDIZATION, I_STACKING, I_LOOP, I_VOLUME
from nascent.graph_sim_nx.utils import prod, WC, rcompl, AttrDict, tupleify
from nascent.graph_sim_nx.looptracking import effectuate_loop_changes




def test_loop_tracking_01_fourway():
    """
#         A       .        B
# 5' ------------' `--------------- 3'
# 3' -----a------. .-------b------- 5'
#                | |                  .
# 5' -----d------' `-------c------- 3'
# 3' ------------. .--------------- 5'
#         D                C

reaction h+*:  s2_c > < s3_C    (s2, s3)
reaction h+*:  s1_A > < s4_a    (s1, s4)
reaction h+*:  s1_B > < s2_b    (5412, 7374)
reaction h+ :  s3_D >_< s4_d    (73404)

    """
    ## Case fourway_junction_1 reaction path 01:
    a = Domain("a", seq="GCTA"*4)
    A = Domain("A", seq="TAGC"*4)
    b = Domain("b", seq="CGTA"*4)
    B = Domain("B", seq="TACG"*4)
    c = Domain("c", seq="TA"*8)
    C = Domain("C", seq="TA"*8)
    d = Domain("d", seq="GCTT"*4)
    D = Domain("D", seq="AAGC"*4)
    # domain_pairs = # No need to specify domain_pairs; will use upper/lower pairing scheme if not given.

    s1 = Strand("s1", [A, B])
    s2 = Strand("s2", [b, c]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [C, D])
    s4 = Strand("s4", [d, a])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ReactionMgr([s1, s2, s3, s4], params={}, volume=1e-12, init_reactions=True)
    # Or use a ReactionMgr if you want to use hybridize_and_process


    # result = mgr.hybridize(c, C)
    reaction_pair, reaction_attr = (c, C), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair) # I guess I'm still using frozensets
    reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    cmplx0 = result['new_complexes'][0]

    # result = mgr.hybridize(a, A)
    reaction_pair, reaction_attr = (a, A), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    cmplx1 = result['new_complexes'][0]

    # result = mgr.hybridize(B, b)
    reaction_pair, reaction_attr = (B, b), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    cmplx = result['changed_complexes'][0]
    assert cmplx == cmplx0 or cmplx == cmplx1

    # result = mgr.hybridize(D, d)
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p

    # result = mgr.stack(C.end3p, c.end5p, d.end3p, D.end5p)
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Unstack C-d, then stack again:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]

    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-b then unstack A-b:
    # result = mgr.stack(A.end3p, a.end5p, b.end3p, B.end5p)
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]

    # Unstack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-b, stack C-d, unstack A-b, unstack C-d:
    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Stack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-b, stack C-d, unstack A-b, unstack C-d, again:
    for reaction_pair, reaction_attr in (
        (((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA),
        (((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA),
        (((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA),
        (((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA),
    ):
        reaction_pair = frozenset(reaction_pair)
        reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
        assert cmplx == result['changed_complexes'][0]

    # To quickly get adj list:
    # pp sorted([(str(src), str(tgt), ekey) for src, tgts in self.ends5p3p_graph.adj.items() for tgt, edges in tgts.items() for ekey in edges if ekey != 'b'])




def form_and_check(reaction_pair, reaction_attr, cmplx, mgr, update_reactions=False):

    reaction_str = reaction_to_str(None, reaction_attr, reaction_pair)
    print("\n\nPerforming form_and_check: ", reaction_str)

    elem1, elem2 = tuple(reaction_pair)
    reaction_pair = frozenset(reaction_pair)
    energy_subtotals_before = cmplx.energy_subtotals.copy() # For numpy arrays we don't need to deepcopy
    activity, loop_effects = mgr.loop_formation_effects(elem1, elem2, reaction_attr.reaction_type)
    if not activity > 0:
        print("ERROR: loop_effects['total_loop_activity_change'] = %s" % loop_effects['total_loop_activity_change'])
        pdb.set_trace()
    # Note: loop_formation_effects can return (0, None) if the shortest-path activity is 0, i.e. reaction is impossible.
    # np.prod will return a prod data-type if input is a generator; make sure you use a list input
    assert loop_effects['total_loop_activity_change'] == activity
    assert activity == (loop_effects['new_loop']['activity'] *
                        np.prod([loop['activity_change_ratio'] for loop in loop_effects['changed_loops_by_hash'].values()]))
    dS_loop = ln(activity)
    reaction_energy = np.zeros((4, 2), dtype=float)  # 4 types of energy (rows), each with two scalar values (columns)
    reaction_energy[I_LOOP, I_DS] = dS_loop
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        dHdS_hyb = mgr.hybridization_energy(elem1, elem2)
        reaction_energy[I_HYBRIDIZATION, :] = dHdS_hyb
    else:
        dHdS_stack = mgr.stacking_energy(elem1, elem2)
        reaction_energy[I_STACKING, :] = dHdS_stack
    energy_subtotals_expected = energy_subtotals_before + reaction_energy
    # Make sure cmplx energy has been re-calculated:
    cmplx.recalculate_complex_energy(mgr.volume_entropy)
    # Always use np.all(ndarr) or ndarr.all(), NOT all(ndarr)
    assert np.all(np.isclose(energy_subtotals_before, cmplx.energy_subtotals))
    # assert np.all(energy_subtotals_before == cmplx.energy_subtotals)
    # I believe it defaults to np.isclose for equality - edit: ALWAY use isclose.


    # # Perform reaction:
    # if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
    #     reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    # else:
    #     reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)

    # EDIT: Doing reactions manually step-by-step.

    # 1. De-hybridize/unstack
    result = (mgr.hybridize(elem1, elem2) if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION else
              mgr.stack(elem1, elem2))
    assert cmplx == result['changed_complexes'][0]

    # 3. effectuate_loop_changes:
    effectuate_loop_changes(cmplx, loop_effects, reaction_attr.is_forming) # returns (0, None) if reaction is impossible.

    # 4. assert_state_change:
    reaction_spec_pair = get_reaction_spec_pair(elem1, elem2, reaction_attr)
    cmplx.assert_state_change(reaction_spec_pair, reaction_attr)
    cmplx.recalculate_complex_energy(mgr.volume_entropy)

    energy_subtotals_after = cmplx.energy_subtotals.copy() # For numpy arrays we don't need to deepcopy
    cmplx.recalculate_complex_energy(mgr.volume_entropy)
    # Always use np.all(ndarr) or ndarr.all(), NOT all(ndarr)
    assert np.all(np.isclose(energy_subtotals_after, cmplx.energy_subtotals))


    # 5. update_possible_reactions (optional, if you want to retain reactionmgr possible_reactions data integrity:
    if update_reactions:
        changed_domains = cmplx.domains()
        mgr.update_possible_hybridization_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)
        mgr.update_possible_stacking_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)

    # Check that reaction energy matches difference in before vs after:
    assert np.isclose(energy_subtotals_after, energy_subtotals_before + reaction_energy).all()
    print("-completed form_and_check: ", reaction_str, '\n')
    return loop_effects


def check_complex_ifnodes(cmplx, mgr=None):
    # Check that ifnode.top_delegates() gives the same result as checking ifnode.delegatee is None:
    set1 = set([end.ifnode for d in cmplx.domains() for end in (d.end5p, d.end3p) if end.ifnode.delegatee is None])
    set2 = set([end.ifnode.top_delegate() for d in cmplx.domains() for end in (d.end5p, d.end3p)])
    assert set1 == set2
    if mgr:
        # Check that complex ifnodes are the same as mgr interface graph ifnodes:
        # (If so, then top_delegates should also be the same...)
        assert set(mgr.interface_graph.nodes()) == set([end.ifnode for d in cmplx.domains() for end in (d.end5p, d.end3p)])
    return set1


def break_and_check(reaction_pair, reaction_attr, cmplx, mgr, update_reactions=False):
    # Note: loop_breakage_effects HAS SIDE-EFFECTS and is NOT an idempotent operation.
    # That is, it WILL update existing paths to make them reflect the current ifnode top delegation.
    # However, the following hybridize will not update the loop paths, not even if you invoke loop_formation_effects(),
    # because this IS designed to be an indempotent function as it is used predictively (for a range of
    # possible reactions of which at most one of them is actually performed).
    # You COULD argue that ifnode top_delegation path updating should be done as a separate operation.
    # Would it be possible to do loop_breakage_effects without altering the paths ifnodes,
    # and just do that as an OPTIONAL operation?
    # Maybe, but it would likely require a complete re-design of loop_breakage_effects().

    reaction_pair = frozenset(reaction_pair)
    elem1, elem2 = tuple(reaction_pair)
    reaction_str = reaction_to_str(None, reaction_attr, reaction_pair)
    print("\n\nPerforming break_and_check: ", reaction_str)

    energy_subtotals_before = cmplx.energy_subtotals.copy() # For numpy arrays we don't need to deepcopy
    # Make sure cmplx energy has been re-calculated:
    cmplx.recalculate_complex_energy(mgr.volume_entropy)
    # Always use np.all(ndarr) or ndarr.all(), NOT all(ndarr)
    assert np.all(np.isclose(energy_subtotals_before, cmplx.energy_subtotals))
    # assert np.all(energy_subtotals_before == cmplx.energy_subtotals)
    # I believe it defaults to np.isclose for equality - edit: ALWAY use isclose.

    # Loop breakage must be calculated *after* the reaction but *before* resetting/asserting state change.
    # It should be OK to simulate this by dehybridizing, then calculating loop effects, then re-hybridizing
    # TODO: It would probably be better to do all steps manually:
    # dehybridize, calculate loop_breakage_effects, effectuate_loop_changes, assert_state_change
    # I believe that should be just about it.
    # TODO: Compare activity from loop_formation_effects vs loop_breakage_effects.

    # 1. De-hybridize/unstack
    result = (mgr.dehybridize(elem1, elem2) if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION else
              mgr.unstack(elem1, elem2))
    assert cmplx == result['changed_complexes'][0]

    # 2. calculate loop_breakage_effects:
    loop_effects = mgr.loop_breakage_effects(elem1, elem2, reaction_attr.reaction_type) # OBS: _cached() returns tuple not dict
    if not loop_effects['total_loop_activity_change'] > 0:
        print("ERROR: loop_effects['total_loop_activity_change'] = %s" % loop_effects['total_loop_activity_change'])
        pdb.set_trace()

    # 3. effectuate_loop_changes:
    effectuate_loop_changes(cmplx, loop_effects, reaction_attr.is_forming) # returns (0, None) if reaction is impossible.

    ## 3b: Calculate reaction energy:
    assert isclose(loop_effects['total_loop_activity_change'], (
        1/loop_effects['del_loop']['activity'] * # inverted value since we DELETE, not form the loop
        np.prod([loop['activity_change_ratio'] for loop in loop_effects['changed_loops_by_hash'].values()])))
    dS_loop = ln(loop_effects['total_loop_activity_change'])
    reaction_energy = np.zeros((4, 2), dtype=float)  # 4 types of energy (rows), each with two scalar values (columns)

    # minus because breaking # edit: it seems the energy is already negated because it is natively describing a "breaking" reaction.
    reaction_energy[I_LOOP, I_DS] = dS_loop
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        dHdS_hyb = mgr.hybridization_energy(elem1, elem2)
        reaction_energy[I_HYBRIDIZATION, :] -= dHdS_hyb
    else:
        dHdS_stack = mgr.stacking_energy(elem1, elem2)
        reaction_energy[I_STACKING, :] -= dHdS_stack

    # 4. assert_state_change:
    reaction_spec_pair = get_reaction_spec_pair(elem1, elem2, reaction_attr)
    cmplx.assert_state_change(reaction_spec_pair, reaction_attr)
    cmplx.recalculate_complex_energy(mgr.volume_entropy)

    # 5. update_possible_reactions (optional, if you want to retain reactionmgr possible_reactions data integrity:
    if update_reactions:
        changed_domains = cmplx.domains()
        mgr.update_possible_hybridization_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)
        mgr.update_possible_stacking_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)


    # Check energy difference (after - before) vs reaction energy
    energy_subtotals_after = cmplx.energy_subtotals.copy() # For numpy arrays we don't need to deepcopy
    cmplx.recalculate_complex_energy(mgr.volume_entropy)
    # Always use np.all(ndarr) or ndarr.all(), NOT all(ndarr)
    assert np.all(np.isclose(energy_subtotals_after, cmplx.energy_subtotals))

    # Check that reaction energy matches difference in before vs after:
    assert np.isclose(energy_subtotals_after, energy_subtotals_before + reaction_energy).all()
    print("Comparison:")
    print(reaction_energy)
    print(energy_subtotals_after - energy_subtotals_before)
    print("Difference:")
    print(energy_subtotals_after - energy_subtotals_before - reaction_energy)

    print("-completed break_and_check: ", reaction_str, '\n')

    return loop_effects


def react(reaction_pair, reaction_attr, cmplx, mgr, loop_effects=None, update_reactions=False):
    """ Perform reaction and assert complex state change incl energy re-calculation. """

    elem1, elem2 = tuple(reaction_pair)
    reaction_str = reaction_to_str(None, reaction_attr, reaction_pair)
    print("\n\nPerforming reaction: ", reaction_str)

    # 1. Form/break interaction - Hybridize/stack/dehybridize/unstack:
    if reaction_attr.is_forming:
        result = (mgr.hybridize(elem1, elem2) if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION else
                  mgr.stack(elem1, elem2))
    else:
        result = (mgr.dehybridize(elem1, elem2) if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION else
                  mgr.unstack(elem1, elem2))

    if cmplx is None:
        cmplx = result['new_complexes'][0] if result['new_complexes'] else result['changed_complexes'][0]

    # 3. Effectuate loop changes (if required) - must be done BEFORE asserting state changes..
    if loop_effects:
        effectuate_loop_changes(cmplx, loop_effects, reaction_attr.is_forming)

    # 4. assert_state_change:
    reaction_spec_pair = get_reaction_spec_pair(elem1, elem2, reaction_attr)
    cmplx.assert_state_change(reaction_spec_pair, reaction_attr)
    cmplx.recalculate_complex_energy(mgr.volume_entropy)

    # 5. update_possible_reactions (optional, if you want to retain reactionmgr possible_reactions data integrity:
    if update_reactions:
        changed_domains = cmplx.domains()
        mgr.update_possible_hybridization_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)
        mgr.update_possible_stacking_reactions(
            changed_domains, reaction_pair, reaction_attr, reaction_spec_pair)

    print("-completed reaction: ", reaction_str, '\n')
    return result


reaction_name_by_reaction_attr = {
    RA_STACK_INTRA: "stacking",
    RA_UNSTACK_INTRA: "unstacking",
    RA_HYB_INTRA: "hybridizing",
    RA_DEHYB_INTRA: "dehybridizing",
}

def get_reaction_desc(reaction_pair, reaction_attr, state_desc):

    rev_reaction_attr = ReactionAttrs(reaction_attr.reaction_type, not reaction_attr.is_forming, reaction_attr.is_intra)
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        pair_desc = "-".join(d.name for d in reaction_pair)
    else:
        assert reaction_attr.reaction_type is STACKING_INTERACTION
        pair_desc = "-".join(dupend[0].domain.name for dupend in reaction_pair)
    desc = " ".join([reaction_name_by_reaction_attr[reaction_attr], pair_desc,
                     state_desc.format("form" if reaction_attr.is_forming else "break")])
    reverse_desc = " ".join([reaction_name_by_reaction_attr[rev_reaction_attr], pair_desc,
                             state_desc.format("form" if reaction_attr.is_forming else "break")])
    return desc, reverse_desc



def test_loop_energy_01_fourway():
    """
    Test that the loop activity (dS) matches the difference in
    Complex loop energy sub-total.

#         A       .        B
# 5' ------------' `--------------- 3'
# 3' -----a------. .-------b------- 5'
#                | |                  .
# 5' -----d------' `-------c------- 3'
# 3' ------------. .--------------- 5'
#         D                C

reaction h+*:  s2_c > < s3_C    (s2, s3)
reaction h+*:  s1_A > < s4_a    (s1, s4)
reaction h+*:  s1_B > < s2_b    (5412, 7374)
reaction h+ :  s3_D >_< s4_d    (73404)

    """
    ## Case fourway_junction_1 reaction path 01:
    a = Domain("a", seq="GCTA"*4)
    A = Domain("A", seq="TAGC"*4)
    b = Domain("b", seq="CGTA"*4)
    B = Domain("B", seq="TACG"*4)
    c = Domain("c", seq="TA"*8)
    C = Domain("C", seq="TA"*8)
    d = Domain("d", seq="GCTT"*4)
    D = Domain("D", seq="AAGC"*4)
    # domain_pairs = # No need to specify domain_pairs; will use upper/lower pairing scheme if not given.

    s1 = Strand("s1", [A, B])
    s2 = Strand("s2", [b, c]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [C, D])
    s4 = Strand("s4", [d, a])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ReactionMgr([s1, s2, s3, s4], params={}, volume=1e-12, init_reactions=True)
    # Or use a ReactionMgr if you want to use hybridize_and_process


    # result = mgr.hybridize(c, C)
    reaction_pair, reaction_attr = (c, C), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair) # I guess I'm still using frozensets
    # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    result = react(reaction_pair, reaction_attr, None, mgr, loop_effects=None, update_reactions=False)
    cmplx0 = result['new_complexes'][0]

    # result = mgr.hybridize(a, A)
    reaction_pair, reaction_attr = (a, A), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair)
    # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    result = react(reaction_pair, reaction_attr, None, mgr, loop_effects=None, update_reactions=False)
    cmplx1 = result['new_complexes'][0]

    # result = mgr.hybridize(B, b)
    reaction_pair, reaction_attr = (B, b), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair)
    # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    result = react(reaction_pair, reaction_attr, None, mgr, loop_effects=None, update_reactions=False)
    cmplx = result['changed_complexes'][0]
    assert cmplx == cmplx0 or cmplx == cmplx1


    loop_effects_cache = {}
    reaction_activities = {}
    desc_cstates = {}   # desc => cstate, reaction_pair, reaction_attr

    def loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc=None):
        """
        # To debug:
        pp reaction_pair, reaction_attr, desc, reverse_desc, activity, 1/activity, reaction_activities[desc], reaction_activities.get(reverse_desc)
        pp reaction_pair, reaction_attr, desc, cstate, desc_cstates[desc], cstate_by_desc.get((cstate, reaction_pair, reaction_attr)), desc_cstates
        """
        print([str(v) for v in (reaction_pair, reaction_attr, desc, reverse_desc)])
        assert desc.split()[0] == reaction_name_by_reaction_attr[reaction_attr]
        rev_reaction_attr = ReactionAttrs(reaction_attr.reaction_type, not reaction_attr.is_forming, reaction_attr.is_intra)
        assert reverse_desc.split()[0] == reaction_name_by_reaction_attr[rev_reaction_attr]
        if reaction_attr.is_forming:
            loop_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
        else:
            loop_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
        activity = loop_effects['total_loop_activity_change']
        cstate_by_desc = {v: k for k, v in desc_cstates.items()}
        cstate = cmplx.state_fingerprint()
        if desc in loop_effects_cache:
            # assert loop_effects_cache[desc] == tupleify(loop_effects)
            assert desc_cstates[desc] == (cstate, reaction_pair, reaction_attr)
            assert reaction_activities[desc] == activity  # exact or just "close"?
        else:
            loop_effects_cache[desc] = tupleify(loop_effects)
            reaction_activities[desc] = activity
            assert desc not in desc_cstates
            desc_cstates[desc] = (cstate, reaction_pair, reaction_attr)
        if reverse_desc:
            if reverse_desc in reaction_activities:
                assert isclose(reaction_activities[reverse_desc], 1/activity)
            else:
                print(reverse_desc, "not in reaction_activities")


    ## Forming the first loop by hybridization:
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    desc, reverse_desc = "hybridizing d-d to form loop", "dehybridizing d-d to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # loop_effects_cache[desc] = loop_formation_effects


    ## BREAKING the first loop by hybridization:
    reaction_pair, reaction_attr = (d, D), RA_DEHYB_INTRA
    desc, reverse_desc = "dehybridizing d-d to break loop", "hybridizing d-d to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # loop_effects_cache["dehybridizing d-d to break loop"] = loop_formation_effects
    # # We are currently inverting the activity for breaking reactions:
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change'])


    ## Re-Forming the first loop by hybridization:
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    desc, reverse_desc = "hybridizing d-d to form loop", "dehybridizing d-d to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                loop_effects_cache["hybridizing d-d to form loop"]['total_loop_activity_change'])



    # STACKING ENDS:
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p



    # STACK C-d:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking C-d with A-b unstacked", "unstacking C-d with A-b unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # loop_effects_cache["stacking C-d with A-b unstacked"] = loop_formation_effects


    # Unstack C-d, then stack again:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking C-d with A-b unstacked", "stacking C-d with A-b unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change'])
    # loop_effects_cache["unstacking C-d with A-b unstacked"] = loop_breakage_effects

    # STACK C-d:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking C-d with A-b unstacked", "unstacking C-d with A-b unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change'])


    # Stack A-b then unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b with C-d stacked", "unstacking A-b with C-d stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # loop_effects_cache["stacking A-b with C-d stacked"] = loop_formation_effects

    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b with C-d stacked", "stacking A-b with C-d stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change'])
    # loop_effects_cache["unstacking A-b with C-d stacked"] = loop_breakage_effects


    # Unstack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking C-d with A-b unstacked", "stacking C-d with A-b unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert loop_breakage_effects == loop_effects_cache["unstacking C-d with A-b unstacked"]
    # Edit: Hard to compare dicts, even if we tupleify them, so just comparing total_loop_activity_change values for now.
    # assert isclose(loop_breakage_effects['total_loop_activity_change'],
    #                loop_effects_cache["unstacking C-d with A-b unstacked"]['total_loop_activity_change'])



    # Stack A-b, stack C-d, unstack A-b, unstack C-d:
    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b with C-d unstacked", "unstacking A-b with C-d unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change']) # A-b should be identical to C-d

    # Stack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking C-d with A-b stacked", "unstacking C-d with A-b stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_formation_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert "stacking C-d with A-b stacked" not in loop_effects_cache
    # loop_effects_cache["stacking C-d with A-b stacked"] = loop_formation_effects
    # Assert symmetry:
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                loop_effects_cache["stacking A-b with C-d stacked"]['total_loop_activity_change'])


    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b with C-d stacked", "stacking A-b with C-d stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_formation_effects['total_loop_activity_change'],
    #                1/loop_breakage_effects['total_loop_activity_change'])
    # assert isclose(loop_breakage_effects['total_loop_activity_change'],
    #                loop_effects_cache["unstacking A-b with C-d stacked"]['total_loop_activity_change'])

    # Unstack C-d
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking C-d with A-b unstacked", "stacking C-d with A-b unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert isclose(loop_breakage_effects['total_loop_activity_change'],
    #                loop_effects_cache["unstacking C-d with A-b unstacked"]['total_loop_activity_change'])



    # TODO:  Maybe use "with all hybridized None stacked" as state_desc? But that wouldn't work for the reversible reaction....
    # Then you would indeed have to describe the state in full, and differently for forward vs reverse reaction.


    ### Stack non-neighboring helices:

    # Stack A-C:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_STACK_INTRA
    desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
    assert (desc, reverse_desc) == ("stacking A-C with b-d unstacked", "unstacking A-C with b-d unstacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-C:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    ## Stack A-C, stack b-d, unstack A-C, unstack b-d

    if False:
        # Stack A-C:
        reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_STACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Stack b-d: -- THIS CURRENTLY GIVES AN ERROR!
        reaction_pair, reaction_attr = ((b.end3p, B.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with A-C stacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Unstack A-C:
        reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_UNSTACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d stacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Unstack b-d:
        reaction_pair, reaction_attr = ((b.end3p, B.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with A-C unstacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    ## Stack A-C, stack b-d, unstack b-d, unstack A-C

    if False:  # disabled, gives an error because loop degeneracy (I think)
        # Stack A-C:
        reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_STACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Stack b-d: -- THIS CURRENTLY GIVES AN ERROR!
        reaction_pair, reaction_attr = ((b.end3p, B.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with A-C stacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Unstack b-d:
        reaction_pair, reaction_attr = ((b.end3p, B.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with A-C stacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

        # Unstack A-C:
        reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_UNSTACK_INTRA
        desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
        loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)




    ## BREAKING the loop by hybridization:
    reaction_pair, reaction_attr = (d, D), RA_DEHYB_INTRA
    desc = "dehybridizing d-d to break loop"
    desc, reverse_desc = "dehybridizing d-d to break loop", "hybridizing d-d to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)
    # loop_breakage_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
    # assert desc not in loop_effects_cache
    # loop_effects_cache["dehybridizing d-d to break loop"] = loop_formation_effects
    # We are currently inverting the activity for breaking reactions:
    # assert isclose(loop_breakage_effects['total_loop_activity_change'],
    #                loop_effects_cache[desc]['total_loop_activity_change'])
    # assert isclose(1/loop_breakage_effects['total_loop_activity_change'],
    #                loop_effects_cache["hybridizing d-d to form loop"]['total_loop_activity_change'])






    # TODO: Continue, forming loops by stacking without all strands being hybridized:
    # You can still stack/unstack: A-b, A-C, b-C in all order permutations


    ## Produce the "pure" stacking loop formation/breaking data:

    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b to form loop", "unstacking A-b to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b to break loop", "stacking A-b to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Stack A-C:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-C to form loop", "unstacking A-C to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-C:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-C to break loop", "stacking A-C to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    # Stack C-b:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking b-C to form loop", "unstacking b-C to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack C-b:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking b-C to break loop", "stacking b-C to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)



    ### Mixed stacking + hybridization loop formation/breakage permuations:

    ## Permutation 1: stack, hybridize, dehybridize, unstack

    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b to form loop", "unstacking A-b to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Hybridize d-d:
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    desc, reverse_desc = "hybridizing d-d with A-b stacked", "dehybridizing d-d with A-b stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Dehybridize d-d:
    reaction_pair, reaction_attr = (d, D), RA_DEHYB_INTRA
    desc, reverse_desc = "dehybridizing d-d with A-b stacked", "hybridizing d-d with A-b stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b to break loop", "stacking A-b to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    ## Permutation 2: stack, hybridize, unstack, dehybridize

    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b to form loop", "unstacking A-b to break loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Hybridize d-d:
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    desc, reverse_desc = "hybridizing d-d with A-b stacked", "dehybridizing d-d with A-b stacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b with C-d unstacked", "stacking A-b with C-d unstacked"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Dehybridize d-d:
    reaction_pair, reaction_attr = (d, D), RA_DEHYB_INTRA
    desc, reverse_desc = "dehybridizing d-d to break loop", "hybridizing d-d to form loop"
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    ## Permutation 3: stack A-b, hybridize D, stack C-d, unstack A-b, unstack C-d, dehybridize D


    # Stack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking A-b to form loop", "unstacking A-b to break loop"
    assert desc, reverse_desc == get_reaction_desc(reaction_pair, reaction_attr, "to {} loop")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Hybridize D:
    reaction_pair, reaction_attr = (d, D), RA_HYB_INTRA
    desc, reverse_desc = "hybridizing d-d with A-b stacked", "dehybridizing d-d with A-b stacked"
    assert desc, reverse_desc == get_reaction_desc(reaction_pair, reaction_attr, "with A-b stacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Stack C-d:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    desc, reverse_desc = "stacking C-d with A-b stacked", "unstacking C-d with A-b stacked"
    assert desc, reverse_desc == get_reaction_desc(reaction_pair, reaction_attr, "with A-b stacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack A-b:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking A-b with C-d stacked", "stacking A-b with C-d stacked"
    assert (desc, reverse_desc) == get_reaction_desc(reaction_pair, reaction_attr, "with C-d stacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Unstack C-d:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    desc, reverse_desc = "unstacking C-d with A-b unstacked", "stacking C-d with A-b unstacked"
    assert desc, reverse_desc == get_reaction_desc(reaction_pair, reaction_attr, "with A-b unstacked")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)

    # Dehybridize d-d:
    reaction_pair, reaction_attr = (d, D), RA_DEHYB_INTRA
    desc, reverse_desc = "dehybridizing d-d to break loop", "hybridizing d-d to form loop"
    assert desc, reverse_desc == get_reaction_desc(reaction_pair, reaction_attr, "to {} loop")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


    # print("\n\nReturning pre-maturely...\n\n")
    # return

    # Now you can do all sorts of cross comparison with the cache...

def test_circfb_01():
    """
    Design from ./examples/single_duplex/circfb_1.txt
####          Circfb design (Circular FoldBack):        ####
#                            E          A          T
# 5' .---------------------------3'5'---------------------.
#    |                 .----------------------.            |
#    |                (                        )           |
#    |                 \_____c___3'5'___b_____/            |
#    |_____________________________________________________|
#            D               C          B
# Domains:
# - A, a	CGACTGGAAAGCGGGC, GCCCGCTTTCCAGTCG
# - T, TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
# - B, b	TGACCATGATTACGAA, TTCGTAATCATGGTCA
# - C, c	TGGCGCCCAATACGCA, TGCGTATTGGGCGCCA
# - D	TGATTTATAAGGGATTTTGCCGATTTCGGAAC
# - E, h2b	GGCACGACAGGTTTCC, GGAAACCTGTCGTGCC

    """
    ## Case fourway_junction_1 reaction path 01:
    a = Domain("a", seq="GCTA"*4)
    A = Domain("A", seq="TAGC"*4)
    T = Domain("T", seq="T"*20)
    b = Domain("b", seq="CGTA"*4)
    B = Domain("B", seq="TACG"*4)
    c = Domain("c", seq="TA"*8)
    C = Domain("C", seq="TA"*8)
    # d = Domain("d", seq="GCTT"*4)
    D = Domain("D", seq="AAGC"*4)
    E = Domain("E", seq="GGCACGACAGGTTTCC")
    e = Domain("e", seq="GGAAACCTGTCGTGCC")
    # domain_pairs = # No need to specify domain_pairs; will use upper/lower pairing scheme if not given.

    s1 = Strand("templ", [A, T, B, C, D, E])
    s2 = Strand("splint", [b, a, e, c]) # Make sure domains are ordered 5'->3' !
    # s3 = Strand("s3", [C, D])
    # s4 = Strand("s4", [d, a])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ReactionMgr([s1, s2], params={}, volume=1e-12, init_reactions=True)
    # Use a ReactionMgr if you want to use hybridize_and_process

    reaction_pair, reaction_attr = (A, a), RA_HYB_INTER
    reaction_pair = frozenset(reaction_pair) # I guess I'm still using frozensets
    # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    result = react(reaction_pair, reaction_attr, None, mgr, loop_effects=None, update_reactions=False)
    cmplx = result['new_complexes'][0]


    # reaction_pair, reaction_attr = (b, B), RA_HYB_INTRA
    # reaction_pair = frozenset(reaction_pair) # I guess I'm still using frozensets
    # # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    # result = react(reaction_pair, reaction_attr, cmplx, mgr, loop_effects=None, update_reactions=False)
    # assert cmplx == result['changed_complexes'][0]
    #
    #
    # reaction_pair, reaction_attr = (c, C), RA_HYB_INTRA
    # reaction_pair = frozenset(reaction_pair) # I guess I'm still using frozensets
    # # reacted_pair, result = mgr.hybridize_and_process(reaction_pair, reaction_attr=reaction_attr)
    # result = react(reaction_pair, reaction_attr, cmplx, mgr, loop_effects=None, update_reactions=False)
    # assert cmplx == result['changed_complexes'][0]

    loop_effects_cache = {}
    reaction_activities = {}
    desc_cstates = {}   # desc => cstate, reaction_pair, reaction_attr

    def loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc=None):
        """
        Usage:
            loop_and_check((domain_A, domain_a),
                            RA_HYB_INTRA,
                            "hybridize A-a to form loop",
                            "dehybridize A-a to break loop")
        # To debug:
        pp reaction_pair, reaction_attr, desc, reverse_desc, activity, 1/activity, reaction_activities[desc], reaction_activities.get(reverse_desc)
        pp reaction_pair, reaction_attr, desc, cstate, desc_cstates[desc], cstate_by_desc.get((cstate, reaction_pair, reaction_attr)), desc_cstates
        """
        print([str(v) for v in (reaction_pair, reaction_attr, desc, reverse_desc)])
        assert desc.split()[0] == reaction_name_by_reaction_attr[reaction_attr]
        rev_reaction_attr = ReactionAttrs(reaction_attr.reaction_type, not reaction_attr.is_forming, reaction_attr.is_intra)
        assert reverse_desc.split()[0] == reaction_name_by_reaction_attr[rev_reaction_attr]
        if reaction_attr.is_forming:
            loop_effects = form_and_check(reaction_pair, reaction_attr, cmplx, mgr)
        else:
            loop_effects = break_and_check(reaction_pair, reaction_attr, cmplx, mgr)
        activity = loop_effects['total_loop_activity_change']
        cstate_by_desc = {v: k for k, v in desc_cstates.items()}
        cstate = cmplx.state_fingerprint()
        if desc in loop_effects_cache:
            # assert loop_effects_cache[desc] == tupleify(loop_effects)
            assert desc_cstates[desc] == (cstate, reaction_pair, reaction_attr)
            assert reaction_activities[desc] == activity  # exact or just "close"?
        else:
            loop_effects_cache[desc] = tupleify(loop_effects)
            reaction_activities[desc] = activity
            assert desc not in desc_cstates
            desc_cstates[desc] = (cstate, reaction_pair, reaction_attr)
        if reverse_desc:
            if reverse_desc in reaction_activities:
                assert isclose(reaction_activities[reverse_desc], 1/activity)
            else:
                print(reverse_desc, "not in reaction_activities")

    # Hybridize B-b
    reaction_pair, reaction_attr = (c, C), RA_HYB_INTRA
    desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "to {} loop")
    loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)




    # Stack A-C:
    # reaction_pair, reaction_attr = ((A.end3p, a.end5p), (C.end3p, c.end5p)), RA_UNSTACK_INTRA
    # desc, reverse_desc = get_reaction_desc(reaction_pair, reaction_attr, "with b-d unstacked")
    # loop_and_check(reaction_pair, reaction_attr, desc, reverse_desc)


def test_03():
    # Other reaction combinations that have yielded loop errors in the past...:
    pass



if __name__ == "__main__":

    # # Test 2:
    # test_calculate_loop_activity_2()
    # # Test 3:
    # test_calculate_loop_activity_3()
    #
    # INTRACOMPLEX_ACTIVITY tests:
    # test_loop_formation_effects_1()
    # test_loop_formation_effects_2()
    # test_loop_formation_effects_3()
    # test_loop_formation_effects_4()

    ## LOOP_FORMATION_EFFECTS tests:

    # test_loop_tracking_01_fourway()
    # test_loop_formation_effects_02()

    ## LOOP ENERGY TESTS:

    # test_loop_energy_01_fourway()

    test_circfb_01()


    print("\n\nAll ad-hoc test runs done!\n\n")
