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
# import pdb
try:
    from numpy import isclose
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
from nascent.graph_sim_nx.constants import STACKING_INTERACTION, HYBRIDIZATION_INTERACTION, RA_UNSTACK_INTRA
from nascent.graph_sim_nx.constants import ReactionAttrs, RA_HYB_INTRA, RA_HYB_INTER, RA_DEHYB_INTRA, RA_STACK_INTRA


WC = dict(zip("ATGC", "TACG"))
def rcompl(seq):
    return "".join(WC[b] for b in seq[::-1])


class AttrDict(dict):
    """ Allows you to use AttrDict.key instead of AttrDict["key"] """
    def __getattribute__(self, key):
        try:
            return dict.__getattribute__(self, key)
        except AttributeError:
            return self[key]

    def __setattribute__(self, key, value):
        self[key] = value



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


    # Unstack C-D, then stack again:
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]

    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-B then unstack A-B:
    # result = mgr.stack(A.end3p, a.end5p, b.end3p, B.end5p)
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack A-B:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]

    # Unstack C-D
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-B, stack C-D, unstack A-B, unstack C-D:
    # Stack A-B:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Stack C-D
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_STACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack A-B:
    reaction_pair, reaction_attr = ((A.end3p, a.end5p), (b.end3p, B.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]
    # Unstack C-D
    reaction_pair, reaction_attr = ((C.end3p, c.end5p), (d.end3p, D.end5p)), RA_UNSTACK_INTRA
    reaction_pair = frozenset(reaction_pair)
    reacted_pair, result = mgr.stack_and_process(reaction_pair, reaction_attr=reaction_attr)
    assert cmplx == result['changed_complexes'][0]


    # Stack A-B, stack C-D, unstack A-B, unstack C-D, again:
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

    test_loop_tracking_01_fourway()
    # test_loop_formation_effects_02()

    print("\n\nAll ad-hoc test runs done!\n\n")
