# -*- coding: utf-8 -*-
##    Copyright 2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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

# pylint: disable=C0103,W0201

"""
Test graph_manager module.

"""

import pytest
import sys
from random import choice
from pprint import pprint
import pdb
try:
    from numpy import isclose
except ImportError:
    def isclose(a, b, rtol=1e-5, atol=1e-8):
        """ Default rtol, atol is actually pretty lenient, only compares to about 5 significant digits..."""
        # How numpy does it:
        return abs(a-b) <= (atol + rtol*abs(b))
        # diff = abs(a-b)
        # return diff < atol and diff/b < rtol
sys.path.insert(0, ".")


from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.constants import STACKING_INTERACTION, HYBRIDIZATION_INTERACTION
from nascent.graph_sim_nx.looptracking import effectuate_loop_changes


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



@pytest.fixture
def two_duplex_loop_case2a():
    ## Case 2(a):
    e1 = Domain("e1", seq="GCTA"*4)
    E1 = Domain("E1", seq="TAGC"*4)
    e2 = Domain("e2", seq="CGTA"*4)
    E2 = Domain("E2", seq="TACG"*4)
    e3 = Domain("e3", seq="TA"*3)
    E3 = Domain("E3", seq="TA"*3)
    e4 = Domain("e4", seq="G"*5)
    E4 = Domain("E4", seq="C"*5)
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t1", seq="T"*20)

    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))

    domains = {name: val for name, val in locals().items() if isinstance(val, Domain)}
    strands = {name: val for name, val in locals().items() if isinstance(val, Strand)}

    return locals()


def test_calculate_loop_activity_2():

    ## Case 2(a):
    e1 = Domain("e1", seq="GCTA"*4)
    E1 = Domain("E1", seq="TAGC"*4)
    e2 = Domain("e2", seq="CGTA"*4)
    E2 = Domain("E2", seq="TACG"*4)
    e3 = Domain("e3", seq="TA"*3)
    E3 = Domain("E3", seq="TA"*3)
    e4 = Domain("e4", seq="G"*5)
    E4 = Domain("E4", seq="C"*5)
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*20)

    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0 = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0]
    loop0_path_rev = [end.ifnode.top_delegate() for end in reversed(loop0)]
    loop0_activity = mgr.calculate_loop_activity(loop0_path)
    loop0_activity_rev = mgr.calculate_loop_activity(loop0_path_rev)
    # print("loop0_activity_before:", loop0_activity_before)
    print("loop0_activity:", loop0_activity)
    # assert loop0_activity_before == mgr.calculate_loop_activity(loop0_path_rev)
    assert isclose(loop0_activity, loop0_activity_rev, rtol=1e-8, atol=1e-10)

    ## THE TWO LOOPS FORMED WHEN SPLITTING THE OUTER LOOP BY STACKING e4-3p and e3-5p:
    # Test simple loop activity by primary path:
    # In the future: You are not testing InterfaceNodes/Graph delegation here, just follow a single helix when possible.
    loop1 = (E4.end5p, E4.end3p, t1.end3p, t1.end5p, E3.end5p, E3.end3p)
    loop1_path = [end.ifnode.top_delegate() for end in loop1]
    loop1_path_rev = [end.ifnode.top_delegate() for end in reversed(loop1)]
    loop1_activity_before = mgr.calculate_loop_activity(loop1_path) # Without simulating stacking of the two duplexes
    # For loop1, the stacked double-helix segment is shorter than t1, which currently means that the activity is
    # calculated the same with and without stack
    loop1_activity = mgr.calculate_loop_activity(loop1_path, simulate_reaction=STACKING_INTERACTION)
    print("loop1_activity_before:", loop1_activity_before)
    print("loop1_activity:", loop1_activity)
    assert loop1_activity_before == mgr.calculate_loop_activity(loop1_path_rev)
    assert loop1_activity == mgr.calculate_loop_activity(loop1_path_rev, simulate_reaction=STACKING_INTERACTION)

    # The other loop formed when stacking e4.end3p+e3.end5p:
    loop2 = (e4.end3p, e1.end3p, e1.end5p, t2.end5p, t2.end3p, e2.end3p, e2.end5p, e3.end5p)
    loop2_path = [end.ifnode.top_delegate() for end in loop2]
    loop2_path_rev = [end.ifnode.top_delegate() for end in reversed(loop2)]
    loop2_activity_before = mgr.calculate_loop_activity(loop2_path)
    # Activity as it would be after stacking the duplex ends: (This is what is really relevant)
    # pdb.set_trace()
    loop2_activity = mgr.calculate_loop_activity(loop2_path, simulate_reaction=STACKING_INTERACTION)
    print("loop2_activity_before:", loop2_activity_before)
    print("loop2_activity:", loop2_activity)
    assert loop2_activity_before == mgr.calculate_loop_activity(loop2_path_rev)
    assert loop2_activity == mgr.calculate_loop_activity(loop2_path_rev, simulate_reaction=STACKING_INTERACTION)
    assert loop1_activity_before > loop2_activity_before > loop2_activity
    assert loop2_activity == 0


    ## THE TWO LOOPS FORMED WHEN SPLITTING THE OUTER LOOP BY STACKING E1-3p and E2-5p:






def test_calculate_loop_activity_3():

    ### Case 2(a):   ###
    e1 = Domain("e1", seq="GCTA"*4)   # 16 bp
    E1 = Domain("E1", seq="TAGC"*4)   # 16 bp
    e2 = Domain("e2", seq="CGTA"*4)   # 16 bp
    E2 = Domain("E2", seq="TACG"*4)   # 16 bp
    e3 = Domain("e3", seq="TA"*3)     #  6 bp
    E3 = Domain("E3", seq="TA"*3)     #  6 bp
    e4 = Domain("e4", seq="G"*5)      #  5 bp
    E4 = Domain("E4", seq="C"*5)      #  5 bp
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*50)
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0_desc = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0_desc]
    loop0_path_rev = [end.ifnode.top_delegate() for end in reversed(loop0_desc)]
    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    loop0_activity_before_rev = mgr.calculate_loop_activity(loop0_path_rev)
    print("loop0_activity_before:", loop0_activity_before)
    print("loop0_activity_before_rev:", loop0_activity_before_rev)
    # print("loop0_activity:", loop0_activity)
    assert isclose(loop0_activity_before, loop0_activity_before_rev, rtol=1e-8, atol=1e-10)

    ## THE TWO LOOPS FORMED WHEN SPLITTING THE OUTER LOOP BY STACKING e4-3p and e3-5p:
    # Test simple loop activity by primary path:
    # In the future: You are not testing InterfaceNodes/Graph delegation here, just follow a single helix when possible.
    loop1_desc = (E4.end5p, E4.end3p, t1.end3p, t1.end5p, E3.end5p, E3.end3p)
    loop1_path = [end.ifnode.top_delegate() for end in loop1_desc]
    loop1_path_rev = [end.ifnode.top_delegate() for end in reversed(loop1_desc)]
    loop1_activity_before = mgr.calculate_loop_activity(loop1_path) # Without simulating stacking of the two duplexes
    # For loop1, the stacked double-helix segment is shorter than t1, which currently means that the activity is
    # calculated the same with and without stack
    loop1_activity = mgr.calculate_loop_activity(loop1_path, simulate_reaction=STACKING_INTERACTION)
    print("loop1_activity_before:", loop1_activity_before)
    print("loop1_activity:", loop1_activity)
    assert loop1_activity_before == mgr.calculate_loop_activity(loop1_path_rev)
    assert loop1_activity == mgr.calculate_loop_activity(loop1_path_rev, simulate_reaction=STACKING_INTERACTION)

    # The other (bigger) loop formed when stacking e4.end3p+e3.end5p:
    loop2_desc = (e4.end3p, e1.end3p, e1.end5p, t2.end5p, t2.end3p, e2.end3p, e2.end5p, e3.end5p)
    loop2_path = [end.ifnode.top_delegate() for end in loop2_desc]
    loop2_path_rev = [end.ifnode.top_delegate() for end in reversed(loop2_desc)]
    loop2_activity_before = mgr.calculate_loop_activity(loop2_path)
    # Activity as it would be after stacking the duplex ends: (This is what is really relevant)
    # pdb.set_trace()
    loop2_activity = mgr.calculate_loop_activity(loop2_path, simulate_reaction=STACKING_INTERACTION)
    print("loop2_activity_before:", loop2_activity_before)
    print("loop2_activity:", loop2_activity)
    assert loop2_activity_before == mgr.calculate_loop_activity(loop2_path_rev)
    assert loop2_activity == mgr.calculate_loop_activity(loop2_path_rev, simulate_reaction=STACKING_INTERACTION)
    assert loop1_activity_before > loop2_activity_before > loop2_activity
    assert loop2_activity > 0


    ## THE TWO LOOPS FORMED WHEN SPLITTING THE OUTER LOOP BY STACKING E1-3p and E2-5p:
    ## This is essentially case 1(b)
    # The big loop (through t1): e1+t1+e2
    alt_stack_loop1 = (e1.end5p, e1.end3p, t1.end3p, t1.end5p, e2.end5p, e2.end3p)
    alt_stack_loop1_path = [end.ifnode.top_delegate() for end in alt_stack_loop1]
    alt_stack_loop1_activity_before = mgr.calculate_loop_activity(alt_stack_loop1_path)
    alt_stack_loop1_activity = mgr.calculate_loop_activity(alt_stack_loop1_path, simulate_reaction=STACKING_INTERACTION)
    print("alt_stack_loop1_activity_before:", alt_stack_loop1_activity_before)
    print("alt_stack_loop1_activity:", alt_stack_loop1_activity)
    assert alt_stack_loop1_activity_before > 0
    assert alt_stack_loop1_activity == 0

    # The other loop that only consists of t2:
    alt_stack_loop2 = (E1.end3p, t2.end5p, t2.end3p, E2.end5p)
    alt_stack_loop2_path = [end.ifnode.top_delegate() for end in alt_stack_loop2]
    alt_stack_loop2_activity_before = mgr.calculate_loop_activity(alt_stack_loop2_path)
    alt_stack_loop2_activity = mgr.calculate_loop_activity(alt_stack_loop2_path, simulate_reaction=STACKING_INTERACTION)
    print("alt_stack_loop2_activity_before:", alt_stack_loop2_activity_before)
    print("alt_stack_loop2_activity:", alt_stack_loop2_activity)
    assert alt_stack_loop2_activity > 0


def test_loop_formation_effects_1():
    """
    ### Case 2(a):   ###
    ### Using long enough T-loop to connect the outer arm nodes that the stack can form
    """
    e1 = Domain("e1", seq="GCTA"*4)   # 16 bp
    E1 = Domain("E1", seq="TAGC"*4)   # 16 bp
    e2 = Domain("e2", seq="CGTA"*4)   # 16 bp
    E2 = Domain("E2", seq="TACG"*4)   # 16 bp
    e3 = Domain("e3", seq="TA"*3)     #  6 bp
    E3 = Domain("E3", seq="TA"*3)     #  6 bp
    e4 = Domain("e4", seq="G"*5)      #  5 bp
    E4 = Domain("E4", seq="C"*5)      #  5 bp
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*50)     # 30 or larger should be enough to reach
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    # We use ComponentMgr because that has what we need to prepare the complex assembly: (hybridize, etc)
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    cmplx = res['changed_complexes'][0]

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0_desc = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0_desc]
    loop0_path_spec = cmplx.calculate_loop_path_spec(loop0_path)
    loop0_hash = cmplx.calculate_loop_hash(loop0_path_spec)
    loop0_state_hash = cmplx.calculate_loop_hash_slow(loop0_path)
    assert loop0_hash == loop0_state_hash
    # quick check that loop hash is independent on start point:
    loop0_desc_alt = (E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p, t2.end5p, t2.end3p, E2.end5p)
    loop0_path_alt = [end.ifnode.top_delegate() for end in loop0_desc_alt]
    assert loop0_state_hash == cmplx.calculate_loop_hash_slow(loop0_path_alt)
    assert loop0_state_hash == cmplx.calculate_loop_hash(cmplx.calculate_loop_path_spec(loop0_path_alt))
    # quick test that another path gives a wrong hash:
    loop0_desc_wrong = (E2.end3p, t1.end5p, E1.end5p, t1.end3p, E1.end3p, t2.end5p, t2.end3p, E2.end5p)
    loop0_path_wrong = [end.ifnode.top_delegate() for end in loop0_desc_wrong]
    assert loop0_state_hash != cmplx.calculate_loop_hash_slow(loop0_path_wrong)
    assert loop0_state_hash != cmplx.calculate_loop_hash(cmplx.calculate_loop_path_spec(loop0_path_wrong))


    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    print("loop0_activity_before:", loop0_activity_before)

    ## WIP: Manually adding the loop to the complex until that is automated:
    loopid = 1
    loop_info = {'id': loopid,
                 'path': loop0_path,
                 'activity': loop0_activity_before,
                 'ifnodes': set(loop0_path)
                }
    cmplx.loops[loopid] = loop_info
    cmplx.ifnode_by_hash = None
    cmplx.loopid_by_hash = None
    cmplx.rebuild_loopid_by_hash_index()
    cmplx.rebuild_ifnode_loopids_index()


    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (e4.end3p, e3.end5p), (E3.end3p, E4.end5p)
    # stack1_activity = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    activity, loop_effects = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    # assert stack1_activity == activity
    print("stack1_activity:", activity)
    # pdb.set_trace()


def test_loop_formation_effects_2():
    """
    ### Case 2(a)-2:   ###
    ### Here, the T-loop connecting the outer arm nodes is too short for the stack can form
    """
    e1 = Domain("e1", seq="GCTA"*4)   # 16 bp
    E1 = Domain("E1", seq="TAGC"*4)   # 16 bp
    e2 = Domain("e2", seq="CGTA"*4)   # 16 bp
    E2 = Domain("E2", seq="TACG"*4)   # 16 bp
    e3 = Domain("e3", seq="TA"*3)     #  6 bp
    E3 = Domain("E3", seq="TA"*3)     #  6 bp
    e4 = Domain("e4", seq="G"*5)      #  5 bp
    E4 = Domain("E4", seq="C"*5)      #  5 bp
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*22)     # > 22 nt to reach (remember the two pb links on each side of the domain)
    # 25 ss nt or more should be enough to reach (incl two pb links on each side).
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 25*0.6 = 15 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    cmplx = res['changed_complexes'][0]

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0_desc = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0_desc]
    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    # print("loop0_activity_before:", loop0_activity_before)

    ## WIP: Manually adding the loop to the complex until that is automated:
    loopid = 1
    loop_info = {'id': loopid,
                 'path': loop0_path,
                 'activity': loop0_activity_before,
                 'ifnodes': set(loop0_path)
                }
    cmplx.loops[loopid] = loop_info
    cmplx.ifnode_by_hash = None
    cmplx.loopid_by_hash = None
    cmplx.rebuild_loopid_by_hash_index()
    cmplx.rebuild_ifnode_loopids_index()


    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (e4.end3p, e3.end5p), (E3.end3p, E4.end5p)
    stack1_activity, effects = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    print("stack1_activity:", stack1_activity)
    assert stack1_activity == 0



def test_loop_formation_effects_3():
    """
    ### Case 1(a):   ###
    * Stacking duplexes are directly connected by pb link;
    * Edge e3 is zero.
    * *No* t1 loop between the stacking duplex ends of e3 and e4.
    * Using long enough T-loop to connect the outer arm nodes that the stack can form
    """
    e1 = Domain("e1", seq="GCTA"*4)   # 16 bp
    E1 = Domain("E1", seq="TAGC"*4)   # 16 bp
    e2 = Domain("e2", seq="CGTA"*4)   # 16 bp
    E2 = Domain("E2", seq="TACG"*4)   # 16 bp
    e3 = Domain("e3", seq="TA"*3)     #  6 bp
    E3 = Domain("E3", seq="TA"*3)     #  6 bp
    e4 = Domain("e4", seq="G"*5)      #  5 bp
    E4 = Domain("E4", seq="C"*5)      #  5 bp
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*80)     # > 22 nt to reach (remember the two pb links on each side of the domain)
    # activity vs t2 length:
    #   400nt: 1.14
    #   200nt: 0.683
    #   100nt: 0.228
    #   80 nt: 0.128
    #   70 nt: 8.84e-2
    #   50 nt: 6.85e-3
    #   30 nt: 8.82e-4
    #   25 nt: 1.68e-4
    #   23 nt: 7.07e-5
    #   22 nt: 0
    #   10 nt: 0
    # Note: If you have an extra interface-node in your path, that can produce weird results!

    # 25 ss nt outer arm loop or more should be enough to reach (incl two pb links on each side).
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 25*0.6 = 15 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e4, e3]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    cmplx = res['changed_complexes'][0]

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    # loop0 = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)  # case 2(a)
    loop0_desc = (t2.end5p, t2.end3p, E2.end5p, e3.end3p, e3.end5p, E4.end5p, E4.end3p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0_desc]
    assert len(loop0_desc) == len(set(loop0_path))
    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    print("loop0_activity_before:", loop0_activity_before)

    ## WIP: Manually adding the loop to the complex until that is automated:
    loopid = 1
    loop_info = {'id': loopid,
                 'path': loop0_path,
                 'activity': loop0_activity_before,
                 'ifnodes': set(loop0_path)
                }
    cmplx.loops[loopid] = loop_info
    cmplx.ifnode_by_hash = None
    cmplx.loopid_by_hash = None
    cmplx.rebuild_loopid_by_hash_index()
    cmplx.rebuild_ifnode_loopids_index()

    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (e4.end3p, e3.end5p), (E3.end3p, E4.end5p)
    stack1_activity, effects = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    print("stack1_activity:", stack1_activity)
    # assert stack1_activity == 0


def test_loop_formation_effects_4():
    """
    ### Case 1(b):   ###
    ### stacking duplexes are connected with a t1 loop between the stacking duplex ends of e3 and e4.
    ### - e3 = 0
    ### Using long enough T-loop to connect the outer arm nodes that the stack can form
    """
    e1 = Domain("e1", seq="GCTA"*4)   # 16 bp
    E1 = Domain("E1", seq="TAGC"*4)   # 16 bp
    e2 = Domain("e2", seq="CGTA"*4)   # 16 bp
    E2 = Domain("E2", seq="TACG"*4)   # 16 bp
    e3 = Domain("e3", seq="TA"*3)     #  6 bp
    E3 = Domain("E3", seq="TA"*3)     #  6 bp
    e4 = Domain("e4", seq="G"*5)      #  5 bp
    E4 = Domain("E4", seq="C"*5)      #  5 bp
    t1 = Domain("t1", seq="T"*10)
    t2 = Domain("t2", seq="T"*22)     # > 22 nt to reach (remember the two pb links on each side of the domain)
    # activity vs t2 length:
    #   400nt: 8.82e-3
    #   200nt: 5.41e-3
    #   100nt: 1.86e-3
    #   80 nt: 1.06e-3
    #   70 nt: 7.01e-4
    #   50 nt: 1.91e-4
    #   30 nt: 7.44e-6
    #   25 nt: 1.42e-6
    #   23 nt: 5.98e-7
    #   22 nt: 0
    #   10 nt: 0
    # Note: If you have an extra interface-node in your path, that can produce weird results!

    # 25 ss nt outer arm loop or more should be enough to reach (incl two pb links on each side).
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 25*0.6 = 15 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e4, t1, e3]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    cmplx = res['changed_complexes'][0]

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    # loop0 = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)  # case 2(a)
    loop0_desc = (t2.end5p, t2.end3p, E2.end5p, e3.end3p, e3.end5p, t1.end3p, t1.end5p, E4.end5p, E4.end3p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0_desc]
    assert len(loop0_desc) == len(set(loop0_path))
    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    print("loop0_activity_before:", loop0_activity_before)

    ## WIP: Manually adding the loop to the complex until that is automated:
    loopid = 1
    loop_info = {'id': loopid,
                 'path': loop0_path,
                 'activity': loop0_activity_before,
                 'ifnodes': set(loop0_path)
                }
    cmplx.loops[loopid] = loop_info
    for ifnode in loop0_path:
        cmplx.ifnode_loopids_index[ifnode].add(loopid) # is a defaultdict(set)

    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (e4.end3p, e3.end5p), (E3.end3p, E4.end5p)
    stack1_activity, effects = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    print("stack1_activity:", stack1_activity)
    # assert stack1_activity == 0






def test_loop_formation_effects_01():
    """
    ### Case 3a-4xT+bp:   ###
    ### All four T-loops plus direct link of A1-3p and B1-5p
    """
    # A dict where you can write d.key instead of d["key"]. (BUT NOT: d.key = val; d[val]
    d = AttrDict()
    for arm in "AB":
        for i in range(1, 4):
            seq = "".join([choice("ATGC") for _ in range(16)])
            rseq = rcompl(seq)
            #for case in (arm.upper(), arm.lower())
            name1 = '%s%s' % (arm.upper(), i)
            name2 = '%s%s' % (arm.lower(), i)
            d[name1] = Domain(name1, seq=seq)
            d[name2] = Domain(name2, seq=rseq)
            if arm == 'B' and i == 3:
                # Maybe it would have been simpler to extend the arms??
                d[name1+'a'] = Domain(name1+'a', seq=seq[:8])
                d[name2+'a'] = Domain(name2+'a', seq=rseq[8:])
                d[name1+'b'] = Domain(name1+'b', seq=seq[8:])
                d[name2+'b'] = Domain(name2+'b', seq=rseq[:8])

    d['t0'] = Domain("t1", seq="T"*10)
    d['t1'] = Domain("t1", seq="T"*20)
    d['t2'] = Domain("t2", seq="T"*50)     # 30 or larger should be enough to reach
    d['t3'] = Domain("t2", seq="T"*60)     # 30 or larger should be enough to reach
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    sA = Strand("sA", [d[n] for n in "A1 A2 A3".split()])
    sB = Strand("sB", [d[n] for n in "B1 B2 B3".split()])
    #sAB = Strand("sB", [d[n] for n in "A1 A2 A3 B1 B2 B3".split])
    s1 = Strand("s1", [d[n] for n in "b1 t0 a1".split()]) # Make sure domains are ordered 5'->3' !
    s2 = Strand("s2", [d[n] for n in "b2 t1 a2".split()])
    s3 = Strand("s3", [d[n] for n in "b3a t2 b3b t3".split()])
    # s4 = Strand("s4", [d[n] for n in "b1 t0 a1".split()])
    # We use ComponentMgr because that has what we need to prepare the complex assembly:



def test_loop_formation_effects_02():
    """
    ### Case 3b-4xT+bp:   ###
    ### All four T-loops plus direct link of A1-3p and B1-5p.
    ### Same as _01, but making the structure a bit more regular.
    """
    d = AttrDict() # A dict where you can write d.key instead of d["key"].
    for arm in "AB":
        for i in range(4):
            seq = "".join([choice("ATGC") for _ in range(10)])
            rseq = rcompl(seq)
            #for case in (arm.upper(), arm.lower())
            name1 = '%s%s' % (arm.upper(), i)
            name2 = '%s%s' % (arm.lower(), i)
            d[name1] = Domain(name1, seq=seq)
            d[name2] = Domain(name2, seq=rseq)

    for i, N in enumerate([10, 20, 50, 80]):
        name = "t%s"%i
        d[name] = Domain(name, seq="T"*N)
    # d.t0 = Domain("t0", seq="T"*10)
    # d.t1 = Domain("t1", seq="T"*20)
    # d.t2 = Domain("t2", seq="T"*50)     # 30 or larger should be enough to reach
    # d.t3 = Domain("t2", seq="T"*60)     # 30 or larger should be enough to reach
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    sA = Strand("sA", [d[n] for n in "A0 A1 A2 A3".split()])
    sB = Strand("sB", [d[n] for n in "B3 B2 B1 B0".split()])
    #sAB = Strand("sB", [d[n] for n in "A1 A2 A3 B1 B2 B3".split])
    s0 = Strand("s0", [d[n] for n in "a0 t0 b0".split()]) # Make sure domains are ordered 5'->3' !
    s1 = Strand("s1", [d[n] for n in "a1 t1 b1".split()]) # Make sure domains are ordered 5'->3' !
    s2 = Strand("s2", [d[n] for n in "a2 t2 b2".split()])
    s3 = Strand("s3", [d[n] for n in "a3 t3 b3".split()])

    strands = [sA, sB, s0, s1, s2, s3]
    mgr = ComponentMgr(strands, params={}, volume=1e-12)
    cmplx = None
    for arm in "AB":
        for i in range(4):
            name1 = '%s%s' % (arm.upper(), i)
            name2 = '%s%s' % (arm.lower(), i)
            if not (arm > "A" and i > 0):  # arm == "A" or i == 0
                print("Hybridizing %s and %s" % (d[name1], d[name2]))
                res = mgr.hybridize(d[name1], d[name2])       # No loops registered yet.
                print(" - result: %s" % (res,))
                if cmplx is None:
                    cmplx = res['new_complexes'][0]
            else:
                activity, effects = mgr.loop_formation_effects(d[name1], d[name2], reaction_type=HYBRIDIZATION_INTERACTION)
                print("Hybridizing %s and %s" % (d[name1], d[name2]))
                res = mgr.hybridize(d[name1], d[name2])       # No loops registered yet.
                print(" - result:", res)
                print(" - Applying loop_effects:", effects)
                assert cmplx == res['changed_complexes'][0]
                effectuate_loop_changes(cmplx, effects, is_forming=True)
            cmplx.reset_and_recalculate()
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    for arm in "A":
        for i in range(3):
            # Stack A0/3p with A1/5p, a1/3p + a0/5p
            # Edit: Again, we order by duplex end: dh1end3p, dh1end5p, dh2end3p, dh2end5p
            # So it's: A0/3p, a0/5p, a1/3p, A1/5p
            ends = (d['%s%s' % (arm.upper(), i)].end3p, d['%s%s' % (arm.lower(), i)].end5p,
                    d['%s%s' % (arm.lower(), i+1)].end3p, d['%s%s' % (arm.upper(), i+1)].end5p)
            activity, effects = mgr.loop_formation_effects(ends[:2], ends[2:], reaction_type=STACKING_INTERACTION) # No loops registered yet
            res = mgr.stack(*ends)
            assert cmplx == res['changed_complexes'][0]
            effectuate_loop_changes(cmplx, effects, is_forming=True)
            cmplx.reset_and_recalculate()
    for arm in "B":
        for i in range(3):
            # Stack B1/3p, b1/5p, b0/3p, B0/5p
            ends = (d['%s%s' % (arm.upper(), (i+1))].end3p, d['%s%s' % (arm.lower(), (i+1))].end5p,
                    d['%s%s' % (arm.lower(), i)].end3p, d['%s%s' % (arm.upper(), i)].end5p)
            activity, effects = mgr.loop_formation_effects(ends[:2], ends[2:], reaction_type=STACKING_INTERACTION) # No loops registered yet
            res = mgr.stack(*ends)
            assert cmplx == res['changed_complexes'][0]
            effectuate_loop_changes(cmplx, effects, is_forming=True)
            cmplx.reset_and_recalculate()

    # cmplx = res['changed_complexes'][0]
    assert cmplx == next(iter(mgr.complexes)) # mgr.complexes has set of all complexes.
    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_0 = cmplx.ifnode_by_hash.copy()


    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    # l = AttrDict()
    loop_paths = [[end.ifnode.top_delegate() for end in (
        d["a%s"%i].end5p, d["a%s"% i].end3p, d["t%s"%i].end5p, d["t%s" % i].end3p,
        d["b%s"%i].end5p, d["b%s"% i].end3p, d["t%s"%(i+1)].end3p, d["t%s"%(i+1)].end5p)]
             for i in range(3)]
    loop_hashes_before = [cmplx.calculate_loop_hash_slow(loop_path) for loop_path in loop_paths]
    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_1 = cmplx.ifnode_by_hash.copy()
    assert ifnode_by_hash_0 == ifnode_by_hash_1

    # Primary (shortest-path) activities, not considering reaction's effect on other loop_paths:
    loop_activities_before = [mgr.calculate_loop_activity(loop_path) for loop_path in loop_paths]
    print("loop_activities_before:", loop_activities_before)
    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_2 = cmplx.ifnode_by_hash.copy()
    assert ifnode_by_hash_0 == ifnode_by_hash_2

    ## WIP: Manually adding the loop to the complex until that is automated:
    # loopid = 1
    for loopid, loop_path in enumerate(loop_paths):
        loop_info = {'id': loopid,
                     'path': loop_path,
                     'activity': loop_activities_before[loopid],
                    }
        cmplx.loops[loopid] = loop_info
        for ifnode in loop_path:
            cmplx.ifnode_loopids_index[ifnode].add(loopid) # is a defaultdict(set)

    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_3 = cmplx.ifnode_by_hash.copy()
    assert ifnode_by_hash_0 == ifnode_by_hash_3

    # pdb.set_trace()
    cmplx.rebuild_loopid_by_hash_index()
    loopid_by_hash_0 = cmplx.loopid_by_hash.copy()

    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_4 = cmplx.ifnode_by_hash.copy()
    assert ifnode_by_hash_0 == ifnode_by_hash_4

    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (d.a0.end3p, d.A0.end5p), (d.B0.end3p, d.b0.end5p)
    activity, loop_effects = mgr.loop_formation_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)
    # assert stack1_activity == activity
    print("\n\nstack1_activity:", activity)
    print("loop_effects:")
    pprint(loop_effects)
    # pdb.set_trace()

    # loopid_by_hash_1 = cmplx.loopid_by_hash.copy()
    assert cmplx.loopid_by_hash == loopid_by_hash_0 != {}
    cmplx.rebuild_loopid_by_hash_index()
    assert cmplx.loopid_by_hash == loopid_by_hash_0 != {}

    cmplx.rebuild_ifnode_by_hash_index()
    ifnode_by_hash_5 = cmplx.ifnode_by_hash.copy()
    assert ifnode_by_hash_0 == ifnode_by_hash_5

    cmplx.rebuild_loopid_by_hash_index()
    assert cmplx.loopid_by_hash == loopid_by_hash_0 != {}


    ## Note: The STACKING reaction here does not form any new loops, it merely changes existing loops!

    print("Complex.ifnode_by_hash:")
    print(cmplx.ifnode_by_hash)

    # Try to effctuate loop_effects:
    effectuate_loop_changes(cmplx, loop_effects, is_forming=True)

    # Try the reverse reaction:
    loop_effects = mgr.loop_breakage_effects(dh1, dh2, reaction_type=STACKING_INTERACTION)

    # Try to effctuate reverse loop_effects:
    effectuate_loop_changes(cmplx, loop_effects, is_forming=False)

    # Check if loop_path hashes are the same before and after:
    loop_hashes_after = [cmplx.calculate_loop_hash_slow(loop_path) for loop_path in loop_paths]
    assert loop_hashes_before == loop_hashes_after

    # Check if loop_path hashes are the same before and after:
    cmplx.rebuild_loopid_by_hash_index()
    loop_hashes_after = [cmplx.calculate_loop_hash(loop_path) for loop_path in loop_paths]
    assert loop_hashes_before == loop_hashes_after





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

    test_loop_formation_effects_01()
    test_loop_formation_effects_02()

    print("\n\nAll ad-hoc test runs done!\n\n")
