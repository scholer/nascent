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

# pylint: disable=C0103

"""
Test graph_manager module.

"""

import pytest
import pdb
import sys
sys.path.insert(0, ".")
try:
    from numpy import isclose
except ImportError:
    def isclose(a, b, rtol=1e-5, atol=1e-8):
        """ Default rtol, atol is actually pretty lenient, only compares to about 5 significant digits..."""
        # How numpy does it:
        return abs(a-b) <= (atol + rtol*abs(b))
        # diff = abs(a-b)
        # return diff < atol and diff/b < rtol


from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.constants import STACKING_INTERACTION




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
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)

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
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)

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
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0 = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0]
    loop0_path_rev = [end.ifnode.top_delegate() for end in reversed(loop0)]
    loop0_activity_before = mgr.calculate_loop_activity(loop0_path)
    loop0_activity_before_rev = mgr.calculate_loop_activity(loop0_path_rev)
    print("loop0_activity_before:", loop0_activity_before)
    print("loop0_activity_before_rev:", loop0_activity_before_rev)
    # print("loop0_activity:", loop0_activity)
    assert isclose(loop0_activity_before, loop0_activity_before_rev, rtol=1e-8, atol=1e-10)

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

    # The other (bigger) loop formed when stacking e4.end3p+e3.end5p:
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
    assert alt_stack_loop1_activity_before > alt_stack_loop1_activity == 0

    # The other loop that only consists of t2:
    alt_stack_loop2 = (E1.end3p, t2.end5p, t2.end3p, E2.end5p)
    alt_stack_loop2_path = [end.ifnode.top_delegate() for end in alt_stack_loop2]
    alt_stack_loop2_activity_before = mgr.calculate_loop_activity(alt_stack_loop2_path)
    alt_stack_loop2_activity = mgr.calculate_loop_activity(alt_stack_loop2_path, simulate_reaction=STACKING_INTERACTION)
    print("alt_stack_loop2_activity_before:", alt_stack_loop2_activity_before)
    print("alt_stack_loop2_activity:", alt_stack_loop2_activity)
    assert alt_stack_loop2_activity > 0


def test_intracomplex_activity_1():

    ### Case 2(a):   ###
    ### Using long enough T-loop to connect the outer arm nodes that
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
    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    mgr.hybridize(E2, e2)
    mgr.hybridize(E3, e3)
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    res = mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)
    cmplx = res['changed_complexes'][0]

    ## ACTIVITY FOR THE OUTER LOOP (calculated as though it was not connected):
    # t2 - E2 - t1 - E1
    loop0 = (t2.end5p, t2.end3p, E2.end5p, E2.end3p, t1.end5p, t1.end3p, E1.end5p, E1.end3p)
    loop0_path = [end.ifnode.top_delegate() for end in loop0]
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
        cmplx.loopids_by_interface[ifnode].add(loopid) # is a defaultdict(set)



    ## Intracomplex activity for forming the left stack:
    dh1, dh2 = (e4.end3p, e3.end5p), (E3.end3p, E4.end5p)
    stack1_activity = mgr.intracomplex_activity(dh1, dh2, reaction_type=STACKING_INTERACTION)
    print("stack1_activity:", stack1_activity)


if __name__ == "__main__":
    ## Test 2:
    test_calculate_loop_activity_2()
    ## Test 3:
    test_calculate_loop_activity_3()


    test_intracomplex_activity_1()
