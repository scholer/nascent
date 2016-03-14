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

# pylint: disable=C0103,C0301,W0511,W0611,W0142

## TODO: Enable all pylint checks

"""

Module for testing InterfaceGraph and InterfaceNode classes and related functions.

"""

import pytest
from pprint import pprint
import pdb
from collections import Counter
import sys

if "." not in sys.path:
    sys.path.insert(0, ".")  # Required to import nascent
from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.constants import STACKING_INTERACTION




def test_ifnode_state_fingerprint_1():
    """
    Test to ensure that ifnode state fingerprints are invariant between instances.
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
    domains = [e1, E1, e2, E2, e3, E3, e4, E4, t1, t2]
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    strands = [s1, s2, s3, s4]
    assert set(domains) == set(d for s in strands for d in s.domains)
    domain_ends = [end for d in domains for end in (d.end5p, d.end3p)]

    state_ifnode_fingerprints, state_ifnode_spec_counts = [], []   # 1 dict with ifnodes and fingerprints for each state
    def update_and_assert(match_idx=None, assert_true=True):
        # ifnode_fingerprints.append(cmplx.state_fingerprint())
        top_ifnodes = {end.ifnode for end in domain_ends if end.ifnode.delegatee is None}
        assert top_ifnodes == {end.ifnode.top_delegate() for end in domain_ends}
        ifnode_fingerprints = {ifnode: ifnode.state_fingerprint() for ifnode in top_ifnodes}
        ifnode_spec_counts = Counter(ifnode_fingerprints.values())
        if match_idx is None:
            state_ifnode_fingerprints.append(ifnode_fingerprints)
            # The top_delegates are unlikely to be the same between instances/runs, but
            # dspec_counts should be the same for different instances in the same state:
            state_ifnode_spec_counts.append(ifnode_spec_counts)
        else:
            # assert ifnode_fingerprints == state_ifnode_fingerprints[match_idx]
            assert (ifnode_spec_counts == state_ifnode_spec_counts[match_idx]) is assert_true

        assert state_ifnode_spec_counts[-1].most_common(1)[0][1] == 1  # Count (value) of most common dspec hash is 1

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    update_and_assert()
    mgr.hybridize(E2, e2)
    update_and_assert()
    mgr.hybridize(E3, e3)
    update_and_assert()
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    update_and_assert()
    res = mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)
    update_and_assert()

    # print("state_ifnode_spec_counts:")
    # pprint(state_ifnode_spec_counts)


    ## Delete all and make new instances:
    del mgr
    del res, domains, strands, domain_ends
    del e1, E1, e2, E2, e3, E3, e4, E4, t1, t2
    del s1, s2, s3, s4
    print("Remaining locals:", locals().keys())

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
    domains = [e1, E1, e2, E2, e3, E3, e4, E4, t1, t2]
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    strands = [s1, s2, s3, s4]
    assert set(domains) == set(d for s in strands for d in s.domains)
    domain_ends = [end for d in domains for end in (d.end5p, d.end3p)]

    ## Re-build, checking that we get the same ifnode state fingerprints:
    i = 0
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    update_and_assert(i)
    i += 1
    mgr.hybridize(E2, e2)
    update_and_assert(i)
    i += 1
    mgr.hybridize(E3, e3)
    update_and_assert(i)
    i += 1
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    update_and_assert(i)
    i += 1
    res = mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)
    update_and_assert(i)
    i += 1


    ## Repeat, but this time alter the process:
    del mgr
    del res, domains, strands, domain_ends
    del e1, E1, e2, E2, e3, E3, e4, E4, t1, t2
    del s1, s2, s3, s4
    print("Remaining locals:", locals().keys())

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
    domains = [e1, E1, e2, E2, e3, E3, e4, E4, t1, t2]
    # Stacked double-helix is 16+6+5+16 bp = 43 bp, 0.34nm/bp*43bp = 14.6 nm.
    # t2 contour length is 50*0.6 = 30 nm; t2 dist_ee_sq is 10*0.6*1.8 = 10.8 nm², dist_ee_nm = 3.3 nm.
    # It can reach, but only by stretching.
    # dist_ee_sq = n_nt * ss_rise_per_nt * ss_kuhn_length
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e3, t1, e4]) # Make sure domains are ordered 5'->3' !
    s3 = Strand("s3", [e1])
    s4 = Strand("s4", [e2])
    strands = [s1, s2, s3, s4]
    assert set(domains) == set(d for s in strands for d in s.domains)
    domain_ends = [end for d in domains for end in (d.end5p, d.end3p)]

    ## Re-build, checking that we get the same ifnode state fingerprints:
    i = 0
    mgr = ComponentMgr([s1, s2, s3, s4], params={})
    mgr.hybridize(E1, e1)
    update_and_assert(i)
    i += 1
    mgr.hybridize(E3, e3)
    update_and_assert(i, assert_true=False)  # E2+e2 and E3+e3 switched around
    i += 1
    mgr.hybridize(E2, e2)
    update_and_assert(i, assert_true=True)  # Should now give the same as before.
    i += 1
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack(E2.end3p, e2.end5p, e3.end3p, E3.end5p)
    update_and_assert(i, assert_true=False)  # Stacking switched around
    i += 1
    mgr.stack(E4.end3p, e4.end5p, e1.end3p, E1.end5p)
    update_and_assert(i)
    i += 1



def test_intracomplex_activity_1():
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



if __name__ == "__main__":
    print("NEW TEST STARTED\n"*30)


    test_ifnode_state_fingerprint_1()


