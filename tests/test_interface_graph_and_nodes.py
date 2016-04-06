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
from pprint import pprint

if "." not in sys.path:
    sys.path.insert(0, ".")  # Required to import nascent
from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.constants import STACKING_INTERACTION
from nascent.graph_sim_nx.system_graphs import InterfaceMultiGraph, InterfaceGraph, InterfaceNode


class DomainEndMock(object):
    """
    """
    def __init__(self, name):
        self.name = str(name)

    def __str__(self):
        return str(self.name)

    def state_fingerprint(self):
        return self.name


def check_multi_graph(g):
    for source in g.adj:
        assert source in g.node
        assert isinstance(g.adj[source], dict)
        for target in g.adj[source]:
            assert target in g.node
            assert isinstance(g.adj[source][target], dict)
            assert source in g.adj[target]
            assert isinstance(g.adj[target][source], dict)
            assert g.adj[target][source] is g.adj[source][target]
            for key in g.adj[source][target]:
                assert isinstance(g.adj[source][target][key], dict)



def test_interfacemultigraph_merge_split_1():
    """
    Test basic InterfaceMultiGraph merge(...) and split(...) functionality

    Test case 1, from InterfaceNode docstring:
    For instance, consider the graph
        0------1---2-------3     0---.           .----3     0---.  1        .---3
                                      `.1---2---:                `.----2---:
        4------5---6-------7     4----´ 5   6    `---7      4----´ 5   6    `---7
    Where 5 was delegated to 1, 6 delegated to 2, and then 1 delegated to 2.
    I.e from 5 we have a tree delegation branch: 5 --> 1 --> 2
    While from 2 down we have the tree:      .- 6
                                       2 <-<ˊ
                                            `·- 1 <-- 5

    """
    ## Setup: ##
    g = InterfaceMultiGraph()
    node = ifnodes = [InterfaceNode(DomainEndMock(name)) for name in range(8)]
    g.add_path(ifnodes[:4])
    g.add_path(ifnodes[5:])

    print("Graph.adj:")
    pprint(g.adj)

    # I no longer overwrite native Networkx.MultiGraph methods (add_edge and friends).
    # Instead, once you have added what you need, invoke this to set ifnode.delegated_edges to match graph:
    g.reset_ifnode_delegation_to_current_graph_representation()

    check_multi_graph(g)


    # Merge 5->1, 6->2, then 2->1:
    for delegator, delegatee in [(ifnodes[i], ifnodes[j]) for i, j in [(5, 1), (6, 2), (2, 1)]]:

        print("\nMerging %s, %s" % (delegator, delegatee))
        g.merge(delegator, delegatee) #

        assert delegator.delegatee is delegatee
        assert delegator.top_delegate() is delegatee
        assert delegatee.delegatee is None
        assert delegatee.top_delegate() is delegatee

        # An ifnode should have itself in self.delegated_edges: (ALWAYS, regardless of delegation)
        assert delegator in delegator.delegated_edges
        assert delegatee in delegatee.delegated_edges

        # Delegation of edges from delegator to delegatee:
        assert delegator in delegatee.delegated_edges
        assert delegatee not in delegator.delegated_edges

        check_multi_graph(g)
        print(" - %s, %s was merged OK!\n" % (delegator, delegatee))



    print("\n\ng.adj after all merges:")
    pprint(g.adj)
    print("(ifnode, ifnode.delegated_edges) for ifnode in ifnodes after all merges:")
    pprint(sorted([(ifnode, ifnode.delegated_edges) for ifnode in ifnodes]))


    # Split 1<-2, 1<-5, 2<-6:
    # (delegator vs delegatee should be determined automatically for split - unlike for undelegate() )
    for delegator, delegatee in [(ifnodes[i], ifnodes[j]) for i, j in reversed([(5, 1), (6, 2), (2, 1)])]:
        print("\nSplitting %s, %s" % (delegator, delegatee))

        g.split(delegator, delegatee) #

        assert delegator.delegatee is None
        assert delegator.top_delegate() is delegator
        assert delegatee.delegatee is None
        assert delegatee.top_delegate() is delegatee

        # An ifnode should have itself in self.delegated_edges: (ALWAYS, regardless of delegation)
        assert delegator in delegator.delegated_edges
        assert delegatee in delegatee.delegated_edges

        # Delegation of edges from delegator to delegatee:
        assert delegator not in delegatee.delegated_edges
        assert delegatee not in delegator.delegated_edges

        check_multi_graph(g)
        print(" - %s, %s was split OK!\n" % (delegator, delegatee))

    print("\n\ng.adj after all splits:")
    pprint(g.adj)
    print("ifnode.delegated_edges after all splits:")
    pprint(sorted([(ifnode, ifnode.delegated_edges) for ifnode in ifnodes]))



def test_interfacemultigraph_delegate_undelegate_1():
    """
    Test basic InterfaceMultiGraph delegate(...) and undelegate(...) functionality
    """
    pass





def test_interfacegraph_merge_split_1():
    """
    Test basic InterfaceGraph merge(...) and split(...) functionality
    """
    pass


def test_interfacegraph_delegate_undelegate_1():
    """
    Test basic InterfaceGraph delegate(...) and undelegate(...) functionality
    """
    pass


def test_ifnode_state_fingerprint_1():
    """
    Test to ensure that ifnode state fingerprints are invariant between instances.
    Using ComponentMgr to provide setup. I know this is not good unit testing. - So consider it integration testing :)
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
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
    mgr.hybridize(E1, e1)
    update_and_assert()
    mgr.hybridize(E2, e2)
    update_and_assert()
    mgr.hybridize(E3, e3)
    update_and_assert()
    mgr.hybridize(E4, e4)
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    update_and_assert()
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
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
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
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
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    update_and_assert(i)
    i += 1
    res = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
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
    mgr = ComponentMgr([s1, s2, s3, s4], params={}, volume=1e-12)
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
    mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    update_and_assert(i, assert_true=False)  # Stacking switched around
    i += 1
    mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
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
    # for ifnode in loop0_path:
    #     cmplx.ifnode_loopids_index[ifnode].add(loopid) # is a defaultdict(set)
    cmplx.rebuild_ifnode_by_hash_index()
    cmplx.rebuild_ifnode_loopids_index()







if __name__ == "__main__":
    print("NEW AD-HOC TEST-RUN STARTED\n"*10)


    # test_ifnode_state_fingerprint_1()

    test_interfacemultigraph_merge_split_1()

    print("\n\nAD-HOC TEST-RUN COMPLETE!\n"*10)

