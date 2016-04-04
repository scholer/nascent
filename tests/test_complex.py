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
from pprint import pprint
import pdb
from collections import Counter
import sys

sys.path.insert(0, ".")  # Required to import nascent
#sys.path.insert(0, ".")
from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager
from nascent.graph_sim_nx.constants import STACKING_INTERACTION

# try:
#     from numpy import isclose
# except ImportError:
#     def isclose(a, b, rtol=1e-5, atol=1e-8):
#         """ Default rtol, atol is actually pretty lenient, only compares to about 5 significant digits..."""
#         # How numpy does it:
#         return abs(a-b) <= (atol + rtol*abs(b))




def test_asymmetric_complex_01():
    """
    Simple case, no duplicate strands or domains, complex is obviously asymmetric.
    Structure is case 2(a) from the graph_manager test module.
    """

    ### Case 2(a) from test_graphmgr:   ###
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
    result = mgr.hybridize(E1, e1)
    cmplx = s1.complex
    cstates, dspecs, dspec_counts = [], [], []   # 1 entry for each expected state

    def update_and_assert():
        cstates.append(cmplx.state_fingerprint())
        domain_states = {d: d.state_fingerprint() for d in cmplx.domains()}
        dspecs.append(frozenset(domain_states.items()))
        # dspec_counts should be the same for different instances in the same state:
        dspec_counts.append(Counter(domain_states.values()))
        assert dspec_counts[-1].most_common(1)[0][1] == 1  # Count (value) of most common dspec hash is 1

    update_and_assert()
    result = mgr.hybridize(E2, e2)
    update_and_assert()
    result = mgr.hybridize(E3, e3)
    update_and_assert()
    result = mgr.hybridize(E4, e4)
    update_and_assert()
    # h1end3p, h2end5p, h2end3p, h1end5p   aka   dh1end3p, dh1end5p, dh2end3p, dh2end5p
    result = mgr.stack((E4.end3p, e4.end5p), (e1.end3p, E1.end5p))
    update_and_assert()
    result = mgr.stack((E2.end3p, e2.end5p), (e3.end3p, E3.end5p))
    update_and_assert()

    assert Counter(cstates).most_common(1)[0][1] == 1
    pprint(cmplx)
    pprint(cstates)
    pprint(dspecs)
    pprint(dspec_counts)
    # pdb.set_trace()



def test_symmetric_complex_01():
    """
    Simple symmetric case:
         A ----- B ----- a
         |       |       |
         a ----- B ----- A
    B is self-complementary and palindromic.
    """

    ### Case sym_duplex_01:   ###
    A1 = Domain("A", seq="GCTA"*4)   # 16 bp
    a1 = Domain("a", seq="TAGC"*4)   # 16 bp
    B1 = Domain("B", seq="ATGCAT"*4)  # 24 bp - palindromic
    A2 = Domain("A", seq="GCTA"*4)   # 16 bp
    a2 = Domain("a", seq="TAGC"*4)   # 16 bp
    B2 = Domain("B", seq="ATGCAT"*4)  # 24 bp - palindromic


    s0 = Strand("s", [A1, B1, a1])
    s1 = Strand("s", [A2, B2, a2]) # Make sure domains are ordered 5'->3' !

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s0, s1], params={})
    print("\n\nmgr.hybridize(B1, B2)")
    result = mgr.hybridize(B1, B2)
    cmplx = s1.complex
    assert cmplx in result["new_complexes"]
    cstates, dspecs, dspec_counts = [], [], []   # 1 entry for each expected state

    def update_and_assert(match_idx=None):
        cmplx.reset_state_fingerprint()
        cstate = cmplx.state_fingerprint()
        domain_states = {d: d.state_fingerprint() for d in cmplx.domains()}
        domain_state_count = Counter(domain_states.values())
        if match_idx is None:
            cstates.append(cstate)
            dspecs.append(frozenset(domain_states.items()))
            dspec_counts.append(domain_state_count)
        else:
            assert cstate == cstates[match_idx]
            assert frozenset(domain_states.items()) == dspecs[match_idx]
            assert domain_state_count == dspec_counts[match_idx]
        assert dspec_counts[-1].most_common(1)[0][1] == 1  # Count (value) of most common dspec hash is 1

    update_and_assert()
    print("\n\nmgr.hybridize(A1, a2)")
    result = mgr.hybridize(A1, a2)
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.hybridize(a1, A2)")
    result = mgr.hybridize(a1, A2)
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.stack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))")
    result = mgr.stack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.stack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))")
    result = mgr.stack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\nComplex test-case test_symmetric_complex_01 completed.")
    assert Counter(cstates).most_common(1)[0][1] == 1

    print("Complex:")
    pprint(cmplx)
    print("Complex state fingerprints (collected during assembly):")
    pprint(cstates)
    print("Complex domain state fingerprints (collected during assembly):")
    pprint(dspecs)
    print("Complex domain state fingerprint counts (collected during assembly):")
    pprint(dspec_counts)
    # pdb.set_trace()

    print("\n\nDe-constructing, same order as constructing...")
    i = -1  # Current fingerprints should match
    print("\n\nmgr.unstack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))")
    result = mgr.unstack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.unstack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))")
    result = mgr.unstack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(a1, A2)")
    result = mgr.dehybridize(a1, A2)
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(A1, a2)")
    result = mgr.dehybridize(A1, a2)
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(B1, B2)")
    result = mgr.dehybridize(B1, B2)
    # State should now be obsolete
    assert cmplx in result["obsolete_complexes"]



def test_symmetric_complex_02():
    """
    Complex is asymmetric but only after C has hybridized to one of the c domains.

         A ----- B ----- a ---- c = C
         |       |       |
    c----a ----- B ----- A

    """

    ### Case asym_duplex_01a:   ###
    A1 = Domain("A", seq="GCTA"*4)   # 16 bp
    a1 = Domain("a", seq="TAGC"*4)   # 16 bp
    B1 = Domain("B", seq="ATGCAT"*4)  # 24 bp - palindromic
    A2 = Domain("A", seq="GCTA"*4)   # 16 bp
    a2 = Domain("a", seq="TAGC"*4)   # 16 bp
    B2 = Domain("B", seq="ATGCAT"*4)  # 24 bp - palindromic
    C1 = Domain("C", seq="TAGT"*4)   # 16 bp
    c1 = Domain("c", seq="ACTA"*4)   # 16 bp
    c2 = Domain("c", seq="ACTA"*4)   # 16 bp

    s0 = Strand("s", [A1, B1, a1, c1])
    s1 = Strand("s", [A2, B2, a2, c2]) # Make sure domains are ordered 5'->3' !
    sC = Strand("s", [C1])

    # We use ComponentMgr because that has what we need to prepare the complex assembly:
    mgr = ComponentMgr([s0, s1, sC], params={})
    print("\n\nmgr.hybridize(B1, B2)")
    result = mgr.hybridize(B1, B2)
    cmplx = s1.complex
    assert cmplx in result["new_complexes"]
    cstates, dspecs, dspec_counts = [], [], []   # 1 entry for each expected state

    def update_and_assert(match_idx=None):
        cmplx.reset_state_fingerprint()
        cstate = cmplx.state_fingerprint()
        domain_states = {d: d.state_fingerprint() for d in cmplx.domains()}
        domain_state_count = Counter(domain_states.values())
        if match_idx is None:
            cstates.append(cstate)
            dspecs.append(frozenset(domain_states.items()))
            dspec_counts.append(domain_state_count)
        else:
            assert cstate == cstates[i]
            assert frozenset(domain_states.items()) == dspecs[i]
            assert domain_state_count == dspec_counts[i]
        assert dspec_counts[-1].most_common(1)[0][1] == 1  # Count (value) of most common dspec hash is 1

    update_and_assert()
    print("\n\nmgr.hybridize(A1, a2)")
    result = mgr.hybridize(A1, a2)
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.hybridize(a1, A2)")
    result = mgr.hybridize(a1, A2)
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.stack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))")
    result = mgr.stack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.stack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))")
    result = mgr.stack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\n\nmgr.hybridize(c1, C1)")
    result = mgr.hybridize(c1, C1)
    assert cmplx in result["changed_complexes"]
    update_and_assert()

    print("\nComplex test-case test_symmetric_complex_01 completed.")
    assert Counter(cstates).most_common(1)[0][1] == 1

    print("Complex:")
    pprint(cmplx)
    print("Complex state fingerprints (collected during assembly):")
    pprint(cstates)
    print("Complex domain state fingerprints (collected during assembly):")
    pprint(dspecs)
    print("Complex domain state fingerprint counts (collected during assembly):")
    pprint(dspec_counts)
    # pdb.set_trace()

    print("\n\nDe-constructing, same order as constructing...")
    i = -1  # Current fingerprints should match

    print("\n\nmgr.dehybridize(c1, C1)")
    result = mgr.dehybridize(c1, C1)
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.unstack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))")
    result = mgr.unstack((A2.end3p, a1.end5p), (B1.end3p, B2.end5p))
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.unstack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))")
    result = mgr.unstack((A1.end3p, a2.end5p), (B2.end3p, B1.end5p))
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(a1, A2)")
    result = mgr.dehybridize(a1, A2)
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(A1, a2)")
    result = mgr.dehybridize(A1, a2)
    assert cmplx in result["changed_complexes"]
    i -= 1
    update_and_assert(i)

    print("\n\nmgr.dehybridize(B1, B2)")
    result = mgr.dehybridize(B1, B2)
    # State should now be obsolete
    assert cmplx in result["obsolete_complexes"]



if __name__ == "__main__":
    print("NEW TEST STARTED\n"*30)

    ## Asymmetric complexes:
    # test_asymmetric_complex_01()

    ## Symmetric complexes:
    test_symmetric_complex_01()

    test_symmetric_complex_02()
