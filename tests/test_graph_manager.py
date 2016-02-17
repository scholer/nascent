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

import sys
sys.path.insert(0, ".")

from nascent.graph_sim_nx.domain import Domain
from nascent.graph_sim_nx.strand import Strand
from nascent.graph_sim_nx.graph_manager import GraphManager
from nascent.graph_sim_nx.componentmgr import ComponentMgr
from nascent.graph_sim_nx import graph_manager



def test_calculate_loop_activity():
    ## Case 2(a):
    e1 = Domain("e1", seq="GCTA"*4)
    E1 = Domain("E1", seq="TAGC"*4)
    e2 = Domain("e2", seq="CGTA"*4)
    E2 = Domain("E2", seq="TACG"*4)
    e3 = Domain("e3", seq="TA"*3)
    E3 = Domain("E3", seq="TA"*3)
    e4 = Domain("e4", seq="G"*5)
    E4 = Domain("E4", seq="C"*5)
    t1 = Domain("t1", seq="T"*15)
    t2 = Domain("t1", seq="T"*15)
    s1 = Strand("s1", [E4, E1, t2, E2, E3])
    s2 = Strand("s2", [e4, t1, e3])
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

    #
    loop1_path = [end.ifnode.top_delegate() for end in (E4.end5p, E4.end3p, E3.end5p, E3.end3p)]
    print(mgr.calculate_loop_activity(loop1_path))






if __name__ == "__main__":
    test_calculate_loop_activity()

