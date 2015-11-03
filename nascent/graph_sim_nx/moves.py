#!/usr/bin/env python
# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# pylint: disable=W0142,C0103,C0301,W0141

"""

Module with functional moves for each type of event,
e.g. hybridize or dehybridize.

"""


#import random
#import math
#from collections import deque, OrderedDict
#from itertools import zip_longest, chain#, accumulate
from collections import deque, OrderedDict
from itertools import zip_longest, chain#, accumulate


from .dom_anneal_models import Complex
from .domain import distance_cached


N_AVOGADRO = 6.022e23






#def merge_complexes(c_major, c_minor, new_connections=None):
#    """
#    Merge c_minor complex into c_major complex.
#    """
#
#    assert c_major != c_minor
#    #if other_complex == self:
#    #    return
#    for strand in c_minor.strands:
#        strand.complex = c_major
#    #self.strands += other_complex.strands
#    c_major.strands |= c_minor.strands
#    c_major.N_strand_changes += c_minor.N_strand_changes
#    #c_major.strands_changes.append(("merge", str(other_complex.strands_changes), len(c_major.strands)))
#    #c_major.strands_history.append(str(sorted([str(s) for s in c_major.strands])))
#    c_major.connections |= c_minor.connections
#    # | is the "union" operator for sets, & is intersection and ^ is symdiff
#    if new_connections:
#        for domain1, domain2 in new_connections:
#            c_major.add_connection(domain1, domain2)
#    # del other_complex   # Deleting in this namespace will have no effect.
