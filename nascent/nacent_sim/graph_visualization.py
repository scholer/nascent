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

Module for visualizing oligos and complexes as graphs.

There are three main types of graph:
    (1) Each strand is a node, connected to other strands via hybridizations.
    (2) Each domain is a node, each strand is a directed subgraph of connected domains (5'->3'),
        and hybridized domains are brought together as clusters.
    (3) As above, but each using two nodes for each domain, one for the 5' end of the domain and another
        node for the 3' end of the domain. That way we can visualize the anti-parallel directionality
        of domain hybridizations.


"""
