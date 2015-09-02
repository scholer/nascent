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

Simple module for loading simulation stats back into a python data structure:
    [(Temperature, N_doms_hybridized, %_doms_hybridized), ...]


"""


def load_statsfile(statsfile):
    with open(statsfile) as fp:
        lines = (line.strip() for line in fp if line.strip())
        #data = [[float(val.strip()) for val in  cell]
        #        for line in fp for cell in line.strip().split(",") if line.strip()]
        data = [[float(cell.strip()) for cell in line.split(",")] for line in lines]

    return data


