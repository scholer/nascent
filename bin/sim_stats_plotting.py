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

Plot stats from dom_anneal_sim.

"""

import sys
import os
import argparse
import random

import matplotlib
matplotlib.use("qt4agg")
from matplotlib import pyplot


LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nacent_sim

from nascent.nacent_sim.stats_parser import load_statsfile
from nascent.nacent_sim.dom_utils import parse_strand_domains_file, check_strands



def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("statsfile")

    argns = parser.parse_args(argv)
    return argns


if __name__ == '__main__':
    args = parse_args()
    print("Loading stats from file ", args.statsfile)
    stats = load_statsfile(args.statsfile)
    #[(Temperature, N_doms_hybridized, %_doms_hybridized), ...]
    #for line in stats:
    #    line.append(line[0]-273)

    print("Pre-processing data...")
    xdata = [line[0]-273 + random.gauss(0, 0.05) for line in stats]
    ydata = [line[1] + random.gauss(0, 0.05) for line in stats]

    print("Plotting data...")
    pyplot.plot(xdata, ydata, ",")
    pyplot.show()



