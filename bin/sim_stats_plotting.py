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

Refs:
* https://github.com/jfrelinger/fcm/blob/master/src/graphics/util.py
** https://pythonhosted.org/fcm/basic.html
* https://bitbucket.org/gorelab/flowcytometrytools/src/


"""

import sys
import os
import argparse
import random
import numpy
import math
from collections import defaultdict
import matplotlib
matplotlib.use("qt4agg")
from matplotlib import pyplot


LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nacent_sim

from nascent.nacent_sim.stats_parser import load_statsfile
from nascent.nacent_sim.dom_utils import parse_strand_domains_file, check_strands
from nascent.plotting import colormaps


def bilinear_interpolate(x, y, bins=None):
    """Returns interpolated density values on points (x, y).
    Ref: http://en.wikipedia.org/wiki/Bilinear_interpolation.
    """
    if bins is None:
        bins = int(numpy.sqrt(len(x)))

    z, unused_xedge, unused_yedge = numpy.histogram2d(
        y, x,
        bins=[bins, bins],
        range=[(numpy.min(y), numpy.max(y)), (numpy.min(x), numpy.max(x))])
    xfrac, xint = numpy.modf((x - numpy.min(x)) /
                             (numpy.max(x) - numpy.min(x)) * (bins - 1))
    yfrac, yint = numpy.modf((y - numpy.min(y)) /
                             (numpy.max(y) - numpy.min(y)) * (bins - 1))

    xint = xint.astype('i')
    yint = yint.astype('i')

    z1 = numpy.zeros(numpy.array(z.shape) + 1)
    z1[:-1, :-1] = z

    # values at corners of square for interpolation
    q11 = z1[yint, xint]
    q12 = z1[yint, xint + 1]
    q21 = z1[yint + 1, xint]
    q22 = z1[yint + 1, xint + 1]

    return q11 * (1 - xfrac) * (1 - yfrac) + q21 * (1 - xfrac) * (yfrac) + \
        q12 * (xfrac) * (1 - yfrac) + q22 * (xfrac) * (yfrac)



def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("statsfile")

    argns = parser.parse_args(argv)
    return argns


if __name__ == '__main__':

    # Plot type:
    plottype = "scatter"    # scatter or regular plot

    # Plotting using only one marker per (x, y) speeds up plotting.
    # If only plotting one point, then scatter color will have no effect, only size.
    one_mark_per_xy = True

    # Truncate data to maximun N number of entries. (Downsampling will be evently spread across the whole dataset.)
    max_lines = 100000
    max_lines = None

    # Marker:
    markershape = "o"   # Generally, use "o" for scatter. For plottype="plot", "," can be nice.
    # Minimum and maximum marker size. Default is 20. Denotes area, not radius/diameter, AFAIK.
    max_marker_size = 100
    min_marker_size = 1
    size_scale = lambda x: x    # size scaling function - e.g. linear or sqrt.
    #size_scale = math.sqrt      # Is normalized linearly between min_ms...max_ms.

    # Color (if not using colormap):
    colors = ["b"]
    use_colormap = True         # Will generate a list of colors to replace the colors above
    colormap_from_sizes = True  # Just use sizes as colormap (does not use histogram)
    # Colormap # magma, plasma, inferno, and viridis
    colormap = colormaps.viridis    # viridis is the blue-green-yellow ("Option D") colormap
                                    # magma, inferno and plasma are rather bright at its top.
    colormap_bins = 10

    # Adding noise to each data point makes it stand out better (if you have one marker per point):
    gaussian_point_noise = False       # If not boolean False, this is also used as the standard deviation
    colormap_after_point_noise = False # It might be appropriate to calculate the colormap after randomization.



    args = parse_args()
    print("Loading stats from file ", args.statsfile)
    stats = load_statsfile(args.statsfile)
    print(len(stats), "lines loaded from file.")


    if max_lines and len(stats) > max_lines:
        sample_freq = int(len(stats)/max_lines)
        stats = [line for i, line in enumerate(stats) if i % sample_freq == 0]
        print("Stats reduced to", len(stats), "lines (sample_freq=%s)." % sample_freq)

    print("Pre-processing data...")
    #[(Temperature, N_doms_hybridized, %_doms_hybridized), ...]
    xdata = [line[0]-273 for line in stats]
    ydata = [line[1] for line in stats]

    # Count data points for each (x, y):
    counts = defaultdict(lambda: defaultdict(int))  # A defaultdict with defaultdict(int) default entries
    for x, y in zip(xdata, ydata):
        counts[x][y] += 1
    maxcount = max(c for xdicts in counts.values() for c in xdicts.values())
    print("Max (T, N) count:", maxcount)

    # If reducing to one marker per (x, y) point, this must be done now, before making the sizes array.
    if one_mark_per_xy:
        # Make a set of (x, y) points and plot those with size.
        xdata, ydata = zip(*set((x, y) for x, y in zip(xdata, ydata)))
        print("Only plotting one marker per (x, y) combination; dataset reduced to", len(xdata), "points.")

    # Make marker sizes array:
    sizes = [int(max_marker_size*size_scale(counts[x][y]/maxcount))+1 for x, y in zip(xdata, ydata)]
    print("Average marker size: %0.1f" % (sum(sizes)/len(sizes)))

    # Interpolating colors:
    if use_colormap:
        print("Making colors by interpolations")
        # Is an array of colors, one color per
        # I think in theory this could also be same as sizes...?
        if colormap_from_sizes:
            colors = sizes
        else:
            colors = bilinear_interpolate(xdata, ydata, bins=colormap_bins)
            print("colors array of length:", len(colors))


    if gaussian_point_noise:
        # Making random noise around each (x, y) can give a nice effect (mostly for plot, not scatter).
        #xdata = [line[0]-273 + random.gauss(0, 0.1) for line in stats]
        #ydata = [line[1] + random.gauss(0, 0.1) for line in stats]

        #xdata = [line[0]-273 + round(random.gauss(0, 0.1), 2)  for line in stats]
        #ydata = [line[1] + round(random.gauss(0, 0.1), 2) for line in stats]

        xdata = [x + random.gauss(0, 0.1) for x in xdata]
        ydata = [y + random.gauss(0, 0.1) for y in ydata]
        if use_colormap and colormap_after_point_noise:
            colors = bilinear_interpolate(xdata, ydata, bins=3)


    print("Plotting data...")
    if plottype == "scatter":
        # Using a scatter plot is good for density-based data (e.g. flow cytometry):
        pyplot.scatter(xdata, ydata,
                       marker='o',
                       s=sizes,
                       c=colors,
                       edgecolors='none',
                       cmap=colormaps.viridis)
    else:
        # Regular "one-marker-per-datapoint" plot
        # Different sizes are not supported (you can call plot repeatedly or pass in multiple sets of
        # xdata, ydata, **specs, xdata2, ydata2, **specs -- but that is not really suitable for us.
        pyplot.plot(xdata, ydata, marker=markershape, ms=min_marker_size)

    print("Showing plot...")
    pyplot.show()
    print("Done!")



