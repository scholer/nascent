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
import yaml
from collections import defaultdict
import matplotlib
matplotlib.use("qt4agg")
from matplotlib import pyplot


LIBPATH = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, LIBPATH)

import nascent.nascent_sim
from nascent.nascent_sim.stats_parser import load_statsfile
from nascent.nascent_sim.dom_utils import parse_strand_domains_file, check_strands
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
    """ Parse command line arguments. """
    parser = argparse.ArgumentParser(description="Plot dataset description", prog="sim_stats_plotting")
    parser.add_argument("statsfiles", nargs="+")

    parser.add_argument("--plottype", default="scatter", nargs="?")

    # The field to plot on the y-axis:
    parser.add_argument("--plotfieldx",
                        help="The field to plot on the x-axis.")
    parser.add_argument("--plotfieldy",
                        help="The field to plot on the y-axis.")
    parser.add_argument("--normalize", action="store_true", default=None,
                        help="Normalize the y-data values.")
    # The argument given is
    parser.add_argument("--normalize_by_field",
                        help="The field use to normalize with.")

    parser.add_argument("--legend", action="store_true", default=None,
                        help="Add a legend to the figure.")

    ## defaut config
    parser.add_argument("--default_config",
                        help="Load default config settings from this file.")
    # config
    parser.add_argument("--config",
                        help="Load config settings from this file. This will override command line argument defaults.")
    ## one_mark_per_xy
    parser.add_argument("--one_mark_per_xy", action="store_true", default=None,
                        help="Reduce the dataset to (x, y) = <count>. This will greatly improve plotting performance.")
    #
    ## show plot
    parser.add_argument("--showplot", action="store_true", default=None,
                        help="Show the plot figure after plotting.")
    parser.add_argument("plot_average", action="store_true", default=None)
    ## save plot
    parser.add_argument("--saveplot", action="store_true", default=None,
                        help="Save the plot figure after plotting.")
    ## Interact with the plotting environment showing the plot:
    parser.add_argument("--interactive", action="store_true", default=None,
                        help="Interact with the plotting environment showing the plot.")
    #
    ## colors
    parser.add_argument("--colors")#, default="bgrcmyk")
    parser.add_argument("--use_colormap", default="auto",
                        help="Use colormap when plotting. Default value auto will use colormap if \
                             plotting only one dataset, but not use colormaps if plotting multiple datasets.")
    parser.add_argument("--colormap_base", nargs="?",
                        help="Color each marker according to it's size (when use_colormap is activated). \
                             If not, calculate by binnning.")
    parser.add_argument("--colormap", default="viridis",
                        help="The colormap to use (to mark e.g. frequency count).")
    parser.add_argument("--colormap_bins", default=10, type=int,
                        help="Use this bin size when calculating colormap.")

    ## max_data_lines
    parser.add_argument("--max_data_lines", type=int,
                        help="Trim the dataset to have at most N number of lines. \
                            Trimming is done symmetrically throughout the dataset.")
    parser.add_argument("--markershape", default="o",
                        help="The marker shape to use for each point when plotting.")
    parser.add_argument("--max_marker_size", default=100, type=int,
                        help="Maximum marker size - when the size is proportional to the count of (x, y).")
    parser.add_argument("--min_marker_size", default=1, type=int,
                        help="Minimum marker size - when the size is proportional to the count of (x, y).")
    parser.add_argument("--size_scaling_method", default="linear",
                        help="How to scale size by the (x, y) count. Accepted values are 'linear' and 'sqrt'.")

    ## Adding noise to each data point makes it stand out better (if you have one marker per point):
    parser.add_argument("--gaussian_point_noise", type=int, default=0,
                        help="Add a random gaussian noise to each point (to give a nice effect). \
                             If not boolean False, this is also used as the standard deviation.")
    parser.add_argument("--colormap_after_point_noise", action="store_true",
                        help="Calculate the colormap after applying noise. It might be appropriate \
                             to calculate the colormap after randomization, but usually it is not.")

    return parser.parse_args(argv)


def process_args(args):
    """ Process arguments (load default and overriding configs). """

    print("Command line args:")
    print("\n".join("%s: %s" % kv for kv in sorted(args.items())))

    if args.get('default_config'):
        default_config = yaml.load(open(args['default_config']))
    else:
        default_config = dict(
            # Plot type:
            plottype="scatter",    # scatter or regular plot
            plotfieldx=0,
            plotfieldy=1,
            normalize=True,
            normalize_by_field=None,
            plot_melting_profile=True,
            saveplot=True,
            showplot=True,
            interactive=False,
            plotfilename=None,
            legend=None,

            # Plotting using only one marker per (x, y) speeds up plotting.
            # If only plotting one point, then scatter color will have no effect, only size.
            one_mark_per_xy=True,

            # Truncate data to maximun N number of entries. (Downsampling will be evently spread across the whole dataset.)
            #max_lines = 100000,
            max_lines=None,

            # Marker:
            markershape="o",   # Generally, use "o" for scatter. For plottype="plot", "," can be nice.
            # Minimum and maximum marker size. Default is 20. Denotes area, not radius/diameter, AFAIK.
            max_marker_size=100,
            min_marker_size=1,

            size_scaling_fun=lambda x: x,    # size scaling function - e.g. linear or sqrt.
            #size_scale = math.sqrt      # Is normalized linearly between min_ms...max_ms.

            # Color (if not using colormap):
            #colors=["b"],
            colors=None, #"brgcmyk",
            use_colormap='auto',         # Will generate a list of colors to replace the colors above
            colormap_base='sizes',  # Just use sizes as colormap (does not use histogram)
            # Colormap # magma, plasma, inferno, and viridis
            colormap=colormaps.viridis,
            # viridis is the blue-green-yellow ("Option D") colormap
            # magma, inferno and plasma are rather bright at its top.
            colormap_bins=10,
            gaussian_point_noise=0,
            colormap_after_point_noise=False
            )
    for k, v in default_config.items():
        if args.get(k) is None:
            args[k] = v

    if args['config']:
        args.update(yaml.load(open(args['config'])))

    if isinstance(args['colormap'], str):
        try:
            args['colormap'] = getattr(colormaps, args['colormap'])
        except AttributeError:
            pass

    size_scaling_methods = {'linear': lambda x: x,
                            'sqrt': math.sqrt}
    if isinstance(args['size_scaling_fun'], str):
        args['size_scaling_fun'] = size_scaling_methods[args['size_scaling_fun']]

    print("\nargs (after processing):")
    print("\n".join("%s: %s" % kv for kv in sorted(args.items())))

    return args


def average_by_T(data, Toffset=0):
    """ Calculate melting curve (average number of hybridized domains vs T) based on data. """
    entries_by_T = defaultdict(list)
    for entry in data:
        entries_by_T[entry[0]].append(entry)
    hyb_by_T = {T-Toffset: sum(entry[1] for entry in entries)/len(entries)
                for T, entries in entries_by_T.items()}
    return hyb_by_T


def plot_statsfile(statsfile, color=None, **kwargs):
    """ Plot simulation dataset from statsfile. """
    print("Loading stats from file ", statsfile)
    stats = load_statsfile(statsfile)
    print("-", len(stats), "lines loaded from file.")


    if kwargs['max_lines'] and len(stats) > kwargs['max_lines']:
        sample_freq = int(len(stats)/kwargs['max_lines'])
        stats = [line for i, line in enumerate(stats) if i % sample_freq == 0]
        print("- Stats reduced to", len(stats), "lines (sample_freq=%s)." % sample_freq)

    print("- Pre-processing data...")
    #[(Temperature, N_doms_hybridized, %_doms_hybridized), ...]
    xdata = [line[0]-273 for line in stats]
    if kwargs['normalize'] or kwargs['normalize_by_field']:
        if kwargs['normalize_by_field']:
            ydata = [line[1]/line[kwargs['normalize_by_field']] for line in stats]
        else:
            # Normalize by the maximum value in stats
            ydata = [line[1] for line in stats]
            maxval = max(ydata)
            ydata = [val/maxval for val in ydata]
    else:
        ydata = [line[1] for line in stats]

    # Save the original xdata, ydata (T, n_hyb) pairs.
    xdata_org, ydata_org = xdata, ydata

    # Count data points for each (x, y):
    counts = defaultdict(lambda: defaultdict(int))  # A defaultdict with defaultdict(int) default entries
    for x, y in zip(xdata, ydata):
        counts[x][y] += 1
    maxcount = max(c for xdicts in counts.values() for c in xdicts.values())
    print("- Max (T, N) count:", maxcount)

    # If reducing to one marker per (x, y) point, this must be done now, before making the sizes array.
    if kwargs.get('one_mark_per_xy'):
        # Make a set of (x, y) points and plot those with size.
        xdata, ydata = zip(*set((x, y) for x, y in zip(xdata, ydata)))
        print("- Only plotting one marker per (x, y) combination; dataset reduced to", len(xdata), "points.")

    # Make marker sizes array:
    sizes = [int(kwargs['max_marker_size']*kwargs['size_scaling_fun'](counts[x][y]/maxcount))+kwargs['min_marker_size']
             for x, y in zip(xdata, ydata)]
    print("- Average marker size: %0.1f" % (sum(sizes)/len(sizes)))

    # Interpolating colors:
    if isinstance(kwargs['use_colormap'], str):
        if kwargs['use_colormap'] == "auto":
            kwargs['use_colormap'] = len(kwargs['statsfiles']) < 2
        elif kwargs['use_colormap'][0].lower() == "n":
            kwargs['use_colormap'] = False

    if kwargs['use_colormap']:
        print("- Making colors by interpolations")
        # Is an array of colors, one color per
        # I think in theory this could also be same as sizes...?
        if kwargs['colormap_base'] == 'sizes':
            colors = sizes
        else:
            colors = bilinear_interpolate(xdata, ydata, bins=kwargs['colormap_bins'])
            print("- colors array of length:", len(colors))
    else:
        colors = color

    if kwargs['gaussian_point_noise']:
        # Making random noise around each (x, y) can give a nice effect (mostly for plot, not scatter).
        #xdata = [line[0]-273 + random.gauss(0, 0.1) for line in stats]
        #ydata = [line[1] + random.gauss(0, 0.1) for line in stats]

        #xdata = [line[0]-273 + round(random.gauss(0, 0.1), 2)  for line in stats]
        #ydata = [line[1] + round(random.gauss(0, 0.1), 2) for line in stats]

        xdata = [x + random.gauss(0, 0.1) for x in xdata]
        ydata = [y + random.gauss(0, 0.1) for y in ydata]
        if kwargs['use_colormap'] and kwargs['colormap_after_point_noise']:
            colors = bilinear_interpolate(xdata, ydata, bins=3)

    label = None
    if kwargs['legend']:
        label = os.path.splitext(os.path.basename(statsfile))[0]

    if kwargs['plottype'] == "scatter":
        # Using a scatter plot is good for density-based data (e.g. flow cytometry):
        print("- Plotting data as scatter plot...")
        pyplot.scatter(xdata, ydata,
                       marker=kwargs["markershape"],
                       s=sizes,
                       c=colors,
                       edgecolors='none',
                       cmap=kwargs["colormap"],
                       label=label)
    else:
        print("- Plotting data as regular plot...")
        # Regular "one-marker-per-datapoint" plot
        # Different sizes are not supported (you can call plot repeatedly or pass in multiple sets of
        # xdata, ydata, **specs, xdata2, ydata2, **specs -- but that is not really suitable for us.
        pyplot.plot(xdata, ydata, marker=kwargs['markershape'], ms=kwargs['min_marker_size'])

    if argns.plot_average:
        print("- Making data average plot...")
        #fhyb_by_T = average_by_T(stats)
        fhyb_by_T = average_by_T(zip(xdata_org, ydata_org))
        xdata, ydata = zip(*sorted(fhyb_by_T.items()))
        pyplot.plot(xdata, ydata, '.-k')


def saveplot(plotfilename=None, **kwargs):
    """ Save the plot to file. """

    if plotfilename is None:
        plotfilename = os.path.splitext(kwargs['statsfiles'][0])[0]
        if len(kwargs['statsfiles']) > 1:
            plotfilename += "and_%s_other_plots" % (len(kwargs['statsfiles'])-1)
        plotfilename += '.png'
    print("Saving plot to file:", plotfilename)
    pyplot.savefig(plotfilename)


def showplot(**kwargs):
    """ Show the plot. """
    interactive = kwargs.get('interactive')
    if interactive is not None:
        matplotlib.interactive(interactive)
    print("Showing plot...")
    pyplot.show()
    if interactive:
        # from code import InteractiveConsole
        import code
        code.interact("Interact with the data...")


if __name__ == '__main__':

    argns = parse_args()
    args = argns.__dict__.copy()

    # Load default and overriding config and produce methods, etc
    process_args(args)
    #sys.exit()

    colors = iter(args.get('colors') or "rgbcmyk")

    if args['legend'] is None:
        args['legend'] = len(args['statsfiles']) > 1

    for statsfile in args['statsfiles']:
        color = next(colors)
        plot_statsfile(statsfile, color=color, **args)

    if args['legend']:
        pyplot.legend()

    if args['saveplot']:
        saveplot(**args)

    if args['showplot']:
        showplot(**args)

    print("Done!")
