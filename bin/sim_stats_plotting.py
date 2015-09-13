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
try:
    import numpy
except ImportError:
    pass
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
    parser.add_argument("statsfiles", nargs="*")

    # Different ways to plot/view the data:
    # * f_hyb_hist_vs_T: default, plots T on x-axis and a histogram of f_hyb on the y-axis.
    # * f_hyb_vs_N_steps: Plots N_steps on x-axis and f_hyb on y-axis.
    parser.add_argument("--plot", default="f_hyb_hist_vs_T")

    parser.add_argument("--plottype", default="scatter")
    parser.add_argument("--T-range", nargs=2, type=float)
    parser.add_argument("--figsize", nargs=2, type=float, default=[16, 10])

    parser.add_argument("--nupack_hist_vs_T",
                        help="Plot NuPack histogram files vs T. Input must be directory containing a folder for each T.")

    # complex_size_hist_vs_T
    parser.add_argument("--complex_size_hist_vs_T", nargs="*",
                        help="Plot complex size histogram data vs T.")

    parser.add_argument("--thermocurves", nargs="*")#, default="bgrcmyk")


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

    parser.add_argument("plot_average", action="store_true", default=None)

    #
    ## show plot
    parser.add_argument("--no-showplot", dest="showplot", action="store_false",
                        help="Do not show the plot figure after plotting.")
    parser.add_argument("--showplot", action="store_true",
                        help="Show the plot figure after plotting.")

    ## save plot
    parser.add_argument("--no-saveplot", dest="saveplot", action="store_false",
                        help="Do not save the plot figure after plotting.")
    parser.add_argument("--saveplot", action="store_true",
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


def plot_thermodynamic_meltingcurve(cumstatsfile, KtoC=True, linespec=':', **kwargs):
    """
    Plot meltingcurve cumstatsfile.
    Format is:
        cum_stats[T] = (hybridized, non_hybridized, total_conc, domains_total)
    """
    if isinstance(cumstatsfile, str):
        with open(cumstatsfile) as fd:
            cumstats = yaml.load(fd)
    else:
        cumstats = cumstatsfile
    # Plot fraction hybridized:
    KtoC = 273.15 if KtoC else 0
    Ts, f_hyb = zip(*[(T-KtoC, stats[0]/stats[2]) for T, stats in sorted(cumstats.items())])
    label = kwargs.pop('label', os.path.splitext(os.path.basename(cumstatsfile))[0])
    if "dom_anneal_stats" in label:
        label = label.replace("dom_anneal_stats", "")
        foldername = os.path.basename(os.path.dirname(os.path.abspath(
            cumstatsfile)))
        label = foldername + label

    pyplot.plot(Ts, f_hyb, linespec, label=label, **kwargs)


def parse_complex_size_hist_data(fp):
    data = ({int(T): {int(size): int(count)
                      for size, count in ((int(v) for v in pair.split(","))
                                          for pair in pairs.strip().split("\t") if pair.strip())
                     }
                     # if pairs.strip() else {int(T): {}} # include empty as {T: {}}
            } for T, pairs in (line.strip().split(":") for line in fp if line.strip()))
    return data

def load_complex_size_hist_file(complex_size_hist_file):

    if isinstance(complex_size_hist_file, str):
        filenameroot, ext = os.path.splitext(complex_size_hist_file)
        if ('yaml' in ext or 'yml' in ext):
            print(" - loading as yaml file")
            raise ValueError("YAML no longer supported.")
        #data = yaml.load(open(complex_size_hist_file))
        # Format is:
        # # T: size,count  size,count  ...
        # output is list with list of dicts:
        # [{T: {N: count of complexes with size N}}]
        print("Opening file...")
        with open(complex_size_hist_file) as fp:
            #complex_size_hist_raw = [{int(T): {int(size), int(count) for pair in pairs
            #                                   for size, count in (int(N) for N in pair.split(","))}}
            #                         for T, pairs in ((T, rest.split("\t"))
            #                            for T, rest in (line.split(":") for line in fp))]
            data = parse_complex_size_hist_data(fp)
            data = list(data)   # Read data while file is still open
    else:
        data = parse_complex_size_hist_data(complex_size_hist_file)
        #for T,
        #(T, (pair.split(",") for pair in rest.split("\t"))
        ## for T, rest in (line.split(":") for line in fp))]
        #data = []
        #for line in fp:
        #    try:
        #        T, pairs = line.strip().split(":")
        #    except ValueError as e:
        #        print("ValueError (a),", e)
        #        print(line)
        #        raise e
        #    entry = {}
        #    try:
        #        for pair in pairs.split("\t"):
        #            if pair.strip():
        #                size, count = (int(v) for v in pair.split(","))
        #                entry[size] = count
        #        data.append(entry)
        #        #data.append({int(T): {int(size): int(count)})
        #        #data.append({int(T): {int(size): int(count)
        #        #                      for size, count in [(int(v) for v in pair.split(","))
        #        #                                          for pair in pairs.split("\t") if pair.strip()]
        #        #                     }
        #        #            })
        #    except ValueError as e:
        #        #input("hej")
        #        print("ValueError,", e)
        #        print(line)
        #        print(T)
        #        print(pairs)
        #        #raise(e)
        #    #hist = dict(map(int, pair.split(",")) for pair in sizes_and_counts.split("\t"))
        #    #hist = {int(size): int(count) for pair in sizes_and_counts.split("\t") for k, v in )
        #for line in fp:
        #    T, sizes_and_counts = line.split(":")
        #    #for pair in sizes_and_counts.split("\t"):
        #    #    size, count = (int(v) for v in pair.split(","))
        #        if pairs.strip()]
    return data



def plot_complex_size_hist_vs_T(complex_size_hist_file, color=None, **kwargs):

    # A list of dicts: [{T: {N: number of complexes of size N}}]
    print("\nLoading size_hist_file", complex_size_hist_file)
    with open(complex_size_hist_file) as fp:
        # We can get a generator back if we give a file object rather than string.
        #complex_size_hist_raw = load_complex_size_hist_file(complex_size_hist_file)
        complex_size_hist_raw = load_complex_size_hist_file(fp)
        # print(" - %s entries loaded!" % len(complex_size_hist_raw))
        print(" - Pre-processing size hist data...")
        # Count data points for each (T, sizeN):
        #counts_list = defaultdict(lambda: defaultdict(list))  # A defaultdict with defaultdict(list) default entries
        counts_agg = defaultdict(lambda: defaultdict(int))  # A defaultdict with defaultdict(int) default entries
        # counts = {T: {c_size: [1,2,5,6,8 list of counts]}}
        #grouped_by_T = defaultdict(list)    # {T: [{size: count}}}
        for l, entry in enumerate(complex_size_hist_raw):
            for T, hist in entry.items():
                #grouped_by_T[T].append(hist)
                for c_size, count in hist.items():
                    #counts_list[T][c_size].append(count)
                    counts_agg[T][c_size] += (count)
                    if l % 100000 == 0:
                        print(" - Loaded %s lines from file and counting..." % l)
    print(" - Done loading pre-processing data!")

    #means = {}  # {T: }
    #for T, t_stats in counts.items():
    #    for c_size, size_counts in

    #avg_counts = {(T, c_size): sum(counts_lst)/len(counts_lst)
    #               for T, counts_at_T in counts_list.items()
    #               for c_size, counts_lst in counts_at_T.items()}

    # The highest value of counts_agg[T][c_size]:
    agg_counts_max = max(total for agg_counts_at_T in counts_agg.values() for total in agg_counts_at_T.values())
    #agg_counts_avg = {(T, c_size): total/agg_counts_max
    #                  for T, agg_counts_at_T in counts_agg.items()
    #                  for c_size, total in agg_counts_at_T.items()}

    #xdata, ydata = zip(*[(T, c_size) for T, count_at_T in counts_list.items() for c_size, counts_lst in count_at_T.items()])
    xdata, ydata = zip(*[(T, c_size) for T, agg_count_at_T in counts_agg.items() for c_size, total in agg_count_at_T.items()])

    print("- Max (T, N) count:", agg_counts_max)
    # Make marker sizes array:
    #sizes = [int(kwargs['max_marker_size']*kwargs['size_scaling_fun'](avg_counts[(x,y)]/agg_counts_max))+kwargs['min_marker_size']+5
    sizes = [int(kwargs['max_marker_size']*kwargs['size_scaling_fun'](counts_agg[x][y]/agg_counts_max))+kwargs['min_marker_size']+5
             for x, y in zip(xdata, ydata)]
    print("- Average marker size: %0.1f" % (sum(sizes)/len(sizes)))

    if isinstance(kwargs['use_colormap'], str):
        if kwargs['use_colormap'] == "auto":
            kwargs['use_colormap'] = len(kwargs['complex_size_hist_vs_T']) < 2
        elif kwargs['use_colormap'][0].lower() == "n":
            kwargs['use_colormap'] = False
    if kwargs['use_colormap']:
        colors = sizes
    else:
        colors = color

    label = os.path.splitext(os.path.basename(complex_size_hist_file))[0]
    if "dom_anneal_stats" in label:
        label = label.replace("dom_anneal_stats", "")
        foldername = os.path.basename(os.path.dirname(os.path.abspath(
            complex_size_hist_file)))
        label = foldername + label

    # convert K to C:
    xdata = [T-273.15 for T in xdata]

    if kwargs['plottype'] == "scatter":
        # Using a scatter plot is good for density-based data (e.g. flow cytometry):
        print("- Plotting data as scatter plot... (color=%s)" % color)
        pyplot.scatter(xdata, ydata,
                       marker=kwargs["markershape"],
                       s=sizes,
                       c=colors,
                       edgecolors='none',
                       cmap=kwargs["colormap"],
                       label=label)

    if kwargs.get('plot_average') and False:
        print("- Making data average plot...")
        avg_by_T = average_by_T(zip(xdata, ydata))
        xdata, ydata = zip(*sorted(avg_by_T.items()))
        pyplot.plot(xdata, ydata, marker=' ', c=color, ls='-', markeredgecolor='k')  # '.-k'
    print(" - finished plotting data from", os.path.basename(complex_size_hist_file))


def load_and_process_statsfile(statsfile, **kwargs):
    """
    The data structure is a list of rows, each row with:
        (T, N_domains_hyb, f_domains_hyb, N_strands_hyb, f_strands_hyb)
    """
    print("Loading stats from file ", statsfile)
    stats = load_statsfile(statsfile)
    print("-", len(stats), "lines loaded from file.")
    KtoC = 273.15
    if kwargs['T_range']:
        minT, maxT = [K+KtoC for K in kwargs['T_range']] # Convert to Kelvin
        #print("type(minT), type(maxT), type(stats[0][0]):", type(minT), type(maxT), type(stats[0][0]))
        stats = [line for line in stats if line[0] >= minT and line[0] <= maxT]

    if kwargs['max_lines'] and len(stats) > kwargs['max_lines']:
        sample_freq = int(len(stats)/kwargs['max_lines'])
        stats = [line for i, line in enumerate(stats) if i % sample_freq == 0]
        print("- Stats reduced to", len(stats), "lines (sample_freq=%s)." % sample_freq)

    return stats


def plot_statsfile(statsfile, color=None, **kwargs):
    """ Plot simulation dataset from statsfile. """

    stats = load_and_process_statsfile(statsfile, **kwargs)

    print("- Pre-processing data...")
    #[(Temperature, N_doms_hybridized, %_doms_hybridized), ...]
    xdata = [line[0]-273 for line in stats]
    if kwargs['normalize'] or kwargs['normalize_by_field']:
        if kwargs['normalize_by_field']:
            # As long as the data rows are tuples/lists and not dicts, index must be an integer:
            normalize_by_field = int(kwargs['normalize_by_field'])
            ydata = [line[1]/line[normalize_by_field] for line in stats]
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
    label = os.path.splitext(os.path.basename(statsfile))[0]
    if "dom_anneal_stats" in label:
        label = label.replace("dom_anneal_stats", "")
        foldername = os.path.basename(os.path.dirname(os.path.abspath(
            statsfile)))
        label = foldername + label


    if kwargs['plottype'] == "scatter":
        # Using a scatter plot is good for density-based data (e.g. flow cytometry):
        print("- Plotting data as scatter plot... (color=%s)" % color)
        pyplot.scatter(xdata, ydata,
                       marker=kwargs["markershape"],
                       s=sizes,
                       c=colors,
                       edgecolors='none',
                       cmap=kwargs["colormap"],
                       label=label)
    else:
        print("- Plotting data as regular plot... (color=%s)" % color)
        # Regular "one-marker-per-datapoint" plot
        # Different sizes are not supported (you can call plot repeatedly or pass in multiple sets of
        # xdata, ydata, **specs, xdata2, ydata2, **specs -- but that is not really suitable for us.
        pyplot.plot(xdata, ydata, c=color, marker=kwargs['markershape'], ms=kwargs['min_marker_size'])

    if kwargs.get('plot_average'):
        print("- Making data average plot...")
        #fhyb_by_T = average_by_T(stats)
        fhyb_by_T = average_by_T(zip(xdata_org, ydata_org))
        xdata, ydata = zip(*sorted(fhyb_by_T.items()))
        pyplot.plot(xdata, ydata, marker=' ', c=color, ls='-', markeredgecolor='k')  # '.-k'


def saveplot(plotfilename=None, **kwargs):
    """ Save the plot to file. """

    if plotfilename is None:
        if kwargs['statsfiles']:
            plotfilename = os.path.splitext(kwargs['statsfiles'][0])[0]
            if len(kwargs['statsfiles']) > 1:
                plotfilename += "and_%s_other_plots" % (len(kwargs['statsfiles'])-1)
        elif kwargs['complex_size_hist_vs_T']:
            plotfilename = os.path.splitext(kwargs['complex_size_hist_vs_T'][0])[0]
            if len(kwargs['complex_size_hist_vs_T']) > 1:
                plotfilename += "and_%s_other_plots" % (len(kwargs['complex_size_hist_vs_T'])-1)
        elif kwargs['thermocurves'] and kwargs['thermocurves'][0] != "auto":
            plotfilename = os.path.splitext(kwargs['thermocurves'][0])[0]
        elif kwargs['nupack_hist_vs_T']:
            plotfilename = os.path.splitext(kwargs['nupack_hist_vs_T'][0])[0]
        else:
            plotfilename = "sim_stats_plot"
        if kwargs['plot'] != 'f_hyb_hist_vs_T':
            plotfilename += "_" + kwargs['plot']
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


def plot_nupack_f_complexed_strands_vs_T(nupack_hist_dir, jobid=None, **kwargs):
    """
    Plot the fraction of complexed strands (strands that are not alone) as a function of T,
    using data from nupack.
    """
    # data structure is:
    # dict with dicts:
    # {T: {ComplexID: (ordered-complex-id, (strands stoichiometry list), deltaG, concentration}}
    hist_vs_T_data = load_nupack_hist_vs_T_dir(nupack_hist_dir, jobid=None)
    f_complexed_vs_T = calc_f_complexed_vs_T(hist_vs_T_data)
    x, y = zip(*((k, tup[0]) for k, tup in sorted(f_complexed_vs_T.items())))
    pyplot.plot(x, y, **kwargs)



def load_nupack_hist_vs_T_dir(nupack_hist_dir, jobid=None):
    """
    return a dict of {T: <hist-data-for-T>}
    where <hist-data-for-T> is:
        {ComplexID: (ordered-complex-id, (strands stoichiometry list), deltaG, concentration)}
    The total number of strands in a complex is:
        sum(strands stoichiometry list)
    """
    if jobid is None:
        head, tail = os.path.split(nupack_hist_dir)
        if not tail:
            # folder given as "/path/to/nupack/job/" with trailing slash
            head, tail = os.path.split(head)
            if not tail:
                raise ValueError("nupack_hist_dir does not have a folder name, and jobid was not given.")
        try:
            jobid = int(tail)
        except ValueError:
            jobid = int(tail.rsplit("_", 1)[-1])
    if str(jobid) in os.listdir(nupack_hist_dir):
        nupack_hist_dir = os.path.join(nupack_hist_dir, str(jobid))
    hist_vs_T_data = {}
    for folder in os.listdir(nupack_hist_dir):
        fpath = os.path.join(nupack_hist_dir, folder)
        if not os.path.isdir(fpath):
            continue
        try:
            T = float(folder)
        except ValueError:
            continue
        hist_filename = "%s_hist.dat" % jobid
        hist_file = os.path.join(fpath, hist_filename)
        hist_vs_T_data[T] = parse_nupack_hist_data_file(hist_file)
    return hist_vs_T_data


def calc_f_complexed_vs_T(hist_vs_T_data):
    """
    Return a dict with
        {T: fraction_of_strands_in_a_complex}
    """
    return {T: calc_f_complexed_strands(hist_data) for T, hist_data in hist_vs_T_data.items()}

def calc_f_complexed_strands(hist_data):
    """
    Calculate the fraction of strands that are part of a complex.
    hist_data structure is:
        {ComplexID: (ordered-complex-id, (strands stoichiometry list), deltaG, concentration)}
    Returns a tuple of,
        c_complexed/(c_complexed+c_singles), c_complexed, c_singles
    """
    c_complexed = 0
    c_singles = 0
    for cid, (oid, strands_stroichiometry, deltaG, conc) in hist_data.items():
        N_strands = sum(strands_stroichiometry)
        if N_strands > 1:
            c_complexed += conc*N_strands   # Should we multiply by N_strands?
        else:
            c_singles += conc
    return c_complexed/(c_complexed+c_singles), c_complexed, c_singles





def parse_nupack_hist_data_file(hist_file):
    """
    Parse a NuPack <job-no>_hist.dat data file.
    Returns a <hist-data-for-T> dict:
        {ComplexID: (ordered-complex-id, (strands stoichiometry list), deltaG, concentration)}
    The total number of strands in a complex is:
        sum(strands stoichiometry list)
    """
    with open(hist_file) as fd:
        rows = [line.strip().split("\t") for line in fd if line and line.strip()[0] != "%"]
    hist_data = {}
    for row in rows:
        assert len(row) > 4
        complexid = int(row[0])
        oid = int(row[1])
        try:
            strand_stoichiometry = [int(val) for val in row[2:-2]]
        except ValueError as e:
            print("ValueError parsing strand_stoichiometry in row")
            print("row:       %s" % row)
            print("row[2:-2]: %s" % row[2:-2])
            raise e
        deltaG = float(row[-2])
        conc = float(row[-1])
        hist_data[complexid] = (oid, strand_stoichiometry, deltaG, conc)
    return hist_data


def plot_f_hyb_vs_N_steps(statsfile, plotargs=None, **kwargs):
    """
    Plot fraction of hybridized domains versus line number
    (because we are currently not storing the N_step in the simulation statsfile)
    """
    if plotargs is None:
        plotargs = {}
    stats = load_and_process_statsfile(statsfile, **kwargs)
    if len(stats) < 1:
        print("No stats returned after processing:")
        print(stats)
        return
    y = [row[2] for row in stats]
    print("Plotting %s f_hyb points, average value is %s" % (
        len(y), sum(y)/len(y)))
    pyplot.plot(y, **plotargs)


def main():

    argns = parse_args()
    args = argns.__dict__.copy()

    # Load default and overriding config and produce methods, etc
    process_args(args)
    #sys.exit()

    colors = iter(args.get('colors') or "rgbcmyk")  # make cyclic?

    if args['legend'] is None:
        args['legend'] = len(args['statsfiles']) > 1

    for i, statsfile in enumerate(args['statsfiles']):
        color = next(colors)
        if args['plot'] == 'f_hyb_vs_N_steps':
            plotargs = {'color': color, 'figsize': args['figsize']}
            print("Plotting f_hyb_vs_N_steps...")
            plot_f_hyb_vs_N_steps(statsfile, plotargs=plotargs, **args)
        else:
            plot_statsfile(statsfile, color=color, **args)
        # We may have specified a "thermocurve" for each of the simulations:
        if args['thermocurves'] and (len(args['thermocurves']) > 1 or args['thermocurves'][0] == 'auto'):
            try:
                # Done: Add 'auto' argument.
                if args['thermocurves'][0] == 'auto':
                    cumstatsfile = os.path.splitext(statsfile)[0].rsplit("_", 1)[0]+".thermo_melting.yaml"
                else:
                    cumstatsfile = args['thermocurves'][i]
            except IndexError as e:
                print("\n- WARNING: IndexError while finding cumstatsfilename. args['thermocurves'] =",
                      args['thermocurves'], ". Error is:", e, "\n")
            else:
                print("- Plotting thermodynamic melting curve from", cumstatsfile)
                try:
                    plot_thermodynamic_meltingcurve(cumstatsfile, c=color)
                except OSError as e:
                    print("\n- WARNING: OSError while printing thermocurve:", e, "\n")

    if args['thermocurves'] and len(args['thermocurves']) == 1 \
        and args['thermocurves'][0] != 'auto':
        if args['plot'] != 'f_hyb_vs_N_steps':
            cumstatsfile = args['thermocurves'][0]
            plot_thermodynamic_meltingcurve(cumstatsfile, c='k')

    if args.get('complex_size_hist_vs_T'):
        color = next(colors)
        for size_hist_file in args['complex_size_hist_vs_T']:
            plot_complex_size_hist_vs_T(size_hist_file, color=color, **args)

    ## TODO: Make f_hyb vs N_steps plot (rather than the current f_hyb-vs-T histogram plot)

    ## TODO: Make NuPack f_complexed_hist-vs-T plot
    if args.get('nupack_hist_vs_T'):
        nupack_hist_dir = args['nupack_hist_vs_T']
        if args['plot'] != 'f_hyb_vs_N_steps':
            plot_nupack_f_complexed_strands_vs_T(nupack_hist_dir, color='orange', label="NuPack % complexed strands")

    if args['figsize']:
        fig = matplotlib.pyplot.gcf()
        print("Setting figure size to:", args['figsize'])
        fig.set_size_inches(*args['figsize'])

    if args['legend']:
        pyplot.legend()

    if args['saveplot']:
        saveplot(**args)

    if args['showplot']:
        showplot(**args)

    print("Done!")



if __name__ == '__main__':
    main()
