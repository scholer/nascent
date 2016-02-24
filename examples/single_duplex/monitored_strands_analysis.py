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

Run a single analysis of a single duplex.


"""

from __future__ import absolute_import, print_function
import os
import sys
import webbrowser
import shutil
import pandas as pd
import pdb

# Run from package home dir or have nascent on your python path:
sys.path.insert(0, ".")
print("os.path.abspath('.'):", os.path.abspath('.'))
scriptdir = os.path.dirname(os.path.abspath(__file__))
examples_dir = os.path.dirname(scriptdir)
LIBPATH = os.path.dirname(examples_dir)
try:
    import nascent
except ImportError:
    sys.path.insert(0, LIBPATH)
    import nascent
from nascent.stat_analysis.plotting import load_pyplot, plot_tot_vs_time
from nascent.stat_analysis.processing import load_multiple_stats, get_datafile #, process_stats
from nascent.stat_analysis.meltingcurve import plot_thermodynamic_meltingcurve


colors = list("rgbcmk") + "silver firebrick darkslategrey orange darkmagenta".split()




def main():
    pyplot = load_pyplot()

    plot_tot_hyb = True
    plot_tot_stacked = False
    plot_melting_curve = False

    structure = "duplex2"
    # structure = "duplex_16bp-d2"
    structure = "duplex_16bp_2d_10b-loop"
    structure = "duplex_20bp_2d_5T-loop"
    structure = "duplex_20bp_2d_4T-loop"
    structure = "duplex_20bp_2d_2T-loop"

    # stats, statsfolders = load_stats()
    runidxs = [-1] #
    # runidxs = [-2] #
    # runidxs = [-3] #
    # runidxs = [-1, -2, -3, -6, -7, -8, -9, -10]
    # runidxs = [-11] #
    # runidxs = [-1, -2]
    # runidxs = [-1, -2, -3]
    # runidxs = [-1, -2, -3, -4]
    # runidxs = [-1, -2, -4]
    # runidxs = [-1, -2, -3, -5]
    # runidxs = [-1, -2, -3, -4, -6]
    # runidxs = [-1, -2, -3, -4, -5, -7]
    # runidxs = [-7]
    # runidxs = [-1, -2, -3, -4, -5, -6, -8]
    # runidxs = [-1, -3]
    # runidxs = [-1, -2, -3, -4, -5]
    # runidxs = [-1, -2, -3, -4, -5, -6]
    # runidxs = [-1, -2, -3, -4, -5, -6, -7, -8, -9]
    # runidxs = [-2, -3, -4, -5]
    #runidxs = [-4, -5]
    # runidxs = [-3, -4, -5]
    # runidxs = [-2, -3, -4]
    #stats, statsfolders = load_multiple_stats(runidxs=runidxs, basedir=scriptdir, structure=structure, process=True)
    statsfiles, statsfolders = zip(*[get_datafile(statsfile="monitored_strands_stats.tsv", runidx=runidx,
                                                  basedir=scriptdir, structure=structure)
                                   for runidx in runidxs])
    stats = [pd.read_table(statsfile) for statsfile in statsfiles] # http://pandas.pydata.org/pandas-docs/stable/io.html
    statsfolder = statsfolders[0]

    ## Process (returns a Pandas DataFrame):
    # http://pandas.pydata.org/pandas-docs/stable/groupby.html

    plotfilename = os.path.join(statsfolder, "strand_state_vs_time.png")
    ax = None
    # Instead of passing plot parameters through via function args, consider using the
    # pd.plot_params.use context manager. EDIT: Currently only works for xaxis.compat option :\
    # with pd.plot_params.use('logx', False):
    # "Specific lines can be excluded from the automatic legend element selection
    # by defining a label starting with an underscore."
    for runstats, runidx, color in zip(stats, runidxs, colors[:len(stats)]):
        # runstats = stats for a single run
        # labels, markers, colors, etc all match up against the equivalent data field.
        pltargs = dict(marker="o", markersize=3, markeredgecolor='None', legend=True, alpha=0.3, linewidth=1)
        runstats['system_time_end'] = runstats['system_time'].shift(-1).dropna()
        # The above should be equivalent to
        # runstats['system_time_end'] = runstats['tau'].cumsum()
        # for name, data in DataFrameGroupBy is same as DataFrameGroupBy.group.items()
        for strand, s_trace in runstats.groupby(['strand_uid'], sort=False):
            s_state_traces = s_trace.groupby(['complex_state'], sort=False)
            for cstate, s_state in s_state_traces:
                s_state['s_state_time_cum'] = s_state['tau'].cumsum()
                s_state['s_state_partition'] = s_state['s_state_time_cum']/s_state['system_time_end']
                s_state.dropna(inplace=True) # Rows with NaN can sometimes prevent proper plotting
                # To show multiple select columns:
                #   s_state[['system_time_end', 'tau', 's_state_time_cum', 's_state_partition']]
                #   Note: df[['a','b']] produces a copy, so only use for debug representation...
                if ax is None:
                    pltargs['figsize'] = (16, 10)
                else:
                    pltargs['ax'] = ax
                pltargs['label'] = "s#%s in cstate %s" % (strand, cstate % 1000)
                ax = s_state.plot(kind='line',
                                  # x='system_time',
                                  x='system_time_end',
                                  y='s_state_partition',
                                  **pltargs)
            # Can just plot all groups together?
            #s_state_traces.plot.scatter(x='system_time', y='s_state_partition', **pltargs)

    pyplot.xlim(xmin=0.0) # xmin=-0.05)
    # pyplot.ylim(ymin=-0.05, ymax=1.05)
    pyplot.ylim(ymin=-0.01, ymax=1.01)
    # pyplot.ylim(ymin=0.0, ymax=1.0)
    pyplot.savefig(plotfilename)
    print("Plot save to file:", plotfilename)
    plotfilename_with_runidxs = os.path.join(
        statsfolder, "strand_state_vs_time" + "_".join(str(i) for i in runidxs) + ".png")
    print("Extra copy of plot:", plotfilename_with_runidxs)
    shutil.copy2(plotfilename, plotfilename_with_runidxs)
    webbrowser.open('file://' + plotfilename_with_runidxs)




if __name__ == "__main__":
    main()
