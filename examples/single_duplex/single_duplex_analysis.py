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

Plotting libraries (2D plotting):
* matplotlib
* seaborn and pandas both uses matplotlib as backend.
* ggplot (R-like plotting)
* bokeh (targets webbrowser visualization)
* pygal - svg-based plotting, simpler and more pythonic than matplotlib but less features.
* Chaco - interactive plotting wtih Qt
* plotly - API to online plotting/visualization service.
* GR - OpenGL visualization, can be used as a matplotlib backend with improved performance.
* Vispy - High-performance GPU 2D/3D OpenGL visualization. Can be used as an experimental matplotlib backend.
** Vispy is a combined project of four initial projects, pyqtgraph, visvis, Galry, and Glumpy.


Plotting refs:
* http://pbpython.com/visualization-tools-1.html


"""

import os
import sys
import webbrowser
import shutil

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
from nascent.stat_analysis.processing import load_multiple_stats #, process_stats
from nascent.stat_analysis.meltingcurve import plot_thermodynamic_meltingcurve


colors = list("rgbcmk") + "silver firebrick darkslategrey orange darkmagenta".split()




def main():

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
    runidxs = [-1, -2, -3]
    runidxs = [-1, -2, -3, -4]
    runidxs = [-1, -2, -4]
    runidxs = [-1, -2, -3, -5]
    runidxs = [-1, -2, -3, -4, -6]
    runidxs = [-1, -2, -3, -4, -5, -7]
    runidxs = [-7]
    # runidxs = [-1, -2, -3, -4, -5, -6, -8]
    # runidxs = [-1, -3]
    # runidxs = [-1, -2, -3, -4, -5]
    # runidxs = [-1, -2, -3, -4, -5, -6]
    # runidxs = [-1, -2, -3, -4, -5, -6, -7, -8, -9]
    # runidxs = [-2, -3, -4, -5]
    #runidxs = [-4, -5]
    # runidxs = [-3, -4, -5]
    # runidxs = [-2, -3, -4]
    stats, statsfolders = load_multiple_stats(runidxs=runidxs, basedir=scriptdir, structure=structure, process=True)
    statsfolder = statsfolders[0]

    ## Process (returns a Pandas DataFrame):
    # shift_tau_for_duration is for older stats where stats were collected *after* changing state.
    # stats = [process_stats(runstats) for runstats in stats]
    # Edit: Done in load_multiple_stats through process=True.

    pyplot = load_pyplot()

    ## Plot fraction of hybridized domains:
    if plot_tot_hyb:
        ax = None
        # Instead of passing plot parameters through via function args, consider using the
        # pd.plot_params.use context manager. EDIT: Currently only works for xaxis.compat option :\
        # with pd.plot_params.use('logx', False):
        # "Specific lines can be excluded from the automatic legend element selection by defining a label
        # starting with an underscore."
        for runstats, runidx, color in zip(stats, runidxs, colors[:len(stats)]):
            plotfilename = os.path.join(statsfolder, "f_hybridized_domains_avg_vs_time.png")
            # labels, markers, colors, etc all match up against the equivalent data field.
            fields=('f_hybridized_domains_avg', 'f_partially_hybridized_strands_avg', 'f_fully_hybridized_strands_avg')
            fields=('f_partially_hybridized_strands_avg', 'f_fully_hybridized_strands_avg')
            fields=('f_partially_hybridized_strands_avg', )
            # fields=('f_fully_hybridized_strands_avg', )      # T-loop domains never hybridize...
            ax = plot_tot_vs_time(runstats, filename=None, #plotfilename,
                                  ax=ax,
                                  x='duration_cum',
                                  # x='system_time_end',
                                  figsize=(16, 10),
                                  #linestyles=("-", "None"),
                                  #colors=color*5,
                                  #kind='line',
                                  #kind='scatter',
                                  # fields=('f_hybridized_domains_avg',),
                                  fields=fields,
                                  # fields=('f_hybridized_domains_rolling', 'f_hybridized_domains'),# , ''),
                                  #labels=("f_hyb_dom run %s" % runidx, ),
                                  #labels=("f_hyb_dom run %s" % runidx, 'f_hybridized_domains'),
                                  #labels=("_", "f_hyb_dom run %s" % runidx),
                                  #markers='o.',
                                  marker="o",
                                  markersize=3,
                                  markeredgecolor='None',
                                  #markeredgewidth=0.2,
                                  #markerfacecolor='r',
                                  legend=True,
                                  alpha=0.3,
                                  linewidth=0.5,
                                 )
        # import pdb
        # pdb.set_trace()
        # handles, labels = zip(*[(hdl, lbl) for hdl, lbl in zip(*ax.get_legend_handles_labels())
        #                         if "avg" not in lbl and lbl[0] != "_"])
        # ax.legend(handles=handles, labels=labels, loc=None) # loc="lower right")
        # ax.xlim(0, 100) # ax has no xlim attribute...
        #pyplot.xlim(0, 200)
        pyplot.xlim(xmin=-0.05)
        #pyplot.ylim(ymin=-0.05)
        # pyplot.ylim(ymin=0, ymax=1.0)
        pyplot.savefig(plotfilename)
        webbrowser.open(plotfilename)
        print("Plot save to file:", plotfilename)
        plotfilename_with_runidxs = os.path.join(
            statsfolder, "f_hybridized_domains_avg_vs_time" + "_".join(str(i) for i in runidxs) + ".png")
        print("Extra copy of plot:", plotfilename_with_runidxs)
        shutil.copy2(plotfilename, plotfilename_with_runidxs)


    ## Plot stacked ends:
    if plot_tot_stacked:
        plotfilename = os.path.join(statsfolder, "f_stacked_ends_avg_vs_time.png")
        ax = plot_tot_vs_time(stats, fields=('f_stacked_ends_avg',), add_average=False, filename=plotfilename)
        webbrowser.open(plotfilename)
        print("Plot save to file:", plotfilename)
        plotfilename_with_runidxs = os.path.join(
            statsfolder, "f_stacked_ends_avg_vs_time" + "_".join(str(i) for i in runidxs) + ".png")
        print("Extra copy of plot:", plotfilename_with_runidxs)
        shutil.copy2(plotfilename, plotfilename_with_runidxs)



    ## Plot melting curve?
    if plot_melting_curve:
        meltingcurvefn = "thermo_melting.yaml"
        meltingcurvefn = os.path.join(statsfolder, meltingcurvefn)
        print("Plotting melting curve from file:", meltingcurvefn)
        plot_thermodynamic_meltingcurve(meltingcurvefn, KtoC=False)





if __name__ == "__main__":
    main()
