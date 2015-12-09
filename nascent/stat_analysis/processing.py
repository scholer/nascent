# -*- coding: utf-8 -*-
##    Copyright 2015 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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

# pylint: disable=C0103,W0212,W0142,W0141

"""

General stats loading and processing module.


"""
import os
import pandas as pd

from nascent.graph_sim_nx.stats_manager import StatsReader



def load_stats(statsfile="time_totstats.txt", structure="duplex2",
               statsfolder=None, data_root_directory=None, basedir=".", runidx=-1):
    """
    Load latest stats for structure :structure:
    If statsfile does not immediately exists in current directory,
    will create data_root_directory = <basedir>/simdata/<structure>/
    and statsfolder = <data_root_directory>/<last-directory-after-sorting-by-name>
    and try to load <statsfolder>/<statsfile>.
    (This works in the default situation...)
    """

    statsreader = StatsReader(config=None)

    if not os.path.exists(statsfile):
        if statsfolder is None:
            if data_root_directory is None:
                #scriptdir = os.path.dirname(os.path.abspath(__file__))
                data_root_directory = os.path.join(basedir, "simdata", structure)
                print("Using data root directory:", data_root_directory)
            cands = (os.path.join(data_root_directory, folder) for folder in os.listdir(data_root_directory))
            cands = list(filter(os.path.isdir, cands))
            cands.sort(reverse=False)
            statsfolder = cands[runidx]
            print("Using stats in folder:", statsfolder)
        statsfile = os.path.join(statsfolder, statsfile)
    else:
        statsfolder = "."
    stats = statsreader.load_total_stats_file(statsfile)
    print("%s stats entries loaded from file %s" % (len(stats), statsfile))
    return stats, statsfolder


def load_multiple_stats(statsfile="time_totstats.txt", structure="duplex2",
                        statsfolder=None, data_root_directory=None, basedir=None, runidxs=(-1,),
                        process=False):
    stats, statsfolders = [], []
    for runidx in runidxs:
        runstats, _statsfolder = load_stats(statsfile=statsfile, structure=structure, statsfolder=statsfolder,
                                           data_root_directory=data_root_directory, basedir=basedir, runidx=runidx)
        if process:
            runstats = process_stats(runstats)
        stats.append(runstats)
        statsfolders.append(_statsfolder)
    return stats, statsfolders


def process_stats(stats, shift_tau_for_duration=False):
    """
    Process a basic list of dicts loaded from file,
    calculate a number of derived values for each column,
    and return a pandas dataframe.
    The following derived columns are calculated:
     - duration - tau, shifted 1 row down.
     - duration_cum - system_time, shifted 1 row down.
     - f_hybridized_domains
     - f_stacked_ends
     - duration-weighted, cummulative and cummulative average of
        n/f_hybridized/stacked_domains/ends. I.e.
            n_hybridized_weighted is n_hybridized*duration,
            n_hybridized_cum is cumulative sum of n_hybridized_weighted, and
            n_hybridized_ave is n_hybridized_cum/duration_cum.
    And PS: Please be an adult about the abbreviation of cummulative.
    NOTE: Maybe collect stats and tau after incrementing sys-time but before saving state,
    so that tau "time that system was in the noted state" rather than being "time before switching to noted state".
    """
    df = pd.DataFrame(stats)
    ## In old versions, I would collect tau and system_time *after* performing the reaction.
    ## This made sense because then system_time is the time when the state *starts*.
    ## Tau, is the duration of the *previous* state.
    duration_key = 'duration' if shift_tau_for_duration else 'tau'
    if shift_tau_for_duration:
        df['duration'] = df['tau'].shift(-1)

    df['f_hybridized_domains'] = df['n_hybridized_domains'] / df['N_domains']
    df['f_stacked_ends'] = df['n_stacked_ends'] / df['N_domains'] / 2
    df['duration_cum'] = df[duration_key].cumsum()
    avg_trace_cols = ('n_hybridized_domains', 'f_hybridized_domains', 'n_stacked_ends', 'f_stacked_ends')
    for col in avg_trace_cols:
        col_weighted = col + "_weighted"
        col_weighted_cum = col_weighted + '_cum'
        col_ave = col + "_avg"
        df[col_weighted] = df[col] * df[duration_key]
        df[col_weighted_cum] = df[col_weighted].cumsum()
        df[col_ave] = df[col_weighted_cum] / df['duration_cum']
    print("Done processing stats!")
    return df
