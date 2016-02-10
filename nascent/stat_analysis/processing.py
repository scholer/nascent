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
from pprint import pprint

from nascent.graph_sim_nx.stats_manager import StatsReader


def get_datafile(statsfile="time_totstats.txt", structure="duplex2",
                 statsfolder=None, data_root_dir=None, basedir=".", runidx=-1):
    if not os.path.exists(statsfile):
        if statsfolder is None:
            if data_root_dir is None:
                #scriptdir = os.path.dirname(os.path.abspath(__file__))
                data_root_dir = os.path.join(basedir, "simdata", structure)
                print("Using data root directory:", data_root_dir)
            cands = (os.path.join(data_root_dir, folder) for folder in os.listdir(data_root_dir))
            cands = list(filter(os.path.isdir, cands))
            cands.sort(reverse=False)
            statsfolder = cands[runidx]
            print("Using stats in folder:", statsfolder)
        statsfile = os.path.join(statsfolder, statsfile)
    else:
        statsfolder = "."
    return statsfile, statsfolder


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
    statsfile, statsfolder = get_datafile(statsfile=statsfile, structure=structure, statsfolder=statsfolder,
                                          data_root_dir=data_root_directory, basedir=basedir, runidx=runidx)
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


def process_stats(stats, shift_tau_for_duration=False, shift_systime_for_end=True, rolling_window=1000):
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
    duration_key = 'duration' if shift_tau_for_duration else 'tau'
    ## Pandas DataFrames are a set of Series. Pandas series can
    df = pd.DataFrame(stats)
    ## In old versions, I would collect tau and system_time *after* performing the reaction.
    ## This made sense because then system_time is the time when the state *starts*.
    ## Tau, is the duration of the *previous* state.
    if shift_tau_for_duration:
        df['duration'] = df['tau'].shift(-1)
    if shift_systime_for_end:
        try:
            df['system_time_end'] = df['system_time'].shift(-1)
        except KeyError as e:
            print(type(e), e)
            print("df.keys():")
            pprint(df.keys())
            print(df.head())
            raise e

    # df['duration_s'] = pd.to_timedelta(df['tau'], unit="s")
    # df.set_index('duration_s', drop=True, append=False, inplace=True)

    df['f_hybridized_domains'] = df['n_hybridized_domains'] / df['N_domains']
    df['f_stacked_ends'] = df['n_stacked_ends'] / df['N_domains'] / 2
    df['f_partially_hybridized_strands'] = df['n_partially_hybridized_strands'] / df['N_strands']
    df['f_fully_hybridized_strands'] = df['n_fully_hybridized_strands'] / df['N_strands']
    df['duration_cum'] = df[duration_key].cumsum() # should equal system_time_end?
    avg_trace_cols = ('n_hybridized_domains', 'f_hybridized_domains',
                      'n_stacked_ends', 'f_stacked_ends',
                      'n_partially_hybridized_strands', 'f_partially_hybridized_strands',
                      'n_fully_hybridized_strands', 'f_fully_hybridized_strands',
                     )
    for col in avg_trace_cols:
        col_weighted = col + "_weighted"
        col_weighted_cum = col_weighted + '_cum'
        col_ave = col + "_avg"
        df[col_weighted] = df[col] * df[duration_key]     # is value * dt
        df[col_weighted_cum] = df[col_weighted].cumsum()
        df[col_ave] = df[col_weighted_cum] / df['duration_cum'] # ave is int(value * dt) / int(dt)

    # Set index, to make use of built-in time functions.
    # Note: You must make the series a time series, check index.is_all_dates Series attribute.
    # drop=True means you cannot access the column except as an index.
    # append=True means you get multiple indices.
    # inplace=True, otherwise it returns a copy and leave current DataFrame as it is.
    # df.set_index('system_time', drop=True, append=False, inplace=True)
    # df.set_index('duration_cum', drop=True, append=False, inplace=True)


    # df['duration_rolling'] = pd.rolling_sum(df[duration_key], rolling_window, center=True, min_periods=1)
    # df['f_hybridized_domains_rolling'] = pd.rolling_sum(
    #     df['f_hybridized_domains_weighted'], rolling_window, center=True, min_periods=1) / df['duration_rolling']
    # df['f_hybridized_domains_rolling'] = pd.rolling_mean(
    #     df['f_hybridized_domains'], rolling_window, center=True, min_periods=1)
    # df['f_hybridized_domains_rolling'] = pd.expanding_mean(df['f_hybridized_domains'])  # Works well, actually.
    ## Using expanding is more or less the same as the cumulative versions I make myself above:
    # df['duration_rolling'] = pd.expanding_sum(df[duration_key])
    # df['f_hybridized_domains_rolling'] = pd.expanding_sum(df['f_hybridized_domains_weighted']) / df['duration_rolling']
    # Exponentially-weighted moving average:
    # df['f_hybridized_domains_rolling'] = pd.ewma(df['f_hybridized_domains'], halflife=50)  #
    # I probably need to make a window manually...

    # Custom windows:
    # df['duration_rolling'] = pd.rolling_window(
    #     df[duration_key], win_type='triang', window=rolling_window, min_periods=0, mean=False)
    # df['f_hybridized_domains_rolling'] = pd.rolling_window(
    #     df['f_hybridized_domains_weighted'], win_type='triang', window=rolling_window, min_periods=0, mean=False)


    print("Done processing stats!")
    return df
