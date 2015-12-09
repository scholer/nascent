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

General stats plotting module.


"""


def load_pyplot(backend='agg'):
    """ Load matplotlib, set backend, and load pyplot. """
    import matplotlib
    if backend:
        matplotlib.use(backend)
    from matplotlib import pyplot
    return pyplot


def plot_tot_vs_time(stats, fields=('f_hybridized_domains',), add_average=True, backend="agg", #run_prefixes=None,
                     ax=None, figsize=None, filename=None, labels=None, linestyles=None, colors=None, markers=None,
                     **kwargs):
    """
    Plot "total" stats versus time. (As opposed to plotting per-domain stats).
    :stats: Pandas DataFrame (old: list of dicts) with total stats vs time.
    Be careful:
    tau and sim_system_time refers to the stats in the *previous row*,
    i.e. tau is the duration of the previous state;
    sim_system_time is the time when the state *ended*.

    """
    pyplot = load_pyplot(backend)
    # if not isinstance(stats, list):
    #     stats = [stats]
    # if isinstance(stats[0], list):
    #     stats = pd.DataFrame(stats)
    if add_average:
        fields = [field for tup in ((f, f+'_avg') for f in fields) for field in tup]
    pltargs = {'figsize': figsize} if ax is None else {'ax': ax}
    # for runstats, runprefix in
    customization = [('label', labels),
                     ('ls', linestyles),
                     ('color', colors),
                     ('marker', markers)]
    for i, field in enumerate(fields):
        for key, arg in customization:
            if arg:
                pltargs[key] = arg[i]
        if field not in stats:
            print("field %s not in stats; existing keys: %s" % (field, stats.keys()))
        if 'label' in pltargs and pltargs['label'] is None:
            del pltargs['label']
        pltargs.update(kwargs)
        print("pltargs:", pltargs)
        ax = stats.plot(x="sim_system_time", y=field, **pltargs)
        pltargs = {'ax': ax}
    pyplot.xlabel(kwargs.pop("xlabel", "Simulation time / s"))
    if filename:
        #pyplot.savefig(os.path.join(statsfolder, "f_stacked_ends_ave_vs_time.png"))
        pyplot.savefig(filename)
        print("Plot with fields %s saved to file" % (fields,), filename)
    return ax
