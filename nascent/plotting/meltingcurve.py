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

Module for plotting melting curves (or, in generaal, some value representing state versus temperature).



"""







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
