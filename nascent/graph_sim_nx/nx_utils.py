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

# pylint: disable=W0142,C0103,C0301,W0141

"""
Module with various NetworkX utility functions.

"""

import networkx as nx


def draw_graph_and_save(g, outputfn, pos=None, layout=nx.spring_layout, **kwargs):
    """ Primitive drawing of a networkx graph. """
    try:
        import matplotlib
        matplotlib.use('agg') # agg for png; 'svg' for svg, 'pdf' for pdf, etc.
        from matplotlib import pyplot
    except ImportError:
        print("Matplotlib not available; cannot plot and save graph...")
        return
    kwargs.setdefault("with_labels", True)
    if layout and pos is None:
        if layout == "graphviz":
            # Use graphviz for layout; requires pygraphviz
            try:
                # pos = nx.graphviz_layout(g)
                pos = nx.drawing.nx_agraph.graphviz_layout(g)
            except ImportError as e:
                print(repr(e))
                # It seems (at least some versions), nx.graphviz_layout is determined automatically.
                # It may be networkx.drawing.nx_pydot.graphviz_layout or networkx.drawing.nx_agraph.graphviz_layout
                try:
                    print("Trying pydot graphviz interface...")
                    pos = nx.drawing.nx_pydot.graphviz_layout(g)
                except (ImportError, AttributeError) as e:
                    # nx_pydot is only present if pydot is "sorta" available.
                    print(repr(e))
                    print("Could not find an interface to graphviz, aborting...")
                    return
        else:
            pos = layout(g)
    # Plot each graph on its own figure:
    # f = pyplot.figure(figsize=(12, 8))
    # Edit: Figures are kept after plotting. You must explicitly close them:
    # pyplot.close(f)
    # Maybe use hold=off instead?
    nx.draw(g, pos=pos, hold=False, **kwargs)
    pyplot.savefig(outputfn)
    print("Graph saved to file:", outputfn)
    pyplot.clf()
