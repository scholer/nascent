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

matplotlib_initialized = False

def init_matplotlib():
    global matplotlib_initialized
    if not matplotlib_initialized:
        import matplotlib
        matplotlib.use('agg') # agg for png; 'svg' for svg, 'pdf' for pdf, etc.
        matplotlib_initialized = True
    return matplotlib_initialized


def layout_graph(g, pos=None, layout=nx.layout.spring_layout, save_pos_as_attr=False, **kwargs):
    """
    Layout graph using networkx :layout: function.
    Arguments:
    :g: Graph to layout.
    :pos: Initial node positions.
    :layout: function or str.
    :save_pos_as_attr: If True, will invoke nx.set_node_attributes(g, 'pos', pos) after run.
    For graphviz layout, to use existing positions as starting positions for layout calculation,
    save node pos as node attribute before invoking layout function.
    nx.spring_layout arguments:
        G, dim=2, k=None, pos=None, fixed=None, iterations=50, weight='weight', scale=1.0, center=None
    """
    if pos is None:
        pos = nx.get_node_attributes(g, 'pos')
        if len(pos) == 0:
            pos = None
    if layout == "graphviz":
        # Use graphviz for layout; requires pygraphviz
        try:
            # pos = nx.graphviz_layout(g)
            pos = nx.drawing.nx_agraph.graphviz_layout(g)  # same as pygraphviz_layout(...)
            # Invokes pygraphviz.agraph.Agraph.layout(prog='neato', args=''), which runs graphviz as subprocess
            # and saves graphviz output as via graphviz.agread(...)
        except ImportError as e:
            print(repr(e))
            # It seems (at least some versions), nx.graphviz_layout is determined automatically.
            # It may be networkx.drawing.nx_pydot.graphviz_layout or networkx.drawing.nx_agraph.graphviz_layout
            try:
                print("Trying pydot graphviz interface...")
                pos = nx.drawing.nx_pydot.graphviz_layout(g) # same as pydot_layout(...)
            except (ImportError, AttributeError) as e:
                # nx_pydot is only present if pydot is "sorta" available.
                print(repr(e))
                print("Could not find an interface to graphviz, aborting...")
                raise e
    else:
        if isinstance(layout, str):
            layout = getattr(nx.layout, layout) # pylint: disable=E1101
        pos = layout(g, pos=pos, **kwargs)
    if save_pos_as_attr:
        for node, v in pos.items():
            # pos can sometimes be numpy.ndarray; make sure we save a regular python list.
            # edit: gefx format assumes list attributes to be "dynamic" data [val, start, end]; use tuple instead.
            pos[node] = tuple(v)
        nx.set_node_attributes(g, 'pos', pos)
    return pos


def draw_graph_and_save(g, outputfn, pos=None, layout=nx.spring_layout, clear_graph=True, hold=False, **kwargs):
    """
    Primitive matplotlib.pyplot drawing of a networkx graph.
    What about drawing directly with graphviz? Export to dot and run?
    - See seq_analysis
    """
    try:
        init_matplotlib()
    except ImportError as e:
        print("matplotlib not available, cannot draw graph...")
        return
    else:
        from matplotlib import pyplot

    kwargs.setdefault("with_labels", True)
    if layout and pos is None:
        try:
            pos = layout_graph(g, layout=layout)
        except ImportError:
            print(e, "- Falling back to spring layout...")
            pos = layout_graph(g, layout="spring_layout")
    # Plot each graph on its own figure:
    # f = pyplot.figure(figsize=(12, 8))
    # Edit: Figures are kept after plotting. You must explicitly close them:
    # pyplot.close(f)
    # Maybe use hold=off instead?
    nx.draw(g, pos=pos, hold=hold, **kwargs)
    if outputfn:
        pyplot.savefig(outputfn)
        print("Graph saved to file:", outputfn)
    if clear_graph:
        pyplot.clf()
    return pos


def draw_with_graphviz(G, filename, export=('dot', 'png'), interface='pydot', prog='neato'):
    """
    :G:         The graph to draw.
    :filename:  Filename basis (without the file extension).
    :export:    A tuple/list of file formats to export to
    :prog:      The graphviz layout program to use, typically 'neato' or 'dot'.
    """
    if filename.endswith(".png"):
        filename = filename[:-4]
    if interface == 'pydot' or G.is_multigraph():
        ## Use the pydot interface to graphviz:
        P = nx.to_pydot(G)
        if 'dot' in export:
            P.write_dot(filename + ".dot") # same as P.write(filename, format='dot', prog=None). Default prog is 'dot'.
        if 'png' in export:
            P.write_png(filename + ".png", prog=prog) # returns True on success.
    else:
        # Note: Gohlke's pygraphviz for Windows doesn't seem to support strict=False.
        ## Use the pygraphviz Agraph interface to graphviz:
        A = nx.to_agraph(G)
        A.write(filename + ".dot")
        try:
            A.draw(path=filename + ".png", format='png', prog=prog)
        except OSError as e:
            err, is_ok = e, False
        else:
            is_ok = True
