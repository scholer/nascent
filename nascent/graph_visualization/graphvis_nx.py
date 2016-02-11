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


from __future__ import absolute_import, print_function
import networkx as nx

# Before importing pyplot, you probably want to make sure you have a proper matplotlib backend set up
# matplotlib.use('qt4agg') 'qt4agg for most cases unless you are in a notebook.
# And, as always, don't use pyplot if you are embedding into an existing (gui) app. Use matplotlib directly.
from matplotlib import pyplot


"""

I'm currently struggling with an issue where pyplot makes my python repl unusably slow:
> python
>>> import matplotlib
>>> matplotlib.get_backend()
'agg'
>>> matplotlib.use('qt4agg')
>>> from matplotlib import pyplot
>>> # Interpreter is still fast

>>> pyplot.plot(range(10), range(10))
[<matplotlib.lines.Line2D object at 0x000000000519CB00>]
>>> # Now interpreter is very sluggish.

Investigating...
* I can get it to work reasonably well if I launch it with
    ipython console --matplotlib=qt
    >>> from matplotlib import pyplot
    >>> pyplot.plot(1)
    >>> # Still runs reasonably fast!

Note on agg backends (anti-grain geometry)
* qt4agg: only one I have available
* wxagg, gtkagg, tkagg: N/A

The 'Agg' backend is NOT an interactive backend, it is good at producing PNG files.
* http://matplotlib.org/faq/usage_faq.html#what-is-a-backend
* http://stackoverflow.com/questions/8955869/why-is-plotting-with-matplotlib-so-slow
    (no, is for speeding up rendering, not related to my problem here)


Possible better choices for high-framerate visualization:
* http://www.pyqtgraph.org/ - https://github.com/pyqtgraph/pyqtgraph
* http://code.enthought.com/projects/chaco/




"""

def init_readline(ns):
    import readline
    import rlcompleter
    comp = rlcompleter.Completer(ns)
    readline.set_completer(comp.complete)
    readline.parse_and_bind('tab: complete')



def simple_matplotlib_viz(simulator, uncomplexed_strands=False,
                          #graph_reps=('domain',),
                          rep='domain',
                          ax=None,
                          **kwargs):
    """
    Simple but expensive visualization of the complexes in simulator.
    """
    ax.hold(False)  # will clear the axes before adding new.
    U = nx.union_all(simulator.complexes)
    pos = nx.nx_agraph.graphviz_layout(U)
    nx.draw_networkx(U, pos=pos, ax=ax)
    # if you need to add more stuff, use ax.hold(True) # additive behaviour
