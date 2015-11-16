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

# pylint: disable=C0103,W0212

"""

Module docstring.

"""

import networkx as nx

class ConnectedMultigraph(nx.MultiGraph):
    """
    A NetworkX multigraph that makes it easy to break and merge
    connected component graphs.
    """

    def break_at_edge(self, source, target):
        """
        Break graph at edge and see if it's still connected.
        Four cases:
        Case 0: Nothing breaks off
        Case 1: Two smaller connected components.
        Case 2: One node breaks off.
        Case 3: Two separate nodes.
        returns
            Case-Int, [list of surviving graphs], [list of free nodes]
        """
        self.remove_edge(source, target)
        if len(self) == 2:
            return 3, None, [source, target]
        if len(self[source]) == 1:
            self.remove_node(source)
            return 2, self, [source]
        if len(self[target]) == 1:
            self.remove_node(target)
            return 2, self, [target]
        subgraphs = list(nx.connected_component_subgraphs(self, copy=False))
        if len(subgraphs) == 1:
            return 0, [self], None
        else:
            assert len(subgraphs) == 2
            return 1, subgraphs, None


