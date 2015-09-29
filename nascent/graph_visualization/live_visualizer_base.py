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

# pylint: disable=C0103

"""

Module with live visualizer base class.


"""

import logging
logger = logging.getLogger(__name__)



class LiveVisualizerBase():
    """
    Base class for all live/online graph visualizer classes.
    """

    def __init__(self, config):
        self.config = config
        self.graph_type = config.get('visualization_graph_type', '5p3p')
        self.directed_graph = config.get('visualization_graph_directed', self.graph_type in ('5p3p',))
        # Note: For purely visualization, it doesn't matter if the graph is directed or not;
        # we are not using the visualized graph for any analysis or calculation.

        self.network = None
        self.node_name_to_suid = {}
        self.node_suid_to_name = {}
        # frozenset(source, target) for undirected graphs, (source, target) for directed.
        # If multigraph, this should be a dict of lists:
        self.edge_names_to_suid = {}
        self.edge_suid_to_names = {}
        self.deleted_name_to_suid = []
        self.deleted_suid_to_name = []
        self.id_key = 'SUID'


    def initialize_graph(self, graph, reset=True):
        """
        Initialize the visualization based on graph.
        If reset=True, the current graph visualization will be reset before adding graph.
        """
        raise NotImplementedError("Must be defined by subclass.")
    #
    #def propagate_change(self, change):
    #    """
    #    Propagate a single state change.
    #    """
    #    raise NotImplementedError("Must be defined by subclass.")
    #
    #def propagate_changes(self, changes):
    #    """
    #    Propagate a list of state changes.
    #    """
    #    raise NotImplementedError("Must be defined by subclass.")




    def register_new_nodes(self, new_node_ids):
        """
        Register new nodes.
        :param new_edge_ids: Should be a list of dicts, e.g. [{'SUID': 5455, 'name': 'a', ...],
        """
        node_suid_to_name = {row[self.id_key]: row['name'] for row in new_node_ids}
        node_name_to_suid = {v: k for k, v in node_suid_to_name.items()}
        self.node_suid_to_name.update(node_suid_to_name)
        self.node_name_to_suid.update(node_name_to_suid)


    def register_new_edges(self, new_edge_ids, directed=None):
        """
        Register new edges.
        :param new_edge_ids: Should be a list of dicts, e.g. [{'SUID': 5455, 'source': 4875, 'target': 4876}, ...],
                             or a pandas dataframe.
        """
        if directed is None:
            directed = self.directed_graph
        #if isinstance(new_edge_ids, DataFrame):
        #    edge_suid_to_names = edge_suid_to_names_mapper(new_edge_ids, directed=self.directed_graph)
        #else:
        if directed is True:
            edge_suid_to_names = {d[self.id_key]: (d['source'], d['target']) for d in new_edge_ids}
        elif directed is False:
            edge_suid_to_names = {d[self.id_key]: frozenset((d['source'], d['target']))
                                  for d in new_edge_ids}
        else:
            # directed might be a list of different directed values, one for each node
            edge_suid_to_names = {d[self.id_key]: (d['source'], d['target']) if node_directed
                                             else frozenset((d['source'], d['target']))
                                  for d, node_directed in zip(new_edge_ids, directed)}
        edge_names_to_suid = {v: k for k, v in edge_suid_to_names.items()}
        self.edge_suid_to_names.update(edge_suid_to_names)
        self.edge_names_to_suid.update(edge_names_to_suid)


    def delete_edge(self, source, target, directed=True):
        """ Delete a single node from the graph. """
        if directed:
            key, fallback = ((source, target), frozenset((source, target)))
        else:
            fallback, key = ((source, target), frozenset((source, target)))
        try:
            edge_id = self.edge_names_to_suid.pop(key)
        except KeyError:
            print("Unable to find expected edge key %s in node_name_to_suid map." % key)
            try:
                edge_id = self.edge_names_to_suid.pop(fallback)
                key = fallback
            except KeyError:
                print("- Also unable to find fallback key %s in node_name_to_suid map." % fallback)
                return
        try:
            self.network.delete_edge(edge_id)
        except Exception as e:
            print("Error deleting node %s: %s" % (edge_id, e))
            # Re-insert node_name => node_id entry
            self.edge_names_to_suid[key] = edge_id
        else:
            test_name = self.edge_suid_to_names.pop(edge_id)
            self.deleted_name_to_suid.append({key: edge_id})
            self.deleted_suid_to_name.append({edge_id: key})
            if test_name != key:
                print("WARNING: Mapping issue: node_suid_to_name[node_id] != node_name")

