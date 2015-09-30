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

Module for interacting with Cytoscape.


"""

from pandas.core.frame import DataFrame
import logging
logger = logging.getLogger(__name__)

#import requests
#import py2cytoscape
from py2cytoscape.data import cyrest_client

# Local/relative imports
from .live_visualizer_base import LiveVisualizerBase

def edge_suid_to_names_mapper(frame, directed=False):
    """ Given a data frame, return {suid: (node-source, node-target)} map dict. """
    if directed:
        return {suid: (source, target)
                for suid, source, target in zip(frame.index, frame['source'], frame['target'])}
    else:
        return {suid: frozenset((source, target))
                for suid, source, target in zip(frame.index, frame['source'], frame['target'])}

def edge_names_to_suid_mapper(frame):
    """ Given a data frame, return {(node-source, node-target): suid} map dict. """
    return {v: k for k, v in edge_suid_to_names_mapper(frame).items()}


class CytoscapeStreamer(LiveVisualizerBase):
    """
    Draw graph live using Cytoscape with cyREST plugin.
    """

    def __init__(self, config):
        super().__init__(config)
        host = config.get('visualization_host', '127.0.0.1')
        port = config.get('visualization_port', 1234)
        self.client = cyrest_client.CyRestClient(ip=host, port=port)
        self.network = None     # CyNetwork object for interacting with the Cytoscape graph (network)
        #self.name_to_suid = {} # edit: for now, we store nodes and edges separately. They could probably be joined.
        # While we are debugging, keep track of deletions:
        self.id_key = 'SUID'


    def initialize_graph(self, graph, reset=True):
        """ Initialize the visualization based on graph. """
        if reset:
            self.client.session.delete()
        name = self.config.get('vizualisation_network_name')
        coll = self.config.get('vizualisation_network_collection')
        if graph is None:
            # create args:
            # suid: network id, name: graph name, collection: (?)
            # data (a dict with form {'data': {...}, 'elements': {'nodes': [...], 'edges': [...]}} )
            self.network = self.client.network.create(name=name, collection=coll)
        else:
            # Create create from networkx graph
            self.network = self.client.network.create_from_networkx(graph, collection='Generated by NetworkX')

        # Register the current nodes:
        self.register_new_nodes(self.network.get_node_table(format='cytoscapejs')) # default format is'dataframe'
        self.register_new_edges(self.network.get_edge_table(format='cytoscapejs'))

    #
    #def add_directed_edge(self, source, target, interaction):
    #    """ Add a single directed edge from source to target with <interaction>. """
    #    # py2cytoscape.data.cynetwork is a little odd:
    #    # .add_edge() returns a 1-element list of [{'SUID': 5455, 'source': 4875, 'target': 4876}]
    #    # .add_edges() returns a pandas dataframe.
    #    new_edge_ids = self.network.add_edge(source=self.node_name_to_suid[source],
    #                                         target=self.node_name_to_suid[target],
    #                                         interaction=interaction, directed=True)
    #    self.register_new_edges(new_edge_ids, directed=True)
    #
    #def add_directed_edge_both_ways(self, source, target, interaction):
    #    """ Add a single directed edge from source to target with <interaction>. """
    #    new_edges = [(self.node_name_to_suid[source], self.node_name_to_suid[target], interaction),
    #                 (self.node_name_to_suid[target], self.node_name_to_suid[source], interaction)]
    #    new_edge_ids = self.network.add_edges(new_edges, dataframe=False) # arg 'dataframe' might change to 'format'
    #    self.register_new_edges(new_edge_ids, directed=True)


    def add_node(self, node_name, attributes=None):
        """ Add a single node. """
        new_node_ids = self.network.add_nodes([node_name], dataframe=False)
        self.register_new_nodes(new_node_ids)
        return new_node_ids


    def add_nodes(self, node_name_list, attributes=None):
        """ Add multiple nodes. """
        new_node_ids = self.network.add_node(node_name_list, dataframe=False)
        self.register_new_nodes(new_node_ids)
        return new_node_ids


    def delete_node(self, node_name):
        """ Delete a single node from the graph. """
        node_id = self.node_name_to_suid.pop(node_name)
        try:
            self.network.delete_node(node_id)
        except Exception as e:
            print("Error deleting node %s: %s" % (node_id, e))
            # Re-insert node_name => node_id entry
            self.node_name_to_suid[node_name] = node_id
        else:
            test_name = self.node_suid_to_name.pop(node_id)
            if test_name != node_name:
                print("WARNING: Mapping issue: node_suid_to_name[node_id] != node_name")


    def delete_nodes(self, node_names_list):
        for node_name in node_names_list:
            self.delete_node(node_name)


    def add_edge(self, source, target, interaction, directed, bidirectional, attributes=None):
        """
        :param source:
        :param target:
        :param interaction:
        :param directed:
        :param bidirectional: If True, and directed is also True, add two directed edges,
                              one from source to target and another from target to source.
        """
        # For Cytoscape, we need to convert
        edges = [{'source': self.node_name_to_suid[source],
                  'target': self.node_name_to_suid[target],
                  'interaction': interaction,
                  'directed': directed}]
        if directed and bidirectional:
            edges += [{'source': self.node_name_to_suid[target],
                       'target': self.node_name_to_suid[source],
                       'interaction': interaction,
                       'directed': directed}]
        new_edge_ids = self.network.add_edges(edges, dataframe=False)
        self.register_new_edges(new_edge_ids, directed=directed)
        return new_edge_ids


    def add_edges(self, edges, directed, attributes=None):
        """
        Perform some translation from nascent directive to py2cytoscape directive.
        edges must be a list of dicts as:
            [{'source': <name>, 'target': <name>, 'interaction': <str>, 'directed': <bool>}, ...]
        This method will take care of mapping names->SUIDs before submitting it
        (this happens in-place, so you may want to make a copy of edges if you need it for anything important)
        """
        for edge in edges:
            edge['source'] = self.node_name_to_suid[edge['source']]
            edge['target'] = self.node_name_to_suid[edge['target']]
        new_edge_ids = self.network.add_edges(edges, dataframe=False)
        self.register_new_edges(new_edge_ids, directed=directed)  # Consider being able to add directed here!
        return new_edge_ids


    def delete_edge(self, source, target, directed=True):
        """
        Delete a single edge from the graph.
        :param source:, :param target: the source and target which edge is connecting.
        """
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
            return False
        else:
            test_name = self.edge_suid_to_names.pop(edge_id)
            self.deleted_name_to_suid.append({key: edge_id})
            self.deleted_suid_to_name.append({edge_id: key})
            if test_name != key:
                print("WARNING: Mapping issue: node_suid_to_name[node_id] != node_name")
            return True


    def delete_edges(self, edges):
        """ Delete all edges in :param edges:. """
        all_ok = 1
        for edge in edges:
            try:
                source, target = edge['source'], edge['target']
            except KeyError:
                source, target = edge['nodes']
            directed = edge['directed']
            all_ok *= self.delete_edge(source, target, directed)
        return bool(all_ok)



    #def edge_names_key(self, source, target, directed=True):
    #    if directed:
    #        return (source, target)
    #    else:
    #        return frozenset((source, target))