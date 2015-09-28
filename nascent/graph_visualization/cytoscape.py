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
        #self.config = config
        #self.graph_type = config.get('visualization_graph_type', '5p3p')
        #self.directed_graph = config.get('visualization_graph_directed', self.graph_type in ('5p3p',))
        # Note: For purely visualization, it doesn't matter if the graph is directed or not;
        # we are not using the visualized graph for any analysis or calculation.
        host = config.get('visualization_host', '127.0.0.1')
        port = config.get('visualization_port', 1234)
        self.client = cyrest_client.CyRestClient(ip=host, port=port)
        self.network = None
        #self.name_to_suid = {} # edit: for now, we store nodes and edges separately. They could probably be joined.
        self.node_name_to_suid = {}
        self.node_suid_to_name = {}
        self.edge_names_to_suid = {} # frozenset(source, target) for undirected graphs, (source, target) for directed.
        self.edge_suid_to_names = {}
        # While we are debugging, keep track of deletions:
        self.deleted_name_to_suid = []
        self.deleted_suid_to_name = []


    def initialize_graph(self, graph, reset=True):
        """ Initialize the visualization based on graph. """
        if reset:
            self.client.session.delete()
        name = self.config.get('vizualisation_network_name')
        coll = self.config.get('vizualisation_network_collection')
        self.network = self.client.network.create(name=name, collection=coll) # also suid, data
        node_tbl = self.network.get_node_table(format='cytoscapejs') # format='dataframe' is default
        #self.node_suid_to_name = dict(zip(node_tbl.index, node_tbl['name']))
        self.node_suid_to_name = {r['SUID']: r['name'] for r in node_tbl} # 'cytoscapejs' yields node_tbl['rows']
        self.node_name_to_suid = {v: k for k, v in self.node_suid_to_name.items()}
        edge_tbl = self.network.get_edge_table(format='cytoscapejs')
        if self.directed_graph:
            self.edge_suid_to_names = {r['SUID']: frozenset((r['source'], r['target'])) for r in node_tbl}
        self.edge_suid_to_names = edge_suid_to_names_mapper(edge_tbl, directed=self.directed_graph)
        self.edge_names_to_suid = {v: k for k, v in self.edge_names_to_suid.items()}
        # We need both maps - SUID->name and name->SUID.

    def add_directed_edge(self, source, target, interaction):
        """ Add a single directed edge from source to target with <interaction>. """
        # py2cytoscape.data.cynetwork is a little odd:
        # .add_edge() returns a 1-element list of [{'SUID': 5455, 'source': 4875, 'target': 4876}]
        # .add_edges() returns a pandas dataframe.
        new_edge_id = self.network.add_edge(source=self.node_name_to_suid[source],
                                            target=self.node_name_to_suid[target],
                                            interaction=interaction, directed=True)

    def add_directed_edge_both_ways(self, source, target, interaction):
        """ Add a single directed edge from source to target with <interaction>. """
        new_edges = [(self.node_name_to_suid[source], self.node_name_to_suid[target], interaction),
                     (self.node_name_to_suid[target], self.node_name_to_suid[source], interaction)]
        new_edge_id = self.network.add_edges(source=self.node_name_to_suid[source],
                                            target=self.node_name_to_suid[target],
                                            interaction=interaction, directed=True)

    def register_new_edges(self, new_edge_ids, directed=None):
        """
        Register new edges.
        Format can be either a list of dicts or a pandas dataframe.
        """
        if directed is None:
            directed = self.directed_graph
        if isinstance(new_edge_ids, DataFrame):
            edge_suid_to_names = edge_suid_to_names_mapper(new_edge_ids)
        else:
            edge_suid_to_names = {d['SUID']: (d['source'], d['target']) if directed
                                             else frozenset((d['source'], d['target']))
                                  for d in new_edge_ids}
        edge_names_to_suid = {v: k for k, v in edge_suid_to_names.items()}
        self.edge_suid_to_names.update(edge_suid_to_names)
        self.edge_names_to_suid.update(edge_names_to_suid)

    def propagate_change(self, change):
        """
        Propagate a single state change.
        - change_type: 0 = Add/remove node, 1=Add/remove edge.
        - forming: 1=forming, 0=removing
        """
        #nodes = change['nodes']
        #dt = change['timedelta']
        if change['change_type'] == 1: # 1 = Add/delete edge
            source = change['nodes'][0]
            target = change['nodes'][1]
            if change['forming'] == 1:
                # Add edge
                interaction = change['interaction'] # backbone, hybridization or stacking
                directed = True if interaction == 'stacking' else False
                new_edge_id = self.network.add_edge(source=self.node_name_to_suid[source],
                                                    target=self.node_name_to_suid[target],
                                                    interaction=interaction,
                                                    directed=directed)
                # is a list of [{'SUID': 5455, 'source': 4875, 'target': 4876}]
                #self.edge_names_to_suid = edge_names_to_suid(new_edge_id)


            elif change['forming'] == 0:
                # Remoe edge:
                for node in change['nodes']:
                    self.network.delete_edge(node)
        elif change['change_type'] == 0: # 0 = Add/delete node(s)
            if change['forming'] == 1:
                # Add nodes
                new_node_name_to_suids = self.network.add_nodes(change['nodes'])
                self.node_name_to_suid.update(new_node_name_to_suids)
            else:
                # Delete nodes
                for node in change['nodes']:
                    suid = self.node_name_to_suid[node]
                    self.network.delete_node(suid)
        else:
            raise ValueError("change type value %s not recognized" % repr(change['change_type']))



    def propagate_changes(self, changes):
        """ Propagate a list of state changes. """
        pass


