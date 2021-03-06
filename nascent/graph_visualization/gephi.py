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

# pylint: disable=C0103,W0142

"""

Module for interacting with Gephi through the graph-streaming plugin.

Gephi graph-streaming python libraries:
* https://github.com/panisson/pygephi_graphstreaming
* https://github.com/totetmatt/GephiStreamer
** Very ugly code, but uses requests. Also imports urllib but doesn't use it...
    - urllib.open(...) was replaced by requests.post(...))
* https://github.com/jsundram/pygephi (Jython)
** This is for making use of the gephi libraries from within Python using Jython


About the HTTP requests communication with Gephi graph-streaming plugin:
* Uses the JSONStream HTTP request interface.
* See pygephi_graphstreaming/pygephi/client.py
* Very similar to Cytoscape's REST interface:
* Send GET (or maybe POST) requests to URI in the format of
    http://127.0.0.1:8080/workspace0?operation=updateGraph
  with json-formatted data

Note: I can figure out if the gephi graph-streaming plugin actually conforms to the NetStream protocol or not.
* The Gephi graph-streaming libraries that I have seen
* The NetStream protocol specifies a "message" as consisting of:
        A 4-byte (32 bit) integer indicating the length of the rest of the message
        A Stream ID
        The event it self, e.g.
* An event consists of:
        The first byte (integer values 0-127) specifies the type of event.
        Stream ID (again, apparently)
        Time ID
* Graph event types and identifiers include:
        EVENT_ADD_NODE [0x10] -
        EVENT_DEL_NODE
        EVENT_ADD_EDGE
        EVENT_DEL_EDGE
        EVENT_STEP              Time step
        EVENT_CLEARED           Clear the graph
        EVENT_ADD_GRAPH_ATTR    Add an attribute to the graph
        EVENT_CHG_GRAPH_ATTR    Change an existing attribute on the graph
        EVENT_DEL_GRAPH_ATTR    Remove an attribute from the graph.
        EVENT_ADD_NODE_ATTR, EVENT_CHG_NODE_ATTR, EVENT_DEL_NODE_ATTR
        EVENT_ADD_EDGE_ATTR, EVENT_CHG_EDGE_ATTR, EVENT_DEL_EDGE_ATTR



Other Graph Streaming operations:
* updateGraph
* getGraph      Get information on the current graph
*

Refs:
* https://gephi.wordpress.com/tag/streaming/
* https://github.com/graphstream/gs-netstream/wiki/NetStream-Manual
* https://github.com/graphstream/gs-netstream/wiki/NetStream-Manual#the-netstream-protocol
* https://github.com/gephi/gephi/wiki/Reader%27s-Circle
* https://gephi.wordpress.com/2012/10/04/gsoc-interconnect-gephi-graph-streaming-api-and-graphstream/

General-purpose NetStream libs:
* https://github.com/hhromic/python-netstream - transport via socket or remote host
* https://github.com/graphstream/gs-netstream - has a python sender
    See e.g. NetStreamProxyGraph in python/gs_netstream/sender.py
* These libs are both a bit more involved, but offers more flexibility and performance
    and offers e.g. direct local socket connection.

"""

import math
import logging
logger = logging.getLogger(__name__)

from pygephi import client

# Local/relative imports
from .live_visualizer_base import LiveVisualizerBase
from nascent.graph_sim_nx.constants import (HYBRIDIZATION_INTERACTION, PHOSPHATEBACKBONE_INTERACTION,
                                            STACKING_INTERACTION)

DEFAULT_NODE_ATTR = {"size": 10, 'r': 1.0, 'g': 0.0, 'b': 0.0, 'x': 1, 'y': 0}
DEFAULT_EDGE_ATTR = {}




class GephiMultiGraphStreamer(LiveVisualizerBase):
    """
    Draw graph live using Gephi with graph-streaming plugin.
    """

    def __init__(self, config):
        super().__init__(config)
        # Client:
        use_websocket = config.get('visualization_use_websocket', False)
        host = config.get('visualization_host', '127.0.0.1')
        port = config.get('visualization_port', 8080)
        scheme = config.get('visualization_uri_scheme', 'ws' if use_websocket else 'http')
        workspace = config.get('visualization_workspace', 'workspace0')
        url = config.get('visualization_url', "%s://%s:%s/%s" % (scheme, host, port, workspace))
        # If client.autoflush is True, changes are submitted at once.
        # Otherwise they are aggregated in client.data and only submitted when .flush() is invoked.
        autoflush = config.get('visualization_autoflush', True)
        if use_websocket:
            # Websocket is good if you send a lot of small messages but requires third party library
            from pygephi import websocket_client
            self.client = websocket_client.GephiWsClient(url, autoflush=autoflush)
        else:
            self.client = client.GephiClient(url, autoflush=autoflush)
        self.network = self.client # For gephi, the network and client is the same
        # Styles:
        self.node_attributes = DEFAULT_NODE_ATTR.copy()
        self.node_attributes.update(config.get('visualization_node_attributes',
                                               {'size': 10, 'r': 1.0}))
        self.edge_attributes = DEFAULT_EDGE_ATTR.copy()
        self.edge_attributes.update(config.get('visualization_edge_attributes', {}))

        # Graph network:
        self.network = self.client # For pygephi, the client is a proxy for the network.
        self.id_key = 'id'   # is 'SUID' for Cytoscape, 'id' for graph-streaming
        self.directed = None


    def initialize_graph(self, graph, reset=True):
        """ Initialize the visualization based on graph. """
        if reset:
            self.client.clean()
        self.init_from_networkx(graph)


    def init_from_networkx(self, graph):
        """
        Build Gephi graph from networkx graph.
        """
        autoflush = self.client.autoflush
        self.client.autoflush = False # disable auto-flush while we build the graph
        for node, data in graph.nodes(data=True):
            self.add_node(node, attributes=data)
        directed = self.directed if self.directed is not None else graph.is_directed()
        # print("Gephi streamer: Initializing with graph, directed=%s" % (directed, ))
        ## TODO: Add explicit support for multi-graphs!
        if graph.is_multigraph():
            for source, target, key, data in graph.edges(data=True, keys=True):
                interaction = data.setdefault('interaction', PHOSPHATEBACKBONE_INTERACTION)
                self.add_edge(source, target, directed=directed, key=key, interaction=interaction, attributes=data)
        else:
            for source, target, data in graph.edges(data=True):
                interaction = data.setdefault('interaction', PHOSPHATEBACKBONE_INTERACTION)
                self.add_edge(source, target, directed=directed, interaction=interaction, attributes=data)
        self.client.flush()
        self.client.autoflush = autoflush
        print("node_name_to_suid:", self.node_name_to_suid)
        print("edge_suid_to_names:", self.edge_suid_to_names)


    def add_node(self, node_name, attributes=None):
        """ Add a single node. """
        if attributes:
            # A hack would be attrs = dict(self.edge_attributes, **attributes), but Guido doesn't like that
            attrs = self.node_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.node_attributes
        self.client.add_node(node_name, **attrs)
        self.node_name_to_suid[node_name] = node_name
        self.node_suid_to_name[node_name] = node_name


    def add_nodes(self, node_names_list, attributes=None):
        """ Add all nodes in node_names_list. """
        if attributes:
            attrs = self.node_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.node_attributes
        autoflush = self.client.autoflush
        self.client.autoflush = False # disable auto-flush while we build the graph
        for node_name in node_names_list:
            self.add_node(node_name, **attrs)
        self.client.flush()
        self.client.autoflush = autoflush


    def delete_node(self, node_name):
        """ Delete a single node from the graph. """
        node_id = node_name
        if node_name not in self.node_name_to_suid:
            print("node name", node_name, "not in node_name_to_suid.")
        try:
            self.network.delete_node(node_id)
        except Exception as e:
            print("Error deleting node %s: %s" % (node_id, e))
        self.node_name_to_suid.pop(node_name, None)
        self.node_suid_to_name.pop(node_name, None)


    def add_edge(self, source, target, directed=True, key=None, interaction=None, bidirectional=None, attributes=None):
        """ Add a single edge. Id is auto-generated as {source}-{key}-{target}. """
        ## TODO: ADD SUPPORT FOR MULTI-GRAPHS
        if attributes:
            attrs = self.edge_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.edge_attributes
        if key is None:
            if interaction is None:
                interaction = attrs.get('interaction')
            # If converting key from interaction, always use a standardized way.
            # Using key = interaction is OK as long as we don't have tiny circular constructs.
            key = interaction
        edge_id = "%s--%s-%s%s" % (source, key, (">" if directed else "-"), target)
        # if False and attributes and 'weight' not in attributes and 'len' in attributes:
        #     # 'len' is typically root squared end-to-end distance in nm
        #     attributes['weight'] = math.sqrt(1/attributes['len'])
        self.client.add_edge(edge_id, source, target, directed, **attrs)
        # register_new_edges takes a list of dicts or a single dict with
        if source not in self.node_name_to_suid or target not in self.node_name_to_suid:
            print("Warning:",
                  ("source '%s' not in node_name_to_suid" % source) if source not in self.node_name_to_suid else "",
                  ("target '%s' not in node_name_to_suid" % target) if target not in self.node_name_to_suid else "")
        if edge_id in self.edge_suid_to_names:
            print("Note: edge id", edge_id, "already in edge_suid_to_names.")
        self.register_new_edges([{self.id_key: edge_id, 'source': source, 'target': target}], directed=directed)
        return edge_id


    def add_edges(self, edges, directed, keys=None, attributes=None):
        """
        Takes a list of edges in the form of
            [{'source': <name>, 'target': <name>, 'key': ..., 'interaction': <str>, 'directed': <bool>}, ...]

        """
        if attributes:
            attrs = self.edge_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.edge_attributes
        print("Attributes:", attrs)
        autoflush = self.client.autoflush
        self.client.autoflush = False # disable auto-flush while we build the graph
        # new_edge_ids = []
        # new_edge_ids = [{self.id_key: self.add_edge(source=edge['source'],
        #                                             target=edge['target'],
        #                                             directed=edge.get('directed'),
        #                                             key=edge.get('key', key),
        #                                             attributes=edge.get('attributes', attrs)),
        #                  'source': edge['source'],
        #                  'target': edge['target'],
        #                  'key': edge.get('key')}
        #                 for edge in edges]
        # self.register_new_edges(new_edge_ids)
        # Edit: add_edge will register the edge...
        for edge in edges:
            self.add_edge(source=edge['source'], target=edge['target'],
                          directed=edge.get('directed'), key=edge.get('key'),
                          interaction=edge.get('interaction'), attributes=edge.get('attributes', attrs))
        self.client.flush()
        self.client.autoflush = autoflush


    #def make_edge_id(self, source, target, directed=True):
    #    if directed:
    #        return source+"-->"+target
    #    else:
    #        return source+"---"+target


class GephiRegularGraphStreamer(GephiMultiGraphStreamer):
    """
    Prevent making multi-graph edges.
    """
    def __init__(self, config):
        super().__init__(config)
        self.is_multigraph = False  # Force regular-not-multi graph
        self.current_edges = defaultdict(set) # (source, target) => key or interaction


    def add_edge(self, source, target, directed=None, key=None, interaction=None, bidirectional=None, attributes=None):
        """ Add a single edge. Id is auto-generated as {source}-{key}-{target}. """
        edge = (source, target) if directed else frozenset((source, target))
        key = key if key is not None else (interaction if interaction is not None else attributes.get('interaction'))
        assert key is not None
        super().add_edge(source, target, directed=directed, key=key,
                         interaction=interaction, bidirectional=bidirectional, attributes=attributes)
        self.current_edges[edge].add(key)

    def delete_edge(self, source, target, directed=None, key=None, interaction=None):
        """ Delete edge, but first: check . """
        edge = (source, target) if directed else frozenset((source, target))
        key = key if key is not None else interaction
        assert key is not None
        super().delete_edge(source, target, directed=directed, key=key, interaction=interaction)
        self.current_edges[edge].remove(key)
        if len(self.current_edges[edge]) == 0:
            del self.current_edges[edge]


class MockClient():
    autoflush = True

    def clean(self, *args, **kwargs):
        pass

    def flush(self, ):
        pass

from collections import defaultdict

class GephiRegularGraphStreamerMock(GephiRegularGraphStreamer):
    """
    Test
    """

    def __init__(self, config):
        self.config = config
        self.client = MockClient()
        self.node_attributes = {}
        self.edge_attributes = {}
        self.directed = None
        self.node_name_to_suid = {}
        self.node_suid_to_name = {}
        self.edge_names_to_suid = {}
        self.edge_suid_to_names = {}
        self.is_multigraph = False  # Force regular-not-multi graph
        self.current_edges = defaultdict(set) # (source, target) => key or interaction
        self.deleted_edges = []

    def __getattr(self, *args, **kwargs):
        print("Uncaught getattr() for %s: %s, %s" % (self, args, kwargs))

    def add_node(self, *args, **kwargs):
        pass

    def add_nodes(self, *args, **kwargs):
        pass

    def add_edge(self, source, target, directed=None, key=None, interaction=None, bidirectional=None, attributes=None):
        """ Add a single edge. Id is auto-generated as {source}-{key}-{target}. """
        edge = (source, target) if directed else frozenset((source, target))
        key = key if key is not None else (interaction if interaction is not None else attributes.get('interaction'))
        assert key is not None
        self.current_edges[edge].add(key)

    def delete_edge(self, source, target, directed=None, key=None, interaction=None):
        """ Delete edge, but first: check . """
        edge = (source, target) if directed else frozenset((source, target))
        key = key if key is not None else interaction
        assert key is not None
        self.current_edges[edge].remove(key)
        if len(self.current_edges[edge]) == 0:
            del self.current_edges[edge]


