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

This module uses websocket protocol to communicate with Gephi.
The Gephi graph-streaming plugin supports two communication protocols:
 1. Standard HTTP requests with JSON-encoded data.
 2. Communication via the websocket protocol using Graph-Stream encoded data.

The websocket-based option has the following advantages:
* Websocket connections are kept alive, unlike standard HTTP requests.
* There is no overhead after the hand-shake is made (which is done with HTTP).

This means that for continious data communication of small events,
websocket is 1-2 orders of magnitude faster than standard HTTP requests.

However, using websocket has two draw-backs:
 1. We don't have any standard python libraries, such as pygephi_graphstreaming.
 2. I'm unfamiliar with using websocket modules
 3. We have to encode the data in graph-strem format rather than just
    making a dict and encoding it to json.
    Edit: It seems it the Gephi graph-streaming plugin actually does
    communicate in JSON-encoded format rather than the condensed GraphStream format.
    In fact, I can probably just modify the existing python libraries to support
    websocket communication.


Blog posts:
* http://matthieu-totet.fr/Koumin/2013/09/16/gephi-and-the-streaming-module/
* http://matthieu-totet.fr/Koumin/2014/06/15/lets-play-gephi-streaming-api-the-hidden-websocket/

Python websocket libraries:
* AutoBahn - https://github.com/tavendo/AutobahnPython
    https://pypi.python.org/pypi/autobahn/
    Twisted and asyncio, CPython, PyPy and Jython
* ws4py - https://github.com/Lawouach/WebSocket-for-Python
    https://pypi.python.org/pypi/ws4py/ - pypi entry not updated since 2014 (all others updated at least summer 2015)
    Uses the CherryPy framework, can also use gevent greenlets or Tornado.iostream
    Lots of examples!
* Websocket-Client - https://github.com/liris/websocket-client
    https://pypi.python.org/pypi/websocket-client/
    Python 2.7 and 3.x
    Has a WebSocketApp for easy setup.
    "This version support only hybi-13" (but I think this is just the standard Websocket RFC specification)
    Uses the standard socket library
    Uses threading module for pinging server (not important)
* Websockets - https://github.com/aaugustin/websockets
    https://pypi.python.org/pypi/websockets/
    Uses asyncio - requires python 3.4 or 3.3 with asyncio module.

Graph-stream protocol data-encoding libraries:



Other Gephi graph-streaming python libraries:
* https://github.com/totetmatt/gephi.stream.graph.websocket (JavaScript)
* https://github.com/totetmatt/graphStreamServer (Not websocket)
** Messages still seems to be in JSON format?

About the graph-stream protocol:
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

import logging
logger = logging.getLogger(__name__)

import pygephi

# Local/relative imports
from .live_visualizer_base import LiveVisualizerBase

DEFAULT_NODE_ATTR = {"size": 10, 'r': 1.0, 'g': 0.0, 'b': 0.0, 'x': 1, 'y': 0}
DEFAULT_EDGE_ATTR = {}

class GephiGraphStreamerOld(LiveVisualizerBase):
    """
    Draw graph live using Gephi with graph-streaming plugin.
    """

    def __init__(self, config):
        super().__init__(config)
        # Client:
        host = config.get('visualization_host', '127.0.0.1')
        port = config.get('visualization_port', 8080)
        scheme = config.get('visualization_uri_scheme', 'http')
        workspace = config.get('visualization_workspace', 'workspace0')
        url = config.get('visualization_url', "%s://%s:%s/%s" % (scheme, host, port, workspace))
        # If client.autoflush is True, changes are submitted at once.
        # Otherwise they are aggregated in client.data and only submitted when .flush() is invoked.
        autoflush = config.get('visualization_autoflush', True)
        self.client = pygephi.client.GephiClient(url, autoflush=autoflush)
        self.network = self.client # For gephi, the network and client is the same
        # Styles:
        self.node_attributes = DEFAULT_NODE_ATTR.copy()
        self.node_attributes.update(config.get('visualization_node_attributes', {}))
        self.edge_attributes = DEFAULT_EDGE_ATTR.copy()
        self.edge_attributes.update(config.get('visualization_edge_attributes', {}))

        # Graph network:
        self.network = self.client # For pygephi, the client is a proxy for the network.
        self.id_key = 'id'   # is 'SUID' for Cytoscape, 'id' for graph-streaming


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
        directed = graph.is_directed()
        ## TODO: Add explicit support for multi-graphs!
        for source, target, data in graph.edges(data=True):
            data.setdefault('interaction', 'b')
            self.add_edge(source, target, directed=directed, attributes=data)
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
        """ Add a single edge. Id is auto-generated as source-target"""
        if attributes:
            attrs = self.edge_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.edge_attributes
        if interaction is None:
            interaction = attrs.get('interaction', '-')
        if key is None:
            key = interaction
        # ID vs KEY: An ID identifies a node/edge uniquely within the entire graph.
        # Key is used in multi-graphs to identify one of several edges connecting the same nodes:
        # (source1, target1, key=1) vs (source1, target1, key=2).
        edge_id = source + "--" + key + ("->" if directed else "--") + target
        # e.g. A--h--a or B--s->b
        self.client.add_edge(edge_id, source, target, directed, **attrs)
        # register_new_edges takes a list of dicts or a single dict with
        self.register_new_edges([{self.id_key: edge_id, 'source': source, 'target': target}], directed=directed)
        return edge_id


    def add_edges(self, edges, directed, keys=None, attributes=None):
        """
        Takes a list of edges in the form of
            [{'source': <name>, 'target': <name>, 'interaction': <str>, 'directed': <bool>}, ...]

        """
        if attributes:
            attrs = self.edge_attributes.copy()
            attrs.update(attributes)
        else:
            attrs = self.edge_attributes
        print("Attributes:", attrs)
        default_interaction = attrs.get('interaction', 'h')
        autoflush = self.client.autoflush
        new_edge_ids = []
        self.client.autoflush = False # disable auto-flush while we build the graph

        for edge in edges():
            interact_str = ("--%s->" if edge.get('directed', directed) else "--%s--") % \
                           edge.get('interaction', default_interaction)
            edge[self.id_key] = edge['source'] + interact_str + edge['target']
            edge.update(attrs) # add default/global attributes
            if 'attributes' in edge:
                edge.update(edge.pop('attributes')) # add edge-specific attributes
            self.add_edge(**edge) # add_edge(id, source, target, directed=True, **attributes)
        self.register_new_edges(new_edge_ids)
        self.client.flush()
        self.client.autoflush = autoflush


    #def make_edge_id(self, source, target, directed=True):
    #    if directed:
    #        return source+"-->"+target
    #    else:
    #        return source+"---"+target
