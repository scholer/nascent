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

Main dispatcher module.


The dispatcher is resposible for:
 (1) Receiving state changes (described by domain objects).
 (2) Translating "domain" object instance nodes to string representation.
 (3) Writing state changes to file.
 (4) Propagating the state change to the live graph, by:
    (a) Translating the state change directive to a graph event,
        according to the live graph view (domain-level vs 5p3p-level).
    (b) Invoking the correct graph-updating methods for the events.

Regarding (4): The graph_visualization objects are agnostic to "domains", "strands", etc.
They only know about how to create graphs and add/remove nodes and edges on their respective
output tool. It is the dispatchers job to re-format graph initiation and state change directives
according to the current graph representation and style.
What this means is that:
 1. The dispatcher receive a state change directive, e.g.
        "domain_A:3 hybridized to domain_a:5", or
        "domain_A:4 stacked with domain_B:2" (the 3p end of A to 5p end of B)
 2. The dispatcher rewrites the directive, e.g. if the current graph repr is '5p3p',
    and edge_styles are: {'backbone': {'directed': True},
                          'stacking': {'directed': True},
                          'hybridization': {'directed': False}}
    then these directives must be translated to:
        add_undirected_edge


NOTE: Gephi doesn't support multi-graphs until v0.9 (TBR Dec 20th, 2015)
What to do until then?
Can we use regular graphs?
 - We could, but it requires a lot of extra work.

"""

from itertools import chain
import logging
logger = logging.getLogger(__name__)

# Nascent imports:
from nascent.graph_sim_nx.domain import Domain
from nascent.graph_visualization.graph_utils import directed_for_all_edges
from nascent.graph_sim_nx.constants import (HYBRIDIZATION_INTERACTION, STACKING_INTERACTION,
                                            PHOSPHATEBACKBONE_INTERACTION)
from .graph_translators import (translate_domain_graph_to_5p3p_str_graph,
                                translate_domain_graph_to_domain_str_graph,
                                translate_domain_graph_to_strand_str_graph,
                                translate_domain_change_to_5p3p_str_graph_event,
                                translate_domain_change_to_domain_str_graph_event,
                                translate_domain_change_to_strand_str_graph_event)

STREAMERS = {}
try:
    from nascent.graph_visualization.cytoscape import CytoscapeStreamer
    STREAMERS['cytoscape'] = CytoscapeStreamer
except ImportError as e:
    print(e, "(Streaming to Cytoscape will not be available)")
    CytoscapeStreamer = None
try:
    # Make sure pygephi is available on your python path
    # Use gephi_ws for websocket communication
    from nascent.graph_visualization.gephi import GephiGraphStreamer
    from nascent.graph_visualization.gephi_ws import GephiGraphStreamerWs
    STREAMERS['gephi_old'] = GephiGraphStreamer
    STREAMERS['gephi_ws'] = GephiGraphStreamerWs
except ImportError as e:
    print(e, "(Streaming to Gephi (graph-streaming plugin) will not be available)")
    GephiGraphStreamer = None


class StateChangeDispatcher():
    """
    State change dispatcher class.
    Receives domain-instance graph events, translates them to model-agnostic (str repr) events,
    and forwards them to a file and/or a live graph streamer.

    The usual use-case for a single state_change directive is:

    dispatcher.dispatch(state_change)
        .outputfn.write(.state_change_str(state_change)) # or use cache
        .self.state_changes_to_graph_events(state_changes)
        --> [translate_method(change), ...] => graph_events
        .propagate_graph_events(graph_events):
        --> method_idx = event['change_type']*2 + event['forming']
            (group graph events, if needed)
            .graph_event_methods[method_idx] => method (or *_group_methods if multiple events)
            method(event), e.g. self.add_node(event):
            --> self.live_streamer.add_node(node_name, event.get('attributes'))

    live_streamer.add/delete_node/edge():
    live_streamer.add_node(node_name, attributes=None), add_nodes(node_names_list, attributes=None)
    live_streamer.add_edge(self, source, target, directed=True, interaction=None, bidirectional=None, attributes=None)
        --> client.add_node(...)

    The state_change received by dispatch() is either a dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: 1=backbone, 2=hybridization, 3=stacking (only for edge changes)
                        Edit: Using constants values, e.g. HYBRIDIZATION_INTERACTION (='h').
        - time: Current system time (might be subjective to rounding errors).
        - tau:  Change in system time between this directive and the former.
        - T:    Current system temperature.
        - multi: If True, interpret "nodes" as a list of pairs (u₁, u₂, ...) for NODE-altering directives,
                    or (u₁, v₁, u₂, v₂, ...) for EDGE-altering directives.
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
    """

    def __init__(self, config):
        self.config = config

        ## File output and format: ##
        self.outputfn = config.get('dispatcher_state_changes_fn')
        self.keep_file_open = config.get('dispatcher_keep_file_open')
        # If just specifying the output fields, they will be written as
        #   ",".join(str(state_change[field]) for field in self.state_change_line_fields)
        self.state_change_line_fields = config.get(
            'dispatcher_state_changes_line_fields',
            ['T', 'time', 'tau', 'change_type', 'forming', 'interaction'])
        # If state_change_line_unpack_nodes is True,
        # add all nodes at the end of the line appended to the fields above.
        self.state_change_line_unpack_nodes = config.get('dispatcher_state_changes_unpack_nodes', True)
        # For more control of the output line format, you can specify a string instead:
        self.state_change_line_fmt = config.get('dispatcher_state_changes_line_fmt')
        # If supporting multi node/edge directive, a single state_change can have multiple nodes/edes

        ## State change directies and graph visualization ##
        self.multi_directive_support = config.get('dispatcher_multi_directive_support', True)
        self.graph_translation = config.get('dispatcher_graph_translation', "domain-to-5p3p")
        # 'livestreamer_graph_representation' is replaced by 'dispatcher_graph_translation' above.
        # The "graph visualization adaptors" should be agnostic to the DNA/strand/domain model
        # self.live_graph_repr = config.get('livestreamer_graph_representation', '5p3p')
        self.auto_apply_layout = config.get('livestreamer_auto_apply_layout', 0)
        self.live_graph_layout = config.get('livestreamer_graph_layout_method', 'force-directed')
        streamer = config.get('dispatcher_livestreamer', 'gephi')

        ## Init file output stream: ##
        if self.keep_file_open:
            self.outputfd = open(self.outputfn, 'w') if self.outputfn else None
            self.state_changes_cache_size = None
            self.state_changes_cache = None
        else:
            self.outputfd = None
            self.state_changes_cache_size = config.get('dispatcher_cache_size', 100)
            self.state_changes_cache = []

        ## Init live streamer ##
        self.live_streamer = None
        if streamer:
            try:
                streamer_cls = STREAMERS[streamer]
                if streamer_cls is None:
                    print("\n[Dispatcher] Sorry, the requested '%s' streamer is not available" % streamer,
                          "(error during module import)")
                else:
                    self.live_streamer = streamer_cls(config)
            except KeyError:
                print("'%s' streamer not recognized" % self.live_streamer)

        # Graph event methods: indexed by event['change_type']*2 + event['forming']
        self.graph_event_methods = [self.delete_node, self.add_node,
                                    self.delete_edge, self.add_edge]
        self.graph_event_group_methods = [self.delete_nodes, self.add_nodes,
                                          self.delete_edges, self.add_edges]


    def init_graph(self, graph, reset=True):
        """
        Initialize the live streamer graph visualization and prepare it for streaming.
        Input argument :graph: is a domain-object networkx graph.
        If :reset: is True, the live streamer's current graph/workspace is reset
        before initializing the new graph.
        TODO: Implement this!
        """
        #if self.live_graph_repr == 'strand':
        if self.graph_translation == 'domain-to-strand':
            print("Initialization of Stand graphs not implemented.")
            translated_graph = translate_domain_graph_to_strand_str_graph(graph)
        #elif self.live_graph_repr == 'domain':
        elif self.graph_translation is None or self.graph_translation == 'domain-to-domain':
            # Still need to convert domain objects to domain strings?
            # Or do we rely on string casting in the livestreamer?
            translated_graph = translate_domain_graph_to_domain_str_graph(graph)
        #elif self.live_graph_repr == '5p3p':
        elif self.graph_translation == 'domain-to-5p3p':
            translated_graph = translate_domain_graph_to_5p3p_str_graph(graph)
        else:
            raise ValueError("self.graph_translation '%s' is not a recognized representation." %
                             self.graph_translation)
        print("Initializing streamer using graph with nodes:", sorted(translated_graph.nodes()))
        self.live_streamer.initialize_graph(translated_graph, reset=reset)


    def state_change_str(self, state_change, eol='\n'):
        """
        state_change is a tuple
        """
        if self.state_change_line_fmt:
            return self.state_change_line_fmt.format(**state_change) + eol
        else:
            values = [str(state_change[field]) for field in self.state_change_line_fields]
            if self.state_change_line_unpack_nodes:
                values += [str(node) for node in state_change['nodes']]
            return ",".join(values) + eol


    def flush_state_changes_cache(self):
        """
        Write all state changes in state_changes_cache to file.
        """
        cache_size = len(self.state_changes_cache)
        lines = (self.state_change_str(state_change) for state_change in self.state_changes_cache)
        with open(self.outputfn) as fp:
            #fp.write("".join(lines))
            fp.writelines(lines) # writelines does not add newline characters
        self.state_changes_cache = []
        logger.info("Wrote %s state changes to file %s", cache_size, self.outputfn)


    def dispatch(self, state_change, directive_is_list=False):
        """
        Forward (dispatch) a state change directive to files or live-streamers.
        Propagate a single state change. A state change directive is a single dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: (only for edge types)
            1=backbone, 2=hybridization, 3=stacking - or use the standard INTERACTION constants.
        - time
        - tau (or should we use 'dt'?)
        - T   (the current system temperature)
        - multi: It might be nice to have the option to propagate multiple events with
            a shared set of T, time, dt, interaction, etc. I'm specifically thinking of
            "stacking change cycles", where we re-calculate the stacking state of all possible stacking interactions.
            Note: This is not the same as the 'directive_is_list' argument given to this method, which
            is whether state_change is a single state-change dict or a list of state-change dicts.
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
            IMPORTANT: For stacking interactions, domains must be given in the order of the stack:
                       nodes = (dom1, dom2) means that dom1.end3p stacks with dom2.end5p !

        "State change" vs "graph event": Unlike a state_change dict, a graph event should not have a 'multi' directive.
         - A graph event only relates to things that the graph software can understand.
        Q: Should the state change nodes be exclusively Domain objects, or can they be e.g. DomainEnd instances?
        A: For now, it can only be Domain instances.
        """
        if directive_is_list is None:
            directive_is_list = isinstance(state_change, (list, tuple))
        state_changes = state_change if directive_is_list else [state_change]

        ## Write state change to file (if this dispatcher has a file configured):
        if self.outputfd:
            # Write state changes continuously:
            if directive_is_list:
                for directive in state_changes:
                    self.outputfd.write(self.state_change_str(directive))
            else:
                self.outputfd.write(self.state_change_str(state_change))
            #logger.debug("Writing state change to open file.")
        else:
            # Add state changes to cache and write when buffer is large:
            if directive_is_list:
                self.state_changes_cache.extend(state_changes)
            else:
                self.state_changes_cache.append(state_change)
            if len(self.state_changes_cache) >= self.state_changes_cache_size:
                self.flush_state_changes_cache()
        if self.config.get('dispatcher_debug_print'):
            if directive_is_list:
                print(*[self.state_change_str(directive) for directive in state_changes], sep="\n")
            else:
                print(self.state_change_str(state_change))

        ## Propagate the graph state change directive to live graph streamer, if one is configured:
        if self.live_streamer:
            # Translate the state_changes list-of-dicts to a list of "pure, domain-agnostic" str-repr graph events:
            graph_events = self.state_changes_to_graph_events(state_changes)
            if len(graph_events) == 0:
                print("Dispatcher: NO GRAPH EVENTS!")
                return
            print("len(graph_events):", len(graph_events))
            if len(graph_events) == 1:
                self.propagate_graph_event(graph_events[0])  # or maybe 'stream_change' ?
            else:
                self.propagate_graph_events(graph_events)
            if self.auto_apply_layout:
                self.live_streamer.apply_layout(self.live_graph_layout)


    def propagate_graph_event(self, event):
        """
        Propagate a single graph :param event:.
        A graph event is a single dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: (only for edge types)
            1=backbone, 2=hybridization, 3=stacking
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
        Note: Unlike a state_change dict, a graph event should not have a 'multi' directive.
        A graph event only relates to things that the graph software can understand.
        """
        #nodes = change['nodes']
        #dt = change['timedelta']
        method_idx = event['change_type']*2 + event['forming']
        try:
            method = self.graph_event_methods[method_idx]
        except KeyError:
            raise ValueError("change type value %s, forming=%s not recognized (idx=%s)" %
                             (repr(event['change_type']), event['forming'], method_idx))
        print("propagate: invoking %s(%s)" % (method, event))
        return method(event)


    def propagate_graph_events(self, events, group=True):
        """
        Propagate a list of state changes.
        :param group: If True, the changes will be grouped before dispatch.
                      If False, the changes will be invoked one by one.
                      If None, the changes will be invoked as a group if event 'change_type' and 'forming' is the same
                      for all events.
        """
        # Grouped approach:
        # 1. Group by method_idx = change['change_type']*2 + change['forming']
        # 2. Translate each group appropriately.
        # 3. For each method and changes in method_idx group, invoke the "group" variant,
        #    E.g for all method_idx=3 (change_type 'edge', forming), invoke
        #    self.live_streamer.add_edges(edges)
        if group:
            event_groups = [[] for _ in range(4)]
            for event in events:
                method_idx = event['change_type']*2 + event['forming']
                event_groups[method_idx].append(event)
            for method, event_group in zip(self.graph_event_group_methods, event_groups):
                # partner methods with the corresponding group and invoke:
                if event_group:
                    method(event_group)
        else:
            if group is None and len(set(event['change_type']*2 + event['forming'] for event in events)) == 1:
                method_idx = events[0]['change_type']*2 + events[0]['forming']
                method = self.graph_event_group_methods[method_idx]
                print("propagate evts: invoking %s(%s)" % (method, events))
                method(events)
            else:
                for event in events:
                    method_idx = event['change_type']*2 + event['forming']
                    method = self.graph_event_methods[method_idx]
                    print("propagate evt: invoking %s(%s)" % (method, event))
                    method(event)
                    # or just use:
                    # self.propagate_event(event)

    def unpack_multi_directive(self, state_change):
        """
        If you are using state change directives that pack several nodes/edges into a single state_change.
        """
        ## TODO: Implement multi-node/edge state change directives.
        if state_change['multi_directive'] != True:
            return state_change
        if state_change['change_type'] in (0, "node event"):
            changes = [dict(state_change, nodes=[node]) for node in state_change['nodes']]
        elif state_change['change_type'] in (1, "edge event"):
            node_iter = iter(state_change['nodes'])
            changes = [dict(state_change, nodes=nodes) for nodes in zip(node_iter, node_iter)]
        else:
            raise ValueError("Unrecognized 'change_type' for state_change %s" % state_change)
        return changes

    def state_changes_to_graph_events(self, state_changes):
        """
        Will translate domain nodes state changes to "strand/domain/end" objects to str representation.

        Note that the lenght of the translated events list might be different
        from the input state_changes list. For instance, translating a domain
        hybridization to 5p3p format produces two "add edge" graph events,
            dom1:5p--dom2.3p and dom1:3p--dom2.5p

        # YOU PROBABLY ALSO WANT TO APPLY SOME STYLES AROUND HERE
        # (If using Gephi; Cytoscape is using styles...)
        """
        # When supporting multi node/edge directives, a single state_change dict can have
        # multiple nodes for node-pairs under the "nodes" key.
        if self.multi_directive_support and \
            any(state_change.get('multi_directive') for state_change in state_changes):
            state_changes = chain(self.unpack_multi_directive(state_changes)
                                  if state_change.get('multi_directive')
                                  else [state_change]
                                  for state_change in state_changes)
        if self.graph_translation == "domain-to-5p3p":
            print("Translating domain change to 5p3p graph change...")
            events = [trans_change for state_change in state_changes
                      for trans_change in translate_domain_change_to_5p3p_str_graph_event(state_change)]
        elif self.graph_translation == "domain-to-strand":
            # Requires multi-graph support, in case two strands hybridize with more than one domain,
            # or we have intra-strand hybridization, stacking, etc.
            print("Translating domain change to strand graph change...")
            events = [translate_domain_change_to_strand_str_graph_event(state_change)
                      for state_change in state_changes]
        else:
            print("Translating domain change to domain graph change (str repr)...")
            events = [translate_domain_change_to_domain_str_graph_event(state_change)
                      for state_change in state_changes]
        return events


    ## Live streamer event propagation methods: ##

    # For all *node methods, event['nodes'] must be a single node, NOT a 1-element list.
    # For all *edge methods, event['nodes'] must be a 2-element list.

    def add_node(self, event, attributes=None):
        """ Add a single node to the live streamer's graph. """
        node_name = event['nodes']
        if isinstance(node_name, (list, tuple)):
            if len(node_name) == 1:
                node_name = node_name[0]
            else:
                # We actually have a list with many node names; call the .add_edges method instead
                return self.live_streamer.add_nodes(node_name, event.get('attributes'))
        return self.live_streamer.add_node(node_name, event.get('attributes'))

    def add_nodes(self, events, attributes=None):
        """ Add multiple nodes to the live streamer's graph. """
        node_names = [event['nodes'] for event in events]
        self.live_streamer.add_nodes(node_names, attributes)

    def delete_node(self, event):
        """ Delete a single node from the live streamer's graph. """
        node = event['nodes']
        self.live_streamer.delete_edge(node)

    def delete_nodes(self, events):
        """ Delete multiple nodes from the live streamer's graph. """
        for node in events['nodes']:
            self.live_streamer.delete_node(node)

    def add_edge(self, event, attributes=None):
        """ Add a single edge to the live streamer's graph. """
        interaction = event['interaction'] # backbone, hybridization or stacking
        directed = event.get('directed', interaction == STACKING_INTERACTION)
        source, target = event['nodes']
        return self.live_streamer.add_edge(source=source, target=target,
                                           directed=directed,
                                           key=event.get('key', event.get('interaction')),
                                           interaction=interaction,
                                           attributes=attributes)

    def add_edges(self, events, attributes=None):
        """
        Add multiple edges to the live streamer's graph.
        :events: is a list of graph-event dicts, each dict describing a single node or a single edge.
        """
        edges = [{'source': event['nodes'][0], 'target': event['nodes'][1],
                  'key': event.get('key'),
                  'interaction': event['interaction'],
                  'directed':  event.get('directed', event['interaction'] == STACKING_INTERACTION)}
                 for event in events]
        directed = directed_for_all_edges(edges)
        # Edges is a dict with 'source', 'target', etc keys.
        return self.live_streamer.add_edges(edges, directed=directed, attributes=attributes)

    def delete_edge(self, event):
        """ Delete a single node from the live streamer's graph. """
        directed = event.get('directed', event['interaction'] == STACKING_INTERACTION)
        source, target = event['nodes']
        return self.live_streamer.delete_edge(source=source, target=target,
                                              directed=directed,
                                              key=event.get('key', event.get('interaction')))

    def delete_edges(self, events):
        """ Delete multiple edges from the live streamer's graph. """
        for event in events:
            self.delete_edge(event)




#######################################
### Live streamer event propagation ###
#######################################

# Q: What argument should the dispatcher pass to the 'singular' livestreamer methods?
#    A single graph 'event' dict?
#    An unpacked graph event?
#    I.e. should you call
#    or self.live_streamer.add_node(node_name, ...)  (unpacked)
#    Should we unpack the graph event? (pros/cons)
#    - No: Unpacking graph events passed from the dispatcher is the reason we have the live streamer adaptors,
#        i.e. cytoscape.CytoscapeStreamer and gephi.GephiGraphStreamer.
#        These objects receive a graph event and is responsible for propagating it.
#
#    - Yes: Another way to look at it is that the adaptors are responsible for the *connection*
#    to the graph visualizer, but not for parsing graph events dicts?
#    This really depends on whether all live_streamer classes can share a single, 'unpacked' graph event API,
#    i.e. if all streamers are OK receiving
#        .add_node(node_name, attributes)
#        .add_nodes(node_names_list, attributes)
#
#        .delete_node(node_name)
#        .delete_nodes(node_names_list)
#
#        .add_edge(source, target, directed, interaction, bidirectional, attributes)
#        .add_edges([{source, target, directed, interaction, attributes}, ...] list of dicts)
#
#        .delete_edge(source, target, directed)
#        .delete_edges([{source, target, directed}, ...] list of dicts)
#
#   Conclusion: If I wanted to just pass one or multiple graph event(s) to the live_streamer,
#   there would be little need to have these functions; I could just pass it directly when I
#   do method(event) from .propagate_graph_events(graph_events).
#   In fact, you could argue that I could just move the propagate_graph_events method to the
#   live_streamer as well. It is really just a matter of *where* I want the code to reside,
#   and the translation to take place. Since the live_streamer classes are poly-morphic,
#   I think I would end up with more code duplication if the live_streamers had to take care of
#   event unpacking. As long as all the streamers can share a common API for
#   add/delete_node(s)/edge(s), then it is probably easiest to have it here.
#   That means that graph events are actually just internal data structures used by the dispatcher,
#   which is probably also a good thing.
#
# Q: What argument should the 'multiple' methods receive?
#    (1) event_groups?
#    (2) A list of the exact same arguments that the singular form receives?
# A: Option (2) - for the same reason as above.
#
# Q: Should the singular/multiple methods add/update the graph event, for any method-specific
#    stuff that the translator does not take care of?
#    - No: As much as possible should be done in the translator, and it already knows about add/delete node/edge.
#
