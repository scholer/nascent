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

# pylint: disable=C0103,W0142

"""

Main dispatcher module.


The dispatcher is resposible for:
 (1) Receiving state changes.
 (2) Writing state changes to file.
 (3) Propagating the state change to the live graph, by:
    (a) Translating the state change directive to a graph event,
        according to the live graph view.
    (b) Invoking the correct graph-updating methods for the events.

Regarding (3): The graph_visualization objects are agnostic to "domains", "strands", etc.
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


"""

from itertools import chain
import logging
logger = logging.getLogger(__name__)

# Nascent imports:
from nascent.graph_sim_nx.domain import Domain
from nascent.graph_visualization.graph_utils import directed_for_all_edges

try:
    from ..graph_visualization.cytoscape import CytoscapeStreamer
except ImportError as e:
    print(e, "(Streaming to Cytoscape will not be available)")
    CytoscapeStreamer = None
try:
    from ..graph_visualization.gephi import GephiGraphStreamer
except ImportError as e:
    print(e, "(Streaming to Gephi (graph-streaming plugin) will not be available)")
    GephiGraphStreamer = None


STREAMERS = {'cytoscape': CytoscapeStreamer,
             'gephi': GephiGraphStreamer}


class StateChangeDispatcher():
    """
    State change dispatcher class.

    The usual use-case for a single state_change directive is:

    dispatcher.dispatch(state_change)
        .outputfn.write(.state_change_str(state_change)) # or use cache
        .self.state_changes_to_graph_events(state_changes)
        --> [translate_method(change), ...] => graph_events
        .propagate_graph_events(graph_events):
        --> method_idx = event['change_type']*2 + event['forming']
            (group graph events, if needed)
            .graph_event_methods[method_idx] => method (or *_group_methods if multiple events)
            method(event), e.g.
            add_node(event):
            -->
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
        self.multi_directive_support = config.get('dispatcher_multi_directive_support')
        self.graph_translation = config.get('dispatcher_graph_translation', "domain-to-5p3p")
        self.live_graph_repr = config.get('livestreamer_graph_representation', '5p3p')
        streamer = config.get('dispatcher_livestreamer', 'cytoscape')

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
                    print("Sorry, the requested '%s' streamer is not available" % streamer,
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
        Initialize the graph visualization and prepare it for streaming.
        TODO: Implement this!
        """
        if self.live_graph_repr == 'strand':
            translated_graph = graph
        elif self.live_graph_repr == 'domain':
            translated_graph = graph
        elif self.live_graph_repr == '5p3p':
            translated_graph = graph
        else:
            raise ValueError("self.live_graph_repr = '%s' is not a recognized representation.")
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
        Forward the dispatch.
        Propagate a single state change. A state change directive is a single dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: (only for edge types)
            1=backbone, 2=hybridization, 3=stacking
        - time
        - tau (or should we use 'dt'?)
        - T   (the current system temperature)
        - multi: It might be nice to have the option to propagate multiple events with
            a shared set of T, time, dt, interaction, etc. I'm specifically thinking of
            "stacking change cycles", where we re-calculate the stacking state of all possible stacking interactions.
            Note: This is not the same as the 'multi' argument given to this method, which
            is whether state_change is a single state-change dict or a list of state-change dicts.
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
        Note: Unlike a state_change dict, a graph event should not have a 'multi' directive.
        A graph event only relates to things that the graph software can understand.
        """
        if directive_is_list is None:
            directive_is_list = isinstance(state_change, (list, tuple))
        state_changes = state_change if directive_is_list else [state_change]
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
        graph_events = self.state_changes_to_graph_events(state_changes)
        if self.live_streamer:
            if len(graph_events) == 1:
                self.propagate_graph_event(graph_events[0])  # or maybe 'stream_change' ?
            else:
                self.propagate_graph_events(graph_events)


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
                method(events)
            else:
                for event in events:
                    method_idx = event['change_type']*2 + event['forming']
                    method = self.graph_event_methods[method_idx]
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
        Will translate domain nodes state changes to 5p3p graph representation.
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
            events = chain(self.translate_domain_change_to_5p3p_graph_events(change)
                           for change in state_changes)
        elif self.graph_translation == "domain-to-strand":
            # Requires multi-graph support, in case two strands hybridize with more than one domain,
            # or we have intra-strand hybridization, stacking, etc.
            events = [self.translate_domain_change_to_strand_graph_event(change) for change in state_changes]
        else:
            events = state_changes
        return events


    def translate_domain_change_to_5p3p_graph_events(self, change):
        """ Translate a single state change to 5p3p graph format. """
        event = change
        # The 'multi' key can be used to create several nodes/edges with a single "change",
        # sharing values for time, tau, T, etc. (Compressed to a line line).
        # This is NOT the same as using dispatch([list of state_change directives], multi=True)
        # which would write multiple lines to the file.
        # multi = change['multi'] # uncomment to enable support for multi-node/edge directives.
        if change['change_type'] in (0, 'node event'):
            # Make a single "add_nodes/delete_nodes" event?
            # Or expand to many add_node/delete_nodes events, one for each node?
            event['nodes'] = [nname for nname in
                              [(str(domain)+":5p", str(domain)+":3p")
                               for domain in change['nodes']]]
            return [event]
            # Add style to event?
        elif change['change_type'] in (1, 'edge event'):
            dom1, dom2 = change['nodes']
            if change['interaction'] in (1, 'backbone', 3, 'stacking'):
                # The 3p of the first (5p-most) domain connects to the 5p end of the second (3p-most) domain
                event['nodes'] = [str(dom1)+":3p", str(dom2)+":5p"]
                return [event]
            else:
                print("Unrecognized interaction for change %s" % change)
            if change['interaction'] in (2, 'hybridization'):
                # The 3p of the first domain connects to the 5p end of the second (3p-most) domain,
                # AND connect the 5p end of the first to the 3p of the second.
                # We have already done the former above, only do the latter:
                event2 = change.copy()
                event['nodes'] = [str(dom1)+":3p", str(dom2)+":5p"]
                event2['nodes'] = [str(dom1)+":5p", str(dom2)+":3p"]
                return [event, event2]


    def translate_domain_change_to_strand_graph_event(self, change):
        """
        Assume the default naming format of
            sname#suid:domain#duid:5p/3p
        """
        if isinstance(change['nodes'], Domain):
            change['nodes'] = [domain.strand.instance_name for domain in change['nodes']]
        else:
            # Assume it is a string that we must re-write:
            change['nodes'] = [node.split(":")[0] for node in change['nodes']]
        return change


    ## Live streamer event propagation methods: ##

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
        directed = event.get('directed', True if interaction == 'stacking' else False)
        source, target = event['nodes']
        return self.live_streamer.add_edge(source=source, target=target,
                                           interaction=interaction,
                                           directed=directed,
                                           attributes=attributes)

    def add_edges(self, events, attributes=None):
        """ Add multiple edges to the live streamer's graph. """
        edges = [{'source': event['nodes'][0], 'target': event['nodes'][1],
                  'interaction': event['interaction'],
                  'directed':  event.get('directed', True if event['interaction'] == 'stacking' else False)}
                 for event in events]
        directed = directed_for_all_edges(edges)
        return self.live_streamer.add_edges(edges, directed, attributes)

    def delete_edge(self, event):
        """ Delete a single node from the live streamer's graph. """
        pass

    def delete_edges(self, events):
        """ Delete multiple edges from the live streamer's graph. """
        pass



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