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
    """

    def __init__(self, config):
        self.config = config
        self.outputfn = self.config.get('state_changes_fn')
        self.keep_file_open = self.config.get('dispatcher_keep_file_open')
        if self.keep_file_open:
            self.outputfd = open(self.outputfn, 'w') if self.outputfn else None
            self.state_changes_cache_size = None
            self.state_changes_cache = None
        else:
            self.outputfd = None
            self.state_changes_cache_size = self.config.get('dispatcher_cache_size')
            self.state_changes_cache = []
        self.live_graph_repr = self.config.get('vizualisation_graph_representation', '5p3p')
        self.live_streamer = None
        streamer = self.config.get('dispatcher_live_streamer')
        if streamer:
            try:
                streamer_cls = STREAMERS[streamer]
                if streamer_cls is None:
                    print("Sorry, the requested '%s' streamer is not available (error during module import)"
                          % streamer)
                else:
                    self.live_streamer = streamer_cls(config)
            except KeyError:
                print("'%s' streamer not recognized" % self.live_streamer)

        self.graph_event_methods = [self.add_node, self.delete_node, self.add_edge, self.delete_edge]
        self.graph_event_group_methods = [self.add_nodes, self.delete_nodes, self.add_edges, self.delete_edges]
        self.graph_translation = "domain-to-5p3p"
        # If just specifying the output fields, they will be written as
        #   ",".join(str(state_change[field]) for field in self.state_change_line_fields)
        self.state_change_line_fields = self.config.get('state_changes_line_fields',
                                                        ['T', 'time', 'tau', 'change_type', 'forming', 'interaction'])
        # If state_change_line_unpack_nodes is True, add all nodes at the end of the line appended to the fields above.
        self.state_change_line_unpack_nodes = self.config.get('state_changes_unpack_nodes', True)
        self.state_change_line_fmt = self.config.get('state_changes_line_fmt')


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


    def dispatch(self, state_change, multi_directive=False):
        """
        Forward the dispatch.
        Propagate a single state change. A state change directive is a single dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: (only for edge types)
            1=backbone, 2=hybridization, 3=stacking
        - multi: It might be nice to have the option to propagate multiple events with
            a shared set of T, time, dt, interaction, etc. I'm specifically thinking of
            "stacking change cycles", where we re-calculate the stacking state of all possible stacking interactions.
            Note: This is not the same as the 'multi' argument given to this method, which
            is whether state_change is a single state-change dict or a list of state-change dicts.
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.
        """
        if multi_directive is None:
            multi_directive = isinstance(state_change, (list, tuple))
        state_changes = state_change if multi_directive else [state_change]
        if self.outputfd:
            # Write state changes continuously:
            if multi_directive:
                for directive in state_changes:
                    self.outputfd.write(self.state_change_str(directive))
            else:
                self.outputfd.write(self.state_change_str(state_change))
            #logger.debug("Writing state change to open file.")
        else:
            # Add state changes to cache and write when buffer is large:
            if multi_directive:
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


    def propagate_graph_event(self, event):
        """
        Propagate a single state change. A state change directive is a single dict with keys:
        - change_type: 0 = Add/remove NODE, 1=Add/remove EDGE.
        - forming: 1=forming, 0=eliminating
        - interaction: (only for edge types)
            1=backbone, 2=hybridization, 3=stacking
        - multi: It might be nice to have the option to propagate multiple events with a single set of T, time, dt, etc.
        - nodes: a two-tuple for edge types, a node name or list of nodes for node types.

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
            for method, event_list in zip(self.graph_event_group_methods, event_groups):
                # partner methods with the corresponding group and invoke:
                if event_list:
                    method(event_list)
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

    def unpack_multi_directive(self, changes):
        """
        If you are using state change directives that pack 
        """
        ## TODO: Implement multi-node/edge state change directives.
        return changes

    def state_changes_to_graph_events(self, changes):
        """
        Will translate domain nodes state changes to 5p3p graph representation.
        Note that the lenght of the translated changes list might be different
        from the input. For instance, translating a domain hybridization to 5p3p
        format produces two "add edge" graph events,
            dom1:5p--dom2.3p and dom1:3p--dom2.5p

        # YOU PROBABLY ALSO WANT TO APPLY SOME STYLES AROUND HERE
        """
        if self.multi_directive_support:
            changes = self.unpack_multi_directive(changes)
        if self.graph_translation == "domain-to-5p3p":
            events = chain(self.translate_domain_change_to_5p3p_graph_events(change)
                           for change in changes)
        elif self.graph_translation == "domain-to-strand":
            # Requires multi-graph support, in case two strands hybridize with more than one domain,
            # or we have intra-strand hybridization, stacking, etc.
            events = [self.translate_domain_change_to_strand_graph_event(change) for change in changes]
        else:
            events = changes
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


    #######################################
    ### Live streamer event propagation ###
    #######################################

    # Q: What argument should the 'singular' methods receive?
    #    A single graph 'event' dict?
    #    An unpacked graph event?
    #    - No: Unpacking graph events from the dispatcher is the reason we have the live streamer adaptors,
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
    # Q: What argument should the 'multiple' methods receive?
    #    - event_groups?
    #    - A list of the exact same arguments that the singular form receives?
    #
    # Q: Should the singular/multiple methods add/update the graph event, for any method-specific
    #    stuff that the translator does not take care of?
    #    - No: As much as possible should be done in the translator, and it already knows about add/delete node/edge.
    #

    def add_node(self, event):
        """ Add a single node to the live streamer's graph. """
        node_name = event['nodes']
        if isinstance(node_name, (list, tuple)):
            return self.live_streamer.add_edges(node_name, event.get('attributes'))
        return self.live_streamer.add_edge(node_name, event.get('attributes'))

    def add_nodes(self, events):
        """ Add multiple nodes to the live streamer's graph. """
        node_names = [event['nodes']]
        self.live_streamer.add_nodes(events, event.get('attributes'))

    def delete_node(self, event):
        """ Delete a single node from the live streamer's graph. """
        node = event['nodes']
        self.live_streamer.delete_edge(node)

    def delete_nodes(self, events):
        """ Delete multiple nodes from the live streamer's graph. """
        for node in event['nodes']:
            self.live_streamer.delete_node(node)

    def add_edge(self, event):
        """ Add a single edge to the live streamer's graph. """
        interaction = event['interaction'] # backbone, hybridization or stacking
        directed = event.get('directed', True if interaction == 'stacking' else False)
        source, target = event['nodes']
        return self.live_streamer.add_edge(source=source, target=target,
                                           interaction=interaction, directed=directed)

    def add_edges(self, events):
        """ Add multiple edges to the live streamer's graph. """
        pass

    def delete_edge(self, event):
        """ Delete a single node from the live streamer's graph. """
        pass

    def delete_edges(self, events):
        """ Delete multiple edges from the live streamer's graph. """
        pass
