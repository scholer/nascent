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
Module for analysing reaction graph.

"""

from __future__ import absolute_import, print_function, division
import math
import os
import io
from datetime import datetime
from collections import defaultdict
import yaml
import json
# import networkx as nx
from networkx.classes.digraph import DiGraph
import numbers
import logging
logger = logging.getLogger(__file__)

from .constants import STACKING_INTERACTION, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION
from .utils import sequential_number_generator


def str_or_numeric(key, value):
    return isinstance(value, str) or isinstance(value, numbers.Number)

def load_reaction_graph(fn):
    with open(fn) as fp:
        graph = yaml.load(fp)
    return graph


class ReactionGraph(DiGraph):
    # Should we subclass DiGraph? Or just encapsulate a DiGraph as an attribute?
    # For now, let's just sub-class it, but in the long term it might be better to encapsulate.

    def __init__(self, data=None, params=None, **attr):
        # self.reaction_graph.graph['node'] = {'fontname': 'Courier new'}
        # self.reaction_graph.graph['edge'] = {'fontname': 'Arial',
        #                                      'fontsize': 10.0, # labelfontsize
        #                                      'len': 4.0}
        DiGraph.__init__(self, data, **attr)
        if params is None:
            params = {}
        self.params = params  # or "config" ?
        if 'reaction_graph_default_attrs' in params:
            self.graph.update(params['reaction_graph_default_attrs'])

        # A reaction can have multiple end-state when a complex is being split up:
        # (edge_attrs should be the same for edges belonging to the same reaction for intracomplex reactions
        self.endstates_by_reaction = {}  # [startstate][(reaction_spec_pair, reaction_attr)] = [list of endstates]
        self.endstates_by_reaction[0] = defaultdict(list) # Also adding the "null" node
        self.reverse_reaction_key = {} # get reaction edge key for the opposite direction.
        self.tau_cum_max = 0
        # self.reverse_reaction[(rx_spec_pair, rx_attr)] = (rev_rx_spec_pair, rev_rx_attr)
        # used to be a dict {edge_key => target} but we can have multiple targets for a single reaction.
        self.reaction_graph_complexes_directory = params.get("reaction_graph_complexes_directory")
        if self.reaction_graph_complexes_directory is not None:
            if not os.path.exists(self.reaction_graph_complexes_directory):
                print("Creating directory for complex files:", self.reaction_graph_complexes_directory)
                os.makedirs(self.reaction_graph_complexes_directory)
            assert os.path.isdir(self.reaction_graph_complexes_directory)

        self.dispatchers = []
        # File with changes to the reaction graph, e.g. new nodes/edges and node/edge updates:
        self.reaction_graph_events_file = params.get('reaction_graph_events_file')
        if self.reaction_graph_events_file is None and self.reaction_graph_complexes_directory is not None:
            self.reaction_graph_events_file = os.path.join(self.reaction_graph_complexes_directory,
                                                           "reaction_graph_eventstream.json")

        if self.reaction_graph_events_file:
            #self.reaction_graph_delta_file = open(self.reaction_graph_delta_file, 'a')
            #self.open_files.append(self.reaction_graph_delta_file)
            gs_file_dispatcher = GraphStreamFileDispatcher(self.reaction_graph_events_file)
            gs_file_dispatcher.writer.write("# New reaction graph initialized at %s\n" % datetime.now())
            print("\n\nWriting reaction_graph event stream to file: %s\n" % self.reaction_graph_events_file)
            self.dispatchers.append(gs_file_dispatcher)
        else:
            # raise ValueError("self.reaction_graph_events_file not given: ", self.reaction_graph_events_file)
            print("self.reaction_graph_events_file (%s) not given: Graph events will not be available." %
                  self.reaction_graph_events_file)


    # Same func spec as nx.Graph except added default dHdS:
    def add_node(self, n, attr_dict=None, dHdS=(0, 0), **attr):
        """ Add state node n to reaction graph. """
        assert n not in self.adj   # edge-attrs dict, {src: {tgt1: {edge1_attrs}, ...}}
        # print("ReactionGraph.add_node(%s, %s, %s, %s)" % (n, attr_dict, dHdS, attr))
        if attr_dict is None:
            attr_dict = attr
        elif attr:
            attr_dict.update(attr)
        # First dispatch
        attr_dict['dHdS'] = dHdS
        attr_dict['encounters'] = 1
        for dispatcher in self.dispatchers:
            dispatcher.add_node(n, attr_dict)
        attr_dict['dHdS_count'] = {dHdS: 1}
        # MultiGraph, with edges keyed by (reacted_spec_pair, reaction_attr):
        # reaction_graph.adj[source][target][(reacted_spec_pair, reaction_attr)] = eattr
        #self.reaction_graph.add_node(target_state, node_attrs)
        DiGraph.add_node(self, n, attr_dict)


    def change_node(self, n, attr_dict=None, dHdS=None, **attr): #, source_state=None, edge_key=None, edge_attrs=None):
        """
        Update a strand state. If arriving to strand state node from a previous complex
        and you would like to update the dH_dS_count for the strand state (which should both be 0 for the free state),
        then you can give a complex as well.
        """
        # self.node[n]['encounters'] += 1  # could also simply be 'size'?
        # For now, we keep statistics of the different complex energies that have been calculated
        # https://treyhunner.com/2015/11/counting-things-in-python/
        # Using Counter is fast, but doesn't work well when we want to serialize the graph, so using dict:
        if dHdS is not None:
            try:
                self.node[n]['dHdS_count'][dHdS] += 1
            except KeyError:
                self.node[n]['dHdS_count'][dHdS] = 1
        if attr_dict is None:
            attr_dict = attr
        elif attr:
            attr_dict.update(attr)
        self.node[n].update(attr_dict)
        for dispatcher in self.dispatchers:
            dispatcher.change_node(n, attr_dict)


    #, source_state=None, edge_key=None, edge_attrs=None):
    def add_edge(self, source, target, attr_dict=None, **attr):
        """
        Update a strand state. If arriving to strand state node from a previous complex
        and you would like to update the dHdS_count for the strand state (which should both be 0 for the free state),
        then you can give a complex as well.
        """
        # self.node[n]['encounters'] += 1  # could also simply be 'size'?
        # For now, we keep statistics of the different complex energies that have been calculated
        # https://treyhunner.com/2015/11/counting-things-in-python/
        # Using Counter is fast, but doesn't work well when we want to serialize the graph, so using dict:
        if attr_dict is None:
            attr_dict = attr
        elif attr:
            attr_dict.update(attr)
        # self.adj[source][target].update(attr_dict)
        super().add_edge(source, target, attr_dict, **attr)
        for dispatcher in self.dispatchers:
            dispatcher.add_edge(None, source, target, attr_dict)

    def change_edge(self, source, target, attr_dict=None, **attr):
        """
        Update a strand state. If arriving to strand state node from a previous complex
        and you would like to update the dHdS_count for the strand state (which should both be 0 for the free state),
        then you can give a complex as well.
        """
        # self.node[n]['encounters'] += 1  # could also simply be 'size'?
        # For now, we keep statistics of the different complex energies that have been calculated
        # https://treyhunner.com/2015/11/counting-things-in-python/
        # Using Counter is fast, but doesn't work well when we want to serialize the graph, so using dict:
        if attr_dict is None:
            attr_dict = attr
        elif attr:
            attr_dict.update(attr)
        self.adj[source][target].update(attr_dict)
        for dispatcher in self.dispatchers:
            dispatcher.change_edge(None, source, target, attr_dict)

    def add_time_step(self, step_time):
        """ Add a step event with time (or other value). """
        for dispatcher in self.dispatchers:
            dispatcher.step(step_time)

    def close_all_dispatchers(self, clear=True):
        """ Loop over all dispatchers and close them. """
        for dispatcher in self.dispatchers:
            dispatcher.close()
        if clear:
            self.dispatchers.clear() # pylint: disable=E1101



def node_attr_repr(val):
    return (val if (isinstance(val, (numbers.Number, bool, tuple, list)))
            else str(val))


class GraphStreamEventDispatcher():
    """
    Class for dispatching events in the Graphstreaming format.
    """
    def __init__(self, writer=None):
        # self.writer = open(path, mode)  # Subject to change for the base class

        if writer is None:
            writer = io.StringIO()
        self.writer = writer
        self.create_sequential_edge_ids = True
        self.next_sequential_edge_id = 0
        self.edges_by_edgeid = {}
        self.edgeid_by_nodes = {}
        self.is_directed = True
        # We may not want to include or exclude some node/edge attributes:
        # Include only these attributes (set to None to include all attributes)
        # Filter first, then convert value, or other way around?
        # Order:
        # - attr_filter_include_keys, attr_filter_exclude_keys
        # - attr_filter_include_func, attr_filter_exclude_func
        # - then convert node(s), attrs
        # If a change_node/edge event does not have any attributes left after filtering,
        # it may be reasonable not to emit the event:
        self.emit_empty_change_events = False
        self.attr_filter_include_keys = None
        self.attr_filter_include_keys = {
            # node attrs:
            'pos', 'offset_relative_to', 'node_pos',
            'size', 'scale', 'n_strands',
            # 'tau_cum', 'encounters', 'count',  # Will yield a change_node event for every reaction
            'dG_first', 'dG_std', 'x', 'y', 'z',
            # edge attrs:
            'is_forming', 'is_intra', 'is_joining', 'is_splitting',
            'reaction_str', 'reaction_invocation_count'
            'dH', 'dS', 'dHdS', 'c_j', 'activity',
            'traversals', 'weight',
            'color', 'alpha', 'len', 'length',
        }
        self.attr_filter_include_func = None # str_or_numeric  # passing attr_key, attr_value
        self.attr_filter_exclude_keys = None
        self.attr_filter_exclude_func = None
        self.node_converter = node_attr_repr
        self.attr_converter = node_attr_repr  # Set to None when done debugging

    def get_combined_attr_filter_func(self):
        """
        Return a combined filtering function based on attr_filter_(in/ex)clude_(keys/func).
        Currently NOT optimized, but encapsulating here so we can optimize the current_attr_filter later.
        """
        filter_attributes = (
            self.attr_filter_include_keys, self.attr_filter_exclude_keys,
            self.attr_filter_include_func, self.attr_filter_exclude_func
        )
        if not any(filter_attributes):
            # No filters in place
            current_attr_filter = None
        elif all(att is None for att in
                 (self.attr_filter_exclude_keys, self.attr_filter_include_func, self.attr_filter_exclude_func)):
            # Optimization for when only attr_filter_include_keys is defined (typical situation)
            def current_attr_filter(key, value):
                return key in self.attr_filter_include_keys
        else:
            def current_attr_filter(key, value):
                return (
                    (self.attr_filter_include_keys is None or key in self.attr_filter_include_keys)
                    and (self.attr_filter_exclude_keys is None or key not in self.attr_filter_exclude_keys)
                    and (self.attr_filter_include_func is None or self.attr_filter_include_func(key, value))
                    and (self.attr_filter_exclude_func is None or not self.attr_filter_exclude_func(key, value))
                )
        return current_attr_filter

    def add_node(self, node, attrs):
        """ Add a single new node. """
        # print("GraphStreamEventDispatcher.add_node(%s, %s)" % (node, attrs))
        self.write('an', node, attrs=attrs)

    def change_node(self, node, attrs):
        """ Change/update an existing node's attributes. """
        self.write('cn', node, attrs=attrs)

    def delete_node(self, node):
        """ Delete an new node. """
        self.write('dn', node)

    def add_edge(self, edge_id, source, target, attrs):
        """
        Q: What's the difference between edge_id and edge_key?
        A: edge_ids must be unique across the whole graph, while a key must only
            be unique for a given (source, target) pair (Graph) or tuple (DiGraph).
           You can have many edges with the same key, and you can even have many
           edges from the same source/target!
        Q: But I have a regular Graph/DiGraph and don't need any edge_ids. Is there an easy way to generate
            edge_ids on the fly?
        A: Yes, simply set self.create_sequential_edge_ids = True and a unique edge id is created for each new edge.
        """
        if self.create_sequential_edge_ids:
            edge_key, edge_id = edge_id, self.next_sequential_edge_id
            self.next_sequential_edge_id += 1
            self.edges_by_edgeid[edge_id] = (source, target, edge_key)
            self.edgeid_by_nodes[(source, target)] = edge_id
            if not self.is_directed:
                self.edgeid_by_nodes[(target, source)] = edge_id
        self.write('ae', edge_id, source, target, attrs)

    def change_edge(self, edge_id=None, source=None, target=None, attrs=None):
        """ Change/update an existing edge's attributes. """
        # Note: The Nascent reaction graph edge key is unique *for a particular reactoin* but
        # a reaction can have multiple edges (e.g. reaction that merges or splits complexes)!
        # It might be better to just generate unique/sequential edge ids in ReactionGraph
        # (perhaps even for nodes as well)
        # However, my vispy graph viewer is regular DiGraph so I'm not using edge_id for anything.
        if self.create_sequential_edge_ids:
            # edge_key, edge_id = edge_id, self.next_sequential_edge_id
            assert source is not None and target is not None
            edge_id = self.edgeid_by_nodes[(source, target)]
        self.write('ce', edge_id, source, target, attrs)

    def delete_edge(self, edge_id=None, source=None, target=None):
        """ Delete an existing edge. """
        if edge_id is None:
            assert source is not None and target is not None
            edge_id = self.edgeid_by_nodes[(source, target)]
        self.write('de', edge_id, source, target)

    def step(self, step_time):
        """ Add a step event with time (or other value). """
        self.write('st', step_time)

    def encode(self, event, identifier, source=None, target=None, attrs=None):
        """
        Create a properly encoded event string, default is:
            {event: {identifier: attrs}}
        param :event: any of:
            'an', 'cn', 'dn' - add/change/delete node
            'ae', 'ce', 'de' - add/change/delete edge
            'cg', 'cl'       - change graph attrs; clear whole graph
            'st'             - Add new step (for e.g. dynamic graphs)
        identifier can be either:
            a node id  (for node operations)
            an edge id (for edge operations)
            a numeric value (for step events)
        For edges, we save source and target in the attrs dict.
        Note that this is a little different from the DGS file format, which explicitly has source and target in ae:
            <event>      ::= ( <an> | <cn> | <dn> | <ae> | <ce> | <de> | <cg> | <st> | <cl> ) ( <comment> | <EOL> )
            <an>         ::= "an" <id> <attributes>
            <cn>         ::= "cn" <id> <attributes>
            <dn>         ::= "dn" <id>
            <ae>         ::= "ae" <id> <id> ( <direction> )? <id> <attributes>  # ae e2 n1 > n2 weight:40
            <ce>         ::= "ce" <id> <attributes>
            <de>         ::= "de" <id>
            <cg>         ::= "cg" <attributes>
            <st>         ::= "st" <real>
            <cl>         ::= "cl"
        More info:
            https://github.com/panisson/pygephi_graphstreaming/blob/master/pygephi/client.py
            http://graphstream-project.org/doc/Advanced-Concepts/The-DGS-File-Format/
            http://igraph.org/python/doc/igraph.drawing.graph.GephiGraphStreamingDrawer-class.html
            http://igraph.org/python/doc/igraph.remote.gephi-module.html

        """
        # Order:
        # - attr_filter_include_keys, attr_filter_exclude_keys
        # - attr_filter_include_func, attr_filter_exclude_func
        # - then convert node(s), attrs
        if attrs is None:
            attrs = {}
        else:
            # attrs = attrs.copy()
            # Filter/conversion order:
            # - attr_filter_include_keys, attr_filter_exclude_keys
            # - attr_filter_include_func, attr_filter_exclude_func
            # - then convert node(s), attrs
            attr_filter_func = self.get_combined_attr_filter_func()
            if attr_filter_func:
                attrs = {k: v for k, v in attrs.items() if attr_filter_func(k, v)}
        if self.attr_converter:
            attrs = {k: self.attr_converter(v) for k, v in attrs.items()}

        if not self.emit_empty_change_events and event[0] == "c" and not attrs:
            return None

        if self.node_converter is None:
            if source:
                attrs['source'] = source
            if target:
                attrs['target'] = target
        else:
            if source is None and target is None:
                if event[1] == "n": # Node event, identifier is a node:
                    identifier = self.node_converter(identifier)
            else:
                if source:
                    attrs['source'] = self.node_converter(source)
                if target:
                    attrs['target'] = self.node_converter(target)
        return json.dumps({event: {identifier: attrs}}) + '\n'


    # Or should we call this method `send()` ?
    def write(self, event, identifier, source=None, target=None, attrs=None):
        """ Write event to default writer using default encoder. """
        # print("GraphStreamEventDispatcher.write(%s, %s, %s, %s, %s)" % (
        #     event, identifier, source, target, attrs))
        data = self.encode(event, identifier, source, target, attrs)
        if data:
            self.writer.write(data)
            return 1
        return 0


    def close(self):
        """ Close writer if open. """
        self.writer.close()



class GraphStreamFileDispatcher(GraphStreamEventDispatcher):

    def __init__(self, path, mode='a', **kwargs):
        writer = open(path, mode)
        super().__init__(writer, **kwargs)






def graph_state_partitions(g, T, update=False, unit='R'):
    """
    Calculate dG and partitions for reaction graph nodes (states).
    params:
    :g:         Graph
    :T:         Temperature in K
    :update:    If True, will update graph nodes with result.
    :unit:      The units of dH and dS. Currently only units of R*K and R are supported.
    """
    dHdS = {node: nattr.get('dHdS_first', (0.0, 0.0)) for node, nattr in g.node.items()}
    dG = {node: dH - T*dS for node, (dH, dS) in dHdS.items()}
    # Universe entropy, dS_uni = dS_surroundings + dS_system, in units of R:
    dSuni_R = {node: dS-dH/T for node, (dH, dS) in dHdS.items()}
    partitions = {node: math.exp(dS) for node, dS in dSuni_R.items()}
    Z = sum(partitions.values())
    partitions = {node: part/Z for node, part in partitions.items()}
    if update:
        for node, nattr in g.node.items():
            nattr['partition'] = partitions[node]
            nattr['dSuni'] = dSuni_R[node]
            nattr['dG'] = dG[node]
    return partitions


def reaction_attr_to_str(reaction_attr):
    """
    Alternative name: reaction_attr_repr ?   ("to_str" vs "repr")
    Would be nice to have a generic, overloadable, Julia-style multiple-dispatch, global "repr" function...
    Return
        h+ for forming  hybridization reactions,
        h- for breaking hybridization reactions,
        s+ for forming  stacking reactions,
        s- for breaking stacking reactions,
    Appends an asterix (*) for inter-molecular reactions.
    """
    return "".join((reaction_attr.reaction_type,
                    "+" if reaction_attr.is_forming else "-",
                    " " if reaction_attr.is_intra else "*"))


def reaction_spec_pair_to_str(reaction_spec_pair):
    """
    Alternative name: reaction_spec_pair_repr ?
    Return

    Appends an asterix (*) for inter-molecular reactions.
    """
    return ", ".join(str(val) for val in sorted(reaction_spec_pair))


def reaction_to_str(reaction_spec_pair, reaction_attr):
    """
    Alternative names: reaction_repr, reaction_as_str, reaction_str_repr ?
    Return
        h+*: s1_A > < s2_a  (0, 0)  for inter-complex hybridization
        h+ : s1_A >_< s2_a          for intra-complex hybridization
        h- : s1_A-<_>-s2_a          for de-hybridization (always intra-complex)
    for hybridization reactions, and
        s+ : s1_A3p >_< s2_B5p / s3_a5p >_< s4_b3p   (fp, fp, fp, fp)   for intra-complex stacking
        s- : s1_A3p:<_>:s2_B5p / s3_a5p:<_>:s4_b3p   (fp, fp, fp, fp)   for intra-complex stacking
    for stacking reactions, and
        b+ : s1_A3p >_< s2_B5p  (fp, fp)   for intra-complex ligation
        b- : s1_A3p.<_>.s2_B5p  (fp, fp)   for intra-complex backbone hydrolysis
    for backbon reactions.
    Some of the information is intentionally redundant to verify that everything is as it should be:
    E.g. the reaction_attr_str before ':' can be extracted from the reaction_spec_pair after the colon.
    Also the current domain hybridization state (indicated by dash after domain name in -<_>-) and the
    current domain-end stacking state (indicated by colon after end in 3p:<_>:) is (or should be)
    redundant with the indication of whether the reaction is forming (indicated by "> <" for hybridization/stacking
    reactions) or breaking (indicated by "< >" for dehybridization/unstacking reactions).
    """
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        d1, d2 = {}, {}
        d1['fp'], d2['fp'] = sorted(tuple(reaction_spec_pair))
        for d in (d1, d2):
            dspecie, is_hyb, cstate, icid = d['fp']
            d['dspecie'] = "_".join(dspecie)     # dspecie = (strand.name, domain.name)
            d['hyb'] = "-" if is_hyb else " "   # whether domain is currently hybridized
            d['cstate'] = str(cstate % 100000) if isinstance(cstate, int) else cstate  # complex state (fingerprint)
            d['icid'] = "" if icid == 0 else str(icid)  # in-complex identifier (for symmetric complexes)
        cstates = []
        for d in (d1, d2):
            if d['cstate'] not in cstates:
                cstates.append(d['cstate'])
        states_str = ", ".join(cstates)
        form_str = (">%s<" if reaction_attr.is_forming else "<%s>") % ("_" if reaction_attr.is_intra else " ")
        fmt = "{ra}:  {d1[dspecie]}{d1[hyb]}{form_str}{d2[hyb]}{d2[dspecie]}    ({states_str})"
        return fmt.format(ra=reaction_attr_to_str(reaction_attr), form_str=form_str, d1=d1, d2=d2, states_str=states_str)
    elif reaction_attr.reaction_type is STACKING_INTERACTION:
        h1end3p, h2end5p, h2end3p, h1end5p = {}, {}, {}, {}
        for d, efp in zip((h1end3p, h2end5p, h2end3p, h1end5p),
                          (efp for pair in tuple(reaction_spec_pair) for efp in pair)):
            # DomainEnd fingerprint = (domain.state_fingerprint, self.end,self.is_stacked())
            dfp, d['end'], is_stacked = efp
            # domain attributes:
            dspecie, is_hyb, cstate, icid = dfp
            d['dspecie'] = "_".join(dspecie)     # dspecie = (strand.name, domain.name)
            d['stacked'] = ":" if is_stacked else " "   # whether domain is currently hybridized
            d['cstate'] = str(cstate % 100000)   # complex state (fingerprint)
            d['icid'] = "" if icid == 0 else str(icid)  # in-complex identifier (for symmetric complexes)
        cstates = []
        for d in (h1end3p, h2end5p, h2end3p, h1end5p):
            if d['cstate'] not in cstates:
                cstates.append(d['cstate'])
        states_str = ", ".join(cstates)
        form_str = (">%s<" if reaction_attr.is_forming else "<%s>") % ("_" if reaction_attr.is_intra else " ")
        #             h1end3p         h1end5p
        # Helix 1   ----------3' : 5'----------
        # Helix 2   ----------5' : 3'----------
        #             h2end5p         h2end3p
        # Note: Biopython's "stack_string" is h1end3p h1end5p/h2end5p h2end3p
        # But the format below is in format   h1end3p h1end5p/h2end3p h2end5p
        fmt = ("{ra}:"
               "{h1end3p[dspecie]}{h1end3p[end]}{h1end3p[stacked]}{form_str}"
               "{h1end5p[stacked]}{h1end5p[dspecie]}{h1end5p[end]} / "
               "{h2end3p[dspecie]}{h2end3p[end]}{h2end3p[stacked]}{form_str}"
               "{h2end5p[stacked]}{h2end5p[dspecie]}{h2end5p[end]}    "
               #"({h1end3p[cstate]}, {h1end5p[cstate]}, {h2end3p[cstate]}, {h2end5p[cstate]})"
               "({states_str})"
              )
        return fmt.format(ra=reaction_attr_to_str(reaction_attr), form_str=form_str, states_str=states_str,
                          h1end3p=h1end3p, h2end5p=h2end5p, h2end3p=h2end3p, h1end5p=h1end5p)


# def reaction_eattrs(reaction_attr, activity, c_j, throttle_factor, dHdS):
#     """
#     """
#     return


# def reaction_edge_label(eattr):
#     """
#     Add a suitable label to reaction edge attribute dict :eattr:.
#     eattr dict should already contain the following entries:
#     * reaction_attr tuple values:
#     * reaction_type, is_forming, is_intra
#     * is_forming_str
#     * is_intra_str
#     * dS
#     * dH
#     * loop type enum (0=no loop, 1, 2)
#     * activity
#     * c_j
#     * throttle_factor
#     """
#     return "{reaction_type}{is_forming_str}{is_intra_str} {c_j:0.1e}"
