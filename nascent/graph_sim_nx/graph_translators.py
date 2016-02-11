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

# pylint: disable=C0103,W0611

"""
Module for conversion of graphs at different levels.
Contains two types of functions:
    1) Convert graphs within the model domain, e.g. domain-level to 5p3p-level (nodes are still objects).
    2) Convert graphs from the model domain to string representation, e.g. Domain(nae="A", duid=5) -> "A#5".



"""


from __future__ import absolute_import, print_function, division
import networkx as nx

from .domain import Domain
from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION,
                        interactions, valid_interactions)

# same as constants.interactions
interaction_idx_to_str = [None, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION]


def translate_domain_graph_to_domain_str_graph(graph, directed=None, node_attrs=None, edge_attrs=None):
    """
    Translate a full domain-object graph to domain str-repr graph.
    """
    if directed is None:
        try:
            directed = graph.is_directed()
        except AttributeError:
            directed = False
    trans_graph = (nx.MultiDiGraph() if directed else nx.MultiGraph()) if graph.is_multigraph() else \
                  (nx.DiGraph() if directed else nx.Graph())
    if node_attrs is None:
        node_attrs = {}
    if edge_attrs is None:
        edge_attrs = {}

    trans_graph.add_nodes_from([(str(domain), attrs) for domain, attrs in graph.nodes_iter(data=True)])
    if graph.is_multigraph():
        trans_graph.add_edges_from([(str(domain1), str(domain2), key, attrs)
                                    for domain1, domain2, key, attrs in graph.edges_iter(data=True, keys=True)])
    else:
        trans_graph.add_edges_from([(str(domain1), str(domain2), attrs)
                                    for domain1, domain2, attrs in graph.edges_iter(data=True)])

    return trans_graph


def translate_domain_graph_to_5p3p_str_graph(graph, directed=None):
    """
    Translate a full domain-object graph to 5p3p str-repr graph.
    Uhm... Why not just initialize from the ends5p3p graph??
    """
    if directed is None:
        try:
            directed = graph.is_directed()
        except AttributeError:
            directed = False
    trans_graph = (nx.MultiDiGraph() if directed else nx.MultiGraph()) if graph.is_multigraph() else \
                  (nx.DiGraph() if directed else nx.Graph())
    for domain, attrs in graph.nodes(data=True):
        trans_graph.add_node(str(domain)+":5p", attr_dict=attrs)
        trans_graph.add_node(str(domain)+":3p", attr_dict=attrs)
        trans_graph.add_edge(str(domain)+":5p", str(domain)+":3p", key=PHOSPHATEBACKBONE_INTERACTION, attr_dict=attrs)

    if graph.is_multigraph():
        edge_tups = graph.edges(keys=True, data=True)
    else:
        edge_tups = [(domain1, domain2, None, attrs) for domain1, domain2, attrs in graph.edges(data=True)]

    for domain1, domain2, key, attrs in edge_tups:
        interaction = attrs.get('interaction')
        if interaction in (PHOSPHATEBACKBONE_INTERACTION, STACKING_INTERACTION):
            # For edges connecting domain backbone, we connect
            # the 3p end of the 5p-domain to the 5p domain of the 3p domain
            # 5'-domainA-3' 5'-domainB-3'
            trans_graph.add_edge(str(domain1)+":3p", str(domain2)+":5p", key=key, attr_dict=attrs)
        elif interaction == HYBRIDIZATION_INTERACTION:
            # Always connecting 3p end to 5p end:
            trans_graph.add_edge(str(domain1)+":3p", str(domain2)+":5p", key=key, attr_dict=attrs)
            trans_graph.add_edge(str(domain1)+":3p", str(domain2)+":5p", key=key, attr_dict=attrs)
        else:
            print("Unrecognized interaction '%s' between '%s' and '%s'" %
                  (attrs.get('interaction'), domain1, domain2))
            trans_graph.add_edge(str(domain1)+":3p", str(domain2)+":5p", key=key, attr_dict=attrs)

    return trans_graph


def translate_domain_graph_to_strand_str_graph(graph, directed=None, node_attrs=None, edge_attrs=None):
    """
    Translate a full domain-object graph to domain str-repr graph.
    """
    if directed is None:
        directed = getattr(graph, 'is_directed', False)
    is_multi = graph.is_multigraph()
    trans_graph = nx.MultiDiGraph() if directed else nx.MultiGraph() # strand graphs are always multigraphs
    if node_attrs is None:
        node_attrs = {}
    if edge_attrs is None:
        edge_attrs = {}

    nodes = {str(domain.strand): attrs for domain, attrs in graph.nodes(data=True)} # strand-name => attrs
    trans_graph.add_nodes_from(nodes.items())
    # Don't worry too much about multiple occurences.
    if is_multi:
        edges = [(str(domain1.strand), str(domain2.strand),
                  "-".join((domain1.instance_name, key, domain2.instance_name)),
                  attrs)
                 for domain1, domain2, key, attrs in graph.edges(data=True, keys=True)
                 if domain1.strand != domain2.strand # don't include loops (src, tgt connected to the same node)
                ]
    else:
        edges = [(str(domain1.strand), str(domain2.strand),
                  "-".join((domain1.instance_name, domain2.instance_name)),
                  attrs)
                 for domain1, domain2, attrs in graph.edges(data=True)
                 if domain1.strand != domain2.strand # don't include loops (src, tgt connected to the same node)
                ]
    trans_graph.add_edges_from(edges)

    return trans_graph


def make_edge_key(source, target, interaction, level="domain"):
    """
    """
    if level == "5p3p":
        key = interaction
    elif level == "domain":
        key = make_domain_edge_key(source, target, interaction)
    elif level == "strand":
        key = make_domain_edge_key(source, target, interaction)
    return key

def make_5p3p_edge_key(source, target, interaction):
    return interaction

def make_domain_edge_key(source, target, interaction):
    if interaction == HYBRIDIZATION_INTERACTION:
        return interaction
    else:
        # Stacking or backbone interaction:
        # For undirected graphs we need to distinguish
        # ({A, B}, "A3p-s-B5p") from ({A, B}, "B3p-s-A5p")
        # Although this distinction is only needed for circular constructs, which are rather rare.
        # You could also use domain1.instance_name, domain2.instance_name, key
        #return "%s:3p-%s-%s5p" % (source, interaction, target)
        return interaction

def make_strand_edge_key(source, target, interaction):
    """
    Alternatives:
     key = "-".join((domain1.instance_name, domain2.instance_name, key))
    """
    if interaction == HYBRIDIZATION_INTERACTION:
        return "%s-%s-%s" % (source, interaction, target)
    else:
        return "%s:3p-%s-%s5p" % (source, interaction, target)


def translate_domain_change_to_domain_str_graph_event(change):
    """ Translate a single domain-object state change to domain str-repr graph format. """
    #print("Translating domain change to domain graph event...")
    # TODO: Support multi-graphs and edge keys.
    #       Maybe just provide change['keys'] for keyed edges, propagate the state change as-is,
    #       and let the livestream adaptor take care of the conversion? (Like it already does in many cases)
    if change['change_type'] in (0, 'node event'):
        dom1, = change['nodes']
        change['nodes'] = [str(dom1)]
    elif change['change_type'] in (1, 'edge event'):
        dom1, dom2 = change['nodes']
        change['nodes'] = [str(dom1), str(dom2)]
        interaction = change['interaction']
        if interaction:
            if isinstance(interaction, int):
                change['interaction'] = interaction_idx_to_str[interaction]
            elif interaction not in valid_interactions:
                print("Weird: Unexpected interaction %s for change %s" % (interaction, change))
        # Suggested key: Specifying the ends is not strictly necessary. Actually, since a key only has to be
        # valid/unique together with a pair of (source, target), specifying the ends would be enough.
        # However, for undirected graphs, we need to distinguish
        # ({A, B}, "A3p-s-B5p") from ({A, B}, "B3p-s-A5p")
        # Although we really only need this distinction for circular constructs, which are rather rare.
        # if interaction == HYBRIDIZATION_INTERACTION:
        #     change['key'] = interaction
        # else:
        #     # stacking and backbone can connect in two ways, so we need to specify what ends are interacting:
        #     # change['key'] = "%s-%s-%s" % (dom1.end3p, interaction, dom2.end5p)
        #     change['key'] = "3p-%s-5p" % (interaction,)
        change['key'] = make_domain_edge_key(dom1, dom2, interaction)
    else:
        raise ValueError("Change type %s not recognized for change %s" % (change['change_type'], change))
    return change


def translate_domain_change_to_5p3p_str_graph_event(change):
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
                          [(str(domain.end5p), str(domain.end3p))
                           for domain in change['nodes']]]
        return [event]
        # Add style to event?
    elif change['change_type'] in (1, 'edge event'):
        dom1, dom2 = change['nodes']
        if isinstance(change['interaction'], int):
            change['interaction'] = interaction_idx_to_str[change['interaction']]
        # Suggested key: For ends5p3p, the interaction should be sufficient.
        change['key'] = change['interaction']
        if change['interaction'] in (1, 'pb', 'p', 'backbone', 3, 'bs', 's', 'stacking'):
            # The 3p of the first (5p-most) domain connects to the 5p end of the second (3p-most) domain
            event['nodes'] = [str(dom1.end3p), str(dom2.end5p)]
            return [event]
        elif change['interaction'] in (2, 'h', 'hh', 'hybridization'):
            # The 3p of the first domain connects to the 5p end of the second (3p-most) domain,
            # AND connect the 5p end of the first to the 3p of the second.
            # We have already done the former above, only do the latter:
            event2 = change.copy()
            event['nodes'] = [str(dom1.end3p), str(dom2.end5p)]
            event2['nodes'] = [str(dom1.end5p), str(dom2.end3p)]
            return [event, event2]
        else:
            print("Unrecognized interaction for change %s" % change)


def translate_domain_change_to_strand_str_graph_event(change):
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
