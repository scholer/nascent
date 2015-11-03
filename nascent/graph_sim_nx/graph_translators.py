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



import networkx as nx
from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION,
                        interactions)


interaction_idx_to_str = [None, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION]




def translate_domain_graph_to_5p3p_graph(graph):
    directed = getattr(graph, 'is_directed', False)
    trans_graph = nx.DiGraph() if directed else nx.Graph()

    for domain, attrs in graph.nodes(data=True):
        trans_graph.add_node(str(domain)+":5p", attrs)
        trans_graph.add_node(str(domain)+":3p", attrs)
        trans_graph.add_edge(str(domain)+":5p", str(domain)+":3p", attrs)

    for source, target, attrs in graph.edges(data=True):
        interaction = attrs.get('interaction')
        if attrs.get('interaction') in (1, 'pb', 'p', 'backbone', 3, 'bs', 's', 'stacking'):
            # For edges connecting domain backbone, we connect
            # the 3p end of the 5p-domain to the 5p domain of the 3p domain
            # 5'-domainA-3' 5'-domainB-3'
            trans_graph.add_edge(str(source)+":3p", str(target)+":5p", attrs)
        elif attrs.get('interaction') in (2, 'h', 'hh', 'hybridization'):
            # Always connecting 3p end to 5p end:
            trans_graph.add_edge(str(target)+":3p", str(source)+":5p", attrs)
            trans_graph.add_edge(str(source)+":3p", str(target)+":5p", attrs)
        else:
            print("Unrecognized interaction '%s' between '%s' and '%s'" %
                  (attrs.get('interaction'), source, target))
            trans_graph.add_edge(str(target)+":3p", str(source)+":5p", attrs)

    return trans_graph


def translate_domain_change_to_domain_graph_event(change):
    """ Translate a single state change to 5p3p graph format. """
    #print("Translating domain change to domain graph event...")
    if change['change_type'] in (0, 'node event'):
        pass
    elif change['change_type'] in (1, 'edge event'):
        dom1, dom2 = change['nodes']
        if change['interaction']:
            if isinstance(change['interaction'], int):
                change['interaction'] = interaction_idx_to_str[change['interaction']]
            else:
                # TODO: Is this really weird? What is the standard? Is the standard an integer?
                print("Weird: Non-integer %s interaction for change %s" %
                      (change['interaction'], change))
    return change


def translate_domain_change_to_5p3p_graph_events(change):
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
        if change['interaction']:
            change['interaction'] = interaction_idx_to_str[change['interaction']]
        if change['interaction'] in (1, 'pb', 'p', 'backbone', 3, 'bs', 's', 'stacking'):
            # The 3p of the first (5p-most) domain connects to the 5p end of the second (3p-most) domain
            event['nodes'] = [str(dom1)+":3p", str(dom2)+":5p"]
            return [event]
        else:
            print("Unrecognized interaction for change %s" % change)
        if change['interaction'] in (2, 'h', 'hh', 'hybridization'):
            # The 3p of the first domain connects to the 5p end of the second (3p-most) domain,
            # AND connect the 5p end of the first to the 3p of the second.
            # We have already done the former above, only do the latter:
            event2 = change.copy()
            event['nodes'] = [str(dom1)+":3p", str(dom2)+":5p"]
            event2['nodes'] = [str(dom1)+":5p", str(dom2)+":3p"]
            return [event, event2]


def translate_domain_change_to_strand_graph_event(change):
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
