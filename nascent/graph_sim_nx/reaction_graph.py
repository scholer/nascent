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
from datetime import datetime
from collections import defaultdict
import yaml
import json
import networkx
nx = networkx
from networkx.classes.digraph import DiGraph

from .constants import STACKING_INTERACTION, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION


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
        self.params = params
        if 'reaction_graph_default_attrs' in params:
            self.graph.update(params['reaction_graph_default_attrs'])
        # File with changes to the reaction graph, e.g. new nodes/edges and node/edge updates:
        self.endstates_by_reaction = {}  # [startstate][(reaction_spec_pair, reaction_attr)] = endstate
        self.endstates_by_reaction[0] = defaultdict(list) # Also adding the "null" node
        self.reverse_reaction_key = {} # get reaction edge key for the opposite direction.
        # self.reverse_reaction[(rx_spec_pair, rx_attr)] = (rev_rx_spec_pair, rev_rx_attr)
        # used to be a dict {edge_key => target} but we can have multiple targets for a single reaction.
        self.reaction_graph_complexes_directory = params.get("reaction_graph_complexes_directory")
        if self.reaction_graph_complexes_directory is not None:
            if not os.path.exists(self.reaction_graph_complexes_directory):
                print("Creating directory for complex files:", self.reaction_graph_complexes_directory)
                os.makedirs(self.reaction_graph_complexes_directory)
            assert os.path.isdir(self.reaction_graph_complexes_directory)

        self.dispatchers = []
        self.reaction_graph_events_file = params.get('reaction_graph_events_file')
        if self.reaction_graph_events_file:
            #self.reaction_graph_delta_file = open(self.reaction_graph_delta_file, 'a')
            #self.open_files.append(self.reaction_graph_delta_file)
            gs_file_dispatcher = GraphStreamFileDispatcher(self.reaction_graph_events_file)
            gs_file_dispatcher.writer.write("# New reaction graph initialized at %s" % datetime.now())
            self.dispatchers.append(gs_file_dispatcher)




    def add_node(self, n, attr_dict=None, **attr):
        """ Consider using the graph stream package spec? """
        DiGraph.add_node(self, n, attr_dict, **attr)


    def close_all_dispatchers(self):
        """ Loop over all dispatchers and close them. """
        for dispatcher in self.dispatchers:
            dispatcher.close()




class GraphStreamEventDispatcher():


    def __init__(self, path, mode='a'):
        self.writer = open(path, mode)  # Subject to change for the base class


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
        Note that this is a little different from the DGS file format:
            <event>      ::= ( <an> | <cn> | <dn> | <ae> | <ce> | <de> | <cg> | <st> | <cl> ) ( <comment> | <EOL> )
            <an>         ::= "an" <id> <attributes>
            <cn>         ::= "cn" <id> <attributes>
            <dn>         ::= "dn" <id>
            <ae>         ::= "ae" <id> <id> ( <direction> )? <id> <attributes>
            <ce>         ::= "ce" <id> <attributes>
            <de>         ::= "de" <id>
        More info:
            https://github.com/panisson/pygephi_graphstreaming/blob/master/pygephi/client.py
            http://graphstream-project.org/doc/Advanced-Concepts/The-DGS-File-Format/
            http://igraph.org/python/doc/igraph.drawing.graph.GephiGraphStreamingDrawer-class.html
            http://igraph.org/python/doc/igraph.remote.gephi-module.html

        """
        if attrs is None:
            attrs = {}
        if source:
            attrs['source'] = source
        if target:
            attrs['target'] = target
        return json.dumps({event: {identifier: attrs}})


    def write(self, event, identifier, source, target, attrs):
        """ Write event to default writer using default encoder. """
        self.writer.write(self.encode(event, identifier, source, target, attrs))


    def close(self):
        """ Close writer if open. """
        self.writer.close()



class GraphStreamFileDispatcher(GraphStreamEventDispatcher):

    def __init__(self, path, mode='a'):
        self.writer = open(path, mode)






def graph_state_partitions(g, T, update=False, unit='R'):
    """
    Calculate dG and partitions for reaction graph nodes (states).
    params:
    :g:         Graph
    :T:         Temperature in K
    :update:    If True, will update graph nodes with result.
    :unit:      The units of dH and dS. Currently only units of R*K and R are supported.
    """
    dHdS = {node: nattr.get('dH_dS_first', (0.0, 0.0)) for node, nattr in g.node.items()}
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
