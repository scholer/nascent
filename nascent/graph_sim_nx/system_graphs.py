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

# pylint: disable=C0103


"""

Module for representing system-level graphs (where each complex form it's own subgraph)
 - domains
 - strands


"""

from collections import defaultdict
import networkx as nx
# import pdb
# from pprint import pprint

from nascent.graph_sim_nx.constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION,
                                            STACKING_INTERACTION)



class InterfaceMultiGraph(nx.MultiGraph):  # Graph or MultiGraph?
    """
    Using a MultiGraph would make delegation/undelegation MUCH easier, since we could
        use the original source delegator as key in MultiGraph.adj.
        Then we wouldn't have to spend as much time determining whether edges should
        be removed or not in undelegate and which graph should be used for representation.
    """
    pass


class InterfaceGraph(nx.Graph):  # Graph or MultiGraph?
    """
    Graph or MultiGraph?
     - Hybridization and stacking interaction edges in ends5p3p graph will always be
        represented by a node merge/delegation; the only edges we have are backbone edges.

     - However, using a MultiGraph would make delegation/undelegation MUCH easier, since we could
        use the original source delegator as key in MultiGraph.adj.
        Then we wouldn't have to spend as much time determining whether edges should be removed or not.
    """

    def check_all_edges(self, do_raise=True):
        """ Quickly check that self.adj[source][target] is equivalent to self.adj[target][source] """
        all_ok = True
        for source in self.adj:
            for target in self.adj[source]:
                try:
                    assert source in self.adj[target]
                except AssertionError as e:
                    all_ok = False
                    print("FAIL: source %s not in self.adj[%s] = %s" % (source, target, self.adj[target]))
                    if do_raise:
                        raise e
                try:
                    assert target in self.adj[source]
                except AssertionError as e:
                    all_ok = False
                    print("FAIL: target %s not in self.adj[%s] = %s" % (target, source, self.adj[source]))
                    if do_raise:
                        raise e
                try:
                    assert self.adj[source][target] == self.adj[target][source]
                except AssertionError as e:
                    all_ok = False
                    print("FAIL: self.adj[%s][%s] != self.adj[%s][%s]" % (source, target, target, source))
                    if do_raise:
                        raise e
                except KeyError:
                    pass
        if all_ok:
            print("All edged checked OK.")


    def delegate(self, delegator, delegatee):
        """
        Delegate graph responsibility from delegator node to delegatee node.
        Both nodes should be contained in this graph and both should be
        an instance of InterfaceNode.
        """
        # Graph.edge[source][target] = edge_attr_dict
        # MultiGraph.edge[source][target][key] = edge_attr_dict

        ## 0. Make sure everything looks right:
        assert delegator.delegatee is None
        assert delegator not in delegatee.delegated_edges
        # assert len(set(delegator.delegated_edges.keys()) & set(delegatee.delegated_edges.keys())) == 0
        # Equivalently, use dict_keys.isdisjoint method:
        assert delegatee.delegated_edges.keys().isdisjoint(delegator.delegated_edges)

        ## 1. Update delegatee.delegated_edges with entries from delegator.
        # Keep delegator.delegated_edges intact.
        delegatee.delegated_edges.update(delegator.delegated_edges)
        # And set delegator's delegatee attribute:
        delegator.delegatee = delegatee

        ## 2. Update graph representation:
        # delegatee.delegated_edges[delegator] = []
        # for target, data in self.edge[delegator].items():
        #     self.remove_edge(delegator, target)
        #     self.add_edge(delegatee, target)
        #     delegatee.delegated_edges[delegator].append((target, data))
        ## Optimized version:
        # target_nodes = self.adj[delegator]  # {target1: {}, target2: {}, ...} dict
        # delegate_edges = delegator.delegated_edges  # {delegator: {target: {}, ...}}
        # #self.adj[delegatee].update(target_nodes)
        # #pdb.set_trace()
        # Update target nodes:
        ## 2a. Move all incoming edges from delegator to delegatee
        #for target, eattr in delegator.delegated_edges.items():
        for target, eattr in self.adj[delegator].items():
            if target not in self.adj[delegatee]:
                assert delegatee not in self.adj[target]
                # print("Making new edge from delegatee %s to target %s" % (delegatee, target))
                self.adj[delegatee][target] = eattr
                self.adj[target][delegatee] = eattr
                if 'edge_count' not in self.adj[target][delegatee]:
                    self.adj[target][delegatee]['edge_count'] = 1
                    # Do you *know* that edge_count is 1??
                    # If the edge already represents a "collapsed" representation of e.g. a duplex, edge_count
                    # would be 2. Thus, we only set it to 1 if it isn't set yet.
            else:
                assert delegatee in self.adj[target]
                # print("delegatee %s existing target/edge: %s/%s" % (delegatee, target, self.adj[target][delegatee]))
                # Determine which edge is the better/shorter representation...
                # Although they should be equal...
                try:
                    self.adj[target][delegatee]['edge_count'] += 1
                except KeyError:
                    self.adj[target][delegatee]['edge_count'] = 2
            del self.adj[target][delegator]
            # del self.adj[delegator][target]  # cannot delete entries from dict while iterating over it.
            # (but, you can just convert the dict.keys() iterator to a list and THEN you can delete entries in-place.)
        # 2b. Remove all outgoing edges from delegator:
        self.adj[delegator] = {} # Do not use self.adj[delegator].clear(); that will clear target_nodes in-place
        # delegatee.delegated_edges[delegator] = target_nodes
        # print("edges to %s delegated from %s to %s" %
        #       (list(delegatee.delegated_edges[delegator].keys()), delegator, delegator.delegatee))


    def undelegate(self, delegator, delegatee):
        """
        Reverse the effect of delegate.
        """
        # print("adj before undelegating %s from %s:" % (delegator, delegatee))
        # pprint(self.adj)
        #pdb.set_trace()

        ## 0. Make sure everything looks right:
        assert delegator.delegatee is delegatee
        assert delegator in delegatee.delegated_edges
        ## Which one of these do we prefer? Probably the latter; we keep delegator node in the graph after all.
        # assert delegator not in self.adj or self.adj[delegator] is None or len(self.adj[delegator]) == 0
        assert len(self.adj[delegator]) == 0

        ## 1. Update delegatee.delegated_edges, removing all delegated entries now belonging to delegator again.
        for k in delegator.delegated_edges:
            del delegatee.delegated_edges[k]
        # And reset delegator's delegatee attribute:
        delegator.delegatee = None

        ## 2. Update graph representation:
        # Two approaches:
        #  (a) adj-based: Look at delegatee's current edge targets in graph.adj and determine if
        #       ( i) there should be an edge from delegator to that target, and
        #       (ii) if the edge from delegatee to target should still be in place.
        #  (b) delegated_edges-based: Re-generate delegator edges from scratch.
        #       For each original source-target edge in
        #       delegator.delegated_edges = {source: {target: {}, ...}},
        #       add an edge from delegator to the top-delegate for all targets.
        #       (It should be ok to assume that delegator is now the top-delegatee for all sources, but check anyways.)
        assert all(source.top_delegate() is delegator for source in delegator.delegated_edges)  ## TODO: Remove check
        # top_delegate_edges = [(source, target.top_delegate(), eattr)
        #                       for source, targets in delegator.delegated_edges.items()
        #                       for target, eattr in targets.items()]
        # top_delegate_edges = sorted(top_delegate_edges, key=lambda tup: (tup[0] == delegator))
        # delegator_targets = {tup[1]: tup for tup in top_delegate_edges}
        delegator_targets = defaultdict(dict)
        for source, targets in delegator.delegated_edges.items():
            for target, eattr in targets.items():
                delegator_targets[target.top_delegate()][source] = eattr
        remaining_delegatee_targets = {target.top_delegate()
                                       for source, targets in delegatee.delegated_edges.items()
                                       for target in targets}

        ## Approach 2a: Move or "copy" edges from self.adj:
        edges = [(target, eattr) for target, eattr in self.adj[delegatee].items() if target in delegator_targets]
        # print("undelegate delegator %s  from  delegatee %s" % (delegator, delegatee))
        # print("delegator_targets:", delegator_targets)
        # print("remaining_delegatee_targets:", remaining_delegatee_targets)
        # print("Re-considering %s edges: [%s]", (delegatee, edges))
        for target, eattr in edges:
            if target in remaining_delegatee_targets:
                # print("Target %s still in remaining delegatee targets, making edge using other eattr..." % target)
                # Make a new edge...
                if delegator in delegator_targets[target]:
                    new_eattr = delegator_targets[target][delegator]
                else:
                    # pick edge with highest edge_count?
                    new_eattr = max(delegator_targets[target].values(), key=lambda attr: attr.get('edge_count', 0))
                self.adj[delegator][target] = new_eattr
                self.adj[target][delegator] = new_eattr
            else:
                # Move edge back from delegatee to delegator:
                # print("Moving edge to target %s from delegatee %s to delegator %s" % (target, delegatee, delegator))
                self.adj[delegator][target] = eattr
                self.adj[target][delegator] = eattr
                # self.adj[target][delegator] = eattr
                del self.adj[delegatee][target]
                if target is not delegatee:
                    # If the edge is a self-loop (target is delegatee) we would get a KeyError for removing twice:
                    del self.adj[target][delegatee]

        ## Approach 2b: Re-generate delegator edges from delegator.delegated_edges:
        # for source, target_delegate, eattr in top_delegate_edges:
        #     assert target_delegate in self.adj[delegatee]
        #     # Remove edge to target_delegate from delegatee (or decrement edge_count)
        #     self.adj[delegatee][target_delegate]['edge_count'] -= 1
        #     if self.adj[delegatee][target_delegate]['edge_count'] > 0:
        #         ## We should still have another edge between delegatee and target_delegate
        #         ## edge_count, target_delegate not in remaining_delegatee_targets and
        #         ## eattr is not self.adj[delegatee][target_delegate] should all be equivalent.
        #         # > 0 because we have already subtracted 1.
        #         assert target_delegate not in remaining_delegatee_targets
        #         # eattr should not be for the edge between delegatee and target_delegate.
        #         # If it is, then we would have to find another eattr dict to replace the one
        #         # we are moving back to delegator.
        #         assert eattr is not self.adj[delegatee][target_delegate]
        #     else:
        #         assert target_delegate in remaining_delegatee_targets
        #         assert eattr is not self.adj[delegatee][target_delegate]
        #         del self.adj[delegatee][target_delegate]
        #
        #     # Add edge from delegator to taget_delegatee (or increment edge count)
        #     if target_delegate in self.adj[delegator]:
        #         self.adj[delegator][target_delegate]['edge_count'] += 1
        #     else:
        #         self.adj[delegator][target_delegate] = eattr
        #         try:
        #             assert self.adj[delegator][target_delegate]['edge_count'] == 1
        #         except KeyError:
        #             self.adj[delegator][target_delegate]['edge_count'] = 1




        ### Old code from when I was using the 'nested' delegated_edges data structure.
        # target_nodes = delegatee.delegated_edges.pop(delegator)
        # # del delegatee.delegated_edges[delegator]
        # remaining_delegatee_targets = {target for delegator, targets in delegatee.delegated_edges.items()
        #                                for target in targets}
        # #self.adj[delegator].update(target_nodes) # done below in the for-loop
        # for target, eattr in target_nodes.items():
        #     assert delegator not in self.adj[target]
        #     assert target not in self.adj[delegator]
        #     self.adj[target][delegator] = eattr
        #     self.adj[delegator][target] = eattr
        #     print("self.adj[%s][%s] = %s" % (delegator, target, eattr))
        #     if target not in remaining_delegatee_targets:
        #         if target not in self.adj[delegatee] or delegatee not in self.adj[target]:
        #             print(target in self.adj[delegatee], delegatee in self.adj[target])
        #             pdb.set_trace()
        #         print("del self.adj[%s][%s]" % (delegatee, target))
        #         del self.adj[delegatee][target]
        #         print("del self.adj[%s][%s]" % (target, delegatee))
        #         del self.adj[target][delegatee]
        #     else:
        #         # We still have an edge from delegatee to target.
        #         pass
        # # Alternative to above: map(self.adj[delegatee].__delitem__, target_nodes)
        # delegator.delegatee = None
        # print("adj *after* undelegating %s from %s:" % (delegator, delegatee))
        # pprint(self.adj)


    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        """
        Add edges as "native" edges between nodes.
        """
        ebunch = list(ebunch) # In case it is a generator...
        super().add_edges_from(ebunch, attr_dict, edge_count=1, **attr)
        #print("Updating IF nodes delegated edges for ebunch %s..." % (ebunch,))
        for e in ebunch:
            u, v = e[:2]
            u.delegated_edges[u][v] = self.adj[u][v]
            #print("%s.delegated_edges: %s" % (u, u.delegated_edges))
            v.delegated_edges[v][u] = self.adj[v][u]
            #print("%s.delegated_edges: %s" % (v, v.delegated_edges))


    def add_edge(self, u, v, attr_dict=None, **attr):
        """
        Add edges as "native" edges between nodes.
        """
        print("Updating delegated edges for IF nodes  %s and %s..." % (u, v))
        super().add_edge(u, v, attr_dict, edge_count=1, **attr)
        ## Q: How much do we want the "representation" eattrs to be tied to eattrs in delegated_edges?
        u.delegated_edges[u][v] = self.adj[u][v]
        v.delegated_edges[v][u] = self.adj[v][u]
        # print("%s.delegated_edges: %s" % (u, u.delegated_edges))
        # print("%s.delegated_edges: %s" % (v, v.delegated_edges))


    def print_delegate_info(self):
        """
        print delegation information for this node.
        """
        for node in sorted(self.nodes()):
            if node.delegatee is not None:
                print(node, "--> edges delegated to", node.delegatee)
                if len(self.adj[node]) > 0:
                    print(" - WARNING: len(self.adj[%s]): %s" % (node, len(self.adj[node])))
            else:
                print(node, "<-- delegated edges (including self):")
                for delegator, target_nodes in node.delegated_edges.items():
                    print(" -%s%s: %s%s" % ("(" if delegator is node else " ", delegator,
                                            target_nodes, ")" if delegator is node else " "))
                    if node is not delegator:
                        all_targets_from_all_sources = {v for u in self.adj for v in self.adj[u]}
                        outgoing = len(self.adj[delegator]) > 0
                        incoming_targets = any(delegator in self.adj[target] for target in target_nodes)
                        incoming_all = delegator in all_targets_from_all_sources
                        if outgoing or incoming_targets or incoming_all:
                            #print("    - adj[%s] is empty: %s" % (delegator, ))
                            print("    - WARNING: delegator %s has edges! " % (delegator, ))
                            print("       - Outgoing edges:   ", self.adj[delegator])
                            print("       - Incoming, targets:", not incoming_targets)
                            print("       - Incoming, all:    ", not incoming_all)
                            print("       -", outgoing, incoming_targets, incoming_all)




class InterfaceNode():
    """
    How should delegate_edges be ordered? Should the delegator be the "original" node, or the
    For instance, consider the graph
        0------1---2-------3        0---.  1        .----3
                                         `.----2---:
        4------5---6-------7        4----´ 5   6    `---7
    Where 5 was delegated to 1 and (1 and 6) was delegated to 2.
    I.e from 5 we have a tree delegation branch: 5 --> 1 --> 2
    What should the datastructure of IN2 delegated_edges be w.r.t. 5 ?
     (a)    5 => {4: {}, 6: {}}, 1 => {0: {}, 2: {}} - i.e. store the *original* {source: targets} edges.
     (b)    1 => {4: {}, 6: {}, 0: {}, 2: {}}        - and then in node 1 have {5 => {4: {}, 6: {}}}
     (c)    1 => {5 => {4: {}, 6: {}}, 1 => {0: {}, 2: {}}}  - i.e. nested, 1 level for each delegation.

    Option (a) seems simplest: Any element in delegate_edges will *always* be an *original* edge,
    corresponding 1:1 to an edge in the ends5p3p graph.
    Usage of (a):
        delegation:   Just add all edges from delegator.delegated_edges to delegatee.delegated_edges,
            and update the graph.adj accordingly.
        undelegation: Since delegator still has all the original edges that have been delegated to it
            available in delegator.delegated_edges, just remove all these from delegatee.delegated_edges,
            and update graph.adj accordingly.
        When adding target edges to graph.adj, make sure to use InterfaceNode.top_delegate as target node.
            This is the only way to ensure that the target is correct and its representation have not
            been delegated to another node.
        Potential issue: What if there is a key in delegator.delegated_edges that is also
            already in delegatee.delegated_edges? When could that happen?
            Since the only original edges we consider are phosphate-backbone connections,
            Any node will have at most 2 original edges: upstream or downstream neighbor DomainEnd.
            Thus, the only way a "target" can be in both delegator and delegatee is if
            a top delegatee or any of the delegators below it shares an edge to the same target.
            One way this could happen is if we have a single-domain circular duplex
            where the 3p end stack with the 5p end of the same domain. We generally do not permit that.

    """
    def __init__(self, domain_end):
        self.domain_end = domain_end
        self.delegatee = None # If this node is represented by another node, this node is stored in this attribute.
        # {delegator => {target_node: target_edge_attr_dict, ... (for all original edges)},
        #  ... for all delegators and their delegated edges}
        # delegator can be self for native nodes.
        self.delegated_edges = {self: {}}


    def top_delegate(self, ):
        """ Find the top delegate node. """
        if self.delegatee is None:
            return self
        else:
            return self.delegatee.top_delegate()

    def __str__(self, ):
        return "I:" + str(self.domain_end)

    def __repr__(self, ):
        return str(self) # + " at " + str(hex(id(self)))

    def __lt__(self, other_node):
        return str(self) < str(other_node)


    def print_delegate_info(self, g=None):
        """
        print delegation information for this node.
        """
        if self.delegatee is not None:
            print(self, "edges delegated to", self.delegatee)
            if g:
                print(" - Empty g.adj[self]:", len(g.adj[self]) == 0)
        else:
            print(self, "delegated edges (including self):")
            for delegator, target_nodes in self.delegated_edges.items():
                print(" - %s: %s" % (delegator, target_nodes))
                if g:
                    print(" -- g.adj[delegator] is empty:", len(g.adj[delegator]) == 0)
                    print(" -- delegator has no incoming edges:",
                          all(delegator not in g.adj[target] for target in target_nodes))






class Domain5p3pGraph(nx.Graph):
    """
    Network graph where strands and connections are represented by the 5' and 3' end of domains.
    """

    def init_from_domain_graph(self, graph):
        is_directed = getattr(graph, 'is_directed', False)

        for domain, attrs in graph.nodes(data=True):
            self.add_node(str(domain)+":5p", attrs)
            self.add_node(str(domain)+":3p", attrs)
            self.add_edge(str(domain)+":5p", str(domain)+":3p", interaction='pb')

        for source, target, attrs in graph.edges(data=True):
            interaction = attrs.get('interaction')
            if interaction in (1, 'pb', 'p', 'backbone', 3, 'bs', 's', 'stacking'):
                # For edges connecting domain backbone, we connect
                # the 3p end of the 5p-domain to the 5p domain of the 3p domain
                # 5'-domainA-3' 5'-domainB-3'
                # If input graph is not directed, how to know which of source or target is the 5p vs 3p end?
                if is_directed:
                    self.add_edge(source+":3p", target+":5p", attrs)
                else:
                    upstream = attrs.get('upstream')
                    if upstream not in (source, target):
                        raise ValueError("Graph edge %s---%s %s is not directed and has no upstream designation." %
                                         (source, target, attrs))
            elif interaction in (2, 'h', 'hh', 'hybridization'):
                # Always connecting 3p of [a] to 5p end of [A] vice versa:
                self.add_edge(target+":3p", source+":5p", attrs) # Uh, this needs to either be target.model3p or str(target)+":3p".
                self.add_edge(source+":3p", target+":5p", attrs)
            else:
                print("Unrecognized interaction '%s' between '%s' and '%s'" %
                      (attrs.get('interaction'), source, target))
                self.add_edge(target+":3p", source+":5p", attrs)

    def init_from_strands(self, strands):
        """ Initialize graph from strands. """
        for strand in strands:
            domains = strand.domains
            self.add_nodes_from([node for d in domains for node in (d.model5p, d.model3p)])
            # Add a5p->a3p connections:
            for d in domains:
                self.add_edge(d.model5p, d.model3p, interaction=PHOSPHATEBACKBONE_INTERACTION, type="5p3p")
            if len(domains) > 1:
                for upstream, downstream in zip(domains[:-1], domains[1:]): # source, target
                    self.add_edge(upstream.model3p, downstream.model5p, interaction='pb', weight=1, upstream=upstream)

    def hybridize(self, dom1, dom2):
        self.add_edge(dom1.model5p, dom2.model3p)
        self.add_edge(dom2.model5p, dom1.model3p)

    def dehybridize(self, dom1, dom2):
        self.remove_edge(dom1.model5p, dom2.model3p)
        self.remove_edge(dom2.model5p, dom1.model3p)


class Domain5p3pDiGraph(nx.DiGraph):
    """
    Network graph where strands and connections are represented by the 5' and 3' end of domains.
    """

    def init_from_strands(self, strands):
        """ Initialize graph from strands. """
        for strand in strands:
            domains = strand.domains
            self.add_nodes_from([node for d in domains for node in (d.model5p, d.model3p)])
            # Add a5p->a3p connections:
            for d in domains:
                self.add_edge(d.model5p, d.model3p, interaction=PHOSPHATEBACKBONE_INTERACTION, type="5p3p")
            if len(domains) > 1:
                for upstream, downstream in zip(domains[:-1], domains[1:]): # source, target
                    self.add_edge(upstream.model3p, downstream.model5p, interaction=PHOSPHATEBACKBONE_INTERACTION)

    def hybridize(self, dom1, dom2):
        self.add_edge(dom1.model5p, dom2.model3p)
        self.add_edge(dom2.model5p, dom1.model3p)

    def dehybridize(self, dom1, dom2):
        self.remove_edge(dom1.model5p, dom2.model3p)
        self.remove_edge(dom2.model5p, dom1.model3p)


class DomainGraph(nx.Graph):
    """
    Network graph where strands and connections are represented by domains.
    """
    def init_from_strands(self, strands):
        """ Initialize graph from strands. """
        for strand in strands:
            domains = strand.domains
            self.add_nodes_from(domains)
            if len(domains) > 1:
                for edge in zip(domains[:-1], domains[1:]): # source, target
                    self.add_edge(*edge, interaction='pb', weight=sum(d.length for d in edge)/2,
                                  upstream=edge[0])

    def hybridize(self, dom1, dom2):
        """ Hybridize domain :dom1: and domain :dom2: """
        self.add_edge(dom1, dom2)

    def dehybridize(self, dom1, dom2):
        """ De-hybridize domain :dom1: and domain :dom2: """
        self.remove_edge(dom1, dom2)


class DomainDiGraph(nx.DiGraph):
    """
    Directed graph where strands and connections are represented by domains.
    A directed graph is well-suited for denoting the 5'->3' directionality of strands.
    However, it is not very well-suited for distance calculations.
    """
    def init_from_strands(self, strands):
        """ Initialize graph from strands. """
        for strand in strands:
            domains = strand.domains
            self.add_nodes_from(domains)
            if len(domains) > 1:
                for edge in zip(domains[:-1], domains[1:]): # source, target
                    self.add_edge(*edge, interaction='pb', weight=sum(d.length for d in edge)/2)

    def hybridize(self, dom1, dom2):
        """ Hybridize domain :dom1: and domain :dom2: """
        self.add_edge(dom1, dom2)
        self.add_edge(dom2, dom1)

    def dehybridize(self, dom1, dom2):
        """ De-hybridize domain :dom1: and domain :dom2: """
        self.remove_edge(dom1, dom2)
        self.remove_edge(dom2, dom1)


class StrandGraph(nx.MultiGraph):
    """
    Network graph where strands and connections are represented simply by the strands.
    """
    def init_from_strands(self, strands):
        """ Initialize graph from strands. """
        self.add_nodes_from(strands)

    def hybridize(self, dom1, dom2):
        """ Hybridize domain :dom1: and domain :dom2: """
        self.add_edge(dom1.strand, dom2.strand)

    def dehybridize(self, dom1, dom2):
        """ De-hybridize domain :dom1: and domain :dom2: """
        self.remove_edge(dom1.strand, dom2.strand)


class HelixGraph(nx.MultiGraph):
    """
    Network graph representing dsDNA helices.
    Each node is a ds helix.
    The edges connecting nodes/helices are phosphate backbone connections between domains on the same strand.
    If two nodes/helices have more than one edge connecting them, then the two helices
    must be parallel (at least in the helix segments between the two connections).
    Note: Make sure not to duplicate functionality of structural_elements.helix.DsHelix
    """

    def join(self, helix1, helix2):
        pass

    def disjoin(self, helix, now_sure_how_to_do_this):
        pass

    def hybridize(self, dom1, dom2):
        """ Hybridize domain :dom1: and domain :dom2: """
        self.add_edge(dom1.strand, dom2.strand)

    def dehybridize(self, dom1, dom2):
        """ De-hybridize domain :dom1: and domain :dom2: """
        self.remove_edge(dom1.strand, dom2.strand)


class HelixEndsGraph(nx.MultiGraph):
    """
    Network graph representing dsDNA helices.
    """
    pass


class Sheet2DGraph(nx.MultiGraph):
    """
    Network graph where strands and connections are represented by the 5' and 3' end of domains.
    """
    pass

class HelixBundleGraph(nx.MultiGraph):
    """
    Network graph where the structure is represented by connected helix bundles.
    """
    pass
