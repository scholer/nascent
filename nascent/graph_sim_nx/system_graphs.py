# -*- coding: utf-8 -*-
##    Copyright 2015-2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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
 - DomainEnds graph
 - Interface graph (where hybridized/stacked DomainEnds are merged/collapsed into a single node.)
 - Domains graph
 - Strands graph


"""

from __future__ import absolute_import, print_function, division
from collections import defaultdict
import networkx as nx
import pdb
from pprint import pprint

# Local imports:
from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)



class InterfaceMultiGraph(nx.MultiGraph):
    """
    MultiGraph version of InterfaceGraph.
    Using a MultiGraph make delegation/undelegation MUCH easier, since we can
    use the original source delegator as key in MultiGraph.adj.
    Then we don't have to spend as much time determining whether edges should
    be removed or not in undelegate and which graph should be used for representation.


    ## What edge key to use? ##

    For domain long edges (connecting 5p of a domain to 3p of the same domain):
    - Use 0 or "", since we want duplexed domain edges be represented as a single edge in the InterfaceGraph.
    Edge key for backbone connections between the 3' of one domain (h1end3p) and the 5' domain of downstream domain (h1end5p):
    - Originally I thought simply a tuple of (h1end3p, h1end5p) would suffice. However this would give edge degeneracies
        for a palindromic duplex:   A---a hybridizing to A---a.
    - So now i'm thinking keying by
    - Obviously this is not cacheable, so for caching the edge key I'm thinking to use:
        tuple((end.name, end.state_fingerprint()) for end in edge_key) if end else 0
      Then if there is multiple edges from source to target ifnode (when recreating a loop path),
       first check end.name part, and if that is still ambigious (RARE EVENT), then check the end.state_fingerprint() part.
      Now, obviously that will require Complex to keep an up-to-date domainend_by_hash index, or figure out a way to
      do a "minimal" update (i.e. predicting if the index will be needed). This is probably now something I want to
      implement right *now*, but at least it's an option, so the solution should be future-proof.

    Note that interface_graph is typically derived from ends5p3p_graph (MultiDiGraph), keyed by interaction type.


    ## How does delegation work for MultiGraphs? ##

    Stacking of a strands with two duplexes each to form a single double-helix:
            0------1---2-------3     0---.           .----3     0---.  1        .---3
                                          `.1---2---:                `.----2---:
            4------5---6-------7     4----´ 5   6    `---7      4----´ 5   6    `---7
        Aka:

            A3----A5---B3-----B5     A3--.           .--B5     A3--.  A5       .--B5
                                          `.A5--B3--:               `.----B3--:
            a5----a3---b5-----b3     a5---´ a3  b5   `--b3     a5---´ a3  b5   `--b3
        Where a3 was delegated to A5, b5 delegated to B3, and then A5 delegated to B3.
        Delegation hierarchy/tree:      .- b5
                                 B3 <-<ˊ
                                       `·- A5 <-- a3

    All nodes retains their delegtated edges even after delegation to a delegatee above. Thus
    the a3.delegated_edges is a subset of A5.delegated_edges, A5.delegated_edges is a subset of B3.delegated_edges, etc:

    Q: Should delegated_edges contain edge_keys?
        (a) Yes: Delegated nodes (for each ifnode): {delegator: {key: {target: {edge_key: edge_attrs}}}}

        where edge_key = (h1end3p, h1end5p) # i.e. 3p-to-5p


    # The structure of delegated_nodes: #
    Attention: While the structure of delegated_edges looks deceptively similar to MultiGraph.adj, the edge groups
    are not the same dict objects. Modifying adj[delegator][target] will not modify the delegator.delegated_edges.
    It may be nice to emphasize this difference by flattening the structure of delegated_edges to be
    delegated_edges[delegator][(target, key)] rather than delegated_edges[delegator][target][key]

        ifnode.delegated_edges = {
            delegator1: {
                target1: {key1: eattrs1},
            }
        }
        b3.delegated_edges = {
            'B3': {
                'B5': {'': {}},
                'A5': {('B3','A5'): {}}
            },
            'b5': {
                'b3': {'': {}},
                'a3': {('a3', 'b5'): {}}
            },
            'A5': {
                'A3': {'': {}},
                'B3': {('B3','A5'): {}}
            },
            'a3': {
                'a5': {'': {}},
                'b5': {('a3', 'b5'): {}}
            }
        }
        b5.delegated_edges = {
            'b5': {
                'b3': {'': {}},
                'a3': {('a3', 'b5'): {}}
            },
        }
        A5.delegated_edges = {
            'A5': {
                'A3': {'': {}},
                'B3': {('B3','A5'): {}}
            },
            'a3': {
                'a5': {'': {}},
                'b5': {('a3', 'b5'): {}}
            }
        }
        a3.delegated_edges = {
            'a3': {
                'a5': {'': {}},
                'b5': {('a3', 'b5'): {}}
            }
        }

    Note that whether edges are a multi-edge dict, {key: eattr}, or just a single edge attr, makes little difference,
    since we never overwrite anything in the delegatee.

    """

    def __init__(self, data=None, **attr):
        # node attr key to increment to reflect the number of DomainEnd nodes represented by a the top_delegate ifnode:
        self.node_delegation_size_key = 'size'
        super().__init__(data, **attr)


    def reset_ifnode_delegation_to_current_graph_representation(self):
        """
        Update all ifnode.delegated_edges to reflect the current state of the graph, that is:
            ifnode.delegated_edges = {self: self.adj[self]}

        """
        print("Resetting ifnode delegation to:")
        pprint(self.adj)
        # Make sure you make COPIES of the edge groups, so we don't risk modifying the delegated edges in-place
        # when updating graph.adj representation during delegate/undelegate()
        for delegatee, target_edgegroups in self.adj.items():
            delegatee.delegated_edges = {
                delegatee: {
                    target: {key: eattr for key, eattr in edge_group.items()}
                    for target, edge_group in target_edgegroups.items()
                }
            }
            print("%s.delegated_edges set to:" % delegatee, delegatee.delegated_edges)


    def merge(self, node1, node2):
        """
        Merge representation of node1 and node2.
        Returns the delegatee.
        Mostly equivalent to calling self.delegate(node1, node2), but will also set node attributes e.g. size.
        Synonyms: merge, fuse, join, unite, connect, link, collect, combine, collapse,
        """
        self.delegate(node1, node2)
        try:
            self.node[node2][self.node_delegation_size_key] += 1
            self.node[node1][self.node_delegation_size_key] -= 1
        except KeyError:
            self.node[node2][self.node_delegation_size_key] = 2
            self.node[node1][self.node_delegation_size_key] = 0
        return node2


    def split(self, node1, node2):
        """
        Undo collapse/merge of node1 and node2 by determining which node is delegator and which is delegatee,
        and then calling undelegate(delegator, delegatee) accordingly.
        Returns the delegatee.
        Synonyms: split, part, separate, divide, disjoin, sever, slice, disunite
        """
        if node1.delegatee is None:
            # Node1 does not have any delegatee, so it should itself be the top delegate:
            delegator, delegatee = node2, node1
        elif node2.delegatee is None:
            # Node2 does not have any delegatee, so it should itself be the top delegate:
            delegator, delegatee = node1, node2
        elif node1.delegatee is node2:
            # We might be exchanging the delegation of nodes that are already delegated. This should be OK, but isn't really supported.
            delegator, delegatee = node1, node2
            assert node2.delegatee is not node1
        elif node2.delegatee is node1:
            delegator, delegatee = node2, node1
        else:
            raise ValueError("Either node1 must be delegatee of node2 or node2 must be delegatee of node1.")
        self.undelegate(delegator, delegatee)
        self.node[delegator][self.node_delegation_size_key] += 1
        self.node[delegatee][self.node_delegation_size_key] -= 1
        return delegatee


    def delegate(self, delegator, delegatee):
        """
        Delegate graph responsibility from delegator node to delegatee node.
        Both nodes should be contained in this graph and both should be InterfaceNode instances.
        This is a little different from InterfaceGraph.delegate since we have a MultiGraph:
        - No longer any need to keep an "edge_count" edge attribute (we retain ALL edges... except domain edges..)
        """
        print("Delegating delegator %s to delegatee %s..." % (delegator, delegatee))
        # Networkx API:
        # Graph.edge[source][target] = edge_attr_dict
        # MultiGraph.edge[source][target][key] = edge_attr # adj = {'src': {'tgt': {'key': {'eattr1': 0, 'eattr2': 1}}}}

        ## 0. Make sure everything looks right:
        assert delegator.delegatee is None
        assert delegator not in delegatee.delegated_edges
        # Make sure delegator and delegatee haven't been delegated edges from the same sub-delegator:
        # assert len(set(delegator.delegated_edges.keys()) & set(delegatee.delegated_edges.keys())) == 0
        assert delegatee.delegated_edges.keys().isdisjoint(delegator.delegated_edges.keys())


        # 1a: set delegator's delegatee attribute:
        delegator.delegatee = delegatee
        ## 1b. Update delegatee.delegated_edges with entries from delegator.
        # Keep delegator.delegated_edges intact.
        delegatee.delegated_edges.update(delegator.delegated_edges)
        # The above is actually updating self.adj[delegator] IN PLACE!
        # Be very carefull when updating delegated_edges[delegator] vs adj[source] - they may be the same dict.


        # Concern: How do we ensure that we don't modify the edge_groups in delegator.delegated_edges in-place?
        # Maybe we should just flatten ifnode.delegated_edges[target][key] = eattr
        # to delegated_edges[(target, key)] = eattr to avoid confusion?
        # Edit: edge_groups in self.adj[source][target] are COPIES of the values in delegated_edges, so it shouldn't matter.

        ## 2. Update graph representation:
        ## 2a. Move all incoming edges from delegator to delegatee
        for target, edges in self.adj[delegator].items():
            # We could just do a self.adj[delegatee].update(edges), but...
            if target not in self.adj[delegatee]: # We can move the full edge_group from delegator to delegatee..
                assert delegatee not in self.adj[target]
                self.adj[delegatee][target] = edges
                self.adj[target][delegatee] = edges
            else:
                # We have to evaluate the move of each edge one by one, depending on whether an edge with same key is
                # already present in adj[delegatee][target]
                assert delegatee in self.adj[target]
                for key, eattr in edges.items():
                    if key not in self.adj[delegatee][target]:
                        # print("Making new edge from delegatee %s to target %s" % (delegatee, target))
                        self.adj[delegatee][target][key] = eattr
                        self.adj[target][delegatee][key] = eattr
                        if 'edge_count' not in self.adj[target][delegatee][key]:
                            self.adj[target][delegatee][key]['edge_count'] = 1
                            # Do you *know* that edge_count is 1??
                            # If the edge already represents a "collapsed" representation of e.g. a duplex, edge_count
                            # would be 2. Thus, we only set it to 1 if it isn't set yet.
                    else:
                        # We are merging two edges with the same key. This should be e.g. None or "" used by
                        # domain-edges to explicitly get just a single edge between duplex ends.
                        # We do not try to attempt to determine which edge is "best" since they should be equal.
                        assert key in self.adj[target][delegatee]
                        # print("delegatee %s existing target/edge: %s/%s" % (delegatee, target, self.adj[target][delegatee]))
                        # assert eattr['len_contour'] == self.adj[delegatee][target]['len_contour']
                        # If the above assertion fails, you probably have to make sure you select the shortest edge.
                        # Thought: Could we update edge len_contour automatically? Overlapping edges (same src, tgt, key)
                        # should only happen for duplex ends. When that happends, assign eattr['len_contour'] = eattr['ds_len'] or similar?
                        try:
                            self.adj[target][delegatee][key]['edge_count'] += 1
                        except KeyError:
                            self.adj[target][delegatee][key]['edge_count'] = 2
            # Concern: Should we delete or just update? Is it safe to update the edge group?
            # The edge group may have been moved fully to delegatee...
            # In fact, how do we ensure that we don't modify the edge_groups in delegator.delegated_edges in-place?
            del self.adj[target][delegator]
        # 2b. Remove all outgoing edges from delegator:
        self.adj[delegator] = {}




    def undelegate(self, delegator, delegatee):
        """
        Reverse the effect of delegate.
        That is, take all the edges that delegator has delegated to delegatee and transfer them back to delegator.

        # Two approaches: - Using strategy (a)
        #  (a) graph adj-based: Look at delegatee's current edge targets in graph.adj and determine if
        #       ( i) there should be an edge from delegator to that target, and
        #       (ii) if the edge from delegatee to target should still be in place.
        #  (b) ifnode delegated_edges-based: Re-generate delegator edges from scratch.
        #       For each original source-target edge in
        #       delegator.delegated_edges = {source: {target: {}, ...}},
        #       add an edge from delegator to the top-delegate for all targets.
        #       (It should be ok to assume that delegator is now the top-delegatee for all sources, but check anyways.)
        """
        print("Undelegating delegator %s from delegatee %s..." % (delegator, delegatee))

        ## 0. Make sure everything looks right from the start:
        assert delegator.delegatee is delegatee # delegatee is indeed the ifnode that delegator has delegated to
        assert delegator in delegatee.delegated_edges # delegatee knows that is is representing edges from delegator
        assert len(self.adj[delegator]) == 0 # delegator should not currently have any edges

        ## 1. Update delegatee.delegated_edges, removing all delegated entries now belonging to delegator again.
        ## delegator.delegated_edges contains all edges that was delegated to delegator (incl itself) before they were redelegated to delegatee.
        for k in delegator.delegated_edges:
            # if an edge_group was delegated to delegator, then delegatee must have obtained it through delegator
            del delegatee.delegated_edges[k]
        # And reset delegator's delegatee attribute (marking that delegator no longer has a parent delegatee).
        delegator.delegatee = None


        assert all(source.top_delegate() is delegator for source in delegator.delegated_edges)  ## TODO: Remove check

        ## 2. Update graph representation:
        ## Things to be aware of:
        ##  a) Targets in delegated_edges may not currently be the top_delegate; make sure when you create edges
        ##      to delegator that they go to the top ifnode.

        # Make a set of the edges that was reclaimed by delegator:
        delegator_reclaimed_edges = {(target.top_delegate(), key)
                                     # edgegroups_by_target = out_edges = {target1: {key1: {}, ..}, ..}
                                     for source, edgegroups_by_target in delegator.delegated_edges.items()
                                     for target, edge_group in edgegroups_by_target.items()
                                     for key in edge_group}
        # Determine which edges (target, key) should also still be in delegatee, i.e. when an edge from delegator
        # has the same (target, key) as another of delegatee's edges.
        remaining_delegatee_edges = {(target.top_delegate(), key)
                                     for source, edgegroups_by_target in delegatee.delegated_edges.items()
                                     for target, edge_group in edgegroups_by_target.items()
                                     for key in edge_group}
        # For MultiGraph strategy there should be NO overlap when we are naming the edges between domains:
        # Edit: There should, for un-keyed edges, i.e. overlapping duplex edges:
        overlapping_edges = delegator_reclaimed_edges & remaining_delegatee_edges
        if len(overlapping_edges) > 0:
            print("After undelegation, the following delegagtor target-key edges also remain valid for delegatee:",
                  overlapping_edges)
        # Edit: What about domain (long, traversing) edges? - target should differ...
        # However, the original assertion may still be valid: checking targets and completely moving full sets
        # of edges: {key1: eattrs2, key2: eattrs2}
        remaining_delegatee_targets = {target.top_delegate()
                                       for source, targets in delegatee.delegated_edges.items()
                                       for target in targets}

        # Structure with delegator's edges to top ifnodes: [target-top-ifnode][source aka self or sub-delegator] = eattr
        # (Reversed order of target and source so we can do a quick "if delegator in delegator_targets[target]" below)
        delegator_targets = defaultdict(dict)  # [target][source] = {key1: eattr1, key2: eattr2}
        for source, targets in delegator.delegated_edges.items():
            for target, edge_group in targets.items():
                delegator_targets[target.top_delegate()][source] = edge_group

        ## Approach 2a: Move or "copy" edges from self.adj
        ## - using the current graph representation rather than edges in delegator.delegated_edges.
        ## This avoids problems having to determine which edges in delegated_edges should actually be visible,
        ## and may be slightly faster than approach 2b since in some cases we can just take the full edge_groups dict
        ## and move it from delegatee to delegator.
        # Strategy concern: Should you try to move whole (source, target) edge groups, or just take each edge one at a time?
        # --> I'll first try to move complete edge groups and if I can't, then I'll consider each edge individually.
        # Edge groups that should possibly be moved:
        edge_groups = [(target, edge_group) for target, edge_group in self.adj[delegatee].items() if target in delegator_targets]
        # print("undelegate delegator %s  from  delegatee %s" % (delegator, delegatee))
        # print("delegator_targets:", delegator_targets)
        # print("remaining_delegatee_targets:", remaining_delegatee_targets)
        # print("Re-considering %s edges: [%s]", (delegatee, edges))
        # Note: We need to support

        # Concern: How do we ensure that we don't modify the edge_groups in delegator.delegated_edges in-place?
        # We could "always copy edge_attrs only" and NOT try to modify edge_groups in place..
        # But that requires self.adj[source][target] edge group to be a COPY of the edge group in delegated_edges.
        # OK, so let's guarantee that.

        for target, edge_group in edge_groups:
            if target in remaining_delegatee_targets:
                # We cannot simply move the edge_group dict from self.adj[delegatee][target] and assign it fully to delegator
                # Instead, go over each edge (target, key) and either move it (if the edge is not in remaining_delegatee_edges
                # or find another edge to use by looking at top ifnodes in delegator.delegated_edges (=delegator_targets).
                print("Target %s still in remaining delegatee (%s) targets, moving individual edges (target, key)..." % (target, delegatee))
                for key, eattr in edge_group.items():
                    # This is basically the "for target, eattr in edges" loop in regular InterfaceGraph..
                    # See if you can move eattr from delegatee to delegator:
                    if (target, key) in remaining_delegatee_edges:
                        # We cannot simply remove the eattr dict from self.adj[delegatee][target]
                        # Find the best possible edge to represent (unusual, given that edges should have distinct keys even after delegation)
                        if delegator in delegator_targets[target]:
                            # If delegator is a legitimate source (and not a sub-delegator), then use the original delegator-to-target edge.
                            other_eattr = delegator_targets[target][delegator][key]
                        else:
                            # pick a random eattr:
                            other_eattr = next(eattr for egroup in delegator_targets[target].values() for eattr in egroup.values())
                        print("(%s, %s) still in remaining_delegatee_edges, using other eattr created from delegated_edges: %s" % (target, key, other_eattr))
                        #self.adj[delegator].setdefault(target, {})[key] = other_eattr
                        # self.adj[target].setdefault(delegator, {})[key] = other_eattr
                        # Make sure [delegator][target] and to use the same group_dict, not just equivalent:
                        try:
                            self.adj[delegator][target][key] = other_eattr
                        except KeyError:
                            new_group_dict = {}
                            self.adj[delegator][target] = self.adj[target][delegator] = new_group_dict
                            self.adj[target][delegator][key] = other_eattr
                            assert other_eattr is self.adj[delegator][target][key]
                        else:
                            self.adj[target][delegator][key] = other_eattr
                    else: # Move edge back from delegatee to delegator:
                        print("Moving edge eattr to target %s from delegatee %s to delegator %s" % (target, delegatee, delegator))
                        self.adj[delegator][target][key] = eattr
                        assert eattr is self.adj[target][delegator][key]
                        if self.adj[delegatee][target] is self.adj[target][delegatee]:
                            # Deleting eattr from one edge_group should also delete it from the other as they are the same:
                            del self.adj[delegatee][target][key]
                            assert key not in self.adj[target][delegatee]
                        else:
                            print("\nWARNING: %s.adj[%s][%s] edge group dict differs from %s.adj[%s][%s] edge group dict!\n" % (delegatee, target, target, delegatee))
                            del self.adj[delegatee][target][key]
                            if target is not delegatee:
                                # If the edge is a self-loop (target is delegatee) we would get a KeyError for removing twice:
                                del self.adj[target][delegatee][key]
            else:
                # We can completely move the whole edge group in adj[delegatee][target] from delegatee to delegator:
                # print("Moving edge to target %s from delegatee %s to delegator %s" % (target, delegatee, delegator))
                self.adj[delegator][target] = edge_group
                self.adj[target][delegator] = edge_group
                # self.adj[target][delegator] = eattr
                del self.adj[delegatee][target] # Do we delete edge_groups or just clear them?
                # If the edges are self-loops (target is delegatee) we would get a KeyError for removing twice,
                # so check that before deleting (or use
                if target is not delegatee:
                    del self.adj[target][delegatee]


        # Alternative approach, revisited:
        # First re-evaluate delegatee remaining nodes, deleting those that shouldn't still be there,
        # then create new edges for delegator.

        ## Approach 2b: Re-generate delegator edges from delegator.delegated_edges:
        ## (This may be slightly simpler, especially when dealing with edge-groups in multi-graphs)
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
        #     # Add edge from delegator to target_delegatee (or increment edge count)
        #     if target_delegate in self.adj[delegator]:
        #         self.adj[delegator][target_delegate]['edge_count'] += 1
        #     else:
        #         self.adj[delegator][target_delegate] = eattr
        #         try:
        #             assert self.adj[delegator][target_delegate]['edge_count'] == 1
        #         except KeyError:
        #             self.adj[delegator][target_delegate]['edge_count'] = 1

    #
    # def add_edges_from(self, ebunch, attr_dict=None, **attr):
    #     """
    #     Add edges as "native" edges between nodes.
    #     """
    #     ebunch = list(ebunch) # In case it is a generator...
    #     print(locals())
    #     pdb.set_trace()
    #     super(InterfaceMultiGraph, self).add_edges_from(ebunch, attr_dict, edge_count=1, **attr)
    #     #print("Updating IF nodes delegated edges for ebunch %s..." % (ebunch,))
    #     for e in ebunch:
    #         u, v = e[:2]
    #         u.delegated_edges[u][v] = self.adj[u][v]
    #         #print("%s.delegated_edges: %s" % (u, u.delegated_edges))
    #         v.delegated_edges[v][u] = self.adj[v][u]
    #         #print("%s.delegated_edges: %s" % (v, v.delegated_edges))
    #
    #
    # def add_edge(self, u, v, key, attr_dict=None, **attr):
    #     """
    #     Add edges as "native" edges between nodes.
    #     Note: This MAY cause big problems if you have already done delegations and you then add_edge(...)
    #     (which will overwrite u,v.delegated_edges
    #     """
    #     print("Updating delegated edges for IF nodes  %s and %s..." % (u, v))
    #     super(InterfaceGraph, self).add_edge(u, v, attr_dict, edge_count=1, **attr)
    #     ## Q: How much do we want the "representation" eattrs to be tied to eattrs in delegated_edges?
    #     u.delegated_edges[u][v] = self.adj[u][v]
    #     v.delegated_edges[v][u] = self.adj[v][u]
    #     # print("%s.delegated_edges: %s" % (u, u.delegated_edges))
    #     # print("%s.delegated_edges: %s" % (v, v.delegated_edges))

    def top_delegate(self, node):
        """
        Find the top delegate node. This graph-level method relies on in-graph node attribute 'delegatee',
        i.e. instead of ifnode.delegatee object attribute, it calls self.node[node]['delegatee'] dict entry.
        """
        try:
            node_delegatee = self.node[node]['delegatee']
        except KeyError:
            self.node[node]['delegatee'] = node_delegatee = None
        if node_delegatee is None:
            return self
        else:
            return self.top_delegate(node_delegatee)


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



class InterfaceGraph(nx.Graph):
    """
    Graph or MultiGraph?
     - Hybridization and stacking interaction edges in ends5p3p graph will always be
        represented by a node merge/delegation; the only edges we have are backbone edges.

     - However, using a MultiGraph would make delegation/undelegation MUCH easier, since we could
        use the original source delegator as key in MultiGraph.adj.
        Then we wouldn't have to spend as much time determining whether edges should be removed or not.

    Pro MultiGraph: Loops will be easier to handle.
     - Not having a multi-edged InterfaceGraph makes it slightly harder to deal with loops of just a few nodes:
        These loops quickly turns degenerate when merging nodes.
         For instance, loops consisting of just two nodes are effectively not a loop but just an edge.

    """
    # TODO, WIP: Convert InterfaceGraph to a MultiGraph - Edit: I retain both

    def merge(self, node1, node2):
        """
        Merge representation of node1 and node2.
        Returns the delegatee.
        Mostly equivalent to calling self.delegate(node1, node2), but will also set node attributes e.g. size.
        Synonyms: merge, fuse, join, unite, connect, link, collect, combine, collapse,
        """
        self.delegate(node1, node2)
        try:
            self.node[node2]['size'] += 1
            self.node[node1]['size'] -= 1
        except KeyError:
            self.node[node2]['size'] = 2
            self.node[node1]['size'] = 0
        return node2

    def split(self, node1, node2):
        """
        Undo collapse/merge of node1 and node2 by determining which node is delegator and which is delegatee,
        and then calling undelegate(delegator, delegatee) accordingly.
        Returns the delegatee.
        Synonyms: split, part, separate, divide, disjoin, sever, slice, disunite
        """
        # if node1.delegatee is not None:
        #     # This is not sufficient check; both node1 and node2 can have a delegatee, which itself has a delegatee
        #     delegator, delegatee = node1, node2
        # else:
        #     delegator, delegatee = node2, node1
        if node1.delegatee is None:
            # This, OTOH, should be a valid check. Node1 does not have any delegatee,
            # so it should itself be the top delegate
            delegator, delegatee = node2, node1
        elif node2.delegatee is None:
            # This, OTOH, should be a valid check.
            delegator, delegatee = node1, node2
        elif node1.delegatee is node2:
            delegator, delegatee = node1, node2
            assert node2.delegatee is not node1
        elif node2.delegatee is node1:
            delegator, delegatee = node2, node1
        else:
            raise ValueError("Either node1 must be delegatee of node2 or node2 must be delegatee of node1.")
        self.undelegate(delegator, delegatee)
        self.node[delegator]['size'] += 1
        self.node[delegatee]['size'] -= 1
        return delegatee


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
        # assert delegatee.delegated_edges.keys().isdisjoint(delegator.delegated_edges)

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
                # Although they should be equal... # Edit: Why is that again? - Asserting to make sure...
                assert eattr['len_contour'] == self.adj[delegatee][target]['len_contour']
                # If the above assertion fails, you probably have to make sure you select the shortest edge.
                try:
                    self.adj[target][delegatee]['edge_count'] += 1
                except KeyError:
                    self.adj[target][delegatee]['edge_count'] = 2
            del self.adj[target][delegator]
            # del self.adj[delegator][target]  # cannot delete entries from dict while iterating over it.
            # (but, you could have converted the dict.keys() iterator to a list and THEN delete entries in-place.)
        # 2b. Remove all outgoing edges from delegator:
        self.adj[delegator] = {} # Do not use self.adj[delegator].clear(); that will clear target_nodes in-place
        # delegatee.delegated_edges[delegator] = target_nodes
        # print("edges to %s delegated from %s to %s" %
        #       (list(delegatee.delegated_edges[delegator].keys()), delegator, delegator.delegatee))


    def undelegate(self, delegator, delegatee):
        """
        Reverse the effect of delegate.
        """

        ## 0. Make sure everything looks right:
        assert delegator.delegatee is delegatee
        assert delegator in delegatee.delegated_edges
        ## Which one of these do we prefer? Probably the latter; we keep delegator node in the graph after all.
        # assert delegator not in self.adj or self.adj[delegator] is None or len(self.adj[delegator]) == 0
        assert len(self.adj[delegator]) == 0

        ## 1. Update delegatee.delegated_edges, removing all delegated entries now belonging to delegator again.
        for k in delegator.delegated_edges:
            # if an edge was delegated to delegator, then delegatee must have obtained it through delegator
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
        # Structure with delegator's edges to top ifnodes: [target-top-ifnode][source aka self or sub-delegator] = eattr
        # (Reversed order of target and source so we can do a quick "if delegator in delegator_targets[target]" below)
        delegator_targets = defaultdict(dict)
        for source, targets in delegator.delegated_edges.items():
            for target, eattr in targets.items():
                delegator_targets[target.top_delegate()][source] = eattr
        remaining_delegatee_targets = {target.top_delegate()
                                       for source, targets in delegatee.delegated_edges.items()
                                       for target in targets}

        ## Approach 2a: Move or "copy" edges from self.adj:
        ## (This may be slightly faster than approach 2b since in some cases we can remove the eattr dict from delegatee to delegator.
        edges = [(target, eattr) for target, eattr in self.adj[delegatee].items() if target in delegator_targets]
        # print("undelegate delegator %s  from  delegatee %s" % (delegator, delegatee))
        # print("delegator_targets:", delegator_targets)
        # print("remaining_delegatee_targets:", remaining_delegatee_targets)
        # print("Re-considering %s edges: [%s]", (delegatee, edges))
        for target, eattr in edges:
            if target in remaining_delegatee_targets:
                # We cannot simply remove the eattr dict from self.adj[delegatee][target]
                # Instead, find another eattr dict in delegator_targets and use that to update self.adj[delegator]
                # print("Target %s still in remaining delegatee targets, making edge using other eattr..." % target)
                if delegator in delegator_targets[target]:
                    # If delegator is a legitimate source (and not a sub-delegator), then use the original delegator-to-target edge.
                    new_eattr = delegator_targets[target][delegator]
                else:
                    # No "original" delegator-to-target edge is available, all edges are adopted from sub-delegators;
                    # Use edge with highest edge_count as representation in the graph:
                    new_eattr = max(delegator_targets[target].values(), key=lambda attr: attr.get('edge_count', 0))
                self.adj[delegator][target] = new_eattr
                self.adj[target][delegator] = new_eattr
            else: # Move edge back from delegatee to delegator:
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
        #     # Add edge from delegator to target_delegatee (or increment edge count)
        #     if target_delegate in self.adj[delegator]:
        #         self.adj[delegator][target_delegate]['edge_count'] += 1
        #     else:
        #         self.adj[delegator][target_delegate] = eattr
        #         try:
        #             assert self.adj[delegator][target_delegate]['edge_count'] == 1
        #         except KeyError:
        #             self.adj[delegator][target_delegate]['edge_count'] = 1


    def add_edges_from(self, ebunch, attr_dict=None, **attr):
        """
        Add edges as "native" edges between nodes.
        """
        ebunch = list(ebunch) # In case it is a generator...
        super(InterfaceGraph, self).add_edges_from(ebunch, attr_dict, edge_count=1, **attr)
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
        super(InterfaceGraph, self).add_edge(u, v, attr_dict, edge_count=1, **attr)
        ## Q: How much do we want the "representation" eattrs to be tied to eattrs in delegated_edges?
        u.delegated_edges[u][v] = self.adj[u][v]
        v.delegated_edges[v][u] = self.adj[v][u]
        # print("%s.delegated_edges: %s" % (u, u.delegated_edges))
        # print("%s.delegated_edges: %s" % (v, v.delegated_edges))

    def top_delegate(self, node):
        """
        Find the top delegate node. This graph-level method relies on in-graph node attribute 'delegatee',
        i.e. instead of ifnode.delegatee object attribute, it calls self.node[node]['delegatee'] dict entry.
        """
        try:
            node_delegatee = self.node[node]['delegatee']
        except KeyError:
            self.node[node]['delegatee'] = node_delegatee = None
        if node_delegatee is None:
            return self
        else:
            return self.top_delegate(node_delegatee)


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




class InterfaceNode(object):
    """
    Node for use by InterfaceGraph.
    Each node represents one or more DomainEnds, which, when stacked or hybridized, is collapsed into a single
    InterfaceNode representation. This is akind to having a single DomainEnd node representing multiple DomainEnds.
    Since stacking and hybridization causes InterfaceNodes to be merged, rather than connecting them by a new edge,
    all edges between InterfaceNodes are phosphate-backbone edges.
    Attributes:
    :domain_end:    The original DomainEnd instance represented by this InterfaceNode.
    :delegatee:     If graph representation of this node has been delegated to another node, this node is stored here.
    :delegated_edges: A dict holding *all* edges represented by this node, including the original set of edges,
                    with the form: {original_delegator_node1: {target_node1: eattr_dict, ...}, ...}.
                    Since delegated_edges includes this node it self, it is guaranteed that delegated_edges will
                    include self as key-element, with the original targets-dict as values. There are a number of
                    alternative datastructure implementations for this model, c.f. discussion below.

    Old discussion: How should delegate_edges be ordered?  [Selected implementation: (a)]
    For instance, consider the graph
        0------1---2-------3     0---.           .----3     0---.  1        .---3
                                      `.1---2---:                `.----2---:
        4------5---6-------7     4----´ 5   6    `---7      4----´ 5   6    `---7
    Where 5 was delegated to 1, 6 delegated to 2, and then 1 delegated to 2.
    I.e from 5 we have a tree delegation branch: 5 --> 1 --> 2
    While from 2 down we have the tree:      .- 6
                                       2 <-<ˊ
                                            `·- 1 <-- 5
    Considerations on the datastructure of delegated_edges, e.g. node 1 be w.r.t. node 5 (1-level merge):
     (a)    {5: {4: {}, 6: {}}, 1: {0: {}, 2: {}}}       - i.e. store the *original* {source: targets} edges.
     (b)    {1: {4: {}, 6: {}, 0: {}, 2: {}}}            - and then in node 1 have {5 => {4: {}, 6: {}}}
     (c)    {1: {5: {4: {}, 6: {}}, 1: {0: {}, 2: {}}}}  - i.e. nested, 1 level for each delegation.

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
    We've selected option (a), so delegated_edges [for the merged example to the right] has structure:
        Node 2 delegated_edges: {2: {1: {}, 3: {}},   # The original phosphate-backbone targets for node 2
                                 6: {5: {}, 7: {}},   # Phosphate-backbone targets delegated from node 6.
                                 5: {4: {}, 6: {}},   # Targets delegated from node 5
                                 1: {0: {}, 2: {}}}   # Delegated from node 1 (which were delegated to 5 and then 2).
        Node 5 delegated_edges: {5: {4: {}, 6: {}},   # Node 5's original backbone targets.
                                 1: {0: {}, 2: {}}}   # Targets delegated from node 1.
        Node 1 delegated_edges: {1: {0: {}, 2: {}}}   # Node 1's original backbone targets.
        Node 6 delegated_edges: {6: {5: {}, 7: {}}}   # Node 6's original backbone targets.
        Node 0 delegated_edges: {0: {1: {}}}          # Node 0's original backbone targets.
        Node 3, 4, 7 are similar to node 0.
        As you can see, the content of an InterfaceNode's delegated_edges does not change when it is being delegated
        to another node. This makes it cheaper to "reverse" a node-merge.
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


    def state_fingerprint(self):
        """ Return state-specific fingerprint hash for this node. Currently NOT cached. """
        # Changed to use always use the fingerprint of the top_delegate. This ensures that src/tgt.state_fingerprint()
        # in Complex.calculate_loop_hash() works more consistently regardless of delegation scheme.
        return hash(frozenset(ifnode.domain_end.state_fingerprint()
                              for ifnode in self.top_delegate().delegated_edges.keys()))


    def __str__(self, ):
        return "I" + self.domain_end.name

    def __repr__(self, ):
        #return "I:" + str(self.domain_end)
        return str(self) # + " at " + str(hex(id(self)))

    def __key__(self):
        """ Unfortunately python doesn't support a constant key function for comparison :( """
        return self.state_fingerprint()

    def __lt__(self, other):
        """ Make InterfaceNodes sortable. """
        # return str(self) < str(other_node)
        # return self.__key__() < other_node.__key__()
        return self.state_fingerprint() < other.state_fingerprint()

    def __gt__(self, other):
        return self.state_fingerprint() > other.state_fingerprint()


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
