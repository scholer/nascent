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

# pylint: disable=C0103


"""

Module for representing system-level graphs (where each complex form it's own subgraph)
 - domains
 - strands


"""

import networkx as nx


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
                self.add_edge(d.model5p, d.model3p, interaction='pb', weight=d.length, type="5p3p")
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
                self.add_edge(d.model5p, d.model3p, interaction='pb', weight=d.length, type="5p3p")
            if len(domains) > 1:
                for upstream, downstream in zip(domains[:-1], domains[1:]): # source, target
                    self.add_edge(upstream.model3p, downstream.model5p, interaction='pb', weight=1)

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
