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

Module for

"""

import itertools
import networkx as nx
import numpy as np

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        interactions)
from .strandcolor import StrandColor


# Module-level constants and variables:
make_sequential_id = sequential_number_generator()
colormgr = StrandColor()
#colormgr = StrandColorRandom()



class Strand(nx.MultiGraph):
    """
    Graph representing a DNA strand/oligo.
    Should it be directed?
    - No:  We have actual connection in both directions.
    - Yes: Makes it easier to draw/visualize.
    """
    def __init__(self, name, domains, start_complex=None):
        # all **kwargs given to nx.Graph is used as attributes.
        # a graph's attributes are simply a dict stored in Graph.graph member
        self.uuid = next(sequential_uuid_gen)   # sequential number unique across all objects
        self.suid = next(make_sequential_id)    # sequential number unique across strands
        super().__init__(name=name,
                         complex=start_complex,
                         suid=self.suid,
                        )
        #self.name = name  # Is a native Graph property, linked to self.graph
        self.instance_name = "%s#%s" % (self.name, self.suid)
        self.ends5p3p_graph = nx.MultiGraph()
        self.hue = colormgr.strand_color(self)
        self.color_rgb = colormgr.domain_rgb(self, frac=0.5)
        # suid = Strand unique id. Can be the same as a domain's id. But two strands will never have the same id.
        if domains:
            self.set_domains(domains)

    ## TODO: Add strand-unique node/edge colors (hue)
    ## All domain in the same strand have the same hue.
    ## The domains then change their saturation/value from the 5p to the 3p end of the strand.
    ## Edit: Hue depends on strand *specie*.

    def set_domains(self, domains):
        """ Set this strands domains. Should only be invoked once. """
        self.domains = domains
        self.length = sum(d.length for d in domains)

        # self is a graph of domains:
        # self.add_nodes_from(domains) # Add all domains as nodes
        n_domains = len(domains)
        # for idx, domain in domains:
        # connect nodes from the 5p end to the 3p end:
        # edge colors: backbone=blue, hyb=red, stack=green/cyan,
        edge_attrs = {'interaction': PHOSPHATEBACKBONE_INTERACTION,
                      'R': 0, 'G': 0, 'B': 1,  # pygephi colors
                      #'type': PHOSPHATEBACKBONE_INTERACTION,
                      #'direction': None
                     }
        #edges = zip(iter(domains), iter(domains[1:]), [edge_attrs]))
        #self.add_edges_from(edges)
        from itertools import zip_longest
        for idx, (domain, domain2) in enumerate(zip_longest(self.domains, self.domains[1:])):
            ## Add domain nodes:
            frac = (idx+1)/(n_domains+1)
            hsv = colormgr.domain_hsv(self, frac)
            hsv_str = "%0.03f %0.03f %0.03f" % hsv
            # Edit: The above works for graphviz, but we generally don't plot with graphviz...
            # networkx.draw* uses matplotlib, which takes rgb tuples:
            # TODO: This should probably be visualizer-specific properties
            rgb_tup = colormgr.domain_rgb(self, frac)
            node_attrs = {'color_rgb': rgb_tup,
                          'color_hsv': hsv_str,
                          'size': len(domain), # pygephi: 200 is a reasonable size.
                          'Label': domain.name,
                          'R': rgb_tup[0],
                          'G': rgb_tup[1],
                          'B': rgb_tup[2],
                         }
            self.add_node(domain, attr_dict=node_attrs)

            ## Add edge between domain and domain2:
            if domain2 is not None:
                #s_edge_key = (frozenset((domain1.universal_name, domain2.universal_name)), PHOSPHATEBACKBONE_INTERACTION)
                s_edge_key = PHOSPHATEBACKBONE_INTERACTION
                avg_length = sum((len(domain), len(domain2)))/2
                # Note: If strand is not a multigraph, then key will be an edge attribute, NOT a proper key
                # That means when you do multigraph.add_edges_from(strand.edges(data=True)), then
                # the multigraph's edges will have an attribute key=akey,
                # BUT THEY WILL NOT ACTUALLY HAVE ANY KEYS!
                self.add_edge(domain, domain2, key=s_edge_key,
                              len=avg_length, weight=10/avg_length, **edge_attrs)

            ## Update ends5p3p_graph:
            self.ends5p3p_graph.add_node(domain.end5p, attr_dict=node_attrs)
            self.ends5p3p_graph.add_node(domain.end3p, attr_dict=node_attrs)

            # The "long" edge from 5p of d to 3p of domain:
            long_edge_attrs = dict(edge_attrs, weight=10/len(domain), len=len(domain),
                                   _color=rgb_tup, style='bold')
            self.ends5p3p_graph.add_edge(domain.end5p, domain.end3p, key=PHOSPHATEBACKBONE_INTERACTION,
                                         attr_dict=long_edge_attrs)
            # The "short" edge connecting 3p of one domain to 5p of the next domain:
            if domain2 is not None:
                short_edge_attrs = dict(edge_attrs,
                                        weight=10,  # Typically used as spring constant. Higher -> shorter edge.
                                        len=4       # With of ds helix ~ 2 nm ~ 5 bp
                                       )
                self.ends5p3p_graph.add_edge(domain.end3p, domain2.end5p, key=PHOSPHATEBACKBONE_INTERACTION,
                                             attr_dict=short_edge_attrs)
                # Update domain 5p3p ends:
                # (Yes, this duplicates the info in strand.ends5p3p_graph and system 5p3p graph, which is bad.)
                domain.end3p.pb_downstream = domain2.end5p
                domain2.end5p.pb_downstream = domain.end3p


        # ends_domain_edges = [(d.end5p, d.end3p, PHOSPHATEBACKBONE_INTERACTION,
        #                       dict(edge_attrs, weight=10/len(d), len=len(d)))
        #                      for d in domains]
        # ## Note: For spring-based layouts, weight is sometimes the spring constants.
        # ## Thus, larger weight -> shorter edge length.
        # edge_attrs['weight'] = 10   # Not the edge betwen 5p-3p in one domain, but 3p of one domain to 5p of the next.
        # # Length from the 3p of domain1 to 5p of domain2:
        # edge_attrs['len'] = 4
        # ends_domain_interfaces = zip([d.end5p for d in domains], [d.end3p for d in domains[1:]],
        #                         itertools.cycle([PHOSPHATEBACKBONE_INTERACTION]), itertools.cycle([edge_attrs]))
        # self.ends5p3p_graph.add_edges_from(itertools.chain(ends_domain_edges, ends_domain_interfaces))

        # Make sure all domains have this strand registered as their parent:
        for domain in domains:
            domain.set_strand(self)

        ## Adjacency and distance matrices ##
        # could be {frozenset((d1, d2)): dist} or a matrix.
        # self.strand_domain_distances = {
        #     frozenset((d1, d2)): 1
        #     for d1, d2 in zip(iter(domains), iter(domains[1:]))}
        # n = len(domains)
        # self.adjacency_matrix = np.ndarray((n, n), dtype=bool)
        # self.distance_matrix = np.ndarray((n, n), dtype=int)
        # # nt_distance_matrix: only domains between d1 and d2 (but not including either)
        # self.nt_distance_matrix = np.ndarray((n, n), dtype=int)
        #
        # for i, j in itertools.combinations_with_replacement(range(n), 2):
        #     # j >= i for combinations_with_replacement:
        #     self.distance_matrix[i, j] = self.distance_matrix[j, i] = np.absolute(i - j)
        #     self.adjacency_matrix[i, j] = self.adjacency_matrix[j, i] = np.absolute(i - j) == 1
        #     self.nt_distance_matrix[i, j] = self.nt_distance_matrix[j, i] = \
        #         sum(dom.length for dom in domains[i+1:j])


    @property
    def complex(self):
        return self.graph['complex']

    @complex.setter
    def complex(self, cmplx):
        self.graph['complex'] = cmplx

    # 'name' is already a built-in property of nx.Graph objects
    # @property
    # def name(self):
    #     return self.graph['name']
    #
    # @name.setter
    # def name(self, name):
    #     self.graph['name'] = name

    def sequence(self, sep=""):
        """ Return the full strand sequence as the sequence of each domain, separated by <sep>. """
        return sep.join(d.sequence for d in self.domains)

    def is_hybridized(self):
        """ Return whether any domain in strand is hybridized. """
        return any(domain.partner is not None for domain in self.domains)

    def fqdn(self):
        """ Return Complex:Strand[Domain] """
        return "%s:%s[%s]" % (self.complex, self.name, self.suid)

    def __repr__(self):
        return self.fqdn()

    def __str__(self):
        #return self.fqdn()
        # string representation must be INVARIANT during the course of the simulation!
        return self.instance_name
