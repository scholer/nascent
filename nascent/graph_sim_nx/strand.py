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

# Module-level constants and variables:
make_sequential_id = sequential_number_generator()





class Strand(nx.Graph):
    """
    Graph representing a DNA strand/oligo.
    Should it be directed?
    - No:  We have actual connection in both directions.
    - Yes: Makes it easier to draw/visualize.
    """
    def __init__(self, name, domains, start_complex=None):
        # all **kwargs given to nx.Graph is used as attributes.
        # a graph's attributes are simply a dict stored in Graph.graph member
        super().__init__(name=name,
                         complex=start_complex,
                         suid=next(make_sequential_id),
                        )
        #self.name = name  # Is a native Graph property, linked to self.graph
        self.uuid = next(sequential_uuid_gen)   # sequential number unique across all objects
        self.suid = next(make_sequential_id)    # sequential number unique across strands
        self.instance_name = "%s#%s" % (self.name, self.uuid)
        self.ends5p3p_graph = nx.Graph()
        # suid = Strand unique id. Can be the same as a domain's id. But two strands will never have the same id.
        if domains:
            self.set_domains(domains)


    def set_domains(self, domains):
        """ Set this strands domains. Should only be invoked once. """
        self.domains = domains
        self.length = sum(d.length for d in domains)

        # self is a graph of domains:
        self.add_nodes_from(domains) # Add all domains as nodes
        # connect nodes from the 5p end to the 3p end:
        edge_attrs = {'weight': 1,  # Not the edge betwen 5p-3p in one domain, but 3p of one domain to 5p of the next.
                      'type': PHOSPHATEBACKBONE_INTERACTION,
                      'direction': None}
        edges = zip(iter(domains), iter(domains[1:]), itertools.cycle([edge_attrs]))
        self.add_edges_from(edges)

        ## Update ends5p3p_graph:
        domain_edges = [(d.end5p, d.end3p, dict(edge_attrs, weight=len(d), len=len(d)))
                        for d in domains]
        domain_interfaces = zip([d.end5p for d in domains], [d.end5p for d in domains[1:]],
                                itertools.cycle([edge_attrs]))
        self.ends5p3p_graph.add_edges_from(itertools.chain(domain_edges, domain_interfaces))

        # Make sure all domains have this strand registered as their parent:
        for domain in domains:
            domain.set_strand(self)

        ## Adjacency and distance matrices ##
        # could be {frozenset((d1, d2)): dist} or a matrix.
        self.strand_domain_distances = {
            frozenset((d1, d2)): 1
            for d1, d2 in zip(iter(domains), iter(domains[1:]))}
        n = len(domains)
        self.adjacency_matrix = np.ndarray((n, n), dtype=bool)
        self.distance_matrix = np.ndarray((n, n), dtype=int)
        # nt_distance_matrix: only domains between d1 and d2 (but not including either)
        self.nt_distance_matrix = np.ndarray((n, n), dtype=int)

        for i, j in itertools.combinations_with_replacement(range(n), 2):
            # j >= i for combinations_with_replacement:
            self.distance_matrix[i, j] = self.distance_matrix[j, i] = np.absolute(i - j)
            self.adjacency_matrix[i, j] = self.adjacency_matrix[j, i] = np.absolute(i - j) == 1
            self.nt_distance_matrix[i, j] = self.nt_distance_matrix[j, i] = \
                sum(dom.length for dom in domains[i+1:j])


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
        return all(domain.partner is not None for domain in self.domains)

    def fqdn(self):
        """ Return Complex:Strand[Domain] """
        return "%s:%s[%s]" % (str(self.graph['complex']), self.graph['name'], self.graph['suid'])

    def __repr__(self):
        return self.fqdn()

    def __str__(self):
        #return self.fqdn()
        # string representation must be INVARIANT during the course of the simulation!
        return self.instance_name
