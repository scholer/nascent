#!/usr/bin/env python
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

# pylint: disable=W0142,C0103,C0301,W0141

"""


TODO: Change all mentions of "strand" to "oligo" to conform with cadnano's nomenclature.

"""


#import random
import math
#from collections import deque, OrderedDict
#from itertools import zip_longest, chain#, accumulate
import itertools
import networkx as nx
import numpy as np

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen,
                    PHOSPHATE_BACKBONE, DUPLEX_HYBRIDIZATION, STACKING_INTERACTION)
from .algorithms import connectivity_rings

# Module-level constants and variables:
make_sequential_id = sequential_number_generator()

N_AVOGADRO = 6.022e23


dont_follow_stacking_interactions = lambda eattrs: eattrs.get('type') != STACKING_INTERACTION



class Domain():
    """
    A domain represents the molecular instance of a strand segment.
    For instance, the molecular instance of a strand "sA" species could
    be composed of domains da, db, dc:
               da       db       dc
        5' -------- --------- -------- 3'
    """
    def __init__(self, name, strand, seq=None, partner=None):
        self.name = name
        self.strand = strand
        self.domain_strand_specie = (strand.name, name)
        self.instance_name = "%s#%s" % (self.name, self.uuid)
        self.universal_name = "%s:%s" % (self.strand.instance_name, self.instance_name)
        self.sequence = seq
        self.partner = partner
        # Domain unique id. Can be the same as a strand's id. But two domains will never have the same id.
        self.duid = next(make_sequential_id)
        self.uuid = next(sequential_uuid_gen)

        # Cached values:
        self._in_complex_identifier = None
        self._specie_state_fingerprint = None



    def in_complex_identifier(self, strandgraph=None):
        """
        Return a hash that can be used to identify the current domain within the complex.
        """
        if self._in_complex_identifier is None:
            s = domain.strand
            c = s.complex
            if not c or len(c.strands_by_name[s.name]) < 2:
                # If only one strand in the complex, no need to calculate strand-domain identifier:
                self._in_complex_identifier = 0
            else:
                # We have a complex with more than one copy of the domain's strand specie:
                # probably just probing the local environment is enough.
                # It won't tell *exactly* which domain we have, but it will usually be good enough.
                # Note: Using hybridization graph, because no reason to vary based on stacking here
                neighbor_rings = connectivity_rings(c, s, 5,
                                                    edge_filter=dont_follow_stacking_interactions)
                self._in_complex_identifier = hash(tuple(d.domain_strand_specie for d in neighbor_rings))
        return self._in_complex_identifier


    def domain_state_fingerprint(domain, strandgraph=None):
        """
        This is a hash of:
            (domain_strand_specie, complex-state, in-complex-identifier)
        """
        if not domain._specie_state_fingerprint:
            c = s.complex
            dspecie = domain.domain_strand_specie,  # e.g. (strandA, domain1)
            # the complex's state:
            c_state = c.state_fingerprint() if c else 0,
            domain._specie_state_fingerprint = hash((dspecie, c_state, self.in_complex_identifier()))
        return domain._specie_state_fingerprint


    def state_change_reset(self):
        self._in_complex_identifier = None
        self._specie_state_fingerprint = None


    def fqdn(self):
        return "%s:%s[%s]" % (self.strand.fqdn(), self.name, self.did)

    def __repr__(self):
        return self.fqdn()

    def __str__(self):
        return self.fqdn()





def print_connection(parents, domain):
    """ """
    print("->".join(str(d) for d in gen_parents_connection(parents, domain)))


def print_domain_distances(distances):
    """
    Since distances {{d1, d2}: dist} can have d1==d2, so just be {d1}: dist
    then it is a little tricky/tedious to print.
    """
    print("\n".join("%s<->%s: %s" % tuple([d.name for d in ds]
                                          +([] if len(ds) > 1 else [d.name for d in ds])
                                          +[dist])
                    for ds, dist in distances.items()))
