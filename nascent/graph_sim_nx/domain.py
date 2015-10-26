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
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION,
                        N_AVOGADRO, ss_kuhn_length, ss_rise_per_nt, ds_rise_per_bp)
from .algorithms import connectivity_rings


# Module-level constants and variables:
make_sequential_id = sequential_number_generator()

#N_AVOGADRO = 6.022e23


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
        self.complex = None
        self.domain_strand_specie = (strand.name, name)
        self.instance_name = "%s#%s" % (self.name, self.uuid)
        self.universal_name = "%s:%s" % (self.strand.instance_name, self.instance_name)
        self.sequence = seq
        self.length = self.n_nt = len(seq)
        # E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (λˢˢ)², N = N_nt∙lˢˢ/λˢˢ
        #       = N_nt * lˢˢ * λˢˢ
        # Could be ss_mean_squared_end2end_distance or ss_msqee_dist or ss_ersq
        self.ss_length_sq = self.n_nt * ss_rise_per_nt * ss_kuhn_length
        self.ss_length_nm = math.sqrt(self.ss_length_sq)
        self.ds_length_nm = self.n_nt * ds_rise_per_bp
        self.ds_length_sq = self.ds_length_nm**2
        self.partner = partner  # duplex hybridization partner
        self.end5p = Domain5pEnd(self)
        self.end3p = Domain3pEnd(self)
        #self.stacked_upstream = None  # Edit, is self.end5p.stack_partner.domain
        #self.stacked_downstream = None
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
                # TODO: It is very important that the in-complex-identifier is unique.
                #       Fortunately, you have the option to VERIFY that it actually is.
                neighbor_rings = connectivity_rings(c, s, 5,
                                                    edge_filter=dont_follow_stacking_interactions)
                # TODO: THe hash should probably include the strand specie as well
                self._in_complex_identifier = hash(tuple(d.domain_strand_specie for d in neighbor_rings))
                specie_instances = c.domains_by_name[self.name]
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

    def __len__(self):
        return len(self.sequence)


class DomainEnd():
    """
    Attributes:
    :domain: parent domain
    :end: string indicating "5p" or "3p" end.
    :hyb_partner: Another DomainEnd, hybridized to this end.
        Typically:
         * if isinstance(self, Domain5pEnd), then isinstance(self.hyb_partner, Domain3pEnd)
         * self.domain.partner == self.hyb_partner.domain
    :pb_upstream: The DomainEnd connected on the phosphate backbone, upstream (in the 5' direction) relative to this domain.
        Typically:
        * if isinstance(self, Domain3pEnd), then self.pb_upstream.domain == self.domain
    :pb_downstream: The DomainEnd connected on the phosphate backbone, downstream
        (in the 3' direction) relative to this domain. Typically:
        * if isinstance(self, Domain5pEnd), then self.pb_downstream.domain == self.domain
    :stack_partner: A 5p3p end that that this end is stacking with.
        This is *not* the stacking between the 5p end and 3p end of a hybridized domain,
        since we can easily determine this dynamically (based on self.end and self.hyb_partner).
        Thus, stack_partner is always an end on a different domain, facing "away" from the center of this domain.
        Note that stack_partner *can* equal pb_upstream (for a 5p end) or pb_downstream (for a 3p end),
        but it doesn't have to. For instance, at a holliday junction, they will be different.
    """
    def __init__(self, domain, end):
        self.domain = domain
        self.end = end
        self.hyb_partner = None
        self.pb_upstream = None
        self.pb_downstream = None
        self.stack_partner = None

class Domain5pEnd(DomainEnd):
    def __init__(self, domain):
        super().__init__(self, domain, end="5p")

    def stacked_upstream(self, ):
        return self.stack_partner

    def stacked_downstream(self, ):
        """ If a domain is hybridized, we assume it is the domain's ends are stacked. """
        if self.hyb_partner is not None:
            return self.pb_downstream


class Domain3pEnd(DomainEnd):
    def __init__(self, domain):
        super().__init__(self, domain, end="3p")

    def stacked_upstream(self, ):
        """ If a domain is hybridized, we assume it is the domain's ends are stacked. """
        if self.hyb_partner is not None:
            return self.pb_upstream

    def stacked_downstream(self, ):
        return self.stack_partner





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
