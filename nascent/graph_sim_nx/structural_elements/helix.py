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


"""

Module representing helices.

"""
# sequential_number_generator returns a generator; sequential_uuid_gen *is* a generator.
from nascent.graph_sim_nx.utils import (sequential_number_generator, sequential_uuid_gen)
from nascent.graph_sim_nx.constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)

# A generator, get next sequential id with next(make_helix_sequential_id)
make_helix_sequential_id = sequential_number_generator()


def domain_stack_upstream(domain, include_start=True):
    if include_start:
        yield domain
    while domain.stacked_upstream:
        domain = domain.stacked_upstream
        yield domain

def domain_stack_downstream(domain, include_start=True):
    if include_start:
        yield domain
    while domain.stacked_downstream:
        domain = domain.stacked_downstream
        yield domain

def ends5p3p_stack_upstream(end, include_start=True):
    if include_start:
        yield end
    end_upstream = end.stacked_upstream()
    while end_upstream:
        yield end
        end_upstream = end.stacked_upstream()

def ends5p3p_stack_downstream(end, include_start=True):
    if include_start:
        yield end
    end_downstream = end.stacked_downstream()
    while end_downstream:
        yield end
        end_downstream = end.stacked_downstream()

def ends5p3p_stack_from_domain(domain, direction, include_start=True):
    if direction == "upstream":
        if include_start:
            return ends5p3p_stack_upstream(domain.end3p, include_start=True)
        else:
            return ends5p3p_stack_upstream(domain.end5p, include_start=False)
    elif direction == "downstream":
        if include_start:
            return ends5p3p_stack_downstream(domain.end5p, include_start=True)
        else:
            return ends5p3p_stack_downstream(domain.end3p, include_start=False)
    raise ValueError("Unrecognized direction '%s'" % direction, direction)


class DoubleHelix():
    """
    Class representing a single, fully-stacked double helix.
    NOTE: Make sure you are not replicating the functionality of system_graphs.HelixGraph

    Also: A new double helix is formed every time a domain hybridizes. Double-helices can expand/merge
    by stacking, but only by merging with an existing double-helix. You probably do not want to instantiate a
    new object for every hybridization reaction.
    If you need to refer to helices uniquely, you can use a module-level sequential_number_generator.
    """
    rise = 0.4  # Rise from one bp to the next in nm.
    bp_per_turn = 10.5
    twist = 1/bp_per_turn # in revolutions = 2pi radians.


    def __init__(self, start_domain):
        if start_domain.partner is None:
            raise ValueError("start_domain must be hybridized to create a DsHelix")
        self.huid = next(make_helix_sequential_id)
        self.direction_vector = [1.0, 0.0, 0.0]  # [x, y, z], could/should be a quarternion.
        partner = start_domain.partner
        # What is the order of the helices domain?
        # Is it: "both-5p-to-3p"
        # helix1: [5p-most-domain, ...., 3p-most-domain]
        # helix2: [5p-most-domain, ...., 3p-most-domain]
        # Or: "paired-domains"
        # helix1: [5p-most-domain, ...., 3p-most-domain]
        # helix2: [3p-most-domain, ...., 5p-most-domain]

        # A double-helix consists of two helices (helix "strands")
        # - We use helix "strand" although they it does not need to be a single strand *molecule*.
        # - Each helix strand in the double-helix is merely composed of domains.
        self.helices_domains = [list(domain_stack_upstream(start_domain, include_start=False))[::-1]\
                                + list(domain_stack_downstream(start_domain)),
                                list(domain_stack_upstream(partner, include_start=False))[::-1]\
                                + list(domain_stack_downstream(partner))]

        self.helix_strand_sets = [set(domains) for domains in self.helices_domains]
        self.helix_strand_parity_by_domain = {}
        for helix_parity, helix_strand in enumerate(self.helices_domains):
            # helix parity is either 0 or 1 and correspond to the index in helices_domains.
            # This also functions to hold a set of all domains in the ds-helix.
            for domain in helix_strand:
                self.helix_strand_parity_by_domain[domain] = helix_parity
        # domain ends (5p3p)
        # Version 1, assuming "both-5p-to-3p" helix domain ordering:
        self.helices_5p3p_ends = [[end for domain in helix for end in [domain.end5p, domain.end3p]]
                                  for helix in self.helices_domains]

        # Version 2, assuming "both-5p-to-3p" helix domain ordering:
        self.helices_5p3p_ends = [list(ends5p3p_stack_upstream(start_domain.end5p))[::-1]\
                                  + list(ends5p3p_stack_downstream(start_domain.end3p)),
                                  list(ends5p3p_stack_upstream(partner.end5p))[::-1]\
                                  + list(ends5p3p_stack_downstream(partner.end3p))]

        # Note: This zipping assumes "paired-domains" helix domain ordering:
        self.domain_pairs = list(zip(self.helices_domains))
        self.domains = [domain for pair in self.domain_pairs for domain in pair]
        self.lengths = [sum(domain.length for domain in helix) for helix in self.helices_domains]
        if self.lengths[0] != self.lengths[1]:
            print("WARNING: self.lengths[0] != self.lengths[1]:", self.lengths)
        self.length_bp = self.lengths[0]
        self.length_nm = self.length_bp*0.34


    def split(self,):
        """ Split the helix into two, e.g. after an unstacking reaction. """
        pass

    def __lt__(self, other_helix):
        return self.huid < other_helix.huid
    def __gt__(self, other_helix):
        return self.huid > other_helix.huid

    def __len__(self, ):
        return sum(len(domain) for domain in self.helices_domains[0])



class CashableDsHelix(DoubleHelix):
    """
    Class representing a single, fully-stacked double helix.
    This subclass only uses text strings rather than objects and can thus be used
    between complexes.
    Not sure if this is really worth the trouble...
    """

    def __init__(self, start_domain):
        super(CashableDsHelix, self).__init__(self, start_domain)
        self.helices_domains = [[domain.domain_strand_specie for domain in helix] for helix in self.helices_domains]
