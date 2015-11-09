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
from nascent.graph_sim_nx.utils import (sequential_number_generator, sequential_uuid_gen)
from nascent.graph_sim_nx.constants import (PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION)


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


class DsHelix():
    """
    Class representing a single, fully-stacked double helix.
    NOTE: Make sure you are not replicating the functionality of system_graphs.HelixGraph
    """
    rise = 0.4  # Rise from one bp to the next in nm.
    bp_per_turn = 10.5
    twist = 1/bp_per_turn # in revolutions = 2pi radians.


    def __init__(self, start_domain):
        if start_domain.partner is None:
            raise ValueError("start_domain must be hybridized to create a DsHelix")
        partner = start_domain.partner
        # What is the order of the helices domain?
        # Is it: "both-5p-to-3p"
        # helix1: [5p-most-domain, ...., 3p-most-domain]
        # helix2: [5p-most-domain, ...., 3p-most-domain]
        # Or: "paired-domains"
        # helix1: [5p-most-domain, ...., 3p-most-domain]
        # helix2: [3p-most-domain, ...., 5p-most-domain]

        self.helices_domains = [list(domain_stack_upstream(start_domain, include_start=False))[::-1]\
                                + list(domain_stack_downstream(start_domain)),
                                list(domain_stack_upstream(partner, include_start=False))[::-1]\
                                + list(domain_stack_downstream(partner))]

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

    def __len__(self, ):
        return sum(len(domain) for domain in self.helices_domains[0])



class CashableDsHelix(DsHelix):
    """
    Class representing a single, fully-stacked double helix.
    This subclass only uses text strings rather than objects and can thus be used
    between complexes.
    """

    def __init__(self, start_domain):
        super().__init__(self, start_domain)
        self.helices_domains = [[domain.domain_strand_specie for domain in helix] for helix in self.helices_domains]
