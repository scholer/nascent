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

Module representing strands. Not currently used.

"""


def ss_upstream(domain, include_start=True):
    if include_start:
        yield domain
    while domain.backbone[0] and domain.partner is None:
        domain = domain.backbone[0]
        yield domain

def ss_downstream(domain, include_start=True):
    if include_start:
        yield domain
    while domain.backbone[1] and domain.partner is None:
        domain = domain.backbone[1]
        yield domain


class SingleStrand():

    rise = 0.4

    def __init__(self, start_domain, domains=None):
        if domains is not None:
            self.domains = domains
        else:
            self.domains = list(ss_upstream(start_domain, include_start=False))[::-1]\
                           + list(ss_downstream(start_domain))
        self.start_domain_idx = self.domains.index(start_domain)
        self.length_bp = sum(domain.length for domain in self.domains)
        self.length_nm = self.length_bp * self.rise

    def len_nm(self):
        """ Return the length of this multi-domain single-strand in nm. """
        return len(self) * self.rise

    def __len__(self):
        return sum(len(domain) for domain in self.domains)
