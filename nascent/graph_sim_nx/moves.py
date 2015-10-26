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

Module with functional moves for each type of event,
e.g. hybridize or dehybridize.

"""


#import random
#import math
#from collections import deque, OrderedDict
#from itertools import zip_longest, chain#, accumulate
from collections import deque, OrderedDict
from itertools import zip_longest, chain#, accumulate


from .dom_anneal_models import Complex
from .domain import distance_cached


N_AVOGADRO = 6.022e23




def hybridize(domain1, domain2):
    """
    Splitting out logic in preparation for Julia implementation.
    returns
        (new_complexes, obsolete_complexes)
    """

    assert domain1 != domain2
    assert domain1.Partner is None
    assert domain2.Partner is None

    dset = frozenset((domain1, domain2))
    #sset = frozenset(domain1.strand, domain2.strand)
    domain1.Partner = domain2
    domain2.Partner = domain1

    strand1 = domain1.strand
    strand2 = domain2.strand
    c1 = strand1.complex
    c2 = strand2.complex

    if strand1 == strand2:
        # If forming an intra-strand connection, no need to make or merge any Complexes
        if strand1.complex:
            strand1.complex.Connections.add(dset)
        return None, None

    new_complexes, obsolete_complexes = None, None

    ## Update complex:
    if c1 and c2:
        if c1 == c2:
            # Intra-complex hybridization
            assert strand1 in c2.strands and strand2 in c1.strands
            c1.Connections.add(dset)
            c_major = c1
        else:
            # Merge the two complexs:
            #obsolete_complexes = merge_complex()
            c_major, c_minor = (c1, c2) if (len(c1.strands) >= len(c2.strands)) else (c2, c1)
            for strand in c_minor.strands:
                strand.complex = c_major
            c_major.strands |= c_minor.strands
            c_major.N_strand_changes += c_minor.N_strand_changes
            c_major.Connections |= c_minor.Connections
            c_major.Connections.add(dset)
            c_minor.strands = []
            c_minor.Connections = {}
            obsolete_complexes = [c_minor]
    elif c1:
        # domain2/strand2 is not in a complex:
        c1.add_strand(strand2)
        c_major = c1
    elif c2:
        c2.add_strand(strand1)
        c_major = c2
    else:
        # Neither strands are in existing complex; create new complex
        new_complex = Complex()
        new_complexes = [new_complex]
        new_complex.strands |= {strand1, strand2}
        strand1.complex = strand2.complex = new_complex
        c_major = new_complex

    c_major.domain_distances = {} # Reset distances
    c_major.Connections.add(dset)

    assert strand1.complex == strand2.complex != None

    return new_complexes, obsolete_complexes



def dehybridize(domain1, domain2):
    """ Dehybridize domain2 from domain1. """

    assert domain1 != domain2
    assert domain1.Partner == domain2 != None
    assert domain2.Partner == domain1 != None

    strand1 = domain1.strand
    strand2 = domain2.strand
    c = strand1.complex
    c.domains_distances = {}    # Reset distances.
    assert c == strand2.complex

    dset = frozenset((domain1, domain2))
    #sset = frozenset(domain1.strand, domain2.strand)
    domain1.Partner = None
    domain2.Partner = None

    if c:
        c.Connections.remove(dset)

    if strand1 == strand2:
        # Wait, just because the domains are on the same strand
        # doesn't mean we dont have to update the complex.
        return None, None, None
    assert c is not None
    new_complex, obsolete_complexes = None, None

    ## Break the complex:
    # Update distance:
    c.domain_distances = {}
    dist = distance_cached(domain1, domain2) # We have just reset, so needs to be re-calculated
    if dist is None:
        dom1_cc_oligos = strand1.connected_oligos()
        dom2_cc_oligos = strand2.connected_oligos()
        assert strand2 not in dom1_cc_oligos
        assert strand1 not in dom2_cc_oligos
        dom1_cc_size = len(dom1_cc_oligos) + 1  # cc = connected component
        dom2_cc_size = len(dom2_cc_oligos) + 1
        print("dom1_cc_size=%s, dom2_cc_size=%s, len(c.strands)=%s" %
              (dom1_cc_size, dom2_cc_size, len(c.strands)))

        assert len(c.strands) == dom1_cc_size + dom2_cc_size -1

        if dom2_cc_size > 1 or dom1_cc_size > 1:
            # Case (a) Two smaller complexes or case (b) one complex and one unhybridized strand
            # Determine which of the complex fragments is the major and which is the minor:
            if dom1_cc_size > dom2_cc_size:
                domain_major, domain_minor = domain1, domain2

                new_complex_oligos = set(dom2_cc_oligos + [strand2])
            else:
                domain_major, domain_minor = domain2, domain2
                # Remember to add the departing domain.strand to the new_complex_oligos list:
                new_complex_oligos = set(dom1_cc_oligos + [strand1])

            if dom2_cc_size > 1 and dom1_cc_size > 1:
                # Case (a) Two smaller complexes - must create a new complex for detached domain:
                c_new = Complex()
                c.strands -= new_complex_oligos
                c_new.strands |= new_complex_oligos
                new_complex_domains = {odomain for oligo in new_complex_oligos
                                       for odomain in oligo.domains}
                new_complex_connections = {edge for edge in c.Connections
                                           if any(dom in new_complex_domains for dom in edge)}
                c.Connections -= new_complex_connections
                c_new.Connections |= new_complex_connections
                new_complexes = [c_new]
            else:
                # case (b) one complex and one unhybridized strand - no need to do anything further
                c.strands.remove(domain_minor.strand)
        else:
            # Case (c) Two unhybridized strands
            c.strands -= {strand1, strand2}
            assert c.strands == set()
            obsolete_complexes = [c]
            strand1.complex, strand2.complex = None, None

    assert domain1.Partner is None
    assert domain2.Partner is None
    return new_complex, obsolete_complexes, dist



#def merge_complexes(c_major, c_minor, new_connections=None):
#    """
#    Merge c_minor complex into c_major complex.
#    """
#
#    assert c_major != c_minor
#    #if other_complex == self:
#    #    return
#    for strand in c_minor.strands:
#        strand.complex = c_major
#    #self.strands += other_complex.strands
#    c_major.strands |= c_minor.strands
#    c_major.N_strand_changes += c_minor.N_strand_changes
#    #c_major.strands_changes.append(("merge", str(other_complex.strands_changes), len(c_major.strands)))
#    #c_major.strands_history.append(str(sorted([str(s) for s in c_major.strands])))
#    c_major.Connections |= c_minor.Connections
#    # | is the "union" operator for sets, & is intersection and ^ is symdiff
#    if new_connections:
#        for domain1, domain2 in new_connections:
#            c_major.add_connection(domain1, domain2)
#    # del other_complex   # Deleting in this namespace will have no effect.



