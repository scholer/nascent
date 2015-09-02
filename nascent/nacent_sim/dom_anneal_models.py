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


import random
import math
from collections import deque
from itertools import zip_longest, chain


N_AVOGADRO = 6.022e23


class Complex():
    """
    Probably not strictly needed; if you have a strand, it should be easy to
    determine if it is part of a complex.
    """
    def __init__(self, strands=None):
        self.Strands = set(strands) if strands else set()
        for strand in self.Strands:
            strand.Complex = self
        self.Connections = set()
        self.ruid = random.randint(0, 2147483647)   # for np.random.randint, max integer is 2**31-1 (signed int)

    def add_strand(self, strand):
        """ Make strand part of this complex. """
        if strand in self.Strands:
            raise ValueError("strand %s already in this complex (%s)" % (strand, self))
        self.Strands.add(strand)
        strand.Complex = self

    def remove_strand(self, strand):
        """ Make strand not part of this complex. """
        if strand not in self.Strands:
            raise ValueError("strand %s is not in this complex (%s)" % (strand, self))
        self.Strands.remove(strand)
        strand.Complex = self


    def attach_domain(self, domain, to_domain=None):
        """
        Attach (hybridize) <domain> to this complex.
            <to_domain> is the domain that is already part of this complex, to which the domain hybridizes.
        """
        if to_domain is None:
            to_domain = domain.Partner
        # self.Strands should be a set, so it doesn't matter if it is already there.
        self.add_strand(domain.Strand)
        self.add_connection(domain, to_domain)

    def detach_domain(self, domain, from_domain=None):
        """
        Detach (dehybridize) domain from complex.

        Returns
            dist, new_complex
        Where dist is the new distance from domain to from_domain (None if the connection was broken)
        and new_complex is the new complex that was generated in case detaching domain broke up the complex.

        Note that the the domain's strand can still be part of the complex if it has another domain still attached.
        Question: How to keep track of the network connectivity?
        """
        if from_domain is None:
            from_domain = domain.Partner
        try:
            self.remove_connection(domain, from_domain)
        except KeyError as e:
            print("\n\nWARNING: KeyError %s for domain, from_domain = %s, %s\n\n" \
                  % (e, domain, from_domain))
        new_complex = None
        delete_this = False
        # How to know if the strand is still part of this complex???
        # Probably have to dig up some graph theory knowledge
        # Essentially find whether the two domains are still part of the same "connected component".
        # Edges = (domain hybridization connections + domain connected by the same strand)
        # Vertices = domains... or strands? Yeah, it might be easier to assume strands.
        # The easiest way might be to do a Breadth-first search from one domain to the other
        from_domain.Partner = None  # Make sure to break this connection first before assessing the distance...
        domain.Partner = None
        dist = from_domain.distance(domain) # Returns None if no connection could be found
        if dist is None:
            # Connection is broken. We have either of the following situations:
            #       (a) Two smaller complexes
            #       (b) One complex and one unhybridized strand (oligo)
            #       (c) Two unhybridized strands
            # Uhm, wait... maybe calculate connected component sizes for the two now separate graphs?
            # There is, after all, also a probability that from_domain is now "no longer part of this complex".
            # It might be better (and more efficient) to do this at the strand level, rather than at the domain level.
            # I.e. a graph where each node is a strand, instead of a graph where each node is a domain.
            dom1_cc_oligos = domain.Strand.connected_oligos()
            dom0_cc_oligos = from_domain.Strand.connected_oligos()
            dom1_cc_size = len(dom1_cc_oligos) + 1  # cc = connected component
            dom0_cc_size = len(dom0_cc_oligos) + 1
            if dom0_cc_size > 1 or dom1_cc_size > 1:
                # Case (a) Two smaller complexes or case (b) one complex and one unhybridized strand
                if dom1_cc_size > dom0_cc_size:
                    # We should revert our thinking, and remove from_domain from this complex:
                    from_domain, domain = domain, from_domain
                    new_complex_oligos = dom0_cc_oligos
                else:
                    new_complex_oligos = dom1_cc_oligos
                # from_domain is now guaranteed to be on the "bigger" complex, which should continue to be this one;
                # Remove domain.Strand (which is on the new, smaller complex)
                self.remove_strand(domain.Strand)
                domain.Strand.Complex = None
                if dom0_cc_size > 1 and dom1_cc_size > 1:
                    # Case (a) Two smaller complexes - must create a new complex for detached domain:
                    # Connections division needs to be calculated...
                    # Make a set of domains for the new complex:
                    new_complex_domains = {odomain for oligo in new_complex_oligos for odomain in oligo.Domains}
                    new_complex_connections = {edge for edge in self.Connections
                                               if any(dom in new_complex_domains for dom in edge)}
                    # | is the "union" operator for sets, & is intersection, ^ is symdiff, - is difference.
                    self.Connections -= new_complex_connections
                    new_complex = Complex(strands=new_complex_oligos)
                    new_complex.Connections = new_complex_connections
                # else:
                #   # case (b) one complex and one unhybridized strand - no need to do anything further
            else:
                # Case (c) Two unhybridized strands
                domain.Strand.Complex = from_domain.Strand.Complex = None
                delete_this = True

        return dist, new_complex, delete_this

    def add_connection(self, domain1, domain2):
        print("Adding connection to Complex: %s<->%s" % (domain1, domain2))
        self.Connections.add(frozenset((domain1, domain2)))

    def remove_connection(self, domain1, domain2):
        print("Removing connection from Complex: %s<->%s" % (domain1, domain2))
        self.Connections.remove(frozenset((domain1, domain2)))


    def merge(self, other_complex, new_connections=None):
        """
        Merge this complex with another.
        You should compare sizes before determining which complex is the "surviving" one.
        """
        for strand in other_complex.Strands:
            strand.Complex = self
        #self.Strands += other_complex.Strands
        self.Strands |= other_complex.Strands
        self.Connections |= other_complex.Connections
        # | is the "union" operator for sets, & is intersection and ^ is symdiff
        if new_connections:
            for domain1, domain2 in new_connections:
                self.add_connection(domain1, domain2)
        del other_complex
        return self




class Strand():
    """
    """
    def __init__(self, name, domains, start_complex=None):
        self.ruid = random.randint(0, 2147483647)   # Random unique integer
        self.Name = name
        self.Complex = start_complex
        self.Domains = domains          # List of strand domains from 5p to 3p.
        for domain in domains:
            domain.Strand = self

    def connected_oligos_gen(self, include_distance=False):
        """
        Return a generator of connected strands using breadth-first search.
        """
        distances = {}
        parents = {}
        Q = deque()
        distances[self] = 0     # distance to self is 0
        Q.append(self)          # We append "to the right"
        while Q:    # is not empty
            oligo = Q.popleft()     # We pop "from the left"
            for hyb_oligo in oligo.hybridized_oligos():
                if hyb_oligo not in distances:
                    dist = distances[hyb_oligo] = distances[oligo] + 1
                    parents[hyb_oligo] = oligo
                    Q.append(hyb_oligo)
                    if include_distance:
                        yield hyb_oligo, dist
                    else:
                        yield hyb_oligo

    def connected_oligos(self):
        """
        Return a list of connected oligos.
        """
        return list(self.connected_oligos_gen())

    def connected_oligos_count(self):
        """
        Return a list of connected oligos.
        """
        return sum(1 for oligo in self.connected_oligos())

    def hybridized_oligos(self):
        """
        Return a generator of hybridized oligos:
        """
        return (domain.Partner.Strand for domain in self.Domains if domain.Partner)

    def is_hybridized(self):
        return any(domain.Partner for domain in self.Domains)

    def distance(self, strand):
        """
        Find distance to strand using Breadth-first search.
        This is the strand-equivalent of the same method for domains.
        Note: It might be very useful to cache these calculations...!
        """
        # First create an iterator strategy for visiting all nodes in order of distance to self:
        ## TODO: Needs to be implemented!
        distances = {}
        parents = {}
        Q = deque()
        distances[self] = 0     # distance to self is 0
        Q.append(self)          # We append "to the right"
        while Q:    # is not empty
            u = Q.popleft()     # We pop "from the left"
            adjacent_domains = None
            for d in adjacent_domains:
                if d == strand:
                    # We've found the shortest (hopefully) distance to domain.
                    return distances[u]
                if d not in distances:
                    distances[d] = distances[u] + 1
                    parents[d] = u
                    Q.append(d)

    def domains5p(self, from_domain):
        """
        Return a slice of domains from from_domain to the 5p end of this strand,
        not including from_domain itself.
        """
        # Strand domains are from 5p to 3p: 5'-[domain1 - domain2 - domain3]-3'
        idx = self.Domains.index(from_domain)
        if idx == 0:
            # If we have the very first element, return an empty list.
            return []
        # We go backwards from idx-1:
        return self.Domains[idx-1::-1]

    def domains3p(self, from_domain):
        """
        Return a slice of domains from from_domain to the 5p end of this strand,
        not including from_domain itself.
        """
        if from_domain is None:
            return self.Domains[-1]
        # Strand domains are from 5p to 3p: 5'-[domain1 - domain2 - domain3]-3'
        idx = self.Domains.index(from_domain)
        # We go forward from idx+1:  (We don't have to worry about IndexErrors when slicing)
        return self.Domains[idx+1:]


    def domain5p(self, from_domain):
        """ Return the first neighbor towards the 5p end on the same strand from this domain. """
        if from_domain is None:
            return self.Domains[0]
        idx = self.Domains.index(from_domain)
        if idx == 0:
            # If we have the very first element, return an empty list.
            return None
        # We go backwards from idx-1:
        return self.Domains[idx-1]

    def domain3p(self, from_domain):
        """ Return first neighbor towards the 3p end on the same strand from this domain. """
        idx = self.Domains.index(from_domain)
        # We go forward from idx+1:  (We don't have to worry about IndexErrors when slicing)
        try:
            return self.Domains[idx+1]
        except IndexError:
            return None



    def __repr__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "%s[..%s] at %s" % (self.Name, self.ruid % 100, hex(id(self)))

    def __str__(self):
        return "%s[..%s]" % (self.Name, self.ruid % 100)


class Domain():
    """
    Is a domain a "generic" encapsulement of a domain, or is it
    a specific *instance* of a domain, i.e. one domain on one strand?
    It could also be a conceptual "hybridization", i.e. a domain
    is the possible interaction between seq A and it reverse complement a.
    """
    def __init__(self, name, strand, seq=None, partner=None):
        self.Name = name
        self.Strand = strand
        self.Sequence = seq
        self.Partner = partner
        self.ruid = random.randint(0, 2147483647)

    @property
    def Complex(self):
        return self.Strand.Complex

    def hybridize(self, domain):
        """
        Hybridize domain to this domain.
        Arguments:
            domain  the domain to hybridize to
        """
        assert self.Partner is None
        self.Partner = domain
        domain.Partner = self

        # Update Complex:
        if self.Strand.Complex:
            if domain.Strand.Complex:
                # Both domains are in a pre-existing complex; merge the two complexes
                self.Strand.Complex.merge(domain.Strand.Complex, new_connections=[(self, domain)])
            else:
                # Add domain to this complex:
                self.Strand.Complex.attach_domain(domain, self)
        elif domain.Strand.Complex:
            # Other domain strand is in a complex
            domain.Strand.Complex.attach_domain(self, domain)
        else:
            # Neither strands are in existing complex; create new complex
            new_complex = Complex(strands={self.Strand})
            new_complex.attach_domain(domain, self)     # attach domain to self within the complex


    def dehybridize(self, domain=None):
        """
        De-hybridize domain from this domain.
        Passing domain is only used to issue a warning if self.Partner != domain
        """
        if domain:
            if self.Partner != domain:
                print("WARNING: Trying to de-hybridize self (%s) from domain %s, but domain != self.Partner (%s) " \
                      % (self, domain, self.Partner))
                raise ValueError("self.Partner != domain -- %s, %s" % (self.Partner, domain))
        else:
            domain = self.Partner

        assert domain is not None

        # Update Complex:
        if self.Strand.Complex:
            # Detach domain from this complex:
            self.Strand.Complex.detach_domain(domain, self)
        elif domain.Strand.Complex:
            domain.Strand.Complex.detach_domain(self, domain)
        else:
            # Neither strands are in existing complex; that shouldn't happen
            print("WARNING: Cannot dehybridize %s from %s: Neither domains are in a complex." % (self, domain))
            raise Exception("WARNING: Cannot dehybridize %s from %s: Neither domains are in a complex." % (self, domain))

        assert self.Partner is None
        assert domain.Partner is None


    def effective_activity(self, other_domain, volume, oversampling=1):
        """
        Returns this domain's effective activity against another domain.
        E.g. if the two domains share the same complex, then the activity is higher because intramolecular reaction.
        """
        if self.Strand.Complex and other_domain.Strand.Complex and self.Strand.Complex == other_domain.Strand.Complex:
            # They are in the same complex
            # Currently we just use the domain-graph node distance as a measure of how close they are.
            return math.sqrt(self.distance(other_domain))
        else:
            return oversampling/N_AVOGADRO/volume

    def distance(self, domain):
        """
        Find distance to domain using Breadth-first search.
        Note: It might be very useful to cache these calculations...!
        """
        # First create an iterator strategy for visiting all nodes in order of distance to self:
        distances = {}
        parents = {}
        Q = deque()
        distances[self] = 0     # distance to self is 0
        Q.append(self)          # We append "to the right"
        while Q:    # is not empty
            u = Q.popleft()     # We pop "from the left"
            adjacent_domains = u.connected_neighbors()
            #print("adjacent_domains:", adjacent_domains)
            for d in adjacent_domains:
                if d == domain:
                    # We've found the shortest (hopefully) distance to domain.
                    return distances[u]
                if d not in distances:
                    distances[d] = distances[u] + 1
                    parents[d] = u
                    Q.append(d)

    def connected_components(self):
        """
        Generator with connected domains using breadth-first search.
        """
        visited = set()

    def immediate_neighbors(self):
        """
        Return list of immediate neighbors a distance-based.
        This list only includes domains that are directly to the 5p or 3p of this domain, or this domain's Partner.
        """
        from itertools import accumulate
        return [domain for domain in (self.domain5p(), self.domain3p(), self.Partner)
                if domain is not None]

    def connected_neighbors(self):
        """
        Return list of immediate neighbors a distance-based.
        Unlike immediate_neighbors, this method returns all domains connected to the same strand as self,
        as well as self.Partner (in third place).
        """
        from itertools import accumulate
        return [domain for domain in chain(*zip_longest(self.domains5p(), self.domains3p(), [self.Partner]))
                if domain is not None]

    def domains5p(self):
        """ Return a slice of all domains towards the 5p end on the same strand from this domain. """
        return self.Strand.domains5p(self)

    def domains3p(self):
        """ Return a slice of all domains towards the 3p end on the same strand from this domain. """
        return self.Strand.domains3p(self)

    def domain5p(self):
        """ Return the first neighbor towards the 5p end on the same strand from this domain. """
        return self.Strand.domain5p(self)

    def domain3p(self):
        """ Return first neighbor towards the 3p end on the same strand from this domain. """
        return self.Strand.domain3p(self)

    def domain5p_is_hybridized(self):
        neighbor = self.domain5p()
        return bool(neighbor and neighbor.Partner)

    def domain3p_is_hybridized(self):
        neighbor = self.domain3p()
        return bool(neighbor and neighbor.Partner)

    def partner(self):
        return self.Partner

    def upper(self):
        return self.Name.upper()

    def lower(self):
        return self.Name.lower()

    def __repr__(self):
        return "%s:%s[..%s] at %s" % (str(self.Strand), self.Name, self.ruid % 100, hex(id(self)))

    def __str__(self):
        return "%s[..%s]" % (self.Name, self.ruid % 100)




class DomainBasket():
    """
    A basket with all domains.
    Could probably be a simple dict data structure, no need for a dedicated class:
        db[<name>] = <list of domains ... or maybe list of strands?>
    """



class Tube():
    """
    """

    def __init__(self, volume, strands=None):
        self.Volume = volume
        self.Strands = strands
