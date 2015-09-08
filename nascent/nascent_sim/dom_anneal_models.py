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
from collections import deque, OrderedDict
from itertools import zip_longest, chain#, accumulate


N_AVOGADRO = 6.022e23


class Complex():
    """
    Probably not strictly needed; if you have a strand, it should be easy to
    determine if it is part of a complex.
    """
    def __init__(self, strands=None, origin="o"):
        self.Strands = set(strands) if strands else set()
        for strand in self.Strands:
            strand.Complex = self
        self.Connections = set()
        self.ruid = random.randint(0, 2147483647)   # for np.random.randint, max integer is 2**31-1 (signed int)
        # Make sure to store immutables, not sets or list or anything like that:
        self.N_strand_changes = 0
        self.Strands_changes = [(self.N_strand_changes, origin,
                                 repr(self.Strands),
                                 len(self.Strands))]
        self.Strands_history = [str(sorted([str(s) for s in self.Strands]))]

    def add_strand(self, strand, origin="+add"):
        """ Make strand part of this complex. """
        if strand in self.Strands:
            raise ValueError("strand %s already in this complex (%s)" % (strand, self))
        self.Strands.add(strand)
        strand.Complex = self
        self.N_strand_changes += 1
        self.Strands_changes.append((self.N_strand_changes, origin, repr(strand), len(self.Strands)))
        self.Strands_history.append(str(sorted([str(s) for s in self.Strands])))

    def remove_strand(self, strand, origin="-rm"):
        """ Make strand not part of this complex. """
        if strand not in self.Strands:
            raise ValueError("strand %s is not in this complex (%s)" % (strand, self))
        self.Strands.remove(strand)
        strand.Complex = None
        self.N_strand_changes -= 1
        self.Strands_changes.append((self.N_strand_changes, origin, repr(strand), len(self.Strands)))
        self.Strands_history.append(str(sorted([str(s) for s in self.Strands])))

    def attach_domain(self, domain, to_domain=None, origin="+att_d"):
        """
        Attach (hybridize) <domain> to this complex.
            <to_domain> is the domain that is already part of this complex, to which the domain hybridizes.
        """
        if to_domain is None:
            to_domain = domain.Partner
        # self.Strands should be a set, so it doesn't matter if it is already there.
        # Note: At this point we have already have set domain.Partner = to_domain and vice versa.
        assert to_domain.Strand in self.Strands

        self.add_strand(domain.Strand, origin=origin)
        self.add_connection(domain, to_domain)

        # After attaching, make sure everything is right:
        try:
            assert self.Strands == set(to_domain.Strand.connected_oligos()) | {to_domain.Strand}
        except AssertionError as e:
            print("AssertionError during attach_domain:", e)
            print("self.Strands:", self.Strands)
            print("to_domain.Strand.connected_oligos():", to_domain.Strand.connected_oligos())
            print("to_domain.Strand:", to_domain.Strand)
            print("domain:", repr(domain), " - to_domain:", repr(to_domain))
            raise e


    def detach_domain(self, domain, from_domain=None, origin="-det_d"):
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
        delete_complex = False
        # How to know if the strand is still part of this complex???
        # Probably have to dig up some graph theory knowledge
        # Essentially find whether the two domains are still part of the same "connected component".
        # Edges = (domain hybridization connections + domain connected by the same strand)
        # Vertices = domains... or strands? Yeah, it might be easier to assume strands.
        # The easiest way might be to do a Breadth-first search from one domain to the other
        from_domain.Partner = None  # Make sure to break this connection first before assessing the distance...
        domain.Partner = None
        dist = from_domain.distance(domain) # Returns None if no (other) connection could be found
        if dist is None:
            #print("No connection between domains", domain, "and", from_domain,
            #      " -- breaking up complex", domain.Complex,
            #      " with strands:", domain.Complex.Strands)
            # There is nolonger any connection between the two domains.
            # We have either of the following situations:
            #       (a) Two smaller complexes
            #       (b) One complex and one unhybridized strand (oligo)
            #       (c) Two unhybridized strands
            # Uhm, wait... maybe calculate connected component sizes for the two now separate graphs?
            # There is, after all, also a probability that from_domain is now "no longer part of this complex".
            # It might be better (and more efficient) to do this at the strand level, rather than at the domain level.
            # I.e. a graph where each node is a strand, instead of a graph where each node is a domain.
            dom1_cc_oligos = domain.Strand.connected_oligos()
            dom0_cc_oligos = from_domain.Strand.connected_oligos()
            # Pay attention: Is domain.Strand itself included in connected_oligos?
            # - No. Then add 1 to lenght. (If complex consists of two strands,
            #   then len(connected_oligos()) will be 1, but the complex size should be 2.
            dom1_cc_size = len(dom1_cc_oligos) + 1  # cc = connected component
            dom0_cc_size = len(dom0_cc_oligos) + 1
            if dom0_cc_size > 1 or dom1_cc_size > 1:
                # Case (a) Two smaller complexes or case (b) one complex and one unhybridized strand
                if dom1_cc_size > dom0_cc_size:
                    # We should revert our thinking, and remove from_domain from this complex:
                    from_domain, domain = domain, from_domain
                    new_complex_oligos = dom0_cc_oligos + [domain.Strand]
                else:
                    # Remember to add the departing domain.Strand to the new_complex_oligos list:
                    new_complex_oligos = dom1_cc_oligos + [domain.Strand]
                # from_domain is now guaranteed to be on the "bigger" complex, which should continue to be this one;
                # Remove domain.Strand (which is on the new, smaller complex)
                # domain.Strand.Complex = None    # is also done by self.remove_strand()
                if dom0_cc_size > 1 and dom1_cc_size > 1:
                    #print("detach: case (a)")
                    # Remove all strands in the new complex from the old:
                    # OMG, forgetting to do this took me 4+ hours to debug!
                    # (I was just doing self.remove_strand(domain.Strand) -- grr)
                    for oligo in new_complex_oligos:
                        self.remove_strand(oligo, origin="-detach case (a)")
                    try:
                        assert domain.Strand in new_complex_oligos
                        assert domain.Strand not in self.Strands
                        assert from_domain.Strand in self.Strands
                        assert from_domain.Strand not in new_complex_oligos
                    except AssertionError as e:
                        print("---")
                        print("domain.Strand:", domain.Strand)
                        print("new_complex_oligos:", new_complex_oligos)
                        print("self.Strands:", self.Strands)
                        print("from_domain.Strand", from_domain.Strand)
                        print("---")
                    # Case (a) Two smaller complexes - must create a new complex for detached domain:
                    # Connections division needs to be calculated...
                    # Make a set of domains for the new complex:
                    new_complex_domains = {odomain for oligo in new_complex_oligos for odomain in oligo.Domains}
                    new_complex_connections = {edge for edge in self.Connections
                                               if any(dom in new_complex_domains for dom in edge)}
                    # | is the "union" operator for sets, & is intersection, ^ is symdiff, - is difference.
                    self.Connections -= new_complex_connections
                    new_complex = Complex(strands=set(new_complex_oligos),
                                          origin="-new via detach case (a)")
                    new_complex.Connections = new_complex_connections
                else:
                    # case (b) one complex and one unhybridized strand - no need to do anything further
                    #print("detach: case (b)")
                    #pass
                    self.remove_strand(domain.Strand, origin="-detach case (b)")
            else:
                # Case (c) Two unhybridized strands
                #print("detach: case (c)")
                domain.Strand.Complex = None
                from_domain.Strand.Complex = None
                #print("domain.Complex:", domain.Complex,
                #      "domain.Strand.Complex:", domain.Strand.Complex)
                #print("from_domain.Complex:", from_domain.Complex,
                #      "from_domain.Strand.Complex:", from_domain.Strand.Complex)

                assert domain.Complex == from_domain.Complex == domain.Strand.Complex \
                    == from_domain.Strand.Complex == None
                self.remove_strand(domain.Strand, origin="-detach case (c1)")
                self.remove_strand(from_domain.Strand, origin="-detach case (c2)")
                delete_complex = True
                #print("detach: case (c): Complex disintegrated - two unhybridized strands.")
        #else:
        #    print("Domains are still in the same complex.")
        return dist, new_complex, delete_complex

    def add_connection(self, domain1, domain2):
        #print("Adding connection to Complex: %s<->%s" % (domain1, domain2))
        # Using frozenset to represent un-directed connection - (a, b) is equivalent to (b, a)
        self.Connections.add(frozenset((domain1, domain2)))

    def remove_connection(self, domain1, domain2):
        #print("Removing connection from Complex: %s<->%s" % (domain1, domain2))
        self.Connections.remove(frozenset((domain1, domain2)))


    def merge(self, other_complex, new_connections=None):
        """
        Merge this complex with another.
        You should compare sizes before determining which complex is the "surviving" one.
        """
        assert other_complex != self
        #if other_complex == self:
        #    return
        for strand in other_complex.Strands:
            strand.Complex = self
        #self.Strands += other_complex.Strands
        self.Strands |= other_complex.Strands
        self.N_strand_changes -= 1
        self.Strands_changes.append(("merge", str(other_complex.Strands_changes), len(self.Strands)))
        self.Strands_history.append(str(sorted([str(s) for s in self.Strands])))
        self.Connections |= other_complex.Connections
        # | is the "union" operator for sets, & is intersection and ^ is symdiff
        if new_connections:
            for domain1, domain2 in new_connections:
                self.add_connection(domain1, domain2)
        # del other_complex   # Deleting in this namespace will have no effect.
        return self


    def FQDN(self):
        return "C[%s]" % (self.ruid % 1000)

    def __repr__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "Complex[%s] at %s" % (self.ruid % 1000, hex(id(self)))

    def __str__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "C[%s]" % (self.ruid % 100)





class Strand():
    """
    Class representing a DNA strand/oligo.
    TODO: Consider changing nomenclature from "strand" to "oligo" to conform to cadnano's notations.
    """
    def __init__(self, name, domains, start_complex=None):
        self.ruid = random.randint(0, 2147483647)   # Random unique integer
        self.Name = name
        self.Complex = start_complex
        self.Domains = domains          # List of strand domains from 5p to 3p.
        for domain in domains:
            domain.Strand = self

    def sequence(self, sep=""):
        return sep.join(d.Sequence for d in self.Domains)


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
        """ Return whether any of the strand's domains are hybridized. """
        return any(domain.Partner for domain in self.Domains)

    def distance(self, strand):
        """
        Find distance to strand using Breadth-first search.
        This is the strand-equivalent of the same method for domains.
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
        if from_domain is None:
            return self.Domains[-1]
        idx = self.Domains.index(from_domain)
        # We go forward from idx+1:  (We don't have to worry about IndexErrors when slicing)
        try:
            return self.Domains[idx+1]
        except IndexError:
            return None



    def FQDN(self):
        return "%s:%s[%s]" % (str(self.Complex), self.Name, self.ruid % 1000)


    def __repr__(self):
        # return "%s[%s]" % (self.Name, self.ruid % 100)
        return "%s:%s[%s] at %s" % (self.Complex, self.Name, self.ruid % 1000, hex(id(self)))

    def __str__(self):
        #return "%s[..%s]" % (self.Name, self.ruid % 100)
        return self.FQDN()


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
        Returns:
            a list of new complexes or None if no new complexes were created.
            a list of obsolete complexes to be deleted.
        """
        assert self.Partner is None
        assert domain.Partner is None

        self.Partner = domain
        domain.Partner = self

        if self.Strand == domain.Strand:
            # If forming an intra-strand connection, no need to make or merge any Complexes
            # Although this will screw up the complex.Connections.
            if self.Strand.Complex:
                self.Strand.Complex.add_connection(self, domain)
            return None, None

        new_complexes, obsolte_complexes = None, None
        # Update Complex:
        if self.Strand.Complex and domain.Strand.Complex:
            # Both domains are in a pre-existing complex; merge the two complexes
            if self.Strand.Complex == domain.Strand.Complex:
                #print("Intra-complex hybridization!")
                self.Strand.Complex.add_connection(self, domain)
            else:
                new_complexes, obsolte_complexes = \
                    self.Strand.Complex.merge(domain.Strand.Complex, new_connections=[(self, domain)])
            assert self.Strand.Complex == domain.Strand.Complex
        elif self.Strand.Complex:
            # self.Strand is in an existing complex; domain.Strand is not.
            self.Strand.Complex.attach_domain(domain, self, origin="+d.hyb att (a)")
        elif domain.Strand.Complex:
            # Other domain strand is in a complex; self.Strand is not.
            domain.Strand.Complex.attach_domain(self, domain, origin="+d.hyb att (b)")
        else:
            # Neither strands are in existing complex; create new complex
            new_complex = Complex(strands={self.Strand},
                                  origin="New complex from hybridizing two individual strands")
            new_complex.attach_domain(domain, self, origin="+d.hyb att (c)")     # attach domain to self within the complex
            new_complexes = [new_complex]
            assert self.Strand.Complex == domain.Strand.Complex != None
            assert self.distance(domain) != None
            assert self.Strand.Complex.Strands == {self.Strand, domain.Strand}
            assert self.Strand.connected_oligos() == [domain.Strand]
            assert domain.Strand.connected_oligos() == [self.Strand]

        return new_complexes, obsolte_complexes


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
        assert self.Partner == domain
        assert domain.Partner == self
        assert self.Strand.Complex == domain.Strand.Complex     # Can be None if hair-pin

        if self.Strand == domain.Strand:
            # If forming an intra-strand connection, no need to make or merge any Complexes
            return

        # Update Complex:
        if self.Strand.Complex:
            # Detach domain from this complex:
            dist, new_complex, delete_this = self.Strand.Complex.detach_domain(domain, self, origin="-d.dehyb det (a)")
        elif domain.Strand.Complex:
            assert not "This should never happen. \
            Two hybridized domains must be in the same Complex, possibly None if hair-pin."
            dist, new_complex, delete_this = domain.Strand.Complex.detach_domain(self, domain, origin="-d.dehyb det (b)")
        else:
            # Neither strands are in existing complex; that shouldn't happen either
            print("WARNING: Cannot dehybridize %s from %s: Neither domains are in a complex." % (self, domain))
            raise Exception("WARNING: Cannot dehybridize %s from %s: Neither domains are in a complex." % (self, domain))

        assert self.Partner is None
        assert domain.Partner is None
        return dist, new_complex, delete_this


    def effective_activity(self, other_domain, volume, oversampling=1):
        """
        Returns this domain's effective activity against another domain.
        E.g. if the two domains share the same complex, then the activity is higher because intramolecular reaction.
        """
        if self.Strand.Complex and other_domain.Strand.Complex and self.Strand.Complex == other_domain.Strand.Complex:
            # They are in the same complex
            # Currently we just use the domain-graph node distance as a measure of how close they are.
            try:
                dist = self.distance(other_domain)
                return math.sqrt(dist)
            except TypeError as e:
                print("WARNING: TypeError %s when calculating effective_activity based on distance." % e)
                print("- self = %s, other_domain = %s" % (self, other_domain))
                print("- self.Partner: %s, other_domain.Partner: %s" % (self.Partner, other_domain.Partner))
                print("- self.Strand.Complex: %s, other_domain.Strand.Complex: %s" \
                      % (self.Strand.Complex, other_domain.Strand.Complex))
                print("self.distance(other_domain): %s" % dist)
                print("self.distance(other_domain): %s" % self.distance(other_domain))
                print("other_domain.distance(self): %s" % other_domain.distance(self))
                print("self.Strand.connected_oligos():", self.Strand.connected_oligos())
                print("other_domain.Strand.connected_oligos():", other_domain.Strand.connected_oligos())
                print("self.Complex.Strands - set(self.Strand.connected_oligos()",
                      self.Complex.Strands - set(self.Strand.connected_oligos()))
                print("other_domain.Complex.Strands - set(other_domain.Strand.connected_oligos()",
                      other_domain.Complex.Strands - set(other_domain.Strand.connected_oligos()))
                if self.Complex == other_domain.Complex:
                    print("self.Complex strands:")
                    print("\n".join("--- %s: %s" % (strand, strand.Domains) for strand in self.Complex.Strands))
                    print("other_domain.Complex strands:")
                    print("\n".join("--- %s: %s" % (strand, strand.Domains) for strand in other_domain.Complex.Strands))
                    print("self.connected_domains:")
                    domains, parents, distances = self.connected_domains()
                    print("--- domains:", domains)
                    print("--- parents:", parents)
                    print("--- distances:", distances)
                    print("other_domain.connected_domains:")
                    domains, parents, distances = other_domain.connected_domains()
                    print("--- domains:", domains)
                    print("--- parents:", parents)
                    print("--- distances:", distances)
                raise e
        # If the two strands are not part of the same complex, return oversampling * concentration
        return oversampling/N_AVOGADRO/volume


    def distance(self, domain):
        """
        Find distance to domain using Breadth-first search.
        Note: It might be very useful to cache these calculations...!
        """
        # First create an iterator strategy for visiting all nodes in order of distance to self:
        distances = {}
        #parents = {}
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
                    #parents[d] = u
                    Q.append(d)

    def connected_domains(self):
        """
        List (not -Generator-) with connected domains using breadth-first search,
        starting with this domain itself as the first entry.
        Returns
            domains, parents, distances
        """
        # First create an iterator strategy for visiting all nodes in order of distance to self:
        distances = OrderedDict()   # {}
        parents = OrderedDict()     # {}
        domains = []
        Q = deque()
        distances[self] = 0     # distance to self is 0
        Q.append(self)          # We append "to the right"
        while Q:    # is not empty
            u = Q.popleft()     # We pop "from the left"
            #yield u, parents, distances        # yield here to include self as first entry
            domains.append(u)                   # returning a list rather than a generator
            # Use self.immediate_neighbors() if you want to count domains on same strand
            # but separated by N domains as N extra distance.
            adjacent_domains = u.connected_neighbors()
            for d in adjacent_domains:
                if d not in distances:
                    distances[d] = distances[u] + 1
                    parents[d] = u
                    #yield d, parents, distances         # yield here to NOT include self as first entry
                    Q.append(d)
        return domains, parents, distances

    def immediate_neighbors(self):
        """
        Return list of immediate neighbors a distance-based.
        This list only includes domains that are directly to the 5p or 3p of this domain, or this domain's Partner.
        """
        return [domain for domain in (self.domain5p(), self.domain3p(), self.Partner)
                if domain is not None]

    def connected_neighbors(self):
        """
        Return list of immediate neighbors a distance-based.
        Unlike immediate_neighbors, this method returns all domains connected to the same strand as self,
        as well as self.Partner (in third place).
        """
        if self.Partner:
            yield self.Partner
        yield from (domain for domain in chain(*zip_longest(self.Strand.domains5p(self), self.Strand.domains3p(self)))
                    if domain is not None)

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

    def FQDN(self):
        return "%s:%s[%s]" % (self.Strand.FQDN(), self.Name, self.ruid % 1000)

    def __repr__(self):
        return "%s:%s[%s] at %s" % (self.Strand.FQDN(), self.Name, self.ruid % 1000, hex(id(self)))

    def __str__(self):
        #return "%s:%s[..%s]" % (self.Strand.FQDN(), self.Name, self.ruid % 100)
        return self.FQDN()




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
