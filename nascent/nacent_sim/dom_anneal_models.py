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





class Complex():
    """
    Probably not strictly needed; if you have a strand, it should be easy to
    determine if it is part of a complex.
    """
    def __init__(self, strands):
        self.Strands = strands
        for strand in strands:
            strand.Complex = self

    def add_strand(self, strand):
        """ Make strand part of this complex. """
        if strand in self.Strands:
            raise ValueError("strand %s already in this complex (%s)" % (strand, self))
        self.Strands.append(strand)
        strand.Complex = self

    def merge(self, other_complex):
        """
        Merge this complex with another.
        You should compare sizes before determining which complex is the "surviving" one.
        """
        for strand in other_complex.Strands:
            strand.Complex = self
        self.Strands += other_complex.Strands
        del other_complex
        return self




class Strand():
    """
    """
    def __init__(self, name, domains, start_complex=None):
        self.Name = name
        self.Domains = domains
        self.Complex = start_complex



class Domain():
    """
    Is a domain a "generic" encapsulement of a domain, or is it
    a specific *instance* of a domain, i.e. one domain on one strand?
    It could also be a conceptual "hybridization", i.e. a domain
    is the possible interaction between seq A and it reverse complement a.
    """
    def __init__(self, name, strand, seq=None, TM=None, partner=None):
        self.Name = name
        self.Strand = strand
        self.Sequence = seq
        self.TM = TM
        self.Partner = None

    def get_energy()


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


