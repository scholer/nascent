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

Utility functions for dom_simulation.

"""

from collections import defaultdict


from .dom_anneal_models import Strand, Domain



def parse_strand_domains_file(filepath, sep1="\t", sep2=",", n_clones_default=1):
    """
    Parse strand domain definition file and return as list of strands.

    Typical file structure is: (case: 4-way /holliday junction)
        Strand	Domains 	Sequence
        strandA	H1A, H1B	...
        strandB	H2B, H2A	...
        strandC	h2a, h1a	...
        strandD	h1b, h2b	...

    Nomenclature is that h1a is complementary to H1A. (Upper vs lower case)

    If a line begins with "#", it is assumed to be a comment and is ignored.
    """
    with open(filepath) as fp:
        _ = next(fp).split()    # header
        lines = [[cell.strip() for cell in line.strip().split(sep1)] for line in fp if line[0] != "#"]
    strands = []
    for i, line in enumerate(lines, 2):
        strand_name = line[0]
        try:
            n_clones = line[3]
        except IndexError:
            n_clones = n_clones_default
        for _ in range(n_clones):
            dom_names = [domain_name.strip() for domain_name in line[1].split(sep2)]
            dom_seqs = [domain_seq.strip() for domain_seq in line[2].split(sep2)]
            if len(dom_names) != len(dom_seqs):
                print("Warning: dom_names (%s) != dom_seqs (%s) for line %s" % \
                      (len(dom_names), len(dom_seqs), i))
            domains = [Domain(dom_name, strand=None, seq=dom_seq)
                       for dom_name, dom_seq in zip(dom_names, dom_seqs)]
            strand = Strand(strand_name, domains=domains)
            strands.append(strand)
    return strands


def check_strands(strands, unpaired_expected=None):
    """
    Check strands for common errors/mistakes, e.g. unpaired domains.
    """
    domains = {domain for strand in strands for domain in strand.Domains}
    domains_by_name = defaultdict(list)
    for d in domains:
        domains_by_name[d.Name].append(d)
    unpaired_domains = [d for d in domains \
                        if not (d.Name.lower() if d.Name.upper() == d.Name else d.Name.upper())
                        in domains_by_name.keys()]
    print("Unpaired domains: (%s of %s = %s %%)\n" % \
          (len(unpaired_domains), len(domains), int(100*len(unpaired_domains)/len(domains))))
    print(", ".join(str(d) for d in unpaired_domains))
    unexpected_unpaired = None
    if unpaired_expected:
        unexpected_unpaired = set(d.Name for d in unpaired_domains) - set(unpaired_expected)
        if unexpected_unpaired:
            print("Unexpected unpaired domains:", ", ".join(name for name in unexpected_unpaired))
    return unpaired_domains, unexpected_unpaired
