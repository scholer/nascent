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

# pylint: disable=W0142,C0103,C0301,W0141

"""


TODO: Change all mentions of "strand" to "oligo" to conform with cadnano's nomenclature.

"""


from __future__ import absolute_import, print_function, division
import random
import math
from collections import deque, OrderedDict
try:
    from itertools import zip_longest, chain
except ImportError:
    from itertools import izip_longest as zip_longest, chain


N_AVOGADRO = 6.022e23


def distance_cached(domain1, domain2, return_parents=False):
    """
    """
    distances = domain1.strand.complex.domain_distances
    parents = {}
    ftset = frozenset((domain1, domain2)) # "from-to" set
    if ftset in distances:
        #return distances[ftset]
        return (distances[ftset], parents) if return_parents else distances[ftset]

    # Complex.domain_distances should incorporate Strand.domain_distances
    # when strand is added,
    if domain2 == domain1:
        print("domain2 == domain1 - shouldn't happen, but OK.")
        distances[ftset] = 0
        return (0, parents) if return_parents else 0

    print("Determining distance from %s to %s in complex %s" %
          (domain1.name, domain2.name, domain1.strand.complex))
    Q = deque()
    dset = frozenset((domain1, domain1)) # = {domain1}, can only have unique elements.
    #assert dset in distances and distances[dset] == 0
    distances.setdefault(dset, 0)     # distance to domain1 is 0
    Q.append(domain1)          # We append "to the right"
    ## TODO: Search could be optimized if you start by using a "strand node graph",
    ## or, equivalently, if you use the domain-domain dist matrix where you search for zeros.
    while Q:    # is not empty
        u = Q.popleft()     # We pop "from the left"
        #print("u node = %s" % u.name)
        uset = frozenset((domain1, u))
        fset = frozenset((domain2, u))
        if fset in distances:
            distances[ftset] = distances[uset] + distances[fset]

        if u.partner is not None:
            #print("Checking u's partner, %s" % u.partner)
            parents[u.partner] = u
            dset = frozenset((domain1, u.partner))
            if dset not in distances:
                #print(" - Adding %s = %s to distances (same as %s)." % (dset, distances[uset], uset))
                distances[dset] = distances[uset]
                Q.append(u.partner)
            if u.partner == domain2:
                return (distances[dset], parents) if return_parents else distances[dset]
        for d in u.domain5p3p():
            # TODO: This is not the most efficient way to iterate; if we are searching in the
            # 3p direction from a domain, domain5p3p() will also return the domain at the 5p.
            # However, the old didn't give the right behaviour for strands with more than 4 domains.
            if d:   # could be None
                parents[d] = u
                #print("At child domain %s" % d.name)
                dset = frozenset((domain1, d))
                if dset not in distances:
                    distances[dset] = distances[uset] + 1
                    #print(" - Adding dist for dset %s = %s (parent = %s)." % (dset, distances[dset], uset))
                    Q.append(d)
                elif distances[uset] + 1 < distances[dset]:
                    print("Found a faster route for dset: %s (%s < %s" %
                          (dset, distances[uset] + 1, distances[dset]))
                    distances[dset] = distances[uset] + 1
                    Q.append(d)
                if d == domain2:
                    #print("Found target domain %s" % domain2.name)
                    return (distances[dset], parents) if return_parents else distances[dset]

def gen_parents_connection(parents, domain):
    """
    parents is a dict with {parent: child} entries
    from a root node outwards.
    this will generate a sequence of nodes from leaf <domain> to root.
    """
    yield domain
    while domain in parents:
        domain = parents[domain]
        yield domain

def check_acyclic_graph(tree, leaf):
    """
    tree is a dict specified by connected nodes as
        {parent: child}
    """
    visited = []
    is_cyclic = False
    # generates leaf it self as the first node
    for p in gen_parents_connection(tree, leaf):
        if p in visited:
            print("Found cyclic graph:")
            print(" - At p = %s" % p)
            print(" - visited = %s" % visited)
            is_cyclic = True
            break
        visited.append(p)
    return is_cyclic

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
