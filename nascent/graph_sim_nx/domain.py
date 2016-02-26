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
#import random
import math
import pdb
from collections import deque #, OrderedDict
#from itertools import zip_longest, chain#, accumulate
# import itertools
# import networkx as nx
# import numpy as np
import inspect

# Relative imports
from .utils import (sequential_number_generator, sequential_uuid_gen)
from .constants import STACKING_INTERACTION #, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION,
# from .constants import N_AVOGADRO
from .constants import ss_kuhn_length, ss_rise_per_nt, ds_rise_per_bp
from .algorithms import connectivity_rings
from .debug import printd, pprintd
from .system_graphs import InterfaceNode


# Module-level constants and variables:
make_sequential_id = sequential_number_generator()

#N_AVOGADRO = 6.022e23


dont_follow_stacking_interactions = lambda eattrs: eattrs.get('type') != STACKING_INTERACTION



class Domain(object):
    """
    A domain represents the molecular instance of a strand segment.
    For instance, the molecular instance of a strand "sA" species could
    be composed of domains da, db, dc:
               da       db       dc
        5' -------- --------- -------- 3'
    """
    def __init__(self, name, strand=None, seq=None, partner=None):
        """
        Arguments:
            :name: Domain name
            :strand: (optional, the strand that this domain is a part of)
        Note that strand is optional. The typical instantiation scheme is:
            First instantiate the domains,
            Then instantiate the strands, giving the domain-list as init argument.
        """
        # Domain unique id (duid). Can be the same as a strand's id. But two domains will never have the same id.
        self.duid = next(make_sequential_id)    # sequential number unique across domains
        self.uuid = next(sequential_uuid_gen)   # Universally unique id. Unique across all objects.
        # duid vs uuid? I prefer d/s/c uid versions. (complex#1, strand#1, domain#1).
        # uuid's are used mostly for debuggin.
        self.name = name
        self.instance_name = "%s#%s" % (self.name, self.duid)
        if strand:
            self.set_strand(strand)
        else:
            self.strand = None
            self.domain_strand_specie = (None, name)
            self.universal_name = "%s_%s" % (None, self.instance_name)
        self.sequence = seq
        ## TODO: Reach a consensus on "length" attributes
        ## Suggestions:
        ## - Always use 'ss' vs 'ds' and 'sq' vs 'nm'.
        ## - For end-to-end distances, use: 'ee' or 'rms' or 'dist' or 'ee_dist' or 'Er' or 'R' or 'rRsq' or 'rsq' or...?
        ## - For contour lengths, use: 'length', 'len', 'contour',
        self.n_nt = len(seq)
        self.length_nt = self.n_nt
        # E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (λˢˢ)², N = N_nt∙lˢˢ/λˢˢ
        #       = N_nt * lˢˢ * λˢˢ
        # Could be ss_mean_squared_end2end_distance or ss_msqee_dist or ss_ersq
        self.ss_dist_ee_sq = self.n_nt * ss_rise_per_nt * ss_kuhn_length
        self.ss_dist_ee_nm = math.sqrt(self.ss_dist_ee_sq)
        self.ds_dist_ee_nm = self.n_nt * ds_rise_per_bp
        self.ds_dist_ee_sq = self.ds_dist_ee_nm**2
        self.ss_len_contour = self.n_nt * ss_rise_per_nt
        self.ds_len_contour = self.n_nt * ds_rise_per_bp
        self._partner = partner  # duplex hybridization partner
        self.end5p = Domain5pEnd(self)
        self.end3p = Domain3pEnd(self)
        self.end5p.pb_downstream = self.end3p
        self.end3p.pb_upstream = self.end5p
        #self.stacked_upstream = None  # Edit, is self.end5p.stack_partner.domain
        #self.stacked_downstream = None

        # Cached values:
        self._in_complex_identifier = None
        # self._in_complex_identifier_symmetric = None  # replaced by symmetric_icids
        self._specie_state_fingerprint = None
        # Attributes for version 2 of in_complex_identifier labelling, for symmetric complexes:
        self.hash_ring = []
        self.icid = None   # Should be self.name or more unique after complex.label_complex_nodes().
        self.symmetric_icids = []  # Can be used to find similar domains from one asymmetric unit to another.
        self._sym_break_label = None # Should be 0 or more unique after labelling. Used to break symmetry
        self._asymmetric_unit_label = None # Should be (0,) or more unique after labelling.
        # _asymmetric_unit_label is used to distinguish pairs from different asymetric parts.
        self.history = deque(maxlen=50)


    @property
    def partner(self):
        """ Making this a property for now to ensure that ends5p3p are properly set. """
        return self._partner
    @partner.setter
    def partner(self, partner):
        """ Making this a property for now to ensure that ends5p3p are properly set. """
        # self.history.append("Setting self.partner = %r..." % partner)
        if partner is None:
            if self._partner is not None:
                # self.history.append(" - partner: Resetting current partner %r and ends..." % self._partner)
                self._partner._partner = None # pylint:disable=W0212
                self._partner.end3p.hyb_partner = None
                self._partner.end5p.hyb_partner = None
            self._partner = None
            self.end5p.hyb_partner = None
            self.end3p.hyb_partner = None
        else:
            self._partner = partner
            partner._partner = self # pylint:disable=W0212
            self.end5p.hyb_partner = partner.end3p
            partner.end3p = self.end5p.hyb_partner
            self.end3p.hyb_partner = partner.end5p
            partner.end5p = self.end3p.hyb_partner


    def set_strand(self, strand):
        """ (re-)set domain's strand. """
        # self.history.append("Setting self.strand = %r..." % strand)
        self.strand = strand
        self.domain_strand_specie = (strand.name, self.name)
        ## Concern: If you add support for strand nicking/ligation, then you cannot use strand to specify
        ## domain universal name (which must be invariant throughout the simulation).
        self.universal_name = "%s_%s" % (self.strand.instance_name, self.instance_name)


    def in_complex_identifier(self):
        """ Returns self.icid, asserting that it is not None. """
        if self.strand.complex is None:
            assert self.icid is None            # default to self.name after complex node labelling.
            assert self._sym_break_label is None     # defaults to 0 after complex node labelling.
            assert self._asymmetric_unit_label is None # defaults to (0,) after complex node labelling.
        else:
            assert self.icid is not None            # default to self.name after complex node labelling.
            assert self._sym_break_label is not None     # defaults to 0 after complex node labelling.
            assert self._asymmetric_unit_label is not None # defaults to (0,) after complex node labelling.
        return self.icid


    ## TODO: state_fingerprint and _specie_state_fingerprint should be named similarly.
    def state_fingerprint(self):
        """
        This is a hash of:
            (domain_strand_specie, complex-state, in-complex-identifier)

        Note: Use domain.symmetric_icids[0] to get domain state fingerprint without symmetry-breaking labels

        TODO: Consider not making duplexes state-dependent but rather just depend on their local stacking state.
              (Because a duplex has the same de-hybridization reaction kinetics regardless of complex state.)
        Discussion: Can we EVER use a fingerprint with non-unique ICID?
        No, I don't think so. We *always* use fingerprints in context of
        a reaction pair; Consider this:      A -- C -- B
        The "symmetric state fingerprints"   |         |
        will be the same for both As, Bs,    B -- C -- A
        and Cs. Yet, the reaction pair {A, B} is very different depending on whether we have
        {A (top left), B (top right)} or {A (top left), B (bottom left)}.
        Cases where values are the same for both symmetries:
        * Steric repulsions, ...?

        """
        #self._specie_state_fingerprint = None
        # Do NOT invoke state_change_reset() without invoking label_complex_nodes
        if self._specie_state_fingerprint is None:
            # printd("(re-)calculating state fingerprint for domain %r" % (self, ))
            dspecie = self.domain_strand_specie  # e.g. (strandA, domain1)
            # the complex's state:
            c_state = self.strand.complex.state_fingerprint() if self.strand.complex is not None else 0
            # self._specie_state_fingerprint = hash((dspecie, c_state, self.in_complex_identifier()))
            ## TODO: Clean up: Currently including is_hybridized state in fingerprint; (has been useful for debugging)
            in_complex_identifier = self.in_complex_identifier()
            self._specie_state_fingerprint = (dspecie, self.partner is not None, c_state, in_complex_identifier)
            # print("Calculated new fingerprint for domain %s: %s" % (self, self._specie_state_fingerprint))
            # self.history.append("Domain.state_fingerprint: Calculated domain specie state fingerprint = %r..." % (self._specie_state_fingerprint,))
        return self._specie_state_fingerprint


    def domain_complex_state_fingerprint(self, ):
        """ Should be inlined. """
        return self.strand.complex.state_fingerprint() if self.strand.complex is not None else 0


    def state_change_reset(self):
        """
        Reset state-dependent attributes (must be invoked after a state change).
        Note: After resetting a domain's state fingerprint, you must now invoke label_complex_nodes()

        """
        if self._specie_state_fingerprint is None:
            from inspect import currentframe, getframeinfo
            frameinfo = getframeinfo(currentframe().f_back)
            # print("Possible excessive domain.state_change_reset() called from %s:%s" %
            #       (frameinfo.filename, frameinfo.lineno))
            # print(" - self._in_complex_identifier=%s, self._specie_state_fingerprint=%s, complex._state_fingerprint=%s"
            #       % (self._in_complex_identifier, self._specie_state_fingerprint,
            #          self.strand.complex._state_fingerprint if self.strand.complex else "N/A"))
        self._in_complex_identifier = None
        # self._in_complex_identifier_symmetric = None
        self._specie_state_fingerprint = None
        # self.history.append(("Domain.state_change_reset: Unsetting specie_state_fingerprint and in_complex_identifier)...")
        self.hash_ring.clear()
        # None indicates *invalid value* for all of these:
        self.icid = None
        self._sym_break_label = None             # Used to break symmetry
        self._asymmetric_unit_label = None  # Must be set to a tuple during Complex.label_complex_nodes()
        # _asymmetric_unit_label is used to distinguish pairs from different asymetric parts.
        # DISABLED support for invoking complex reset_state_fingerprint from domain_state_reset
        # if reset_complex and self.strand.complex is not None:
        #     self.strand.complex.reset_state_fingerprint(reset_domains=False, **kwargs)


    def print_history(self, history=None, level=0, indent_str="    ", search_str=None, entrylimit=20):
        if history is None:
            history = self.history
        for (level, entry) in self.gen_history_records(history=history, level=level, search_str=search_str):
            print(indent_str*level + entry)

    def gen_history(self, history=None, level=0, indent_str="    ", search_str=None, sep="\n", entrylimit=20):
        return sep.join((indent_str*level + entry) for level, entry
            in self.gen_history_records(history=history, level=level, search_str=search_str))

    def gen_history_records(self, history=None, level=0, search_str=None, entrylimit=20):
        if history is None:
            history = self.history
        if entrylimit and entrylimit < len(history):
            org_length = len(history)
            history = history[-entrylimit:] # Makes a slice copy
            history[0] = "(...history truncated to %s of %s entries...)" % (len(history), org_length)
        for entry in history:
            if isinstance(entry, str):
                if search_str is None or search_str in entry:
                    yield (level, entry)
            else:
                nextgen = self.gen_history_records(history=entry, level=level+1, search_str=search_str)
                #yield from nextgen
                for val in nextgen:
                    yield val


    def fqdn(self):
        """ Return a "fully qualified" name, typically [complex][strand][domains]. """
        # return "%s:%s[%s]" % (self.strand.fqdn(), self.name, self.duid)
        return "%s:%s#%s" % (self.strand.fqdn(), self.name, self.duid)

    def __repr__(self):
        #frameinfo = inspect.getframeinfo(inspect.currentframe().f_back)
        #print("Domain repr called from:", frameinfo.filename, frameinfo.lineno)
        return "%s-%s#%s" % (self.strand, self.name, self.duid)
        #return self.fqdn()

    def __str__(self):
        #return self.fqdn()
        # fqdn returns in a complex; we are using string representation as nodes
        # which MUST BE INVARIANT throughout the simulation.
        return self.universal_name

    def __len__(self):
        return len(self.sequence)




class DomainEnd(object):
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
        self.name = domain.name+end
        self.instance_name = domain.instance_name+end
        self.base = domain.sequence[0 if end == "5p" else -1] if domain.sequence is not None else None
        self.hyb_partner = None     # TODO: Set this with domain.set_hyb_partner or equivalent.
        self.pb_upstream = None     # end connected by phosphate backbone on same strand
        self.pb_downstream = None   # end connected by phosphate backbone on same strand
        self.stack_partner = None   # stacking partner
        self.stack_string = None    # Stacking string, e.g. "CA/GT" or maybe frozenset("CG", "AT")
        self.ifnode = InterfaceNode(self)
        ## TODO: Many of these attributes are essentially duplicated in the ends5p3p system graph.
        ## Consider consolidating these in some way.

    @property
    def strand(self):
        """ return self.domain.strand """
        return self.domain.strand

    def state_fingerprint(self):
        """ return (self.domain.state_fingerprint(), self.end, self.stack_partner is not None) """
        return (self.domain.state_fingerprint(), self.end, self.stack_partner is not None)

    def upstream_stacking_partners(self, start_list=None):
        """
        This could also be a function and/or return a generator...
        Two recursive approaches:
            (1) Pass a start_list the the next upstream and let it append.
            (2) Start by going all the way to the top, then build+return the list as you fall back down the recursion.
        This can also be done by:
            (3) A for-loop.
            (4) A generator.
        Q: Add self to the list?
        """
        ## Method (1): Receive start_list, add self, pass list to next.
        if start_list is None:
            start_list = []
        start_list.append(self)
        upstream = self.stacked_upstream()
        if upstream is None:
            return start_list
        else:
            return upstream.upstream_stacking_partners(start_list)

    def upstream_stacking_partners_reversed(self, ):
        """ Recursive method (2): Go all the way to the top, then build list as you fall back down. """
        upstream = self.stacked_upstream()
        if upstream is None:
            return []
        upstream_list = upstream.upstream_stacking_partners2()
        upstream_list.append(self)
        return upstream_list

    def upstream_stacking_partners_list(self, start_list=None):
        """ For loop (1):  """
        upstream_list = [] if start_list is None else start_list
        end = self.stacked_upstream()
        while end is not None:
            upstream_list.append(end)
            end = end.stacked_upstream()
        return upstream_list

    def upstream_stacking_partners_generator(self):
        """ Generator (1): """
        end = self.stacked_upstream()
        while end is not None:
            yield end
            end = end.stacked_upstream()

    def downstream_stacking_partners_generator(self):
        """ Generator (1): """
        end = self.stacked_upstream()
        while end is not None:
            yield end
            end = end.stacked_downstream()


    ## IF-nodes: ##
    def upstream_stacked_top_ifnodes_generator(self):
        """ Generator (1): """
        top_ifnode = self.ifnode.top_delegate()
        end = self.stacked_upstream()
        while end is not None:
            if end.ifnode not in top_ifnode.delegated_edges: # (keys)
                top_ifnode = end.ifnode.top_delegate()
                yield top_ifnode
            end = end.stacked_upstream()

    def upstream_stacked_top_ifnodes_list(self, start_list=None):
        """ For loop (1):  """
        nodes_list = [] if start_list is None else start_list
        top_ifnode = self.ifnode.top_delegate()
        end = self.stacked_upstream()
        while end is not None:
            if end.ifnode not in top_ifnode.delegated_edges: # (keys)
                top_ifnode = end.ifnode.top_delegate()
                nodes_list.append(top_ifnode)
            end = end.stacked_upstream()
        return nodes_list

    def downstream_stacked_top_ifnodes_generator(self):
        """ Generator (1): """
        top_ifnode = self.ifnode.top_delegate()
        end = self.stacked_downstream()
        while end is not None:
            if end.ifnode not in top_ifnode.delegated_edges: # (keys)
                top_ifnode = end.ifnode.top_delegate()
                yield top_ifnode
            end = end.stacked_downstream()

    def downstream_stacked_top_ifnodes_list(self, start_list=None):
        """ For loop (1):  """
        nodes_list = [] if start_list is None else start_list
        top_ifnode = self.ifnode.top_delegate()
        end = self.stacked_upstream()
        while end is not None:
            if end.ifnode not in top_ifnode.delegated_edges: # (keys)
                top_ifnode = end.ifnode.top_delegate()
                nodes_list.append(top_ifnode)
            end = end.stacked_upstream()
        return nodes_list


    def stacked_upstream(self):
        """ Override in sub-class. """
        raise NotImplementedError()

    def stacked_downstream(self):
        """ Override in sub-class. """
        raise NotImplementedError()

    def __str__(self):
        return str(self.domain)+"_"+self.end

    def __repr__(self):
        return repr(self.domain)+"-"+self.end # + " at " id(self)


class Domain5pEnd(DomainEnd):
    """ Class representing the 5' end of a domain. """
    def __init__(self, domain):
        #super().__init__(domain, end="5p")
        DomainEnd.__init__(self, domain, end="5p")

    def stacked_upstream(self, ):
        return self.stack_partner

    def stacked_downstream(self, ):
        """
        Convention: If a domain is hybridized, we the domain's ends are stacked.
        This is usually true, unless the duplex is buckled.
        This just makes it easier to test if a helix is rigid all the way from one end to the other.
        I.e. if all 5p3p domain ends are stacked in the path between the 5p end and 3p end, then
        the helix is assumed to be a fully-stacked, rigid/semi-rigid helix all the way.
        """
        if self.hyb_partner is not None:
            return self.pb_downstream


class Domain3pEnd(DomainEnd):
    """ Class representing the 3' end of a domain. """
    def __init__(self, domain):
        #super().__init__(domain, end="3p")
        DomainEnd.__init__(self, domain, end="3p")

    def stacked_upstream(self, ):
        """ If a domain is hybridized, we assume it is the domain's ends are stacked. """
        if self.hyb_partner is not None:
            return self.pb_upstream

    def stacked_downstream(self, ):
        return self.stack_partner





def print_connection(parents, domain):
    """ """
    #print("->".join(str(d) for d in gen_parents_connection(parents, domain)))
    pass


def print_domain_distances(distances):
    """
    Since distances {{d1, d2}: dist} can have d1==d2, so just be {d1}: dist
    then it is a little tricky/tedious to print.
    """
    print("\n".join("%s<->%s: %s" % tuple([d.name for d in ds]
                                          +([] if len(ds) > 1 else [d.name for d in ds])
                                          +[dist])
                    for ds, dist in distances.items()))
