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

# pylint: disable=C0103,W0142,W0212,R0902

"""

Module for managing the whole system.

I think there might be a lot of code involved in managing the system (graphs and structure),
and having all that in the simulator makes it a bit hard to read.

Splitting the "system state" code out to a separate module allows the simulator to be focused
on just the "stochastic simulation" part and be mostly agnostic on the system details.

Question: Who manages "reaction propensity" etc?
 - I'm thinking the system manager takes care of everything except the simulation steps / temperature control.


The "management system" is composed of three classes:
* Graph Manager - Takes care of system-level graphs
* Component Manager - Takes care of "components", i.e. domains, strands, complexes - as well as hybridize/dehybridize.
* Reaction Manager - Everything related to reactions, i.e. energies, c_j, propensity functions, etc.



## This manager will take care of: ##
System state:
 * Strands,
 * Domains,
 * Complexes
 * System graphs
 * Structural elements
Reactions: (this could be a separate object, but for now system state and reactions are integrated)
 * Possible hybridization reactions
 * Energies, energy model
 * Hybridization and dehybridization rates
 * Propensity functions
 *

"""

from __future__ import absolute_import, print_function, division
import os
#import random
from collections import defaultdict, namedtuple, Counter, deque
try:
    from itertools import zip_longest, chain
except ImportError:
    from itertools import izip_longest as zip_longest, chain # pylint: disable=E0611
# Using namedtuple is quite slow in CPython, but faster than dicts in PyPy.
# (Attribute accessing namedtuple.a is faster than dict['a']; not sure about instantiation.)
#import math
from math import exp, log, sqrt #, log as ln
ln = log
#from datetime import datetime
from pprint import pprint
import networkx as nx
# from networkx.algorithms.components import connected_components, connected_component_subgraphs
import numpy as np
import pdb
import random

from nascent.energymodels.biopython import DNA_NN4_R, hybridization_dH_dS
# DNA_NN4_R = energy_tables_in_units_of_R['DNA_NN4']
from .constants import R, N_AVOGADRO # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K:
# from .constants import AVOGADRO_VOLUME_NM3
from .constants import PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION
from .constants import REACTION_NAMES
from .constants import ReactionAttrs, RA_HYB_INTRA, RA_HYB_INTER, RA_DEHYB_INTRA
# from .complex import Complex
from .componentmgr import ComponentMgr
from .nx_utils import draw_graph_and_save, layout_graph, draw_with_graphviz
from .debug import printd, pprintd
from .domain import Domain, DomainEnd
from .utils import sequential_number_generator
from .reaction_graph import ReactionGraph, reaction_attr_to_str, reaction_to_str, reaction_spec_pair_to_str


reaction_graph_sequantial_number_generator = sequential_number_generator()  # Used for reaction_graph filenames.
cmplx_state_sequential_number_generator = sequential_number_generator()     # Used for enumerating complex states
# Using sequential numbers instead of hashes might be a little easier on the eyes but otherwise have no effect.
cmplx_state_hashes = []         # number => cstate hash
cmplx_state_enum_by_hash = {}   # cstate_hash => enum (runtime specific)




debug_test_throttles = False
test_throttles = False
# test_throttles = {True: 0}     # Use a boolean True dict to effectiely freeze throttles at 1.0



class ReactionMgr(ComponentMgr):
    """
    ReactionMgr -> ComponentMgr -> GraphManager

    ReactionMgr class is in charge of calculating reactions and making sure reactions are performed properly.
    - But it does not control *which* state changes are invoked, that is done by the simulation system.
    ReactionMgr currently also holds the system state (volume and temperature), but only because we don't need a
    separate class for that.

    ComponentMgr: In charge of manipulating strands/domains/domain-ends and complexes through the
    hybridize/dehybridize/stack/unstack and join_/break_complex methods.
    ComponentMgr does not know anything about reactions or kinetics, it just performs the required
    manipulations of domais, complexes and system graphs when asked to
    The most important product delivered by ComponentMgr is the "result" dict returned by its methods,
    which is fed to ReactionMgr.post_reaction_processing method.

    GraphManager provides the underlying system graphs (strand-, domain- and DomainEnd-graphs)
    and everything related to structural analysis using those graphs:
    Loops and intra-complex activity calculation, etc.
    """

    def __init__(self, strands, params, volume=None, domain_pairs=None, init_reactions=False):
        ComponentMgr.__init__(self, strands=strands, params=params, volume=volume, domain_pairs=domain_pairs)

        self.include_steric_repulsion = False # True
        self.N_invoked_reactions_count = 0  # Keep count of how many reactions; See also simulator.N_steps

        ## Reaction settings/parameters:
        # When a complex is changed, update hybridization reaction for all pairs,
        # including hybridization to other complexes:
        self.reaction_update_pair_against_all = params.get("reaction_update_pair_against_all", True)

        ## Hybridization reaction parameters:
        self.hybridization_rate_constant_per_nt = params.get("hybridization_rate_constant_per_nt", 5e4)
        self.reaction_dehybridization_include_stacking_energy = params.get(
            "reaction_dehybridization_include_stacking_energy", True)

        ## Stacking reaction parameters:
        self.stacking_rate_constant = params.get('stacking_rate_constant', 1e4)
        self.enable_intercomplex_stacking = params.get("enable_intercomplex_stacking", False)


        ## Reaction graph: ##
        # reaction_graph is directed, but should it be a MultiGraph? I can't imagine why..
        # Can we have multiple edges from a source to a target?
        # There really should be only ONE reaction that takes a complex from state A to state B.
        # However, just to make sure we detect that this is indeed true, for now we use a MultiDiGraph,
        # keyed by (reacted_spec_pair, reaction_attr):
        # On second thought, it's probably better to detect this as early as possible,
        # rather than through later analysis. Let reaction_spec_pair and reaction_attr be a entries in edge_attr,
        # and then check that these are the same whenever we traverse an existing edge.
        # Fixed: Reaction_graph was initially a MultiDiGraph with edges keyed by (reacted_spec_pair, reaction_attr).
        self.reaction_graph = ReactionGraph(params=params) # nx.DiGraph()
        # Set graph-level attributes incl default node and edge attributes:
        self.reaction_graph.graph['node'] = {'fontname': 'Courier new'}
        self.reaction_graph.graph['edge'] = {'fontname': 'Arial',
                                             'fontsize': 10.0, # labelfontsize
                                             'len': 4.0}
        # self.reaction_graph.add_node(0) # The "null" node from which all new complexes emerge.
        self.endstates_by_reaction = {}  # [startstate][(reaction_spec_pair, reaction_attr)] = endstate
        self.endstates_by_reaction[0] = defaultdict(list) # Also adding the "null" node
        self.reverse_reaction_key = {} # get reaction edge key for the opposite direction.
        # self.reverse_reaction[(rx_spec_pair, rx_attr)] = (rev_rx_spec_pair, rev_rx_attr)
        # used to be a dict {edge_key => target} but we can have multiple targets for a single reaction.
        self.reaction_graph_complexes_directory = params.get("reaction_graph_complexes_directory")
        if self.reaction_graph_complexes_directory is not None:
            if not os.path.exists(self.reaction_graph_complexes_directory):
                print("Creating directory for complex files:", self.reaction_graph_complexes_directory)
                os.makedirs(self.reaction_graph_complexes_directory)
            assert os.path.isdir(self.reaction_graph_complexes_directory)


        ## Reaction cycle detection: (currently not used)
        self.reaction_microcycles_slice_size = params.get("reaction_microcycles_slice_size")
        self.known_reaction_cycles = set()          # each element is a tuple of (reaction_pair, reaction_pair, ...)
        self.reaction_cycles_count = Counter()
        self.reaction_cycles_by_pair = defaultdict(set)   # reaction_pair: [list of micro-cycles]
        # Keys include: reaction_spec_pair, reaction_str, reaction_attr_str, activity, c_j(_throttled),
        # throttle_factor, reaction_invocation_count, etc.
        self.reaction_graph_edge_label_fmt = ("{reaction_str}  \n" # "{reaction_attr_str}  \n"
                                              "{activity:0.02g}  {c_j:0.03g}  {reaction_invocation_count}  "
                                              "{throttle_factor:0.02g}  {c_j_throttled:0.03g}")

        # Prepare existing reaction_graph nodes:
        if strands:
            # Add un-complexed strand's names as nodes in the reaction graph:
            # See self.record_new_complex_state
            self.state_counter = Counter([strand.name for strand in self.strands if strand.complex is None])
            for strand in strands:
                x, y, z = (-3.0, 1.0 - 2*random.random(), 0.0)
                # Some graph layout algorithms (e.g. graphviz) can only handle two-valued (x,y) positions:
                # self.update_reaction_graph_state_node(strand.name, n_strands=1, x=x, y=y, z=z) #, pos=[x,y])
                self.reaction_graph.add_node(strand.name, n_strands=1, pos=[x, y, z])
                self.endstates_by_reaction[strand.name] = defaultdict(list)
            existing_complexes = {strand.complex for strand in strands if strand.complex is not None}
            if existing_complexes:
                print("\nReactionMgr: One or more strand came attached to existing complex: %s" % existing_complexes)
                self.complexes |= existing_complexes
                self.state_counter += Counter([cmplx.state_fingerprint() for cmplx in existing_complexes])
                added_state_nodes = set()
                for cmplx in existing_complexes:
                    if cmplx.state_fingerprint() in added_state_nodes: continue
                    x, y, z = (-3.0+len(cmplx.strands), 1.0 - 2*random.random(), 0.0)
                    # Some graph layout algorithms (e.g. graphviz) can only handle two-valued (x,y) positions:
                    # self.update_reaction_graph_state_node(strand.name, n_strands=1, x=x, y=y, z=z) #, pos=[x,y])
                    self.reaction_graph.add_node(strand.name, n_strands=len(cmplx.strands), pos=[x, y, z])
                    self.endstates_by_reaction[strand.name] = defaultdict(list)
                    added_state_nodes.add(cmplx.state_fingerprint())
                del added_state_nodes


        ## Reaction throttle:  (doesn't work, yet)
        self.reaction_throttle = params.get("reaction_throttle", False)

        # Variables for calculating throttle vs Number of reaction invocations:
        self.reaction_throttle_fun = params.get("reaction_throttle_fun")
        self.reaction_throttle_offset = params.get("reaction_throttle_offset", 0)
        self.reaction_throttle_rolling_fraction = params.get("reaction_throttle_rolling_fraction", False)

        # Variables for using relative cache instead of calculating new based on absolute Nric
        self.reaction_throttle_use_cache = params.get("reaction_throttle_use_cache", self.reaction_throttle)
        self.reaction_throttle_factor_base = params.get("reaction_throttle_factor_base", 0.99)
        # Nric = normalized_reaction_invocation_count is calculated as:
        # sysmgr.reaction_invocation_count[reaction_spec]/sum(len(sysmgr.domains_by_name[d.name] for d in (d1, d2)))
        # Throttle reaction c_j based on per-complex reaction counts rather that system-wide count.
        self.reaction_throttle_per_complex = params.get("reaction_throttle_per_complex", True)

        if self.reaction_throttle:
            if self.reaction_throttle_fun is None:
                self.reaction_throttle_fun = lambda Nric: exp(-Nric/10)
            elif isinstance(self.reaction_throttle_fun, str):
                self.reaction_throttle_fun = eval(self.reaction_throttle_fun)  # pylint: disable=W0123
        self.reaction_throttle_cache = {}  # defaultdict(lambda: 1.0)
        self.reaction_throttle_freeze = False

        ## Stats:
        # Every time a complex changes, record the new fingerprint here:
        self.complex_state_encounters = Counter()  # keyed by <c state>, value incremented at every occurrence.
        # How many times has a reaction been fired, by (reaction_pair, reaction_attr)
        self.reaction_invocation_count = Counter() # keyed by <pair>, value incremented for every reaction invocation
        # Keep track of the last N reactions:
        self.reaction_invocation_deque = deque(maxlen=self.params.get("reaction_throttle_rolling_window_size", 40))
        # This is slightly more light-weight than the "complex states at time t" collected by StatsManager.

        ## Symbol nomenclature:
        ## Sᵢ, S₁, S₂ - domain species (domain name). Domains with the same sequence are of the same specie.
        ## Fᵢ, F₁, F₂ - domain state fingerprint - domains of the same specie and in the same strand/complex state.
        ## Dᵢ, D₁, D₂ - domain instances, individual domain molecules.

        ### Caches: ###
        # Cache initialized in ComponentMgr

        ## Relative activity - used to moderate the activity of domains based on their accessibility in a complex.
        # Note: Not sure whether to use this or use a "loop energy" type of activation energy.
        # Relative activities, indexed as:
        #  - [({domain1-specie, domain2-specie}, complex-state-fingerprint)] for intra-complex reactions
        #  - [(domain1-specie, domain1-complex-state-fingerprint)] for inter-complex reactions.
        # relative activity is 1 for free strands
        # For inter-complex reactions (including free strands), the combined relative activity is simply the product:
        #   rel_activity = dom1_rel_activity * dom2_rel_activity

        # Activities and c_j caches, indexed by {d₁.F, d₂.F} or {(h1e3p.F, h2e5p.F), (h2e3p.F, h1e5p.F)}
        # self.cache['intracomplex_activity'] = {}  # Moved to ComponentMgr
        self.cache['stochastic_rate_constant'] = {}
        # stochastic_rate_constant_attrs: For debugging, I believe.
        # Has a lot of information on each reaction in cache['stochastic_rate_constant']
        self.cache['stochastic_rate_constant_attrs'] = defaultdict(list)
        # For surface/steric/electrostatic reduction in k_on.
        # TODO: Check surface-hybridization kinetics papers whether reduction in k_on yields same reduction in k_off.
        #       I.e. if the surface reduces Tm, hybridization energy or equivalent.
        #       For electrostatic and possibly steric, I would expect a reuction in Tm but not for surfaces,
        #       so it is a question whether these can in fact be treated equally.
        self.cache['domain_accessibility'] = {}     # indexed by F₁
        # Note: F₁, F₂ is domain-state-fingerprints. It is supposed to be invariant for domain instances
        # of the same domain specie in the same complex state. E.g. two free "domain A" should be in the same state.
        # Alternative name: dom_specie_spec

        #self.relative_activity_cache = defaultdict(dict)

        self.possible_hybridization_reactions = {}  # {d1, d2} => c_j
        self.reaction_attrs = {}                    # {d1, d2} => ReactionAttrs(reaction_type, is_forming, is_intra)
        self.reaction_spec_pairs = {} # reaction_pair (objects) => reaction_spec_pair
        # reaction_pair (pair of objects) => reaction_spec_pair (state dependent hash, shared between obj instances)

        # TODO: Change to generic (reaction_type, is_forming, is_intra)
        self.reaction_pairs_by_domain = defaultdict(set) # domain => {d1, d2} domain_pair of valid reactions.

        # self.stacking_interactions = {}         # {(end5p, end3p), (end5p, end3p)}
        self.possible_stacking_reactions = {}   # {(end5p, end3p), (end5p, end3p)} => c_j
        ## What stacking interactions to keep track of?
        ## end5p:end3p
        ## Only ends of hybridized domains (duplexes)
        ## Do we include ALL possible stacking interactions?
        ## or only to nearby ends?
        ## Strategies:
        ## a) Calculate equivalently to domain hybridization?
        ## b) (For intra-complex stacking): Go over all stacking interactions and
        ##      shuffle stacking state according to partition function.
        ## Will single-stranded neighbours produce steric interferrence? Will double-stranded?
        ## I need stacking state ends finger-prints. perhaps just FE = (FD, "5p"/"3p")
        ## where FD is domain state fingerprint.

        # When we are not grouping reactions by domain state species, then propensity functions are just the same
        # as propensity constants: a_j == X₁ X₂ c_j == c_j, since  X₁==X₂==X₁₂==1 for ungrouped reactions.
        self.hybridization_propensity_functions = self.possible_hybridization_reactions
        self.stacking_propensity_functions = self.possible_stacking_reactions

        if init_reactions and strands is not None:
            self.init_possible_reactions()

        ## Invoked reactions file, for re-playing the simulation:
        self.invoked_reactions_file = None
        if params.get('save_invoked_reactions_to_file'):
            fn = os.path.join(params.get('working_directory', '.'), "invoked_reactions.py") \
                if params['save_invoked_reactions_to_file'] is True else params['save_invoked_reactions_to_file']
            self.invoked_reactions_file = open(fn, 'w')

        print("Simulation system manager initiated at V=%s with %s strands spanning %s domains." \
              % (self.volume, len(self.strands), len(self.domains)))


    def reset_temperature(self, T):
        """ Set system temperature to T and and reset specified caches. """
        self.temperature = T
        if self.params.get("reaction_throttle_reset_on_temperature_change"):
            self.reaction_invocation_count.clear()
        self.cache['intracomplex_activity'].clear()
        self.cache['stochastic_rate_constant'].clear()
        self.init_possible_reactions()  # Update reaction propensities whenever temperature is changed.
        # pdb.set_trace()


    def init_possible_reactions(self):
        """
        Reactions have:
            reaction_spec => propensity_constant c_j, state_change_vector v_j
        However, I do:
            ({domain1, domain2}, is_forming, is_intracomplex) => propensity_constant.
        Where domspec is a domain-state fingerprint, and
            {domspec1, domspec2} = domspec_pair = {F₁, F₂}
            F₁, F₂ are symbols for state-dependent "finger-prints" aka "state-species" for the domains:
                hash((dspecie, c_state, self.in_complex_identifier()))

        We need to find possible reactions. Possible ways to do that:
         1. Use self.domain_pairs - domain species names.
         2. Go over all domains.

        Although I expect all domains to be in the same start at the onset of the simulation,
        I figure it is best to go through them all, determining their state-specie,
        just to make sure I have it right.

        I am grouping domains into sub-populations according to their state.
        I do this to reduce the memory required for caching and increasing the chance of cache hits.

        I am basically just using a "domain subspecie state hash" (domspec) to denote a set of
        particular domains all in the same state.
        E.g. a domspec 928457 could be denote any/all "domain1 on strandA that are free (not in a complex)"
        To get a list of all actual domains matching a certain domain subpopulation domspec, use:
            self.domain_state_subspecies[domspec] - will you all domains matching <domspec>
        Maybe a better name would be domain_subpopulation_by_state?
        """
        # printd("\nInitializing possible reactions (UNGROUPED)...")
        self.update_possible_hybridization_reactions(changed_domains=None, recalculate_hybridized=True)
        self.update_possible_stacking_reactions(changed_domains=None, recalculate_stacked=True)
        print(len(self.possible_hybridization_reactions), "possible hybridization reactions (re-)initialized.")


    def hybridization_rate_constant(self, d1, d2):
        """ Calculate k_on hybridization rate constant between domains :d1: and :d2:. """
        #T = self.temperature
        # At T >> Tm (ΔG° ~ 0), this depends on ΔG°.
        # Also depends on:
        # * salt / ionic strength
        # * length of d1, d2 (and the full length of their strands)
        # But for now, just return default rate constant of 1e5 /s/M
        # Another question: What about intra-complex hybridization?
        #   Should this be accounted for here at the hybridization rate constant,
        #   or at the propensity constant?
        ## TODO: Make hybridization k_on depend on domain length.
        ## This ensures that if we break a domain in two, the equilibrium remais the same:
        ## K1 = k1_on/k1_off = exp(ΔG1°/RT)
        ## K2 = k2_on/k2_off = exp((ΔG2a°+ΔG2a°)/RT)
        ## K2 = k2_on/k2_off * k2_on/k2_off = exp(ΔG2a°/RT) * exp(ΔG2b°/RT)
        ## But it is not an AND combination, it's an OR: We only need either 2A or 2B to be hybridized.
        ##                  /-- k1b on/off ->  (B:b, C, c)  -- k2c on/off --\
        ## (B, C) + (b, c) <                                                 >  (B:b, C:c)
        ##                  \-- k1c on/off ->  (C:c, B, b)  -- k2b on/off --/
        ## That combination gives us:
        ## K2 = ([B_h, C_u] + [B_u, C_h]) / ([B_u]*[C_u]) = k1b_on/k1b_off + k1c_on/k1c_off
        ## K1 == K2
        ## k1_on = k1_off*exp(dS1-dH1/T), k2_on = k2_off*exp(dS2-dH2/T)
        ## Essentially, the off rate for case 2 (unhybridize B:b AND C:c) is the same as
        ## the off rate for case 1 (unhybridize A:a).
        ## However, since we have two reaction pathways (and the same number of each starting domain),
        ## if k_on was the same, we would have twice as much reacting in case (2) in any time interval.
        return self.hybridization_rate_constant_per_nt*min((len(d1), len(d2)))      # A 20 nt duplex will have k_on = 1e6 /M/s


    def dehybridization_rate_constant(self, d1, d2):
        """
        Calculate dehybridization rate constant:
            K = exp(-ΔG°/RT) = exp(-(ΔH°-TΔS°)/RT) = exp(ΔS°/R-ΔH°/R/T) = exp(ΔS°-ΔH°/T) # ΔH°, ΔS° in units of R, R/K
            k_off = k_on/K,
                  = k_on * exp(+ΔG°/RT)
        If T < Tm at molar concentrations (ΔG° < 0), hybridization rate, k_on, is approximately constant at 1e5-1e6.
        If T > Tm at molar concentrations (ΔG° > 0), hybridization rate goes down.
        In practice, we don't have to worry too much about.
        """
        T = self.temperature
        # dH, dS, dG at standard conditions
        dH, dS = self.dHdS_from_state_cache(d1, d2) # Units of R, R/K
        #dG = dH - T*dS # No need to perform extra operation, just use dH and dS:

        # Consideration: Do we really need to have a state-specific energy cache?
        #   a) Yes: Because we might have bending energies, etc which have complex state dependence
        #   b) No:  It should be pleanty fast to add additional energy contributions on the fly.

        ## Add additional energy contributions:
        ## - Stacking energies
        ## - Other? Bending?

        assert all(end.stack_partner is None for end in (d1.end3p, d1.end5p, d2.end3p, d2.end5p))
        # if self.reaction_dehybridization_include_stacking_energy:
        #     stack_dH, stack_dS = self.duplex_stacking_energy(d1, d2)
        #     if any((stack_dH, stack_dS)):
        #         print("\nAdding stack_dH, stack_dS = %s, %s to dehybridization energy (%s and %s)" %
        #               (stack_dH, stack_dS, d1, d2))
        #     dH += stack_dH
        #     dS += stack_dS

        K = exp(dS - dH/T)
        k_on = self.hybridization_rate_constant(d1, d2)  # 1e5  # Is only valid for if K > 1
        k_off = k_on/K   # = k_on * exp(dH/T-dS)
        # printd("Hybridization energy for domain %s with sequence %s" % (d1, d1.sequence))
        # printd("  dS=%0.03g, dH=%.03g, T=%s, dS-dH/T=%.03g K=%0.02g, k_off=%0.02g" % (dS, dH, T, dS - dH/T, K, k_off))
        return k_off


    def unstacking_rate_constant(self, h1end3p, h2end3p, T=None):
        """
        Rate for breaking stacking interaction between duplex ends
        (giving the domain 3p end of involved duplexes).
        Note: stacking_rate_constant is a constant attribute, not a function.
        """
        if T is None:
            T = self.temperature
        h1end5p = h2end3p.hyb_partner
        h2end5p = h1end3p.hyb_partner
        ## TODO: Remove assertions
        assert h1end3p.stack_string == h2end3p.stack_string
        assert h1end3p.stack_string == h2end3p.stack_string == h1end5p.stack_string == h2end5p.stack_string
        assert h1end3p.stack_string is not None
        stack_dH, stack_dS = DNA_NN4_R[h1end3p.stack_string] # Why does this have None as key?
        ## TODO: Remove stack_string equality assertion
        # K = exp(stack_dS - stack_dH/T)
        # k_off = self.stacking_rate_constant/K
        k_off = self.stacking_rate_constant * exp(stack_dH/T - stack_dS)
        return k_off


    def calculate_c_j(self, elem1, elem2, reaction_attr, reaction_spec_pair):
        """
        General method for calculation of c_j stochastic rate constant.
        Unlike the specialized calculate_hybridization_c_j and calculate_stacking_c_j,
        this method also takes care of caching and throttling.
        Cached values does *not* include throttle factor.
        :elem1:, :elem2: Either two domains (if hybridization reaction) or
                         two two-tuples of (h-end3p, h-end5p) pairs (if stacking interaction).
        :reaction_attr: A ReactionAttr namedtuple with (reaction_type, is_forming, is_intra) values.
        :reaction_spec_pair:    A frozenset of state fingerprints for each instance in elem1, elem2.

        If dehybridization and unstacking are state-independent, then we should use a
        state-independent c_j_cache_key to save c_j for dehybridization and unstacking reactions.

        """
        assert reaction_spec_pair != frozenset((elem1, elem2))
        assert reaction_attr.reaction_type != HYBRIDIZATION_INTERACTION or \
            reaction_attr.is_forming or reaction_attr.is_intra
        reaction_type, is_forming, is_intra = reaction_attr

        # c_j_cache_key = (reaction_spec_pair, reaction_attr)
        # If dehybridization and unstacking are state-independent, then we should use a
        # state-independent c_j_cache_key to save c_j for dehybridization and unstacking reactions.
        if reaction_attr.is_forming:
            c_j_cache_key = (reaction_spec_pair, reaction_attr)
        else:
            if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                c_j_cache_key = (frozenset((elem1.name, elem2.name)), reaction_attr)
            elif reaction_attr.reaction_type is STACKING_INTERACTION:
                assert elem1[0].stack_string == elem1[1].stack_string == elem2[0].stack_string == elem2[1].stack_string
                c_j_cache_key = (elem1[0].stack_string, reaction_attr)
            else:
                raise ValueError("Unknown reaction type %s" % reaction_attr.reaction_type, reaction_attr)
        if c_j_cache_key in self.cache['stochastic_rate_constant']:
            assert reaction_spec_pair is not None # or reaction_attr.is_forming is False
            # TOOD: Make sure that this includes stacking state!
            # TODO: We actually have two kinds of species fingerprints: One that allows symmetry and one that don't
            # TODO: We need to differentiate between intra-complex and inter-complex reactions in the cache,
            #       otherwise we can't distinguish between intra-complex reaction and reaction between two identical
            #       domains in two complexes with the same state (although that is somewhat unlikely to happen).
            #       Probably just cache by (reaction_spec_pair, reaction_attr); reaction_attr = (type, forming, intra)
            c_j = self.cache['stochastic_rate_constant'][c_j_cache_key]
            #assert ra == reaction_attr
        else:
            # Calculate c_j and save result to cache
            if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                if reaction_spec_pair is None:
                    reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
                c_j = self.calculate_hybridization_c_j(elem1, elem2, is_forming=is_forming, is_intra=is_intra,
                                                       reaction_spec_pair=reaction_spec_pair)
            elif reaction_attr.reaction_type is STACKING_INTERACTION:
                if reaction_spec_pair is None:
                    reaction_spec_pair = frozenset(((elem1[0].state_fingerprint(), elem1[0].state_fingerprint()),
                                                    (elem2[0].state_fingerprint(), elem2[0].state_fingerprint())))
                c_j = self.calculate_stacking_c_j(elem1[0], elem1[1], elem2[0], elem2[1],
                                                  is_forming=is_forming, is_intra=is_intra,
                                                  reaction_spec_pair=reaction_spec_pair)
            else:
                raise ValueError("Unknown reaction type %r for reaction_attr %s between %s and %s" %
                                 (reaction_attr.reaction_type, reaction_attr, elem1, elem2))
            # Make sure to save c_j before applying throttle:
            self.cache['stochastic_rate_constant'][c_j_cache_key] = c_j
            print("c_j = %0.02e for reaction %s" % (c_j, reaction_to_str(reaction_spec_pair, reaction_attr)))

        if c_j == 0: # or perhaps something like c_j < 1e-12?
            # Short-circuit; do not attempt to apply throttle.
            return c_j

        ## Apply throttle: ##
        if self.reaction_throttle:
            assert reaction_spec_pair is not None
            #print("THROTTLING!")
            #if self.reaction_throttle_use_cache:
            # Per-reaction throttle is adjusted in post_reaction_processing with each reaction invocation.
            # This should be more efficient than calculating some exponential based on reaction invocation count.
            if self.reaction_throttle_per_complex:
                if reaction_attr.is_intra:
                    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                        cmplx = elem1.strand.complex
                    elif reaction_attr.reaction_type is STACKING_INTERACTION:
                        cmplx = elem1[0].domain.strand.complex
                    if reaction_spec_pair in cmplx.reaction_throttle_cache:
                        throttle_factor = cmplx.reaction_throttle_cache[reaction_spec_pair]
                        # if throttle_factor < 1 or True:
                        #     print("throttle_factor %0.02f for reaction" % throttle_factor,
                        #           elem1, elem2, reaction_attr, reaction_spec_pair)
                    else:
                        # print("reaction_spec_pair not in %s.reaction_throttle: %s, %s" %
                        #       (cmplx, reaction_spec_pair, reaction_attr))
                        #pdb.set_trace()
                        throttle_factor = 1
                else:
                    throttle_factor = 1
            else:
                # Global throttle:
                # Include reaction_attr in cache key to distinguish intra vs inter complex reactions:
                try:
                    throttle_factor = self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)][0]
                except KeyError:
                    if debug_test_throttles:
                        print("\nUn-throttled reaction:\n  ", reaction_to_str(reaction_spec_pair, reaction_attr))
                        print("Reaction throttle cache:")
                        print("\n".join("  %s: %s" % (reaction_to_str(*k), v) for k, v in self.reaction_throttle_cache.items()))
                        pdb.set_trace()
                    # It seems we are calculating c_j for reaction
                    # s+ :s2_A3p >_< s2_B5p / s1_b3p >_< s1_a5p    (67283)
                    # However, state "67283" is the already-stacked state. It should be for the unstacked state 19968.
                    throttle_factor = 1.0
            # else:
            #     raise NotImplementedError("Obsolete settings")
            #     # Calculate reaction throttle factor using function, typically exponentially decreasing with Nric.
            #     throttle_factor = self.calculate_throttle_factor(elem1, elem2, reaction_attr, reaction_spec_pair)
            c_j *= throttle_factor
        # else:
        #     # print("NOT throttling reaction:", elem1, elem2, reaction_attr)
        #     pass
        ## TODO: DEBUGGING code, remove this!
        # if abs(c_j - 48.1) < 0.2:
        #     assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
        #     assert reaction_attr.is_forming is False
        #     assert reaction_attr.is_intra is True
        #     assert elem1.name.lower() == 'b'
        #     cmplx = elem1.strand.complex
        #     assert cmplx == elem2.strand.complex
        #     d1fp = elem1.state_fingerprint()
        #     d2fp = elem2.state_fingerprint()
        #     cmplx.reset_state_fingerprint()
        #     assert d1fp == elem1.state_fingerprint()
        #     assert d2fp == elem2.state_fingerprint()
        #     # s1	a,t1,b	AAGGCACGAC,TT,AGGTTTCCAA
        #     # s2	B,t2,A	TTGGAAACCT,TT,GTCGTGCCTT
        #     b = elem1 if elem1.name == 'b' else elem2
        #     assert b.name == 'b'
        #     assert b.partner.name == 'B'
        #     #      b5p ----- t3p ------- t5p ------- a3p
        #     a = b.end5p.pb_upstream.pb_upstream.pb_upstream.domain
        #     assert a.name == 'a'
        #     # Correct one is where domain a is still hybridized:
        #     assert a.partner is not None
        #     assert a.partner.name == 'A'

        return c_j



    def calculate_hybridization_c_j(self, d1, d2, is_forming, is_intra, reaction_type=HYBRIDIZATION_INTERACTION,
                                    reaction_spec_pair=None):
        """
        :is_intra: whether the reaction is intra-molecular; that is, either intra-complex or intra-strand.
        Calculate propensity constant, c_j, for hybridization or dehybridization of d1 and d2.
        Edit: Changed name from "propensity constant" to just "c_j", to avoid confusion with propensity *function*, a_j.
        Alternative names:
        * specific/molecular/instance/individual reaction frequency
        * specific: specific to a single pair of molecules from S₁, S₂.
        * molecular/instance/individual: as opposed to the total propensity for the entire species population, a_j.
        * "frequency" because of unit s⁻¹.
        ** Could also be "propensity", but could be confused with propensity *function*, a_j.
        * "stochastic rate constant"  ⁽¹,²⁾
        * "specific probability rate constant" (³)
        Regarding "probability" vs "frequency":
        * c_j has units of s⁻¹, so probability (which is unitless) wouldn't be correct.
        * The product "c_j∙dt" gives the probability that a selected pair of instances
            will undergo reaction R_j in timespan dt.

        For bi-molecular reactions, the propensity constant is k_on rate constant times concentration:
            c = k_on * [domain]   # c is propensity, not concentration.
        For uni-molecular reactions, the propensity constant is simply the k_off rate constant:
            c = k_off

        Differences between rate constant, k_j, and c_j:
        * c_j is stochastic, k_j is *mass-action* rate constant.⁽¹⁾
        * c_j = k_j / (N_A ∙ volume)ᴺ⁻¹   # N is the number of reactants, e.g. 2 for bimolecular reactions.
        * v_j = k_j [S₁] [S₁]   # reaction rate, has unit of M/s. (Note: not state change vector ν_j, nu)
        * a_j = c_j x₁ x₂       # propensity function, has unit of s⁻¹

        Question: Where would it be appropriate to include effects of intra-complex reactions:
        * Effects related to effective volume / effective/relative activity should be included
            when calculating c_j. c_j is volume dependent, k_j is independent on volume.
        * Effects related to hybridization energy should be included in k_j,
            since k_j *is* dependent on hybridization energy. This could include:
            * Bending strain energy (decreasing k_on) - although bending could also be interpreted as influencing
                the spatial probability distribution function and thus be kind of volumetric in nature.
            * Zippering strain (increasing k_off),
            * Electrostatic interactions.
        Not sure about steric interactions. That is probably volumetric, affecting spatial probability
        districution function, to be exact. The same could be said about electrostatic...
        I guess if an interaction affects k_off, then it is energetic in nature, while if
        it only affects k_on, then it is volumetric in nature (spatial pdf).
        For a first approximation, maybe only consider the effects that can be interpreted in terms
        of effective concentration / relative activity. Then you can assume constant rate constants,
        k_on/k_off (except for stacking affecting k_off), but otherwise have all effect apply only
        to c_j and never k_j.
        These effect could even be isolated to a single factor, relative_activity(d1, d2).

        Question: What is the best way to capture volume pdf effects?

        Question: What should relative activity entail?
            Should rel act be 1 for two free strands - or 1/volume? -- Thus basically the stochastic activity..
            Maybe start by implementing it here and not worry about "relative" vs actual activity.

        TODO: Account for INTRA-complex reactions.
            Question: Should this be done here or in hybridization_rate_constant(d1, d2) ?
        TODO: Account for steric interferrence in INTER-complex reactions.

        Refs:
            [1]: http://people.uleth.ca/~roussel/C4000biochemnet/slides/14CME.pdf
            [2]: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1303974/ "stochastic rate constant cμ"
            [3]: http://www.async.ece.utah.edu/BioBook/lec4.pdf
        """
        ## Done: Make c_j a cached value [done through calculate_c_j generic method]

        assert is_forming or is_intra   # If we are not forming, we are de-hybridizing; is_intra should be True
        reaction_attr = ReactionAttrs(reaction_type, is_forming, is_intra)

        if is_forming:
            # if d1.strand.complex == d2.strand.complex != None:
            # Edit: We might have a single strand interacting with itself
            if is_intra:
                assert d1.strand.complex == d2.strand.complex
                assert d1.strand.complex is not None or d1.strand == d2.strand
                # Intra-complex reaction:
                # try:
                stochastic_activity = self.intracomplex_activity(d1, d2, reaction_type=reaction_type,
                                                                 reaction_spec_pair=reaction_spec_pair)
                # except nx.NetworkXNoPath as e:
                #     c = d1.strand.complex
                #     print(("ERROR: No path between d1 %s and d2 %s "
                #            "which are both in complex %s. ") % (d1, d2, c))
                #     print("complex.nodes():")
                #     pprint(c.nodes())
                #     print("complex.edges():")
                #     pprint(c.edges())
                #     workdir = self.params.get("working_directory", os.path.expanduser("~"))
                #     draw_graph_and_save(c, outputfn=os.path.join(workdir, "complex_graph.png"))
                #     draw_graph_and_save(self.strand_graph, outputfn=os.path.join(workdir, "strand_graph.png"))
                #     draw_graph_and_save(self.domain_graph, outputfn=os.path.join(workdir, "domain_graph.png"))
                #     draw_graph_and_save(self.ends5p3p_graph, outputfn=os.path.join(workdir, "ends5p3p_graph.png"))
                #     raise e
                if stochastic_activity == 0:
                    return 0
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
            if self.include_steric_repulsion:
                #
                # steric_activity_factor = np.prod([1 if d.strand.complex is None else self.steric_activity_factor(d)
                #                                   for d in (d1, d2)])
                steric_activity_factor = ((1 if d1.strand.complex is None else self.steric_activity_factor(d1)) *
                                          (1 if d2.strand.complex is None else self.steric_activity_factor(d2)))
                stochastic_activity *= steric_activity_factor

            k_on = self.hybridization_rate_constant(d1, d2)
            # hybridization rates of free strands can be approximated as being constant
            # (somewhere between 1e5 /M/s and 1e6 /M/s). This is because the activation energy is fairly negligible.
            # Hybridization rate constant, k, is in unit of /M/s = L/mol/s.
            # Activity is always in units of M - so resultant
            c_j = k_on * stochastic_activity # should be e.g. 0.1 /s

            ## Should we include stacking interactions in k_on? (That would increase k_on...)
            # a) Yes: We include them in k_off, so they must also always be included in k_on (because...)
            # b1) No: Stacking and unstacking is separate from hybridization/dehybridization:
            #           (unhybridized, unstacked) ←→ (hybridized, unstacked) ←→ (hybridized, stacked)
            #           To de-hybridize, a domain must first un-stack.
            #           AND WE ONLY CONSIDER DE-HYBRIDIZATION REACTION FOR UN-STACKED DOMAINS
            #           (No.. because stacking reactions might be much faster than hybridization)
            # b2) No: We have a separate reaction for stacking/unstacking.
            #           If the domain is stacked, that is included in k_off, but the domain can spontaneously unstack.
            # c) From a thermodynamic/kinetic perspective: Not sure.
            #           On one hand, I wouldn't expect stacking to increase hybridization rate.
            #           On second thought:
            # d) Thermodynamically: If you include stacking interactions in k_off, you CANNOT also include them in k_on,
            #           as this would effectively include it twice:
            #           ∆G/RT = ln(K) = ln(k_on/k_off) = ln(k_on) - ln(k_off)
            #                 = ln(k_on0*k_stacking) - ln(k_off0/k_stacking) = ln(k_on) - ln(k_off) + 2*ln(k_stacking)
        else:
            # De-hybridization reaction;
            # k_off depends on ΔG°  (at least when we are at T < Tm at molar concentrations where ΔG° < 0)
            # TODO: Ensure that k_off depends on end stacking interactions!
            k_off = self.dehybridization_rate_constant(d1, d2)
            # For uni-molecular reactions, c_j = k_j
            c_j = k_off
        return c_j #, is_forming


    def calculate_stacking_c_j(self, h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra,
                               reaction_type=STACKING_INTERACTION, reaction_spec_pair=None):
        """
        Calculate stochastic rate constant c_j for forming/breaking a stacking interaction.
        Ends annotation:
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p

        Note: Domains on the same helix may or may not be also connected by their phosphate backbone.
        E.g. you could have a hinge, where one helix is backbone-connected and the other one not.
        This is probably the most common case, e.g. in N-way junctions.
        Nomenclature:
        :reaction_spec_pair: (For state species - previously "reaction_pair_fingerprint")
         - For stacking interactions, this is called "stackspec_pair", for domains its called "domspec_pair".
         - The general case its called "reaction_spec_pair".

        Q: Where do we apply throttling? Here? Or under the "undate_*_reactions" methods?
        A: Maybe you could have a "calculate_c_j" wrapper which takes care of caching and throttling?
        """
        # Caching moved to common calculate_c_j method:
        # if reaction_spec_pair in self.cache['stochastic_rate_constant']:
        #     return self.cache['stochastic_rate_constant'][reaction_spec_pair]

        # Comment out in production:
        # Reset complex (and domain) state fingerprints and re-check:
        # complexes = {c for c in (end.domain.strand.complex for end in (h1end3p, h2end5p, h2end3p, h1end5p))
        #              if c is not None}
        # for c in complexes:
        #     c.reset_state_fingerprint()
        # assert reaction_spec_pair == frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
        #                                         (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        d1 = h1end3p.domain
        d2 = h2end3p.domain
        # printd("Re-setting complex state and re-calculating state fingerprint:")
        # d1.state_change_reset()
        # d2.state_change_reset()
        # d1.strand.complex.reset_state_fingerprint(reset_domains=True)
        # d2.strand.complex.reset_state_fingerprint(reset_domains=True)
        # reaction_spec_pair_fresh = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
        #                                       (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        # try:
        #     assert reaction_spec_pair == reaction_spec_pair_fresh
        # except AssertionError as e:
        #     print("Got different reaction_spec_pair fingerprints after resetting domains/complexes!")
        #     raise e
        assert h1end3p.hyb_partner == h2end5p
        assert h1end3p == h2end5p.hyb_partner
        assert h2end3p.hyb_partner == h1end5p
        assert h2end3p == h1end5p.hyb_partner
        # printd("calculate_stacking_c_j invoked with h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra:")
        # printd(", ".join(str(p) for p in (h1end3p, h2end5p, h2end3p, h1end5p, is_forming, is_intra)))
        if is_forming:
            # if d1.strand.complex == d2.strand.complex != None:
            # Edit: We might have a single strand interacting with itself
            if is_intra:
                assert d1.strand.complex == d2.strand.complex
                assert d1.strand.complex is not None or d1.strand == d2.strand
                # Intra-complex reaction:
                # This should really be the average of dist(h1end3p, h1end5p), dist(h2end3p, h2end5p).
                # However, except for the case where the helices are already backbone connected,
                # we can approximate this by just one distance, dist(h1end3p, h2end3p): (TODO: Re-check this...)
                if h1end3p.pb_downstream == h1end5p and h2end3p.pb_downstream == h2end5p:
                    stochastic_activity = 1
                # (h1end3p, h2end5p), (h2end3p, h1end5p)
                stochastic_activity = self.intracomplex_activity((h1end3p, h2end5p), (h2end3p, h1end5p),
                                                                 reaction_type=STACKING_INTERACTION,
                                                                 reaction_spec_pair=reaction_spec_pair)
                # printd("reaction_spec_pair:")
                # pprintd(reaction_spec_pair)
                # printd(" - stacking stochastic_activity (intra reaction) = %s" % stochastic_activity)
                if stochastic_activity == 0:
                    return 0
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
                # printd(" - stacking stochastic_activity (inter reaction) = %s" % stochastic_activity)
            if self.include_steric_repulsion:
                steric_activity_factor = ((1 if d1.strand.complex is None else self.steric_activity_factor(d1)) *
                                          (1 if d2.strand.complex is None else self.steric_activity_factor(d2)))
                stochastic_activity *= steric_activity_factor

            k_on = self.stacking_rate_constant # (h1end3p, h2end3p)    # is fairly constant as long as ΔG° < 0
            c_j = k_on * stochastic_activity
            # printd(" - stacking stochastic rate constant c_j = %s * %s = %s" % (k_on, stochastic_activity, c_j))
        else:
            # De-hybridization reaction;
            k_off = self.unstacking_rate_constant(h1end3p, h2end3p)  # k_off depends on ΔG°
            c_j = k_off # For uni-molecular reactions, c_j = k_j
            # printd(" - un-stacking stochastic rate constant c_j = k_off = %s" % (c_j,))

        return c_j


    def update_possible_hybridization_reactions(self, changed_domains, reacted_pair=None,
                                                reaction_attr=None, reacted_spec_pair=None,
                                                pair_against_all=True, recalculate_hybridized=None,
                                                recalculate_these=None):
        """
        First, process the reacted_pair: if pair was formed (hybridized), remove each domain
        from other possible reactions involving either of the two domains.

        Second, process changed_domains. If :changed_domains: is None, it defaults to ALL domains.
        We only consider un-hybridized domains, unless :recalculate_hybridized: is set to True.

        Arguments:
        :changed_domains:   A list/set of domains that were affected by a recent reaction.
            If None, re-calculate *all* domains.
        :reacted_pair:      The domains that were just reacted, prompting this re-calculation.
        :reaction_attr:     ReactionAttr namedtuple for the reaction that was just triggered.
        :reacted_spec_pair: Should be frozenset(d for d in reacted_pair), can be provided if you already have it.

        :pair_against_all: If True, consider/recalculate reactions against all other un-hybridized domains.
            Otherwise, only consider/recalculate reactions on within the same complex.
        :recalculate_hybridized: If True, re-calculate all changed domains, even hybridized domains,
            which are by default not re-calculated.
        :recalculate_these: (list/set/tuple) If given, domains in this set are always re-calculated, regardless of
            whether the domain is hybridized or not.

        Variable nomenclature:
            d1, d2: The domains that were just reacted.
            domain, domain2: Loop variables for re-calculating hybridization reactions for changed domain-domain pairs.

        New: Domains can only dehybridize if they are unstacked.
        """
        if recalculate_these is None:
            recalculate_these = set()
        if changed_domains is None:
            changed_domains = self.domains
            if recalculate_hybridized is None:
                recalculate_hybridized = True
        if self.reaction_throttle:
            # If reactions are throttled, we cannot assume that all dehyb/unstack reactions have state-independent c_j.
            recalculate_hybridized = True
        # if reacted_spec_pair is None:
        #     reacted_spec_pair = frozenset((d.state_fingerprint() for d in reacted_pair))
        ## TODO: Consolidate this and init_possible_reactions into a single method.
        # printd(("\nupdate_possible_hybridization_reactions invoked with reacted_pair=%s, reaction_attr=%s, "
        #         "reacted_spec_pair=%s, changed_domains:") % (reacted_pair, reaction_attr, reacted_spec_pair))
        # pprintd(changed_domains)

        # n_possible_start = len(self.possible_hybridization_reactions)
        # n_species_start = len(self.domain_state_subspecies)

        # old_d1d2_doms_specs = frozenset((d1._specie_state_fingerprint, d2._specie_state_fingerprint))
        updated_reactions = set()
        # old_domspecs = {domain: domain._specie_state_fingerprint
        #                 for domain in changed_domains}

        ## PROCESS THE REACTED PAIR:
        if reacted_pair:
            is_forming = reaction_attr.is_forming  # (was forming - For the reaction that was just processed.)

            ## Check if stacking reaction or hybridization reaction.
            if reaction_attr.reaction_type is STACKING_INTERACTION:
                #(h1end3p, h1end5p), (h2end3p, h2end5p) = tuple(reacted_pair)    # Alternative DomainEnd pairing
                (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(reacted_pair)
                if is_forming:
                    # We have just formed a new stacking interaction; stacked domains cannot dehybridize:
                    for d1, d2 in ((h1end3p.domain, h2end5p.domain), (h2end3p.domain, h1end5p.domain)):
                        ## Add a stacked_domains_by_name?  unhybridized <-> hybridized <-> stacked
                        # if is_forming:
                        #     self.hybridized_domains_by_name[d1.name].remove(d1)
                        #     self.hybridized_domains_by_name[d2.name].remove(d2)
                        #     self.stacked_domains_by_name[d1.name].add(d1)
                        #     self.stacked_domains_by_name[d2.name].add(d2)
                        # else:  # reaction
                        #     self.stacked_domains_by_name[d1.name].remove(d1)
                        #     self.stacked_domains_by_name[d2.name].remove(d2)
                        #     self.hybridized_domains_by_name[d1.name].add(d1)
                        #     self.hybridized_domains_by_name[d2.name].add(d2)
                        # IF the duplex is already stacked in the other end, then the *dehybridization* reaction (for
                        # the two domains of that duplex) is already removed from possible_hybridization_reactions.
                        domain_pair = frozenset((d1, d2))
                        try:
                            del self.possible_hybridization_reactions[domain_pair]
                        except KeyError:
                            ## TODO: Check that the other end of the duplex is indeed stacked.
                            pass
                        else:
                            del self.reaction_attrs[domain_pair]
                            del self.reaction_spec_pairs[domain_pair]
                        self.reaction_pairs_by_domain[d1].clear()
                        self.reaction_pairs_by_domain[d2].clear()
                else:
                    # We have just broken a stacking interaction;
                    # if any domain duplex pair is unstacked in both ends,
                    # we should add domain pair to possible_hybridization_reactions.
                    # We do this below when we recalculate changed domains:
                    recalculate_these.add(h1end3p.domain)
                    recalculate_these.add(h2end3p.domain)
            elif reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                d1, d2 = tuple(reacted_pair)

                if is_forming:
                    ## TODO: Add a stacked_domains_by_name?  unhybridized <-> hybridized <-> stacked
                    self.unhybridized_domains_by_name[d1.name].remove(d1)
                    self.unhybridized_domains_by_name[d2.name].remove(d2)
                    self.hybridized_domains_by_name[d1.name].add(d1)
                    self.hybridized_domains_by_name[d2.name].add(d2)
                else:
                    self.hybridized_domains_by_name[d1.name].remove(d1)
                    self.hybridized_domains_by_name[d2.name].remove(d2)
                    self.unhybridized_domains_by_name[d1.name].add(d1)
                    self.unhybridized_domains_by_name[d2.name].add(d2)

                # If d1, d2 is hybridizing, then we need to eliminate reactions involving d1, d2 from possible_rxs.
                if is_forming:
                    obsolete_reactions = [domain_pair for domain_pair in self.possible_hybridization_reactions
                                          if d1 in domain_pair or d2 in domain_pair]
                    expected_reactions = self.reaction_pairs_by_domain[d1] | self.reaction_pairs_by_domain[d2]
                    # the tracked "expected" can easily include reactions that have been removed.
                    # For instance: We add {A#1, a#1} to the set. Then A#1 is hybridized to a#2.
                    # We delete reaction_pairs_by_domain[A#1] and [a#2], but a#1 is still there.
                    # Thus, reaction_pairs_by_domain is only a "suggested" set of rxs.
                    if len(set(obsolete_reactions) - expected_reactions) > 0:
                        print("\nlen(set(obsolete_reactions) - expected_reactions) > 0:")
                        print("set(obsolete_reactions):")
                        pprint(set(obsolete_reactions))
                        print("expected_reactions:")
                        pprint(expected_reactions)
                        print("set(obsolete_reactions)-expected_reactions:")
                        pprint(set(obsolete_reactions)-expected_reactions)
                        print("expected_reactions-set(obsolete_reactions):")
                        pprint(expected_reactions-set(obsolete_reactions))
                    else:
                        # printd("\nobsolete_reactions MATCHES expected  (%s elements)" % len(expected_reactions))
                        pass

                    self.reaction_pairs_by_domain[d1].clear()
                    self.reaction_pairs_by_domain[d2].clear()
                    for domain_pair in obsolete_reactions:
                        del self.possible_hybridization_reactions[domain_pair]
                        del self.reaction_attrs[domain_pair]
                        del self.reaction_spec_pairs[domain_pair]
                    # You can keep track of domain hybridization reactions, but then
                    # please note that the reaction can have become obsolete by the other domain being hybridized.
                    # It is thus only a set of possible reactions that can be deleted.
                    # Add de-hybridization reaction for d1, d2: (will be skipped when processing all changed domains)
                    domain_pair = reacted_pair # frozenset((d1, d2))
                    domain_spec_pair = frozenset((d1.state_fingerprint(), d2.state_fingerprint()))
                    # Calling calculate_hybridization_c_j or any other can produce a call to domain.state_fingerprint()
                    ra = ReactionAttrs(reaction_type=HYBRIDIZATION_INTERACTION, is_forming=False, is_intra=True)
                    # Use calculate_c_j(elem1, elem2, reaction_attr, reaction_spec_pair) to enable throttle and caching.
                    # cj = self.calculate_hybridization_c_j(d1, d2, is_forming=False, is_intra=True)
                    c_j = self.calculate_c_j(d1, d2, reaction_attr=ra, reaction_spec_pair=domain_spec_pair)
                    if c_j > 0:
                        # if abs(c_j - 48.1) < 0.2: pdb.set_trace()  ## TODO: Debug code, remove this!
                        self.possible_hybridization_reactions[domain_pair] = c_j
                        self.reaction_attrs[domain_pair] = ra
                        self.reaction_spec_pairs[domain_pair] = domain_spec_pair
                        self.reaction_pairs_by_domain[d1].add(domain_pair)
                        self.reaction_pairs_by_domain[d2].add(domain_pair)
                else:
                    # d1 and d2 have just been de-hybridized;
                    # Find new possible partners for d1 and d2:
                    # d1 and d2 are processed together with the other changed domains
                    pass
                del d1
                del d2
            else:
                raise ValueError("Unknown reaction type", reaction_attr.reaction_type)
            # end processing reacted_pair

        ## PROCESS CHANGED DOMAINS (all domains in all affected complexes):
        for domain in changed_domains:
            # IMPORTANT: changed_domains must not contain any... what?
            # TODO: Checking old vs new domspec is only needed if reacted_pair is not None
            # old_domspec = old_domspecs[domain]
            # For d1, d2, old_domspec yields the new domspec instead!
            # printd("Re-setting and re-calculating state_fingerprint for %s - old is: %s"
            #       % (domain, domain._specie_state_fingerprint))
            # domain.state_change_reset() # TODO: Reduce complex state reset and checking.
            # new_domspec = domain.state_fingerprint()
            ## TODO: Re-enable domspec change check
            # if new_domspec == old_domspec and reacted_pair is not None:
            #     c = domain.strand.complex
            #     print("\nWeird: new_domspec == old_domspec for changed domain %s: %s == %s" %
            #           (domain, new_domspec, old_domspec))
            #     print(" - complex._state_fingerprint =", c._state_fingerprint)
            #     print(" - complex.strands_fingerprint() =", c.strands_fingerprint())
            #     print(" - complex.hybridization_fingerprint() =", c.hybridization_fingerprint())
            #     print(" - complex.stacking_fingerprint() =", c.stacking_fingerprint())
            #     print(" - complex.reset_state_fingerprint() ...")
            #     print(" - complex.state_fingerprint() =", c.state_fingerprint())

            if domain.partner is None:
                ## No partner (but is has potential partners at the specie level specified in domain_pairs)
                ## consider hybridization reactions with all other unhybridized complementary domains.
                ## Or maybe just re-calculate for complementary domains in changed domains?
                # for d2 in [d for cname in self.domain_pairs[d1.name]
                #            for d in self.unhybridized_domains_by_name[cname]]:
                #pdb.set_trace()
                if domain.name in self.domain_pairs:
                    # is_forming = True
                    for cname in self.domain_pairs[domain.name]:
                        for domain2 in self.unhybridized_domains_by_name[cname]:
                            # Remove old reaction propensity:
                            # TODO: Consider only updating/removing intra-complex hybridization reactions:
                            # Inter-complex reactions are only affected if steric/electrostatic repulsion is included.
                            if (pair_against_all or domain2 in changed_domains) and domain is not domain2:
                                domain_pair = frozenset((domain, domain2))
                                if domain_pair in updated_reactions:
                                    continue
                                domain_spec_pair = frozenset((domain.state_fingerprint(), domain2.state_fingerprint()))
                                # is_intra: intra-complex OR intra-strand reaction:
                                is_intra = (domain.strand.complex is not None and \
                                            domain.strand.complex == domain2.strand.complex) \
                                            or (domain.strand == domain2.strand)
                                ra = ReactionAttrs(
                                    reaction_type=HYBRIDIZATION_INTERACTION, is_forming=True, is_intra=is_intra)
                                assert ra == (RA_HYB_INTRA if is_intra else RA_HYB_INTER)
                                c_j = self.calculate_c_j(domain, domain2, reaction_attr=ra,
                                                         reaction_spec_pair=domain_spec_pair)
                                if c_j > 0:
                                    self.possible_hybridization_reactions[domain_pair] = c_j
                                    self.reaction_attrs[domain_pair] = ra
                                    self.reaction_spec_pairs[domain_pair] = domain_spec_pair
                                    self.reaction_pairs_by_domain[domain].add(domain_pair)
                                    self.reaction_pairs_by_domain[domain2].add(domain_pair)
                                #if domain_pair not in updated_reactions:
                                updated_reactions.add(domain_pair)
            ## Domain has partner:
            elif recalculate_hybridized or domain in recalculate_these:
                # Currently, we don't re-calculate hybridized domains that are already hybridized (other than d1, d2)
                # - unless explicitly told to via recalculate_hybridized argument.
                # If we have a particular pair of domains that we would like to re-calculate (e.g. because of recent
                # stacking reaction), then these can be specified with recalculate_these.
                # EDIT: If throttling is enabled, then we MUST re-calculate
                # d2_domspec = d2.state_fingerprint()
                domain2 = domain.partner
                domain_pair = frozenset((domain, domain2))
                if domain.end5p.stack_partner is not None or domain.end3p.stack_partner is not None:
                    # Domains cannot de-hybridize if they are stacked:
                    assert domain_pair not in self.possible_hybridization_reactions
                    continue
                else:
                    domain_spec_pair = frozenset((domain.state_fingerprint(), domain2.state_fingerprint()))
                    # Only calculate domain_spec_pair if needed, in calculate_c_j.
                    # Edit: We always need it to get throttle_factor.
                    if domain_pair in updated_reactions:
                        continue
                    ra = ReactionAttrs(
                        reaction_type=HYBRIDIZATION_INTERACTION, is_forming=False, is_intra=True)
                    assert ra == RA_DEHYB_INTRA
                    c_j = self.calculate_c_j(domain, domain2, reaction_attr=ra, reaction_spec_pair=domain_spec_pair)
                    assert c_j > 0
                    self.possible_hybridization_reactions[domain_pair] = c_j
                    self.reaction_attrs[domain_pair] = ra
                    self.reaction_spec_pairs[domain_pair] = None # Calculate as-needed.
                    self.reaction_pairs_by_domain[domain].add(domain_pair)
                    self.reaction_pairs_by_domain[domain2].add(domain_pair)

                    updated_reactions.add(domain_pair)

        # end for domain in changed_domains


    def update_possible_stacking_reactions(self, changed_domains=None, reacted_pair=None, reaction_attr=None,
                                           # Edit/new: domains can only dehybridize if already unstacked;
                                           # use reacted_pair=(domain1, domain2) to detect obsolete reactions.
                                           #dehybridized_ends=None,
                                           recalculate_stacked=False,
                                           reacted_spec_pair=None):
        """
        The stacking equivalent to update_possible_hybridization_reactions.
        If changed_domains is None, consider ALL domains.
        (This makes it equivalent to init_possible_stacking_reactions)
        :reacted_pair: If provided, should be a frozenset of {(h1end3p, h2end5p), (h2end3p, h1end5p)}.
        Only provide reacted_pair when processing a recent stacking/unstacking reaction.
        DO NOT provide reacted_pair when processing a hybridization/dehybridization reaction.
        (The hybridized/dehybridized domains are listed in "changed_domains" and processed accordingly.)

        :dehybridized_duplex_ends: is a list of ends [d1end3p, d2end5p, ...]
            for domains that have just been dehybridized (and thus cannot undergo stacking reactions).
        """
        # printd("Updating possible stacking reactions for domains: (None=All)")
        # pprintd(changed_domains)
        # printd("Reacted pairs:", reacted_pairs)
        # printd("Reaction attrs:", reaction_attrs)
        if changed_domains is None:
            # Update possible stacking reactions for *all* domains:
            changed_domains = list(self.domains)
        else:
            assert reacted_pair is not None

        ## NOTE: The reaction could have been either stacking/unstacking OR HYBRIDIZATION/DEHYBRIDIZATION.

        ## First, determine reactions that are obsolete by stacking/un-stacking:
        if reacted_pair is not None:
            # reacted_pair MUST be duplex ends for a recent stacking/unstacking reaction; it should not be
            # given when updating after hybridization/dehybridization reactions, except if a de-hybridization
            # resulted in broken stacking interactions:
            if reaction_attr is None:
                reaction_attr = self.reaction_attrs[reacted_pair]
            # We are only considering newly formed stacking pairs; newly broken stacking pairs
            # should be in changed_domains and treated together with the other changed domains:
            # Stacking interaction between duplex ends was broken and is available for new rections.
            # This is treated together with the changed domains below:
            if reaction_attr.reaction_type is STACKING_INTERACTION:
                # reaction_pair grouped as {duplexend1, duplexend2} = {(h1end3p, h2end5p), (h2end3p, h1end5p)}
                (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(reacted_pair) # for stacking/unstacking reaction
                reacted_ends = (h1end3p, h2end5p, h2end3p, h1end5p)
                # We just stacked duplexend1:duplexend2; Make sure it is actually stacked...
                # assert reaction_attr.reaction_type == STACKING_INTERACTION
                assert all(isinstance(tup, tuple) for tup in reacted_pair)
                # Note: We don't process for if reaction_attr.reaction_type == HYBRIDIZATION_INTERACTION:
                # This case, where we have formed a new duplex, is processed under "for h1end3p in stacked_ends" below.
                ## Note: Being able to check len(reacted_pair & other_pair) > 0 is one reason for having stacking_pairs
                ## as {(h1end3p, h2end5p), (h2end3p, h1end5p)} instead of {(h1end3p, h1end5p), (h2end3p, h2end5p)}
                ## i.e. we have pairs of {duplexend1, duplexend2}.
                if reaction_attr.is_forming:
                    # We have just performed stacking reaction, so the ends should be stacked by now:
                    assert all(end.stack_partner is not None for end in reacted_ends)
                    ## We just formed a new stack, remove other stacking reactions involving these duplex ends:
                    obsolete_reactions = [other_pair for other_pair in self.possible_stacking_reactions
                                          if len(reacted_pair & other_pair) > 0]
                    # "if len(reacted_pair & other_pair) > 0" checks if the two stacking_pairs have any ends in common.
                    # set.intersection(other_set) is slower than '&' for small sets:
                    for reaction_pair in obsolete_reactions:
                        # printd("Removing from possible_stacking_reactions obsolete stacking_pair:", reaction_pair)
                        del self.possible_stacking_reactions[reaction_pair]
                        del self.reaction_attrs[reaction_pair]
                    # check:
                    for pair, ra in self.reaction_attrs.items():
                        if not ra.reaction_type is STACKING_INTERACTION:
                            continue
                        assert pair in self.possible_stacking_reactions
                        if ra.is_forming:
                            assert all(end.stack_partner is None for tup in tuple(pair) for end in tup)  # Error in check_system()
                        else:
                            assert all(end.stack_partner is not None for tup in tuple(pair) for end in tup)  # Error in check_system()
                else:
                    # We have just performed un-stacking reaction, so the ends should be stacked by now:
                    assert all(end.stack_partner is None for end in reacted_ends)
                ## TODO: Check to make sure that all reacted_pairs domains are in changed_domains
                # {(h1end3p, h1end5p), (h2end3p, h2nd5p)}
                for end in reacted_ends: # (end for tup in reacted_pair for end in tup):
                    if not end.domain in changed_domains:
                        # printd(" - update_possible_stacking_reactions: Adding %s's domain to changed_domains." % end)
                        #changed_domains.add(end.domain)
                        changed_domains.append(end.domain)
            else:
                ## We just hybridized/dehybridized domains in reacted_pair and are now updating stacking_reactions
                assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
                d1, d2 = tuple(reacted_pair)
                assert isinstance(d1, Domain) and isinstance(d2, Domain)
                # changed_domains.append(end.domain)
                assert d1 in changed_domains and d2 in changed_domains
                if reaction_attr.is_forming:
                    # We just formed a new duplex (hybridization reaction),
                    # Detecting new stacking reactions should be done below when running through changed_domains below
                    # But we might want to check that the now-hybridized domains are actually hybridized:
                    assert all(d1.partner is not None for d in (d1, d2))
                else:
                    # We just eliminated a duplex (de-hybridization reaction).
                    # Assert that the domains are actually de-hybridized:
                    assert all(d1.partner is None for d in (d1, d2))
                    # remove forming stacking reactions involving any of the de-hybridized domains:
                    dehybridized_duplex_ends = {(d1.end3p, d2.end5p), (d2.end3p, d1.end5p)}
                    obsolete_reactions = [other_pair for other_pair in self.possible_stacking_reactions
                                          if len(dehybridized_duplex_ends & other_pair) > 0]
                    for stacking_pair in obsolete_reactions:
                        del self.possible_stacking_reactions[stacking_pair]
                        del self.reaction_attrs[stacking_pair]


        # if dehybridized_ends is not None:
        #     obsolete_reactions = [stacking_pair for stacking_pair in self.possible_stacking_reactions
        #                           if any(end in ends_tup for end in dehybridized_ends
        #                                  for ends_tup in stacking_pair)]
            # printd("Removing reactions obsolete by dehybridized_duplex_ends:", dehybridized_ends)
            # pprintd(obsolete_reactions)
            ## Edit/new: domains cannot dehybridize if stacked; len(obsolete_reactions) should be 0
            # assert len(obsolete_reactions) == 0
            # for stacking_pair in obsolete_reactions:
            #     del self.possible_stacking_reactions[stacking_pair]


        # Only consider hybridized domains:
        hybridized_domains = [domain for domain in changed_domains if domain.partner is not None]
        stacked_ends = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is not None]
        # printd("stacked_ends:")
        # pprintd(stacked_ends)
        unstacked_duplex_ends3p = [domain.end3p for domain in hybridized_domains if domain.end3p.stack_partner is None]
        if self.enable_intercomplex_stacking:
            # Match against *all* unstacked domains (not just intra-complex)
            all_unstacked_duplex_ends3p = [domain.end3p for domain in self.domains
                                           if domain.partner is not None and domain.end3p.stack_partner is None]
        else:
            # Only consider stacking within changed complexes:
            all_unstacked_duplex_ends3p = unstacked_duplex_ends3p
        #unstacked_ends5p = [domain.end5p for domain in hybridized_domains if domain.end5p.stack_partner is None]
        # for domain in hybridized_domains:
        #     if domain.end3p.stack_partner is None:
        updated_reactions = set()
        # stacked_ends is a list of stacked DomainEnd3p on hybridized domains in changed_domains:
        for h1end3p in stacked_ends:
            ## TODO/Consideration: Will unstacking (off) rates ever change? Is there any point to re-calculating these?
            h1end5p = h1end3p.stack_partner
            h2end5p = h1end3p.hyb_partner
            h2end3p = h2end5p.stack_partner
            ## Assertions:
            assert h2end3p == h1end5p.hyb_partner
            # Multi-graph NetworkX interface is graph[node1][node2][key] => edge_attrs
            assert STACKING_INTERACTION in self.ends5p3p_graph[h1end5p][h1end3p]
            assert STACKING_INTERACTION in self.ends5p3p_graph[h1end3p][h1end5p] # Should be a given.
            assert STACKING_INTERACTION in self.ends5p3p_graph[h2end5p][h2end3p]
            assert STACKING_INTERACTION in self.ends5p3p_graph[h2end3p][h2end5p]

            # if h1end3p.stack_partner == h1end3p.pb_downstream and h2end3p.stack_partner == h2end3p.pb_downstream:
            #     # Performance hack for multiple domains on the same strand on the same duplex. Don't unstack.
            #     # Edit/new: All reactions must be reversible, even this. (No more cyclic reaction graphs).
            #     continue

            dupend1, dupend2 = (h1end3p, h2end5p), (h2end3p, h1end5p)  # duplex ends
            stacking_pair = frozenset((dupend1, dupend2))
            stacking_spec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                            (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
            if stacking_pair in updated_reactions:
                continue
            # For currently-stacked ends, we only need to update if its not in possible_stacking_reactions;
            # Once a stacking_pair is in possible_stacking_reactions, it doesn't change except when unstacked.
            # What if stacking_pair is already in possible_stacking_reactions, but the complex has changed
            # enough to make the first c_j obsolete? Uhm, see above...
            # And for the case looking at forming stacking reaction: the complex has changed,
            # so we *MUST* to re-calculate c_j.
            # stackspec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
            #                             (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
            # if stacking_pair in self.reaction_attrs:
            #     if stackspec_pair == self.reaction_attrs[stacking_pair].stackspec_pair:
            #         # No change since last:
            #         continue
            if stacking_pair not in self.possible_stacking_reactions or \
                (recalculate_stacked and stacking_pair not in updated_reactions):
                ## TODO: Reaction throttle should probably be inserted here. Edit: Use calculate_c_j() to add throttle.
                ra = ReactionAttrs(reaction_type=STACKING_INTERACTION, is_forming=False, is_intra=True)
                # Use calculate_c_j to add caching and throttling:
                c_j = self.calculate_c_j(dupend1, dupend2, reaction_attr=ra, reaction_spec_pair=stacking_spec_pair)
                # self.possible_stacking_reactions[stacking_pair] = \
                #     self.calculate_stacking_c_j(h1end3p, h2end5p, h2end3p, h1end5p, is_forming=False, is_intra=True)
                if c_j > 0:
                    self.possible_stacking_reactions[stacking_pair] = c_j
                    self.reaction_attrs[stacking_pair] = ra
                updated_reactions.add(stacking_pair)
        # unstacked_ends is a list of unstacked DomainEnd3p on hybridized domains *in changed_domains*:
        # all_unstacked_ends is a list of unstacked DomainEnd3p on hybridized domains for *all* duplex ends:
        for h1end3p in unstacked_duplex_ends3p:
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            # reaction_pair grouped as {duplexend1, duplexend2} = {(h1end3p, h2end5p), (h2end3p, h1end5p)}
            for h2end3p in all_unstacked_duplex_ends3p:
                if h2end3p is h1end3p:
                    continue # Ends cannot stack to them self

                # unpack variables:
                h2end5p = h1end3p.hyb_partner
                assert h1end3p.domain.partner == h2end5p.domain
                h1s1 = h1end3p.domain.strand
                h1c1 = h1s1.complex
                h1end5p = h2end3p.hyb_partner
                h1s2 = h1end5p.domain.strand
                h1c2 = h1s2.complex
                is_intra = h1s1 is h1s2 or (h1c1 is not None and h1c1 is h1c2)
                if not self.enable_intercomplex_stacking and (h1c1 != h1c2 if h1c1 is not None else h1s1 != h1s2):
                    # Interaction between different complexes (or strands):
                    assert is_intra is False
                    continue
                if h2end3p.hyb_partner == h1end3p.stacked_upstream() and not self.enable_helix_bending:
                    # Don't allow duplex to hybridize to itself (opposite ends).
                    #                 h1end5p    h1end3p
                    # This end --> 5'--------------------3' <-- cannot hybridize to this.
                    # Helix 2   -> 3'--------------------5' <-
                    #                h2end3p    h2end5p
                    # printd("update_possible_stacking_reactions: Ignoring stacking  %r vs %r" % (h1end3p, h2end3p))
                    continue
                ## Assertions:
                assert all(end is not None for end in (h1end5p, h1end3p, h2end3p, h2end5p))
                if h1end3p in self.ends5p3p_graph[h1end5p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h1end5p][h1end3p]
                if h1end5p in self.ends5p3p_graph[h1end3p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h1end3p][h1end5p]
                if h2end3p in self.ends5p3p_graph[h2end5p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h2end5p][h2end3p]
                if h2end5p in self.ends5p3p_graph[h2end3p]:
                    assert STACKING_INTERACTION not in self.ends5p3p_graph[h2end3p][h2end5p]

                dupend1, dupend2 = (h1end3p, h2end5p), (h2end3p, h1end5p)  # duplex ends
                stacking_pair = frozenset((dupend1, dupend2))
                stacking_spec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                                (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))

                if stacking_pair not in updated_reactions: # self.possible_stacking_reactions:
                    # print("  %s -: :- %s \n  %s -: :- %s  " % (h1end3p, h1end5p, h2end5p, h2end3p))
                    # Use calculate_c_j to add caching and throttling:
                    assert all(end.stack_partner is None for end in (h1end5p, h1end3p, h2end3p, h2end5p))
                    ra = ReactionAttrs(reaction_type=STACKING_INTERACTION, is_forming=True, is_intra=is_intra)
                    c_j = self.calculate_c_j(dupend1, dupend2, reaction_attr=ra,
                                             reaction_spec_pair=stacking_spec_pair)
                    if c_j > 0:
                        self.possible_stacking_reactions[stacking_pair] = c_j
                        self.reaction_attrs[stacking_pair] = ra

                    updated_reactions.add(stacking_pair)
        # pprint(self.reaction_attrs)
        # pprint(locals())
        # pdb.set_trace()


    def hybridize_and_process(self, domain_pair, reaction_attr, reaction_spec_pair=None):
        """
        Will select a random pair of domain instances from the domain species pair
        for hybridization reaction, or a random duplex in case of dehybridization reactions.
        Reaction specie consists of:
            ({domspec1, domspec2}, is_forming, is_intracomplex)
        """
        ## TODO: Merge this with stack_and_process
        # printd("\nreact_and_process invoked with args: domain_pair = %s, reaction_attr = %s" % (domain_pair, reaction_attr))
        # printd("domain domspecs/fingerprints (before reaction):", [d.state_fingerprint() for d in domain_pair])
        if self.invoked_reactions_file:
            print("domain_pair = frozenset((domains_by_duid[%s], domains_by_duid[%s]))" %
                  tuple([d.duid for d in domain_pair]),
                  file=self.invoked_reactions_file)
            print("sysmgr.hybridize_and_process(domain_pair, %s)" % (reaction_attr, ),
                  file=self.invoked_reactions_file)


        d1, d2 = tuple(domain_pair)
        if reaction_spec_pair is None:
            ## TODO: Add a global reaction_spec_pair by reaction_pair cache (making a frozenset is a bit slow).
            reaction_spec_pair = frozenset((d1.state_fingerprint(), d2.state_fingerprint()))
        else: ## TODO: Remove assertion when done debugging.
            assert reaction_spec_pair == frozenset((d1.state_fingerprint(), d2.state_fingerprint()))


        assert reaction_attr.reaction_type == HYBRIDIZATION_INTERACTION
        if reaction_attr.is_intra:
            assert d1.strand.complex == d2.strand.complex
            if d1.strand.complex is None:
                assert d1.strand == d2.strand
        else:
            assert (d1.strand.complex != d2.strand.complex) or (d1.strand.complex is None or d2.strand.complex is None)

        if reaction_attr.is_forming:
            # printd("Hybridizing domain %s %s and %s %s..." % (repr(d1), d1._specie_state_fingerprint, repr(d2), d2._specie_state_fingerprint))
            assert d1.partner is None and d2.partner is None
            result = self.hybridize(d1, d2)
            # printd("Completed hybridization of domain %s and %s..." % (d1, d2))
        else:
            # printd("De-hybridizing domain %s %s and %s %s..." % (d1, d1._specie_state_fingerprint, d2, d2._specie_state_fingerprint))
            assert d1.partner == d2 and d2.partner == d1
            result = self.dehybridize(d1, d2)
            # printd("Completed de-hybridization of domain %s and %s..." % (d1, d2))

        # 4: Update/re-calculate possible_hybridization_reactions and hybridization_propensity_functions
        # - domain_state_subspecies  - this is basically x̄ ← x̄ + νj
        # - possible_hybridization_reactions
        # - hybridization_propensity_functions
        # Note: If evaluating whether an object is boolean False, the steps include:
        # - Does it have a __len__ attribute? - Yes? Return bool(len(obj))
        # Whereas evaluating whether "obj is None" is about 10 times faster.
        # Fixed: Excluding domains that have no partners -- these will never undergo hybridization reaction.

        # printd("result: (is_forming: %s)" % is_forming)
        # pprintd(result)
        # if result['free_strands']:
            # printd("free strands domains:")
            # pprintd([s.domains for s in result['free_strands']])

        changed_domains = []
        if result['changed_complexes']:
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()
                                # if domain.partner is None and domain.name in self.domain_pairs
                                # edit: we need to have all, not just unhybridized - for stacking
                               ]
            # printd("Changed complexes domains:")
            # pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['new_complexes']:
            # self.complexes |= set(result['new_complexes'])
            # printd("Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            # pprintd(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()
                                 # if domain.partner is None and domain.name in self.domain_pairs
                                 # edit: we need to have all, not just unhybridized - for stacking
                                ]
            changed_domains += new_cmplx_domains
            # printd("New complexes domains:")
            # pprintd(new_cmplx_domains)
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains
                               #if domain.partner is None and domain.name in self.domain_pairs
                               # edit: we need to have all, not just unhybridized - for stacking
                              ]
            changed_domains += free_st_domains
        # if result['obsolete_complexes']:
            # self.complexes -= set(result['obsolete_complexes'])
            # self.removed_complexes += result['obsolete_complexes']
            # printd("Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            # pprintd(self.complexes)
            # printd("sysmgr.removed_complexes:")
            # pprintd(self.removed_complexes)

        # if reaction_attr.is_forming:
        #     # d1, d2 are partnered and not included in changed_domains (which currently excludes hybridized domains)
        #     # edit: we now include *all* domains in changed_domains, hybridized and un-hybridized.
        #     changed_domains += [d1, d2]

        if len(changed_domains) != len(set(changed_domains)):
            print("\nWARNING: changed_domains contains duplicates!! THIS SHOULD NOT HAPPEN!\n")
            print("changed_domains:")
            pprint(changed_domains)
            print("result['changed_complexes']:")
            pprint(result['changed_complexes'])
            print("changed_complexes domains: (using complex.nodes())")
            pprint([domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()])
            print("changed_complexes domains: (using complex.strands)")
            pprint([domain for cmplx in result['changed_complexes']
                    for s in cmplx.strands for domain in s.domains])
            print("free_strands:")
            pprint(result['free_strands'])
            print("free_strands domains:")
            pprint([d for s in result['free_strands'] for d in s.domains])

            print("( %s reaction between %s (%s) and %s (%s) )" %
                  ("Hybridization" if reaction_attr.is_forming else "De-hybridization", d1, repr(d1), d2, repr(d2)))
            print("-----------------")


        # Invoke post_reaction_processing after reaction, before updating reactions/propensity constants/functions:
        # Note that we give reaction_spec_pair, not reaction_pair/stacking_pair:
        self.post_reaction_processing(reaction_spec_pair, reaction_attr,
                                      reacted_pair=domain_pair, reaction_result=result)


        # DEBUGGING: Resetting complex fingerprint. Should be done in post_reaction_processing...
        # TODO: Move this hybridization/dehybridization methods and apply conditionally.
        # Also, doing this here would be too late: We've already updated possible reactions, so they would be all wrong.
        # Also, remember to reset domain state fingerprint explicitly if there is no complex!!
        if d1.strand.complex is None:
            d1.state_change_reset()
        else:
            d1.strand.complex.reset_state_fingerprint()
        if d2.strand.complex is None:
            d2.state_change_reset()
        elif d2.strand.complex is not d1.strand.complex:
            d2.strand.complex.reset_state_fingerprint()


        ## Update reactions for changed domains:
        self.update_possible_hybridization_reactions(changed_domains,
                                                     reacted_pair=domain_pair,
                                                     reaction_attr=reaction_attr,
                                                     reacted_spec_pair=reaction_spec_pair)
        # Only provide reacted_pairs if a reaction really has occured:
        # HEY, IF UNHYBRIDIZE() SETS d1.partner = None, then that currently also sets d1.end3p.hyb_partner = None, etc!
        # Then, if update_possible_stacking_reactions tries to "assert h2end5p == h1end3p.hyb_partner", it will fail.
        # Actually, it won't fail until it tries to stack the ends that have just been de-hybridized.
        # if 'unstacking_results' in result:
        #     reacted_stacking_pairs = result['unstacking_results']
        #     # unstacking_results is a dict of: {(h1end3p, h2end5p, h2end3p, h1end5p): result), ....}
        # else:
        #     reacted_stacking_pairs = None
        # Using dehybridized_duplex_ends is much better than trying to pretend that a (stacking) reaction has occoured:
        #dehybridized_ends = [d1.end5p, d1.end3p, d2.end5p, d2.end3p] if not reaction_attr.is_forming else None
        # Edit/new: domains can only dehybridize if they are are unstacked
        # (and if two domains hybridize, they should currently be unstacked; only duplexes can stack---for now.)
        assert all(end.stack_partner is None for end in (d1.end5p, d1.end3p, d2.end5p, d2.end3p))
        ## TODO: Just pass in reacted_pair as usual instead of having a special dehybridized_ends argument.
        self.update_possible_stacking_reactions(changed_domains=changed_domains,
                                                reacted_pair=domain_pair,
                                                reaction_attr=reaction_attr,
                                                reacted_spec_pair=reaction_spec_pair)
        # Note: We still use 'unstacking_results' for forwarding these side-effect unstacking reactions to dispatcher.

        # printd("changed_domains _specie_state_fingerprint:")
        # pprintd([(d, d._specie_state_fingerprint) for d in changed_domains])
        # printd("changed_domains specie_state_fingerprint():")
        # pprintd([(d, d.state_fingerprint()) for d in changed_domains])

        assert set(self.hybridization_propensity_functions.keys()) == set(self.possible_hybridization_reactions.keys())

        if self.invoked_reactions_file:
            print("print('- hybridize_and_process complete. (execfile)')", file=self.invoked_reactions_file)
        # Return the selected domains that were actually hybridized/dehybridized
        return (d1, d2), result



    def stack_and_process(self, stacking_pair, reaction_spec_pair=None):
        """
        The stacking equivalent to hybridize_and_process.
        Perform a stacking reaction, determine which domains have changed,
        and forward that info to reaction updater method.
        stacking_pair = frozenset(((h1end3p, h2end5p), (h2end3p, h1end5p)))
        Ends annotation:
                    h1end3p         h1end5p
        Helix 1   ----------3' : 5'----------
        Helix 2   ----------5' : 3'----------
                    h2end5p         h2end3p
        Returns stacking_pair, result
        """
        ## TODO: Consolidate stacking and hybridization

        # printd("stack_and_process invoked with stacking_pair:")
        # pprintd(stacking_pair)
        (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(stacking_pair)
        if self.invoked_reactions_file:
            print(("stacking_pair = frozenset(("
                   "(getattr(domains_by_duid[%s], 'end%s'), getattr(domains_by_duid[%s], 'end%s'))"
                   "(getattr(domains_by_duid[%s], 'end%s'), getattr(domains_by_duid[%s], 'end%s'))))") %
                  (h1end3p.domain.duid, h1end3p.end,
                   h1end5p.domain.duid, h1end5p.end,
                   h2end3p.domain.duid, h2end3p.end,
                   h2end5p.domain.duid, h2end5p.end),
                  file=self.invoked_reactions_file)
            print("sysmgr.stack_and_process(stacking_pair)", file=self.invoked_reactions_file)

        reaction_attr = self.reaction_attrs[stacking_pair]
        if reaction_spec_pair is None:
            reaction_spec_pair = frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                            (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        else:
            assert reaction_spec_pair == frozenset(((h1end3p.state_fingerprint(), h2end5p.state_fingerprint()),
                                                    (h2end3p.state_fingerprint(), h1end5p.state_fingerprint())))
        # reaction_str = ("" if reaction_attr.is_forming else "de-") + reaction_attr.reaction_type


        # d1 = h1end3p.domain
        # d2 = h2end3p.domain
        # Assertions:
        assert h2end5p == h1end3p.hyb_partner
        assert h2end3p == h1end5p.hyb_partner
        if reaction_attr.is_forming:
            # assert d1.partner is None and d2.partner is None
            assert h1end5p.stack_partner is None and  h1end3p.stack_partner is None
            assert h2end3p.stack_partner is None and  h2end5p.stack_partner is None
        else:
            assert h1end5p == h1end3p.stack_partner
            assert h2end3p == h2end5p.stack_partner
            # assert d1.partner == d2 and d2.partner == d1
        # printd("stack_and_process: %s h1end3p, h2end5p, h2end3p, h1end5p - %s %s and %s %s..." % (reaction_str, h1end3p, h2end5p, h2end3p, h1end5p))
        assert reaction_attr.reaction_type == STACKING_INTERACTION
        if reaction_attr.is_intra:
            assert len(set(e.domain.strand.complex for e in (h1end3p, h2end5p, h2end3p, h1end5p))) == 1
            if h1end3p.domain.strand.complex is None:
                assert len(set(e.domain.strand for e in (h1end3p, h2end5p, h2end3p, h1end5p))) == 1

        if reaction_attr.is_forming:
            assert all(e.stack_partner is None for e in (h1end3p, h2end5p, h2end3p, h1end5p))
            result = self.stack(h1end3p, h2end5p, h2end3p, h1end5p)
        else:
            assert all(e.stack_partner is not None for e in (h1end3p, h2end5p, h2end3p, h1end5p))
            result = self.unstack(h1end3p, h2end5p, h2end3p, h1end5p)
        # printd("stack_and_process: Completed %s of h1end3p, h2end5p, h2end3p, h1end5p - %s %s and %s %s" % (reaction_str, h1end3p, h2end5p, h2end3p, h1end5p))

        # printd("stack_and_process: %s result:" % reaction_str)
        # pprintd(result)

        changed_domains = []
        if result['changed_complexes']:
            # Edit: For stacking, we want duplexed domains, regardless of domain_pairs!!
            # No longer filtering with "if domain.partner is None and domain.name in self.domain_pairs"
            # Instead, do appropriate filtering in update_possible_hybridization_reactions.
            ch_cmplx_domains = [domain for cmplx in result['changed_complexes'] for domain in cmplx.nodes()]
            # printd("stack_and_process: Changed complexes domains:")
            # pprintd(ch_cmplx_domains)
            changed_domains += ch_cmplx_domains
        if result['free_strands']:
            free_st_domains = [domain for strand in result['free_strands'] for domain in strand.domains]
            changed_domains += free_st_domains
        if result['new_complexes']:
            # self.complexes |= set(result['new_complexes'])
            # printd("stack_and_process: Adding new complexes %s to sysmgr.complexes:" % result['new_complexes'])
            # pprintd(self.complexes)
            new_cmplx_domains = [domain for cmplx in result['new_complexes'] for domain in cmplx.nodes()]
            changed_domains += new_cmplx_domains
            # printd("stack_and_process: New complexes domains:")
            # pprintd(new_cmplx_domains)
        # if result['obsolete_complexes']:
            # self.complexes -= set(result['obsolete_complexes'])
            # self.removed_complexes += result['obsolete_complexes']
            # printd("stack_and_process: Removing obsolete complexes %s from sysmgr.complexes:" % result['obsolete_complexes'])
            # pprintd(self.complexes)
            # printd("stack_and_process: sysmgr.removed_complexes:")
            # pprintd(self.removed_complexes)

        assert all(isinstance(domain, Domain) for domain in changed_domains)

        # Invoke post_reaction_processing after reaction, before updating reactions/propensity constants/functions:
        # Note that we give reaction_spec_pair, not reaction_pair/stacking_pair:
        self.post_reaction_processing(reaction_spec_pair, reaction_attr,
                                      reacted_pair=stacking_pair, reaction_result=result)

        # TODO: Consolidate with init_possible_reactions.. Also, consider renaming:
        #   init_possible_reactions -> update_reactions, and
        #   update_possible_hybridization_reactions -> update_hybridization_reactions
        # Update hybridization reactions:
        self.update_possible_hybridization_reactions(
            changed_domains, reacted_pair=stacking_pair, reaction_attr=reaction_attr,
            # edit/new: stacked domains cannot dehybridize; providing reacted_pair=stacking_pair should remove
            # the hybridized duplex domain pair from possible_hybridization_reactions
            # recalculate_these={e.domain for e in (h1end3p, h2end5p, h2end3p, h1end5p)}
            )
        # Update stacking reactions:
        self.update_possible_stacking_reactions(changed_domains, reacted_pair=stacking_pair,
                                                reaction_attr=reaction_attr)
        # DEBUGGING: Resetting complex fingerprint. THIS MAY NOT BE CALLED FOR STACKING REACTIONS,
        # BECAUSE ComponentMgr.join/break_complex_at IS NOT CALLED WHEN INTER-COMPLEX STACKING IS DISABLED.
        # However, we are resetting complex and free strand's domain state fingerprints in post_reaction_processing.
        # TODO: Move this to hybridization/dehybridization methods and apply conditionally.
        # if h1end3p.domain.strand.complex:
        #     h1end3p.domain.strand.complex.reset_state_fingerprint()
        # if h2end3p.domain.strand.complex:
        #     h2end3p.domain.strand.complex.reset_state_fingerprint()

        return stacking_pair, result


    def update_state_times(self, tau):
        """ Update tau_cum for all state nodes in reaction graph and add a time step. """
        # Q: Does it matter if we have one complex that stays in a state for 2 secs,
        # or two complexes staying for 1 s? -- No.
        # But, we cannot be certain that tau represents the full time that a changed complex have been in the state.
        # Actually, we should go over all complexes in self.complexes and increase tau... And we should
        # do this *before* asserting state change for new/changed complexes...

        # Unary operator '+' on Counters will return all elements with count > 0:
        for state, count in (+self.state_counter).items():
            # try:
            #     self.reaction_graph.node[state]['tau_cum'] += tau
            # except KeyError:
            #     self.reaction_graph.node[state]['tau_cum'] = tau
            tau_cum = self.reaction_graph.node[state].get('tau_cum', 0) + tau*count
            if tau_cum > self.reaction_graph.tau_cum_max:
                self.reaction_graph.tau_cum_max = tau_cum
            self.reaction_graph.change_node(state, tau_cum=tau_cum)
        self.reaction_graph.add_time_step(tau)
        ## TODO: This should also be done at the end of the simulation!
        ## Perhaps better to have a dedicated "update_state_times"?


    def post_reaction_processing(self, reacted_spec_pair, reaction_attr, reacted_pair, reaction_result): #, tau=0):
        """
        post_reaction_processing is invoked *just after* a reaction has been performed,
        but *before* updating possible reactions.

        :reacted_pair:  A frozenset of either {domain1, domain2} or {(h1end3p, h2end5p), (h2end3p, h1end5p)}.
        :reaction_attr: A namedtuple with (reaction_type, is_forming, is_intra).
        :reaction_result: A result dict with keys 'changed_complexes', 'new_complexes', 'obsolete_complexes',
                        'free_strands', 'case', etc.

        Note: reaction_result can be None, e.g. for stacking/unstackig where (if configured)
              we don't consider inter-complex reactions (so unstacking will never "result" in new complexes).
        ...but: We still use the results dict to determine *changed* complexes.

        Regarding reaction_result['case']:
            If case == -1 then no join/break_complex action was attempted.
            If case < 2 then we have an intra-molecular reaction.
        Hybridization:
            Case 0: intRA-STRAND hybridization.
            Case 1: IntRA-complex hybridization
            Case 2: IntER-complex hybridization between two complexes. Merge the two complexs:
            Case 3: One complex + one strand.
            Case 4: Neither strands are in existing complex; create new complex
        Breaking:
            CASE 0: NO COMPLEX PRESENT, intra-strand reaction.
            Case 1: The two strands are still connected: No need to do anything further
            Case 2: Two smaller complexes - must create a new complex for detached domain:
            Case 3: One complex and one unhybridized strand - no need to do much further
            Case 4: Two unhybridized strands
        Thus:
            Case 0-1: Inter/uni-molecular reaction; No change in volume energy.
            Case 2-5: Inter/bi-molecular reaction; Change in volume energy.

        """
        ## TODO: Throttle, consider alternatively have reaction_invocation_count on a per-domain basis.

        ## TODO: Use reaction graph to predict known end-state (fingerprints) for the complexes;
        #        and forward that info to assert_state_change.

        ## TODO: Implement a "checking" child-class (for both complex and reactionmgr);
        ##       The "checking" class has all the check to ensuring integrity,
        ##       The base-class can be used as-is for performance.

        edge_key = (reacted_spec_pair, reaction_attr)
        reaction_spec_tuple = tuple(reacted_spec_pair)
        self.reaction_invocation_count[edge_key] += 1
        throttle_factor = self.reaction_throttle_cache[edge_key][0] if edge_key in self.reaction_throttle_cache else 1.0
        reaction_attr_str = (reaction_attr.reaction_type + ("+" if reaction_attr.is_forming else "-")
                             + (" " if reaction_attr.is_intra else "*"))
        reaction_is_joining = reaction_attr.is_forming and not reaction_attr.is_intra
        reaction_is_splitting = not reaction_attr.is_forming and reaction_result['case'] > 1

        # We have not yet updated possible reactions, so we can get the old c_j value as:
        if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
            c_j = self.possible_hybridization_reactions[reacted_pair]
        elif reaction_attr.reaction_type is STACKING_INTERACTION:
            c_j = self.possible_stacking_reactions[reacted_pair]
        else:
            raise ValueError("Unknown reaction type %s" % reaction_attr.reaction_type, reaction_attr)
        # Perhaps compare with self.cache['stochastic_rate_constant']?

        elem1, elem2 = tuple(reacted_pair)
        if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
            domain1, domain2 = tuple(reacted_pair)
            # domain state fingerprint = (dspecie, self.partner is not None, c_state, in_complex_identifier)
            d1fp, d2fp = reaction_spec_tuple
            c1state, c2state = d1fp[2], d2fp[2]
            # Note: It *is* possible to have two separate with identical state merging (or splitting)!
            if reaction_is_joining:
                reaction_spec_source_states_list = [d_fp[2] for d_fp in reaction_spec_tuple] # cstate or strand.name
                reaction_spec_source_states = set(reaction_spec_source_states_list)
            else:
                reaction_spec_source_states = {d_fp[2] for d_fp in reaction_spec_tuple} # cstate or strand.name
                reaction_spec_source_states_list = list(reaction_spec_source_states)
                assert len(reaction_spec_source_states_list) == 1

        elif reaction_attr.reaction_type is STACKING_INTERACTION:
            #             h1end3p         h1end5p
            # Helix 1   ----------3' : 5'----------
            # Helix 2   ----------5' : 3'----------
            #             h2end5p         h2end3p
            ## Regarding order of DomainEnd stacking pairing/tuples (check that they are correct)
            # Stacking is downstream: 3'::5'. But we bundle by duplex ends, so it becomes symmetric, and then intead
            # we organize the ends for each dupel as (h1end3p, h2end5p), even though that may not make much sense.
            # We've had this discussion before about the order of DomainEnds in stacking reaction pairs,
            # ordering either by duplex-end, or by helix:
            #   {(h1end3p, h2end5p), (h2end3p, h1end5p)}   or   {(h1end3p, h1end5p), (h2end3p, h2end5p)}
            # Remember that stacking of duplex ends is a symmetric pair (end1 stack end2 is same as end2 stack end1),
            # so the order shouldn't matter when storing the pair (use frozenset, not tuple).
            # (h1end5p, h2end3p), (h2end5p, h1end3p) = tuple(reacted_pair)  # Wrong
            (h1end3p, h2end5p), (h2end3p, h1end5p) = tuple(reacted_pair)
            assert all(end.end == "5p" for end in (h1end5p, h2end5p))  ## TODO: Remove assertions
            assert all(end.end == "3p" for end in (h1end3p, h2end3p))
            if reaction_is_joining:
                # reaction_spec_tuple is frozenset of tuples of DomainEnd fingerprints:
                # DomainEnd fingerprints: (self.domain.state_fingerprint(), self.end, self.stack_partner is not None)
                reaction_spec_source_states_list = [tup[0][0][2] for tup in reaction_spec_tuple] # cstate
                reaction_spec_source_states = set(reaction_spec_source_states_list)
            else:
                reaction_spec_source_states = {e_fp[0][2] for tup in reaction_spec_tuple for e_fp in tup} # cstate
                reaction_spec_source_states_list = list(reaction_spec_source_states)
                assert len(reaction_spec_source_states_list) == 1
            #
        elif reaction_attr.reaction_type is PHOSPHATEBACKBONE_INTERACTION:
            h1end3p, h1end5p = tuple(reacted_pair)
            assert h1end3p.end == "3p" and h1end5p.end == "5p"  ## TODO: Remove assertion
            reaction_spec_source_states = {d_fp[2] for d_fp in reaction_spec_tuple} # cstate or strand.name
        else:
            raise ValueError("Unknown reaction type '%s'" % (reaction_attr.reaction_type,), reaction_attr.reaction_type)

        ## List of all changed complexes (incl new complexes):
        new_or_changed_complexes = [cmplx
                                    for lst in (reaction_result['changed_complexes'], reaction_result['new_complexes'])
                                    if lst is not None and len(lst) > 0
                                    for cmplx in lst]
        # new_complexes_set = set(reaction_result['new_complexes']) if reaction_result['new_complexes'] else set()


        ### 0. ASSERT STATE CHANGE FOR ALL CHANGED COMPLEXES: ###
        source_states = set()
        asserted_target_states = dict()
        asserted_target_states_list = []
        for cmplx in new_or_changed_complexes:
            source_state = cmplx._historic_fingerprints[-1]  # Can be 0 for new complexes!
            ## TODO: How about when merging, you make fingerprint = frozenset of the two merging complexes?
            assert source_state is 0 or source_state in reaction_spec_source_states  # 0 for new complexes
            source_states.add(source_state)
            expected_state_fingerprint = None
            expected_state_fingerprints = None
            assert source_state in self.endstates_by_reaction # The "before reaction" state should be present.
            if edge_key in self.endstates_by_reaction[source_state]:
                # If a reaction produces multiple new/changed complexes, then that reaction
                # effectively has two outgoing edges and
                #   endstates_by_reaction[source][edge_key] should list MULTIPLE targets.
                expected_state_fingerprints = self.endstates_by_reaction[source_state][edge_key]
                if len(expected_state_fingerprints) == 1:
                    expected_state_fingerprint = expected_state_fingerprints[0]
            ## Assert complex state change
            ## (will update complex state fingerprint - but not reset it first unless you pass reset=True):
            target_state, _ = cmplx.assert_state_change(
                reacted_spec_pair, reaction_attr, expected_state_fingerprint=expected_state_fingerprint)
            asserted_target_states[target_state] = cmplx
            asserted_target_states_list.append(target_state)
            if expected_state_fingerprints is not None:
                assert target_state in expected_state_fingerprints  # List of target states via edge with edge_key.

        ### 0(b): Reset domain state fingerprint and icid for all free strand's domains:
        ### (This is also done in componentmgr.join/break_complex_at, but better safe than sorry)
        if reaction_result['free_strands']:
            for strand in reaction_result['free_strands']:
                for domain in strand.domains:
                    domain.state_change_reset()
                asserted_target_states[strand.name] = strand
                asserted_target_states_list.append(strand.name)

        if reaction_is_splitting:
            assert len(asserted_target_states_list) == 2
        else:
            assert len(asserted_target_states_list) == 1

        ### 0(c) Update reactionmgr attributes:
        new_reaction_graph_nodes = asserted_target_states.keys() - self.reaction_graph.node
        self.state_counter.update(asserted_target_states_list)
        self.state_counter.subtract(reaction_spec_source_states_list)

        ## 1. Form or break loops: ##

        ## ERROR: At this point we have already performed the reaction and closed the loop by forming the bond.
        ## If we try to calculate ifnode fingerprints, they will no longer match the hashes in loop path_spec.
        ## Thus, we should either update the loop_effects dict *before* doing the reaction, injecting the actual
        ## ifnode instances and loopid, or we should register the new loop (and corresponding changes) BEFORE
        ## doing the actual reaction.
        ##
        # We do this here because we need the loop energies for updating complexes and reaction graph energies
        # By definition, if we are forming or breaking a loop, then only one complex is changed.
        # We get info about forming loops from the loop_effects dict obtained in intracomplex_activity
        if False and reaction_attr.is_forming:
            # See if we have any loop_effects dicts registered for this reaction:
            if reacted_spec_pair in self.reaction_loop_effects:
                loop_effects = self.reaction_loop_effects[reacted_spec_pair]
                ## TODO: (Optimization) Check if the old ifnodes path is still valid (if may well be...)
                ## Re-create ifnodes from
                assert len(reaction_result['changed_complexes']) == 1
                assert reaction_result['new_complexes'] is None
                assert reaction_result['free_strands'] is None
                cmplx = reaction_result['changed_complexes'][0]
                # The shortest-path loop should certainly be created:
                # TODO: Consolidate all this with a Complex method register_new_loop
                shortest_loop_path_spec = loop_effects['shortest_path_spec']  # primary loop path
                # loop_path_spec, loop_path, loop_activity, replacing_loop_spec=None
                cmplx.register_new_loop(
                    shortest_loop_path_spec,
                    loop_effects['shortest_path_activity'])

                # Then consider splitted-loops:
                # What about loops that are just adjusted? E.g. by stacking backbone-linked duplexes?
                # old_loop_spec => [replacement_path1, replacement_path2]
                for loop0_hash, replacement_loops in loop_effects['changed_loops'].items():
                    # loop0_hash => list of one or two new loops. (Should always just be 1, right?)
                    # A loop is split in either 1 or two new loops.
                    # if len(replacement_loops) == 1:
                    #     # Perhaps just adjust existing loop?
                    # else:
                    for replacement_loop in replacement_loops:
                        cmplx.update_changed_loop(loop0_hash, replacement_loop)

                # Consider other affected loops, not being split but still affected?
                # if loop_effects['loops_affected']:
                #     for loop0 in loop_effects['loops_affected']:
                #         # Should re-calculate loop0 energy:
                #         if new_loop_paths:
                            # Loop is not




        ### 2. Calculate state energy change used to update complex energies and reaction graph: ###

        ## Done: Fix issue with energies when merging or splitting complexes!
        ## (Fixed by tracking all energy contributions for a complex)
        ## TODO: If reaction edges already exist in reaction_graph, then use energies from that.
        # In theory we could get reaction energy using reaction_spec_pair (source state) and
        # edge_key = (reacted_spec_pair, reaction_attr).

        ## TODO: Completely re-work this so that instead of using reaction dH, dS to calculate complex energy,
        ## we calculate complex energy using individual contributions (hybridizations, stackings, loops, Nstrands),
        ## and then use reaction dH, dS to check the *change* in complex energy (before vs after).

        ### 2a: Calculate reaction energy dH, dS totals:
        dH, dS = 0, 0  # Total reaction enthalpy and entropy
        dH_hyb, dS_hyb, dH_stack, dS_stack, dS_shape, dS_volume = None, None, None, None, None, None # Starting values
        dHdS_contribs = {'hybridization': [0., 0.],
                         'stacking': [0., 0.],
                         'shape': [0., 0.],
                         'volume': [0., 0.],
                        }
        if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
            dH_hyb, dS_hyb = self.dHdS_from_state_cache(domain1, domain2, state_spec_pair=reacted_spec_pair)
            # The dHdS in the cache are energy/entropy of FORMATION. Negate values if we are dehybridizing:
            if not reaction_attr.is_forming:
                dH_hyb, dS_hyb = -dH_hyb, -dS_hyb
            dH += dH_hyb
            dS += dS_hyb
            dHdS_contribs['hybridization'] = [dH_hyb, dS_hyb]
        elif reaction_attr.reaction_type is STACKING_INTERACTION:
            # h1end3p.stack_string is only set when the DomainEnd is actually stacked.
            stack_string = "%s%s/%s%s" % (h1end3p.base, h1end5p.base, h2end5p.base, h2end3p.base)
            if stack_string not in DNA_NN4_R:
                stack_string = stack_string[::-1]
            dH_stack, dS_stack = DNA_NN4_R[stack_string]
            if not reaction_attr.is_forming:
                dH_stack, dS_stack = -dH_stack, -dS_stack
            dH += dH_stack
            dS += dS_stack
            dHdS_contribs['hybridization'] = [dH_stack, dS_stack]
        else:
            raise NotImplementedError("Only HYBRIDIZATION and STACKING interactions implemented ATM.")
        # Note: reaction_attr.is_intra is *always* true for dehybridize/unstack reactions;
        # use result['case'] to determine if the volume energy of the complex is changed.
        # reacted_spec_pair will only occour in self._statedependent_dH_dS when dehybridization_rate_constant
        # has been called. Also: Is this really state dependent? Can't we just let de-hybridization and unstacking
        # be independent of complex state?
        #dH, dS = self._statedependent_dH_dS[reacted_spec_pair]
        if reaction_attr.is_forming:
            # Case 0/1:   IntRA-STRAND/COMPLEX hybridization.
            # Case 2/3/4: IntER-complex hybridization between two complexes/strands.
            if reaction_result['case'] <= 1:
                # IntRA-complex reaction
                assert reaction_attr.is_intra
                activity = self.cache['intracomplex_activity'][reacted_spec_pair]
                if activity == 0:
                    print(("\n\nActivity %s for %s reaction between %s and %s is <= 0; "
                           "shape/loop energy is infinite; reaction should revert.\n\n") %
                          (activity, reaction_attr.reaction_type, elem1, elem2))
                    dS_shape = 100000
                else:
                    ## Activity is in implicit unit of Molar, so loop entropy is just -R*ln(activity)
                    ## However, the "-R*ln(activity)" change in entropy is when releasing the loop.
                    ## When we are forming the loop, the entropy is reduced; activity is << 1, so ln(activity) < 0.
                    dS_shape = ln(activity)
                dS += dS_shape
                dHdS_contribs['shape'][0] += dS_shape
            else:
                # IntER-strand/complex reaction; Two strands/complexes coming together.
                assert reaction_attr.is_intra is False
                # self.volume_entropy is the (positive) increase in entropy (in units of R) when a
                # volume restriction is lifted. Here it's negative because the restriction is formed.
                activity = self.specific_bimolecular_activity # set activity for later use
                dS_volume = -self.volume_entropy  # Minus because we are forming; Multiply by -R*T to get dG.
                dS += dS_volume
                dHdS_contribs['volume'][0] += dS_volume
        else:
            ## Dehybridize/unstack reaction: (is always "intra-molecular" with activity=1)
            assert reaction_attr.is_intra is True
            activity = 1  # Used for edge attrs

            ## TODO: Re-work this so that dS_shape is calculated using complex shape entropy difference:
            ##          dS_shape = S_shape_before - S_shape_after
            ##       (state entropies rather than reaction activity...)
            ## Note: This is only really needed for annotating reaction graph edges.

            if reaction_result['case'] <= 1:
                if reaction_attr.reaction_type is STACKING_INTERACTION:
                    reverse_reaction_spec_pair = frozenset(
                        ((elem1[0].state_fingerprint(), elem1[1].state_fingerprint()),
                         (elem2[0].state_fingerprint(), elem2[1].state_fingerprint())))
                else:
                    assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
                    reverse_reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
                reverse_activity = self.intracomplex_activity(elem1, elem2, reaction_attr.reaction_type,
                                                              reverse_reaction_spec_pair)
                if reverse_activity <= 0:
                    print(("\n\nActivity for %s reaction between %s and %s is <= 0, reaction_attr %s"
                           "shape/loop energy is infinite; reaction should revert.\n\n") %
                          (reverse_activity, elem1, elem2, reaction_attr))
                    dS_shape = 1000
                else:
                    # The loop restriction is lifted, so shape entropy should increase.
                    # activity is << 1, so -ln(activity) (with the minus) is positive.
                    dS_shape = -ln(reverse_activity)
                dS += dS_shape
                dHdS_contribs['shape'][0] += dS_shape
            else:
                # Case > 1: complex is broken up into two complexes/molecules.
                # Increase in entropy (in units of R) when a volume restriction is lifted.
                dS_volume = self.volume_entropy  # Multiply by -R*T to get dG.
                dS += dS_volume
                dHdS_contribs['volume'][0] += dS_volume
        assert (dS_shape is None) != (dS_volume is None)

        ### 2b. Check energies for new/changed complexes:
        # (I've temporarily separated this from the "update reaction graph for all new/changed complexes" to
        # improve readability and clarity...)
        ## TODO/EDIT: Complex energy should be updated by ComponentMgr; you just need to verify that:
        ##      dH, dS = sum(all complexes' dHdS after reaction) - sum(all complexes' dHdS before reaction)
        ## and optionally also for individual contributions (hyb, stack, shape, volume)..
        for cmplx in new_or_changed_complexes:
            ## TODO: This for-loop is pretty long, consider shortening it or moving code to dedicated methods.
            ## TODO: ^^ Perhaps just move energy updates to componentmgr.hybridize (etc) ?
            ### 3(a). Update complex energy: ###
            # cmplx.update_complex_energy(dH_hyb, dS_hyb, dH_stack, dS_stack, dS_shape, dS_volume,
            #                             reaction_attr=reaction_attr, reacted_pair=reacted_pair)
            # EDIT: Complex energies are updated by ComponentMgr when the reaction is done.
            cmplx.recalculate_complex_energy(self.volume_entropy, dS_shape)
            cmplx.check_complex(self.volume_entropy)
            assert cmplx._state_fingerprint == cmplx._historic_fingerprints[-1]
            assert cmplx._historic_fingerprints[-1] in asserted_target_states
            ## Make sure that cmplx.assert_state_change has been invoked before using _historic_fingerprints:
            #source, target = cmplx._historic_fingerprints[-2], cmplx._historic_fingerprints[-1]
            # Uh, source and target state was defined when looping over changed complexes to assert state changes.
            # If we are e.g. merging two complexes, then source_state will match src_state of the last complex
            # in changed_complexes.

            # 2c: Check if complex energy matches the expected complex energy and correct if not.
            if cmplx._state_fingerprint in self.reaction_graph.node:
                state_node = self.reaction_graph.node[cmplx._state_fingerprint]
                expected_dHdS = state_node['dHdS']
                expected_dHdS_contributions = state_node['dHdS_contributions']
                actual_dHdS = cmplx.energy_total_dHdS
                # if actual_dHdS != expected_dHdS:
                if not all(np.isclose(actual_dHdS, expected_dHdS)):
                    print("Complex %s energy %s does not equal expected energy %s from reaction_graph" %
                          (cmplx, actual_dHdS, expected_dHdS))
                    print("Contributions:")
                    print(" - Expected: %s" % expected_dHdS_contributions)
                    print(" - Current : %s" % cmplx.energies_dHdS)
                    print(" - Current energy contributions: %s" % cmplx.energy_contributions)
                    # Find recent reaction cycle (reaction graph cycle):
                    # if cmplx._state_fingerprint in cmplx._historic_fingerprints[::-1]:
                    # Note: list.index(elem) is slower than "elem in list", but worth it if elem is likely in list.
                    try:
                        # n_reactions_since_last_state_encounter = right index:
                        ridx = cmplx._historic_fingerprints[-2::-1].index(
                            cmplx._state_fingerprint) + 1
                    except ValueError:
                        # Complex has not previously been in this state..
                        pass
                    else:
                        # reaction, source, target:
                        reactions_and_states = zip_longest(
                            list(cmplx.reaction_deque)[-2:-ridx-1:-1],
                            cmplx._historic_fingerprints[-3:-ridx-2:-1],
                            cmplx._historic_fingerprints[-2:-ridx-1:-1])
                        reactions_and_states = list(reactions_and_states)[::-1]
                        # Note: cmplx._historic_fingerprints[0] == 0 is not in reaction graph!
                        # Also, no edge between cmplx._historic_fingerprints[-2] and cmplx._historic_fingerprints[-1]
                        # in reaction_graph yet!
                        print("Reactions since we were last in this state (not including this reaction):")
                        print("\n".join(
                            " %s %s -(%s %s)-> %s %s" % (
                                source, self.reaction_graph.node[source]['dHdS'],
                                self.reaction_graph.adj[source][target]['reaction_str'],
                                (self.reaction_graph.adj[source][target]['dH'], self.reaction_graph.adj[source][target]['dS']),
                                target, self.reaction_graph.node[target].get('dHdS'))
                            for reaction, source, target in reactions_and_states))
                        # Find simple cycle (no repeated nodes or edges):
                        state_idx = {} # state: idx
                        simple_reaction_cycle = [] # No repeated nodes or edges
                        for reaction, state, endstate in reactions_and_states:
                            if endstate in state_idx:
                                # We've previously been here, revert to simple_reaction_cycle
                                # as it was when we were here last.
                                simple_reaction_cycle = simple_reaction_cycle[:state_idx[endstate]+1]
                                # print(simple_reaction_cycle)
                            else:
                                simple_reaction_cycle.append((reaction, state, endstate))
                                state_idx[endstate] = len(simple_reaction_cycle) - 1
                                # print(simple_reaction_cycle)
                        reaction, source, target = simple_reaction_cycle[0]
                        src_attrs = self.reaction_graph.node[source]
                        print("%s %s" % (source, src_attrs['dHdS']))
                        for reaction, source, target in simple_reaction_cycle:
                            src_attrs = self.reaction_graph.node[source] if source is not None else {}
                            tgt_attrs = self.reaction_graph.node[target] if target is not None else {}
                            edge_attrs = (self.reaction_graph.adj[source][target]
                                          if source is not None and target is not None else
                                          {'reaction_str': "(no edge attrs)", 'dH': '?', 'dS': '?'})
                            # TODO: Rename dHdS for states to HS to make it easier to distinguish
                            # reaction energies from complex state energies.
                            print("  -(%s %s)-> %s %s" % (
                                edge_attrs['reaction_str'], (edge_attrs['dH'], edge_attrs['dS']),
                                target, tgt_attrs.get('dHdS')))
                    # Done printing reaction cycle info...
                    print("Resetting complex %s energy to %s." % (cmplx, expected_dHdS))
                    cmplx.energies_dHdS = expected_dHdS_contributions
                    cmplx.energy_total_dHdS = expected_dHdS
                # end if actual_dHdS != expected_dHdS
            # end if cmplx._state_fingerprint in self.reaction_graph.node
        # end for cmplx in new_or_changed_complexes






        ### 3. Update reaction graph: ###
        asserted_target_states_list = list(asserted_target_states)
        reaction_spec_source_states_list = list(reaction_spec_source_states)
        new_nodes_added, new_edges_added = 0, 0
        edge_attrs = {
            'edge_key': edge_key,
            'reaction_spec_pair': reacted_spec_pair,
            'reaction_spec_source_states': reaction_spec_source_states_list,
            'asserted_target_states': asserted_target_states_list,
            'reaction_attr': tuple(reaction_attr),
            'reaction_type': reaction_attr.reaction_type,
            'is_forming': reaction_attr.is_forming,
            'is_intra': reaction_attr.is_intra,
            'reaction_attr_str': reaction_attr_str,
            'reaction_spec_pair_str': reaction_spec_pair_to_str(reacted_spec_pair),
            'reaction_str': reaction_to_str(reacted_spec_pair, reaction_attr),
            'reaction_invocation_count': self.reaction_invocation_count[edge_key],
            #'dHdS': dHdS,
            'dH': dH, #dHdS[0],
            'dS': dS, #dHdS[1],   # includes loop/volume entropy
            'activity': activity,
            'c_j': c_j,
            'c_j_throttled': c_j,
            # 'len': 4.0,     # optimal edge length, default = 1.0
            # Above attrs are (should be) constant; below attrs change during simulation:
            'throttle_factor': throttle_factor,
            #'throttle_factor': self.reaction_throttle_cache[edge_key][0],
            # 'traversals': 1, # Number of edge traversals = N_reaction_invocation_count
            # 'weight': 1,    # higher weight = more inclined to have actual edge length close to "len".
        }
        edge_attrs['label'] = self.reaction_graph_edge_label_fmt.format(**edge_attrs)
        # Note: We do this for both changed/existing complexes and *new* complexes;
        # Although we don't need to check that the state-fingerprint has changed,
        # we do need to add the new complex to the reaction_graph, and make the reaction_spec_pair edge.

        ## 3a. Update target state nodes for all new or changed complexes:
        added_reaction_graph_nodes = set()
        for cmplx in new_or_changed_complexes:
            ## TODO: This for-loop is pretty long, consider shortening it or moving code to dedicated methods.
            ## TODO: This would also prevent the bug where dS is altered multiple times!

            target_state = cmplx._state_fingerprint

            # new_nodes_added += self.update_reaction_graph_changed_complex(cmplx)
            # Regarding z-coordinate: Some software uses x, y, z attribute names, others uses a single pos=(x,y,z).
            # Adding a z coordinate may cause display issues in e.g. Gephi.
            # offset_relative_to = [[offset-x, (node1, node2)], [offset-y, (node1, node3)], [offset-z, (node4, node5
            # node_attrs = {
            #     "dHdS": cmplx.energy_total_dHdS,
            # }
            count = self.state_counter[target_state]

            if target_state in new_reaction_graph_nodes:
                ## NEW REACTION GRAPH NODE: ##
                # Add node centered between the source states with a random offset:

                n_strands = len(cmplx.strands)
                # Instead of trying to specify a single particular position, give a "offset_relative_to" attribute.
                dG_first = cmplx.energy_total_dHdS[0] - self.temperature*cmplx.energy_total_dHdS[1]
                dG_std = cmplx.energy_total_dHdS[0] - 298.15*cmplx.energy_total_dHdS[1]
                # dH/dG are in units of R*K, I believe...
                x_scene_offset = -4
                pos = [x_scene_offset + 2*n_strands**1.2, 0, dG_first/1e3]
                grid_pos = pos.copy() # [x_scene_offset + 2*n_strands**1.2, 0, 0] # or use dG_first/1e5 for grid_pos[z]?
                # node_attrs["pos"] = [x, y, z]
                # node_attrs["grid_pos"] = [x, y, z]
                offset_relative_to = [
                    (0, None), (random.random(), reaction_spec_source_states_list), (0, None)]
                # For each new node, the initial node position is calculated as:
                # start_pos + [offset+<average source nodes coord> if source_nodes else 0
                #              for offset, source_nodes in offset_relative_to]
                assert count == 1
                self.reaction_graph.add_node(
                    target_state,
                    dHdS=tuple(cmplx.energy_total_dHdS),
                    dHdS_contributions=cmplx.energies_dHdS,
                    n_strands=n_strands,
                    size=sqrt(n_strands),
                    pos=pos,
                    offset_relative_to=offset_relative_to,
                    encounters=1,
                    count=count
                    # x=x, y=y, z=z
                )
                new_nodes_added += 1
                assert target_state not in self.endstates_by_reaction
                # endstates_by_reaction[target_state][edge_key] = [list of targets for reaction]
                self.endstates_by_reaction[target_state] = defaultdict(list)
                assert target_state in self.reaction_graph
                assert target_state in self.endstates_by_reaction
                # record_new_complex_state just saves the reaction_graph to file and draws it using graphviz..
                # self.record_new_complex_state(cmplx, target_state)  # Move to reaction_graph.
                new_reaction_graph_nodes.remove(target_state)
                added_reaction_graph_nodes.add(target_state)
            # Some graph layout algorithms (e.g. graphviz) can only handle two-valued (x,y) positions:
            # pos = [x, y]  # Even this doesn't work with pydot(plus) graphviz.. not sure how to give pos..
            else:
                # Update existing state node:
                encounters = self.reaction_graph.node[target_state].get('encounters', 1)
                self.reaction_graph.change_node(
                    target_state,
                    dHdS=tuple(cmplx.energy_total_dHdS),
                    encounters=encounters,
                    count=count,
                    # n_strands=n_strands, size=10*sqrt(n_strands),
                    # pos=pos,
                    # x=x, y=y, z=z
                )
                ## Check that target is in self.endstates_by_reaction:
                assert target_state in self.endstates_by_reaction
                assert target_state in self.reaction_graph

            assert target_state in self.reaction_graph
            for source_state in reaction_spec_source_states:
                assert source_state in self.reaction_graph
                new_edges_added += self.update_reaction_graph_edge(source_state, target_state, edge_key, edge_attrs)
                assert source_state in self.reaction_graph
                assert target_state in self.endstates_by_reaction[source_state][edge_key]
        # end for cmplx in changed_complexes:

        if reaction_result['obsolete_complexes']:
            for cmplx in reaction_result['obsolete_complexes']:
                # Should be empty; released strands are available in 'free_strands'.
                assert len(cmplx.strands) == 0
                # self.state_counter is updated using asserted_target_states_list and reaction_spec_source_states_list


        ### 3b. For all freed strands: Update reaction graph: ###
        if reaction_result['free_strands']:
            """
            Wait, is this proper? I mean, what if we are in a state with two strands bound to the scaffold and
            one of them fall off. We cannot distinguish one strand falling off from the other strand falling off
            in the reaction graph edge connecting to the "0" node. If we try to use that edge to obtain information
            on the reaction where the other strand fall off, we would get wrong information.
            """
            ## TODO: Fix the case above by not using a "0" node for all strands but having one node for each strand.
            # free strands can only result from breaking (dehybr/unstack) reactions
            assert reaction_attr.is_forming is False
            assert reaction_attr.is_intra
            assert c1state == c2state  # The states in reaction_spec_pair
            # source_state = c1state
            # We should only have free strands if de-hybridizing, not unstacking
            for strand in reaction_result['free_strands']:
                for source_state in reaction_spec_source_states:
                    self.update_reaction_graph_edge(source_state, strand.name, edge_key, edge_attrs)
                del source_state # Debugging, prevent bugs from variable reuse.

        ## 3c: Update source states (with tau, if given): ##
        for source_state in reaction_spec_source_states:
            count = self.state_counter[source_state]
            assert source_state in self.reaction_graph.node
            # Update: count (the number of complexes currently in that state),
            self.reaction_graph.change_node(source_state, count=count)


        ### 3d: Save reaction graph to file: ###
        if (new_edges_added or new_nodes_added) and self.params.get('reaction_graph_save_when_updated'):
            pass
            # self.save_reaction_graph()  # Disabled, performance. [post_reaction_processing, free strands]




        ### 5. Adjust throttle: ###

        # if self.reaction_throttle_use_cache:
        #     self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)] = \
        #         0.95 * self.reaction_throttle_cache[(reaction_spec_pair, reaction_attr)]
        ## What if the complex breaks up?? We don't want to throttle those, do we?
        ## Unless we also are throttling inter-complex reactions? That might actually be easier to get right.
        ## But, of course, then we don't have any way to keep track of reaction counts.

        if self.reaction_throttle and self.reaction_throttle_use_cache and not self.reaction_throttle_freeze:
            if self.reaction_throttle_per_complex:
                print("Per-complex throttle...")
                ## TODO: This really should be done when looping over changed complexes
                ## (but global throttle seems better at the moment)
                # throttle_cache_key = (reacted_spec_pair, reaction_attr)
                ## Per-complex reaction_throttle_cache
                ## We *could* do this in complex.assert_state_change; but then we would have to pass
                ## in reaction_result as argument (if we want to distinguish intRA vs intER molecular reactions...)
                if reaction_attr.reaction_type is STACKING_INTERACTION:
                    cmplx = h1end3p.domain.strand.complex
                else:
                    # Update when we introduce backbone ligation/nicking reactions:
                    assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
                    cmplx = domain1.strand.complex
                # Save current throttle_factor for reaction_graph edge attr:
                throttle_factor = cmplx.reaction_throttle_cache.get(reacted_spec_pair, 1.0)
                if reaction_attr.is_intra and reaction_result['case'] < 2:
                    # result case < 2 is for intRA complex/strand reactions
                    if reacted_spec_pair not in cmplx.reaction_throttle_cache:
                        print("Adding new reacted_spec_pair to %s.reaction_throttle_cache: %s\n" %
                              (cmplx, reacted_spec_pair))
                        cmplx.reaction_throttle_cache[reacted_spec_pair] = 0.9
                    else:
                        cmplx.reaction_throttle_cache[reacted_spec_pair] *= 0.9
                        # print("%s reaction_throttle %0.02f for reacted_spec_pair %s" %
                        #       (cmplx, cmplx.reaction_throttle_cache[reacted_spec_pair], reacted_spec_pair))
                    # print("Reaction throttle for %s reduced to %s" %
                    #       (reacted_spec_pair, cmplx.reaction_throttle_cache[reacted_spec_pair]))
                else:
                    ## Reaction is not bi-directional intRA-molecular. (i.e. intra-molecular in both ways)
                    # print("\n\nreaction_attr.is_intra=%s and reaction_result['case']=%s" %
                    #       (reaction_attr.is_intra, reaction_result['case']))
                    # pdb.set_trace()
                    pass
            else:
                ## Global reaction_throttle_cache (not per-complex):
                ## TODO: Account for free strands (also add edges for dehybridization to free strands)
                ## Save current throttle_factor for reaction_graph edge attr: (maybe not desired?)
                # throttle_factor = self.reaction_throttle_cache.get(edge_key, 1.0)
                if edge_key not in self.reaction_throttle_cache:
                    # New: Make sure to use the same throttle in both directions.
                    # Check whether we have an edge in the other direction before creating a new:
                    # TODO: Re-implement as try...except KeyError
                    # # if edge_key in self.endstates_by_reaction[source_state]:
                    #     # Uhm... how does this work when splitting a complex in two, i.e. we have two
                    #     # outgoing edges with the same (reaction_spec_pair, reaction_attr) edge key...?
                    #     # TODO: Re-implement endstate(S)_by_reaction to hold a LIST of end-states
                    #     # Also, this gives us the target_state, but not reacted_spec_pair, reaction_attr,
                    #     # we have to get them using:
                    #     # reverse_reaction_source_states = self.endstates_by_reaction[source_state][edge_key]
                    #     rev_reaction_source = self.endstates_by_reaction[source_state][edge_key][0]
                    #     # Do we have an edge from rev_reaction_source (=one of this reaction's targets) ?
                    #     if source_state in self.reaction_graph[rev_reaction_source]:
                    #         # reverse_edge_attrs = self.reaction_graph[rev_reaction_source][source_state]
                    #         # rev_edge_key = (reverse_edge_attrs['reaction_spec_pair'],
                    #         #                 reverse_edge_attrs['reactionn_atttr'])
                    #         rev_edge_key = self.reaction_graph[rev_reaction_source][source_state]['edge_key']
                    #         throttle_factor = self.reaction_throttle_cache[rev_edge_key]
                    # # Edit: Use self.reverse_reaction_key instead.
                    if edge_key in self.reverse_reaction_key:
                        rev_edge_key = self.reverse_reaction_key[edge_key]
                        assert rev_edge_key in self.reaction_throttle_cache
                        throttle_factor = self.reaction_throttle_cache[rev_edge_key]
                    else:
                        throttle_factor = [self.reaction_throttle_factor_base]
                    # else if edge_key not in self.reverse_reaction_key, then we risk
                    # creating two throttle_factor encapsulating lists:
                    self.reaction_throttle_cache[edge_key] = throttle_factor
                elif test_throttles:
                    if reaction_attr_str in test_throttles:
                        if reaction_attr.reaction_type is STACKING_INTERACTION:
                            domain1, domain2 = h1end3p.domain, h2end5p.domain
                        if domain1.name in test_throttles[reaction_attr_str]:
                            self.reaction_throttle_cache[edge_key][0] = test_throttles[reaction_attr_str][domain1.name]
                        elif domain2.name in test_throttles[reaction_attr_str]:
                            self.reaction_throttle_cache[edge_key][0] = test_throttles[reaction_attr_str][domain2.name]
                        elif len(test_throttles) > 1 and debug_test_throttles:
                            pdb.set_trace()
                    elif len(test_throttles) > 1 and debug_test_throttles:
                        if edge_key not in self.reaction_throttle_cache:
                            pdb.set_trace()
                else:
                    self.reaction_throttle_cache[edge_key][0] *= self.reaction_throttle_factor_base
                    if edge_key in self.reverse_reaction_key:
                        rev_edge_key = self.reverse_reaction_key[edge_key]
                        assert self.reaction_throttle_cache[edge_key] is self.reaction_throttle_cache[rev_edge_key]
                    # else:
                    #     print("\nWARNING: edge_key\n\t%s\nnot in self.reverse_reaction_key:"
                    #           % reaction_to_str(*edge_key))
                    #     print("\n".join("%s => %s" % (reaction_to_str(*k), reaction_to_str(*v))
                    #                     for k, v in self.reverse_reaction_key.items()))
                    #     pdb.set_trace()

                # print("Global throttle %0.02f for reaction (%s, %s)" %
                #       (self.reaction_throttle_cache[(reacted_spec_pair, reaction_attr)],
                #        reacted_spec_pair, reaction_attr))

        # print("\n".join("\n%s: %s\n%s:%s" % (
        #     reaction_to_str(*k), v, reaction_to_str(*self.reverse_reaction_key[k]),
        #     self.reaction_throttle_cache[self.reverse_reaction_key[k]])
        #                 for k, v in sorted(self.reaction_throttle_cache.items())))
        # print("\n".join("\n%-52s: %s\n%-52s: %s" % (reaction_to_str(*k), v, reaction_to_str(*self.reverse_reaction_key[k]), self.reaction_throttle_cache[self.reverse_reaction_key[k]]) for k, v in sorted(self.reaction_throttle_cache.items())))


    def update_reaction_graph_edge(self, source_state, target_state, edge_key, edge_attrs):
        """
        Used to handle either:
        (a) A new complex, or
        (b) A complex being split into two free strands
        """
        # if reaction_attr.is_forming:
        #     # A new complex is formed from two free strands:
        #     source_states = [strand.name for strand in cmplx.strands]
        #     target_states = (cmplx._state_fingerprint,)
        #     assert 1 <= len(source_states) <= 2
        # else:
        #     # Complex is broken down into two free strands:
        #     source_states = (cmplx._state_fingerprint)
        #     target_states = [strand.name for strand in cmplx.strands]
        #     assert 1 <= len(target_states) <= 2
        # # if edge_key in self.endstates_by_reaction[source_state]:  # [source][key] = [list of targets]
        # # if target_state in self.reaction_graph[source_state]:     # [source][target][key] = eattrs
        #     # self.endstates_by_reaction[source_state][edge_key] is a list of possible end-states.
        #     # I don't think there is a guarantee that target_state "0" is in end state.
        #     # In fact, we might have created the edge_key entry just above.
        #     # assert target_state in self.endstates_by_reaction[source_state][edge_key]
        # for strand in reaction_result['free_strands']:
        # edit: using the same node "0" for all free strands
        # No need to check that the expected target state is in the graph; it should be
        # for source_state in source_states:
        #     for target_state in target_states:
        if target_state not in self.reaction_graph[source_state]:
            ## There is currently no edge from source to target, so add new edge:
            edge_attrs['traversals'] = 1
            edge_attrs['weight'] = 0.5
            self.reaction_graph.add_edge(source_state, target_state, edge_attrs)
            # self.reaction_graph[source_state][target_state]['traversals'] = 1
            # self.reaction_graph[source_state][target_state]['weight'] = 0.5

            # Update self.endstates_by_reaction[source_state]
            # assert edge_key not in self.endstates_by_reaction[source_state]  # Again, could JUST be added above.
            # Did you mean to write assert edge_key in self.endstates_by_reaction[source_state]?
            self.endstates_by_reaction[source_state][edge_key].append(target_state)
            if edge_key not in self.reverse_reaction_key:
                if source_state in self.reaction_graph.adj[target_state]:
                    rev_edge_key = self.reaction_graph.adj[target_state][source_state]['edge_key']
                    self.reverse_reaction_key[edge_key] = rev_edge_key
                    self.reverse_reaction_key[rev_edge_key] = edge_key
            return True
        else:
            ## UPDATE EDGE ATTRS:
            # target_state (strand.name) is in reaction_graph[source_state]
            # There is already an edge from source to target: Check that the edge is correct
            # Again, there might be issues with multiple edges when complexes break down to form a new complex.
            ## TODO: Fix all this.
            existing_edge_attrs = self.reaction_graph[source_state][target_state]
            assert existing_edge_attrs['reaction_attr'] == edge_attrs['reaction_attr']
            assert existing_edge_attrs['reaction_spec_pair'] == edge_attrs['reaction_spec_pair']
            assert existing_edge_attrs['edge_key'] == edge_attrs['edge_key']
            existing_edge_attrs['traversals'] += 1
            existing_edge_attrs['weight'] += 0.1/edge_attrs['reaction_invocation_count'] # weight log vs traversals
            edge_attrs['traversals'] = existing_edge_attrs['traversals']
            edge_attrs['weight'] = existing_edge_attrs['weight']
            existing_edge_attrs['throttle_factor'] = edge_attrs['throttle_factor'] # = throttle_factor
            existing_edge_attrs['c_j_throttled'] = edge_attrs['c_j_throttled'] # = c_j
            existing_edge_attrs['label'] = edge_attrs['label'] # = edge_label_fmt.format(**edge_attrs)
            # Check self.endstates_by_reaction[source_state]
            self.reaction_graph.change_edge(source=source_state, target=target_state, attrs=edge_attrs)
            assert edge_key in self.endstates_by_reaction[source_state]
            assert target_state in self.endstates_by_reaction[source_state][edge_key]
            if edge_key in self.reverse_reaction_key:
                rev_edge_key = self.reverse_reaction_key[edge_key]
                assert self.reaction_throttle_cache[edge_key] is self.reaction_throttle_cache[rev_edge_key]
            return False



    def record_new_complex_state(self, cmplx, state_fingerprint):
        """
        :cmplx: A complex in a state that have not been seen before.
        # nx.write_yaml(cmplx, path) # records *all* objects in the yaml file. Too much.
        # alternatives:
        # write_(multiline_)edgelist,
        gexf - doesn't work, has problem with python types (tuples, etc) node/edge attrs,
                only works with keys in xml_type: int, float, bool, list, dict.
                (Edit: Modified the networkx source, adding tuple to be the same as list.)
        gml - doesn't work with objects as attr values.
        graphml - doesn't support tuples (again, change xml_type)
        GIS Shapefile - write_shp requires OGR: http://www.gdal.org/
        In order of file size:
        * adjlist
        * edgelist
        * multiline_adjlist
        * pajek
        * gexf
        * yaml
        doesnt_work = ("gml", "graphml", "shp", )
        """
        if state_fingerprint is None:
            state_fingerprint = cmplx.state_fingerprint()
        else:
            assert state_fingerprint == cmplx.state_fingerprint()
        # assert state_fingerprint not in self.reaction_graph.node  # node-attrs dict, {node: {node attrs}}
        # assert state_fingerprint not in self.reaction_graph.adj   # edge-attrs dict, {src: {tgt1: {edge1_attrs}, ...}}
        #
        # dHdS = tuple(cmplx.energy_total_dHdS)
        # node_attrs = {'dHdS_first': dHdS,
        #               'dHdS_count': {dHdS: 1},
        #               'encounters': 1,
        #              }
        # # MultiGraph, with edges keyed by (reacted_spec_pair, reaction_attr):
        # # reaction_graph.adj[source][target][(reacted_spec_pair, reaction_attr)] = eattr
        # self.reaction_graph.add_node(state_fingerprint, encounters=1, dHdS_first=tuple(cmplx.energy_total_dHdS))
        # # self.reaction_graph.node[state_fingerprint]['dHdS_first'] = dHdS
        # self.reaction_graph.node[state_fingerprint]['dHdS_count'] = {} # Counter()
        # self.reaction_graph.node[state_fingerprint]['dHdS_count'][dHdS] = 1
        # #self.reaction_graph.node[state_fingerprint] = node_attrs
        # #self.reaction_graph.adj[state_fingerprint] = {}

        ## TODO: Provide suggested x-coordinate:
        ## Should essentially be proportional to the number of strands, BUT ALSO provide room
        ## depending on the previous state (small adjustments within the same area).


        if self.reaction_graph_complexes_directory is not None:
            # Save complex to file:
            # https://networkx.github.io/documentation/latest/reference/readwrite.html
            # adjlist: not useful; each line simply has: str(source) str(target1) str(target2) ...
            # edgelist: Better, each line describes an edge: source target edge_attr
            # multiline_adjlist: edgelist meets adjlist; N_targets+1 lines for each source.
            # pajek: Has node+edge attrs but no multi-graph edge keys.
            # gefx: xml-based; nice, but very verbose. Has all node and edge attrs AND edge keys and node object attrs.
            # yaml: Dumps *everything*, Domains, DomainEnds, Strands, Complexes, etc.
            # for method in ("edgelist", "adjlist", "multiline_adjlist", "gexf", "pajek", "yaml"):
            for method in ("edgelist", "gexf"):
                path = os.path.join(self.reaction_graph_complexes_directory,
                                    "%s.%s" % (state_fingerprint, method))
                write_function = getattr(nx, "write_" + method)
                write_function(cmplx, path)
            # create a graph visualization:
            # start_pos = cmplx.graph['pos'] # edit: saving 'pos' directly in node attr dict:
            # start_pos = nx.get_node_attributes(cmplx, 'pos')
            path = os.path.join(self.reaction_graph_complexes_directory,
                                "%s.%s" % (state_fingerprint, "png"))
            #pos = layout_graph(cmplx, save_pos_as_attr=True)
            # cmplx.graph['pos'] = pos  # Nope, pos is an attr for each node.
            # nx.set_node_attributes(cmplx, 'pos', pos)

            ## Draw complex graph. Disabled for now; instead, process the dot files *after* simulation.
            #draw_graph_and_save(cmplx, path, pos=pos)
            draw_with_graphviz(cmplx, path)

            ## I stopped using per-complex ends5p3p graphs quite some time ago, it was only used for vizualization.
            # path = os.path.join(self.reaction_graph_complexes_directory,
            #                     "%s_ends5p3p.%s" % (state_fingerprint, "png"))
            # draw_with_graphviz(cmplx.graph_5p3p, path)
            complexes_overview_fn = os.path.join(self.reaction_graph_complexes_directory, "complexes_overview.txt")
            with open(complexes_overview_fn, 'a') as fd:
                fd.write("%s\t%s\n" % (self.system_time, state_fingerprint))
            printd("New complex state %s save to file. (%s)" % (state_fingerprint, cmplx))


    def save_reaction_graph(self, fnpostfix=""):
        """
        Save reaction graph to file...
        See also StatsManager.save_reaction_graph (!)
        """

        if self.reaction_graph_complexes_directory is None:
            print("reaction_graph_complexes_directory is not set; cannot save reaction_graph/system graphs...")
            return

        #self.reaction_graph.node[0]['pos'] = (0, 0)

        ## Save system graphs:
        graphs_to_save = (#"reaction_graph",
                          "ends5p3p_graph",
                          "interface_graph")
        fn_number = next(reaction_graph_sequantial_number_generator)
        for graph_name in graphs_to_save:
            #graph_fnfmt = "%s_%s.png" % (graph_name, fn_number)
            path = os.path.join(self.reaction_graph_complexes_directory,
                                "%s_%s%s.png" % (graph_name, fn_number, fnpostfix))
            if os.path.exists(path):
                print("WARNING: Graph output file already exists:", path)
                return
            graph = getattr(self, graph_name)
            #pos = layout_graph(graph, save_pos_as_attr=True, k=3, weight=None, scale=4.0)
            # (defaults: k=sqrt(N_nodes), scale=1.0, weight='weight')
            #draw_graph_and_save(graph, path, pos=pos) # draw graph with matplotlib
            # draw_with_graphviz(graph, path)
            # TODO: Save yaml/gefx format as well (for better import to Gephi)
            printd("%s saved to file %s" % (graph_name, path))
        reaction_graph_complexes_fn = os.path.join(
            self.reaction_graph_complexes_directory, "reaction_graph_complexes.txt")
        # Make sure that for each reaction graph we save, we can go back and determine which complexes were present:
        with open(reaction_graph_complexes_fn, 'a') as fd:
            fd.write("%s\t%s\t%s\n" % (fn_number, self.system_time, ",".join(
                [str(cmplx.state_fingerprint()) for cmplx in self.complexes])))



    def check_system(self):
        """
        Run a complete systems-check to assert that everything is as expected.
        This is only for testing and should not not be part of the "pure performance" class.
        """
        is_good = True

        ### Check self.(un)hybridized_domains_by_name lookup tables: ###

        unhybridized_domains_by_name = defaultdict(set)
        hybridized_domains_by_name = defaultdict(set)
        for name in self.unhybridized_domains_by_name:
            _ = unhybridized_domains_by_name[name]
        for name in self.hybridized_domains_by_name:
            _ = hybridized_domains_by_name[name]
        for name, domains in self.domains_by_name.items():
            for d in domains:
                if d.partner is None:
                    unhybridized_domains_by_name[name].add(d)
                else:
                    hybridized_domains_by_name[name].add(d)
        if self.unhybridized_domains_by_name != unhybridized_domains_by_name:
            print("\nSystem-check: self.unhybridized_domains_by_name != unhybridized_domains_by_name")
            print("self.unhybridized_domains_by_name:")
            pprint(self.unhybridized_domains_by_name)
            print("unhybridized_domains_by_name (real):")
            pprint(unhybridized_domains_by_name)
            is_good = False
        if self.hybridized_domains_by_name != hybridized_domains_by_name:
            print("\nSystem-check: self.hybridized_domains_by_name != hybridized_domains_by_name")
            print("self.hybridized_domains_by_name:")
            pprint(self.hybridized_domains_by_name)
            print("hybridized_domains_by_name (real):")
            pprint(hybridized_domains_by_name)
            is_good = False

        ### Check complexes: ###
        # Check state_counter (only positive counts/values), dict-like (key+value) comparison
        assert +self.state_counter == (Counter([cmplx._state_fingerprint for cmplx in self.complexes])
                                       + Counter([strand.name for strand in self.strands if strand.complex is None]))
        for cmplx in self.complexes:
            # Check complex strands:
            strands_by_name = defaultdict(set)
            for name in cmplx.strands_by_name:
                _ = strands_by_name[name]
            for strand in cmplx.strands:
                strands_by_name[strand.name].add(strand)
            if cmplx.strands_by_name != strands_by_name:
                print("\nSystem-check: Problem with complex %r" % cmplx)
                print("System-check: cmplx.strands_by_name != strands_by_name")
                print("cmplx.strands_by_name:")
                pprint(cmplx.strands_by_name)
                print("strands_by_name (from cmplx.strands):")
                pprint(strands_by_name)
                is_good = False

            # Check complex domains:
            strand_domains = {d for strand in cmplx.strands for d in strand.domains}
            if strand_domains != set(cmplx.nodes()):
                is_good = False
                print("\nSystem-check: Problem with complex %r" % cmplx)
                print("Domains from complex.strands:")
                pprint(strand_domains)
                print("complex.nodes():")
                pprint(cmplx.nodes())


        ### Check self.reaction_pairs_by_domain lookup tables: ###


        ### Check self.possible_hybridization_reactions: ###
        for domain_pair in self.possible_hybridization_reactions:
            ra = self.reaction_attrs[domain_pair]
            d1, d2 = tuple(domain_pair)
            if ra.is_forming:
                assert d1.partner is None and d2.partner is None
            else:
                assert d1.partner is d2 and d2.partner is d1
            # stacked domains cannot dehybridize; unhybridized domains can't be stacked
            # i.e. if it is stacked, the only reaction a domain can undergo is unstack.
            assert all(end.stack_partner is None for end in (d1.end5p, d1.end3p, d2.end5p, d2.end3p))

        ### Check self.possible_stacking_reactions: ###
        for stacking_pair in self.possible_stacking_reactions:
            ra = self.reaction_attrs[stacking_pair]
            dup_end1, dup_end2 = tuple(stacking_pair)
            (h1end3p, h2end5p), (h2end3p, h1end5p) = dup_end1, dup_end2
            # stacked domains must be hybridized (duplexes):
            assert all(end.domain.partner is not None for end in (h1end3p, h2end5p, h2end3p, h1end5p))
            assert h1end3p.hyb_partner is h2end5p and h2end5p.hyb_partner is h1end3p
            assert h2end3p.hyb_partner is h1end5p and h1end5p.hyb_partner is h2end3p
            if ra.is_forming:
                assert all(end.stack_partner is None for end in (h1end3p, h2end5p, h2end3p, h1end5p))
            else:
                assert h1end3p.stack_partner is h1end5p and h1end5p.stack_partner is h1end3p
                assert h2end3p.stack_partner is h2end5p and h2end5p.stack_partner is h2end3p

        ### Check self.reaction_attrs:
        for reaction_pair, reaction_attr in self.reaction_attrs.items():
            elem1, elem2 = tuple(reaction_pair)
            if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
                assert isinstance(elem1, Domain)
                assert isinstance(elem2, Domain)
                assert reaction_pair in self.possible_hybridization_reactions
                # assert reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION
                if reaction_attr.is_forming:
                    assert elem1.partner is None and elem2.partner is None
                else:
                    assert elem1.partner is elem2 and elem2.partner is elem1
            else:
                assert reaction_attr.reaction_type is STACKING_INTERACTION
                assert isinstance(elem1, tuple)
                assert isinstance(elem2, tuple)
                assert reaction_pair in self.possible_stacking_reactions
                (h1end3p, h2end5p), (h2end3p, h1end5p) = elem1, elem2
                if reaction_attr.is_forming:
                    assert all(end.stack_partner is None for end in (h1end3p, h2end5p, h2end3p, h1end5p))
                else:
                    assert all(end.stack_partner is not None for end in (h1end3p, h2end5p, h2end3p, h1end5p))
                    assert h1end3p.stack_partner is h1end5p and h1end5p.stack_partner is h1end3p
                    assert h2end3p.stack_partner is h2end5p and h2end5p.stack_partner is h2end3p

        ### Check self.propensity functions (grouped sysmgr): ###

        if not is_good:
            pdb.set_trace()


    def get_complexes_stats_strs(self, ):
        """ Get a string with information on the current complexes. """
        lines = ["\t".join(str(hash(itm) if isinstance(itm, frozenset) else itm) for itm in (
            cmplx, cmplx._state_fingerprint, cmplx._strands_fingerprint,
            cmplx._hybridization_fingerprint, cmplx._stacking_fingerprint))
                 for cmplx in self.complexes]
        return lines

    def print_complexes_stats(self, ):
        """ Print information on the current complexes. """
        print("Current complexes and fingerprints:")
        print("Complex   \tstate_fingerprint\tstrands_fingerprint\t"
              "hybridization_fingerprint\tstacking_fingerprint")
        print("\n".join(self.get_complexes_stats_strs()))


    def print_reaction_stats(self):
        """ Print information on reaction pathways. """
        print("Reaction pathways (ungrouped reactions):")
        print("\n".join(self.reaction_stats_strs()))


    def reaction_stats_strs(self, ):
        """ Generate a string with reaction stats. """
        reaction_stats_strs = []
        #keyfun = lambda kv: sorted([d.instance_name for d in kv[0]])
        # We sort the domains to make the stats easier to read:
        # Inner sorted sorts the domain pair (d1, d2 vs d2, d1), outer sorts the reactions by domain
        for reactant_pair_str, c_j, reactant_pair in sorted((sorted([repr(d) for d in k]), v, k) \
            for k, v in self.possible_hybridization_reactions.items()):
            #domspec1, domspec2 = tuple(reactant_pair[0])
            d1, d2 = tuple(reactant_pair)
            domain1, domain2 = reactant_pair_str
            is_forming = d1.partner is None
            is_intra = (d1.strand.complex is not None and d1.strand.complex == d2.strand.complex) or \
                       (d1.strand == d2.strand)  # intra-complex OR intra-strand reaction
            reaction_attr = self.reaction_attrs[reactant_pair]
            #reaction_desc = REACTION_NAMES[reaction_attr.is_forming][reaction_attr.reaction_type]
            # Use %10s or %18r for domain str/repr
            stat = (
                "%18s %9s %s %18s" % (domain1, "hybrid" if is_forming else "de-hyb",
                                      "intra" if is_intra else "inter", domain2),
                ((": %0.03e x 1 x 1 = %03e" if reaction_attr.is_forming else ": %0.03e   x   1 = %03e")
                 % (c_j, self.hybridization_propensity_functions[reactant_pair]))
            )
            reaction_stats_strs.append("".join(stat))
        # Inner sorted sorts the domain pair (d1, d2 vs d2, d1), outer sorts the reactions by domain
        for reactant_pair_str, c_j, reactant_pair in sorted((sorted([repr(d) for d in k]), v, k) \
            for k, v in self.possible_stacking_reactions.items()):
            #domspec1, domspec2 = tuple(reactant_pair[0])
            d1, d2 = tuple(reactant_pair)
            domain1, domain2 = reactant_pair_str
            reaction_attr = self.reaction_attrs[reactant_pair]
            reaction_desc = REACTION_NAMES[reaction_attr.is_forming][reaction_attr.reaction_type]
            # Use %10s or %18r for domain str/repr
            stat = ("%42s %9s %s %42s" % (domain1, reaction_desc[:8], #reaction_attr.is_intra, domain2),
                                          "intra" if reaction_attr.is_intra else "inter", domain2),
                    ((": %0.03e x 1 x 1 = %03e" if reaction_attr.is_forming else ": %0.03e   x   1 = %03e")
                     % (c_j, self.stacking_propensity_functions[reactant_pair]))
                   )
            reaction_stats_strs.append("".join(stat))
        return reaction_stats_strs


    def map_reaction_specs_by_domspec(self):
        """
        Generate an equivalent to self.hybridization_reactions_by_domspec
        from self.possible_hybridization_reactions, i.e. a map of
            domspec => {set of reaction_specs involving this domspec}
        """
        reactions_by_domspec = defaultdict(set)
        for domain_pair in self.possible_hybridization_reactions:
            for domspec in domain_pair:
                reactions_by_domspec[domspec].add(domain_pair)
        return reactions_by_domspec


    def get_reactions_by_domspec(self, domspec):
        """ Filter possible_hybridization_reactions, returning only reactions involving :domspec:. """
        return {reaction_spec: v
                for reaction_spec, v in self.possible_hybridization_reactions.items()
                if domspec in reaction_spec[0]}
