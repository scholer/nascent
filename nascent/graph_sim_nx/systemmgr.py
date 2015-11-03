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

# pylint: disable=C0103

"""

Module for managing the whole system.

I think there might be a lot of code involved in managing the system (graphs and structure),
and having all that in the simulator makes it a bit hard to read.

Splitting the "system state" code out to a separate module allows the simulator to be focused
on just the "stochastic simulation" part and be mostly agnostic on the system details.

Question: Who manages "reaction propensity" etc?
 - I'm thinking the system manager takes care of everything except the simulation steps / temperature control.

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
 * Propencity functions
 *

"""

#import os
#import random
from collections import defaultdict
#import math
from math import exp #, log as ln
#from datetime import datetime

import networkx as nx
from networkx.algorithms.components import connected_components, connected_component_subgraphs
import numpy as np


from nascent.energymodels.biopython import DNA_NN4, hybridization_dH_dS
from .constants import N_AVOGADRO #, R # N_AVOGADRO in /mol, R universal Gas constant in cal/mol/K
from .constants import (PHOSPHATEBACKBONE_INTERACTION,
                        HYBRIDIZATION_INTERACTION,
                        STACKING_INTERACTION,
                        N_AVOGADRO, AVOGADRO_VOLUME_NM3)
from .complex import Complex


class SystemMgr():
    """
    Simulator class to hold everything required for a single simulation.
    """

    def __init__(self, volume, strands, params, domain_pairs=None):

        self.params = params
        self.temperature = params.get('temperature', 300)
        self.volume = volume or params.get('volume')
        # Base single-molecule activity of two free single reactants in volume.
        self.specific_bimolecular_activity = 1/self.volume/N_AVOGADRO # × M
        self.include_steric_repulsion = True
        self.complexes = []
        self.removed_complexes = []
        self.strands = strands
        self.strands_by_name = defaultdict(list)
        for strand in strands:
            self.strands_by_name[strand.Name].append(strand)
        print("Strands in self.strands_by_name:")
        print("\n".join("- %10s: %s species" % (sname, len(strands))
                        for sname, strands in self.strands_by_name.items()))
        self.domains_list = [domain for strand in strands for domain in strand.domains]
        self.domains = set(self.domains_list)  # doesn't change

        # Stats - counts
        self.N_domains = len(self.domains)
        self.N_strands = len(self.strands)
        self.N_domains_hybridized = sum(1 for domain in self.domains_list if domain.partner)
        self.N_strands_hybridized = sum(1 for oligo in self.strands if oligo.is_hybridized())
        # Keep track of domain state depletions.
        # If the (grouped) domain species count occilates between 0 and 1, then the grouped approach might be
        # disadvantageous. Instead, propensities should be evaluated for each domain instance.
        # Or you can have a mixed model where hybridized/duplexed domains and domains on free strands are grouped
        # while unhybridized domains within a complex are considered individually.
        # Or you could consider intra-complex reactions individually while group other reactions.
        self.N_state_depletions = 0
        self.N_state_repletions = 0

        self.domains_by_name = defaultdict(list)
        self.unhybridized_domains_by_name = defaultdict(set)
        self.hybridized_domains_by_name = defaultdict(set)

        for d in self.domains:
            self.domains_by_name[d.Name].append(d)
        print("Domains in self.domains_by_name:")
        print("\n".join("- %10s: %s species" % (dname, len(domains))
                        for dname, domains in self.domains_by_name.items()))
        if domain_pairs is None:
            # mapping: dom_a -> dom_A, dom_A -> dom_a
            # TODO: This could perhaps be a list, if you want to have different types of domains interacting,
            # E.g. dom_a could be perfect match for dom_A, while dom_ax has 1 mismatch:
            # domain_pairs[dom_A.Name] = [dom_a.Name, dom_ax.Name]
            # Or it could be a set of sets: {{da, dA}, {dA, dax}} and then generate partners by:
            # partners_species = set(chain(pair for pair in domain_pairs if dA in pair)) - {dA}
            # However, might as well only do this once and save the list!
            # Also, if you change domain_pairs mapping, remember to adjust Domain_dHdS cache as well.
            domain_pairs = {d.Name: d.name.lower() if d.name == d.name.upper() else d.name.upper()
                            for d in self.domains_list}
            # remove pairs without partner:
            domain_pairs = {d1name: d2name for d1name, d2name in domain_pairs.items()
                            if d2name in self.domains_by_name}
        # allow self-complementarity?
        assert not any(k == v for k, v in domain_pairs.items())
        self.domain_pairs = domain_pairs
        #self.upper_case_domains = [d for d in self.domains if d == d.upper()]
        #self.lower_case_domains = [d for d in self.domains if d == d.lower()]

        print("Simulation system manager initiated at V=%s with %s strands spanning %s domains." \
              % (self.volume, len(self.strands), len(self.domains)))


        ### System graphs ###
        # Not sure whether we should have system-level graphs or a graph for each complex.
        # Domain-level graphs can be created for each graph, in cases where we need cheap reduced graphs.
        # Otherwise, system-level graphs are probably cheaper and easier.
        # Why have graphs at multiple levels?
        # - Some operations, e.g. connected_strands() are much cheaper if we have strand-level graph.
        # MultiGraph vs Graph?
        # - Should be a multigraph to support stacking
        # - domains and ends can connect to the same partner via both backbone and stacking interactions
        # - hairpins without loops has hybridization to same domain as phosphate backbone.
        # - Not sure what the performance penalty of multi graphs are vs regular graphs?
        # - We could have a separate graphs without stacking, but not sure when those would be useful?
        self.ends_graph = nx.MultiGraph()
        self.domain_graph = nx.MultiGraph()
        self.strand_graph = nx.MultiGraph()


        ## Symbol nomenclature:
        ## Sᵢ, S₁, S₂ - domain species (domain name). Domains with the same sequence are of the same specie.
        ## Fᵢ, F₁, F₂ - domain state fingerprint - domains of the same specie and in the same strand/complex state.
        ## Dᵢ, D₁, D₂ - domain instances, individual domain molecules.

        ### Caches: ###
        self.cache = {} # defaultdict(dict)
        # Standard enthalpy and entropy of hybridization,
        # indexed as [frosenset({d1.Name, d2.Name})][0 for enthalpy, 1 for entropy]
        # - Note: index with frozenset((a,b)) or just cache[a, b] = cache[b, a] = value? # d[1,2] same as d[(1,2)]
        # --> Creating a set() or frozenset() takes about 10x longer than to make tuple.
        # --> dict assignment with frozenset is 0.4/0.5 us vs 0.17/0.05 for the "double tuple entry" (python/pypy),
        #       so if you have the memory for it, go for the faster tuples which takes 2x memory.
        # Note: Whenever I use "name", I'm refering to the name of the specie - domain or strand specie.
        self.cache['domain_hybridization_energy'] = self.domain_dHdS = {} # indexed simply as Sᵢ or {S₁, S₂} ?

        ## Relative activity - used to moderate the activity of domains based on their accessibility in a complex.
        # Note: Not sure whether to use this or use a "loop energy" type of activation energy.
        # Relative activities, indexed as:
        #  - [({domain1-specie, domain2-specie}, complex-state-fingerprint)] for intra-complex reactions
        #  - [(domain1-specie, domain1-complex-state-fingerprint)] for inter-complex reactions.
        # relative activity is 1 for free strands
        # For inter-complex reactions (including free strands), the combined relative activity is simply the product:
        #   rel_activity = dom1_rel_activity * dom2_rel_activity

        self.cache['intracomplex_activity'] = {}    # indexed by {F₁, F₂}
        self.cache['stochastic_rate_constant'] = {} # indexed by {F₁, F₂}
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

        ## State-dependent hybridization energy cache
        # I think this was for hybridization that required e.g. bending or zippering...
        self._statedependent_dH_dS = {}  # indexed as {Fᵢ}


        # mapping with the different domains in different states
        # indexed as: Fᵢ => list of domains in this conformation
        # For simple complexes, Fᵢ is just (domain-strand-specie, complex-state-fingerprint)
        # However, if the complex has multiple copies of that strand, the above is insufficient to specifying
        # which strand-domain within the complex we are talking about.
        # Using Fᵢ = domain.domain_state_fingerprint() will give us a proper hash.
        self.domain_state_subspecies = defaultdict(set)
        for domain in self.domains_list:
            self.domain_state_subspecies[domain.domain_state_fingerprint()].append(domain)


        # Up-to-date list of hybridizations that are currently possible:
        # contains [dom_specie_spec, c_j, is_hybridizing] tuples
        # {(domain1, cstate), (domain2, cstate)} => c_j, is_hybridizing
        self.possible_hybridization_reactions = {}  # Rj, j=0..M - however, I index by doms_specs
        self.depleted_hybridization_reactions = {}
        # Quick lookup with {F₁: [{F₁, F₂}, ...]}; equivalent to
        # {domspec for k in possible_hybridization_reactions if domspec in k}
        self.hybridization_reactions_by_domspec = defaultdict(set)

        # Propencity functions:  aj(x)
        # {(domain1, cstate), (domain2, cstate)} => a
        self.propensity_functions = {} # Indexed by indexed by {F₁, F₂}



    def init_possible_reactions(self):
        """
        Reactions have:
            doms_specs => propensity_constant c_j, state_change_vector v_j
            {(domain1-specie, cstate), (domain2-specie, cstate)} => c_j, is_hybridizing
            Edit: {F₁, F₂} => c_j, is_hybridizing
            Where F₁, F₂ are state-dependent "finger-prints" aka "state-species" for the domains:
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
        Rxs = {}
        for d1 in self.domains:
            d1_statespecie = d1.domain_state_fingerprint()
            if d1.partner:
                # If d1 is hybridized, then there is only one possible reaction channel: dehybridizing
                d2 = d1.partner
                d2_statespecie = d2.domain_state_fingerprint()
                doms_specs = frozenset((d1_statespecie, d2_statespecie))
                if doms_specs not in Rxs:
                    # tuple values are (propensity constant c_j, is_hybridizing)
                    Rxs[doms_specs] = self.calculate_c_j(d1, d2, is_hybridizing=False)
            else:
                # find all possible reaction channels for domain d1
                # self.domain_pairs[dname] = [list of complementary domain names]
                # all possible hybridization partners/candidates:
                for d2 in (d for cname in self.domain_pairs[d1.name]
                           for d in self.unhybridized_domains_by_name[cname]):
                    d2_statespecie = d2.domain_state_fingerprint()
                    doms_specs = frozenset((d1_statespecie, d2_statespecie))
                    if doms_specs not in Rxs:
                        # R_j = (c_j, v_j) - propensity constant for reaction j
                        Rxs[doms_specs] = self.calculate_c_j(d1, d2, is_hybridizing=True)
        self.possible_hybridization_reactions = Rxs


    def update_possible_reactions(self, changed_domains, d1, d2, is_hybridizing):
        """
        Maybe it is better to just re-calculate a_j for all changed domains,
        and then filter out reactions with a_j = 0.

        Wow, this is getting really complex.
        I'm almost contemplating going back to the old approach where we do not try to group and "count"
        domain states, but instead we just consider them all independent.
        Pros for simulating without grouping domains by domain state:
        - We don't have to keep track of possible_hybridization_reactions, etc.
        - Each domain instance has its own propensity function; (just beware of (d1, d2) and (d2, d1) duplicates).
        - At least for small copy numbers, grouping and keeping track of domain state count is likely inefficient.
        - Because of the large number of possible states, even for large copy numbers, the state count might occilate
            between 0 and 1. That is particularly unfavorable for the "grouped" strategy.
        - We can still use "domain state fingerprints" for caching.
        Cons:
        - If we have large copy numbers, a_j becomes very large. Well... It becomes the number of domains.
            Although, even if we have 1000 domain instances, it still only takes 0.1 ms to find a suitable a_j.
        - If we have 10 domain "a" copies and 10 complementary domain "A" copies, then there are 10*10=100 combinations.
            If
        - Grouping has other advantages, e.g. caching.

        Instead, I guess we could have
        """
        # Add changed domains to new species
        depleted_domspecs = set()
        changed_domspecs = set()
        new_domspecs = set()
        Rxs = self.possible_hybridization_reactions
        old_d1d2_doms_specs = frozenset(d1._specie_state_fingerprint, d2._specie_state_fingerprint)

        if is_hybridizing:
            self.unhybridized_domains_by_name[d1.name].remove(d1)
            self.unhybridized_domains_by_name[d2.name].remove(d2)
            self.hybridized_domains_by_name[d1.name].add(d1)
            self.hybridized_domains_by_name[d2.name].add(d2)
        else:
            self.hybridized_domains_by_name[d1.name].remove(d1)
            self.hybridized_domains_by_name[d2.name].remove(d2)
            self.unhybridized_domains_by_name[d1.name].add(d1)
            self.unhybridized_domains_by_name[d2.name].add(d2)


        for domain in changed_domains:
            # Set vs list:
            # - Lists are slightly faster to iterate over;
            # - Sets are significantly faster for lookups and removing arbitrary elements (about 80 times faster)
            old_domspec = domain._specie_state_fingerprint
            try:
                # This is really expensive for lists so domain_species should be a set
                self.domain_state_subspecies[old_domspec].remove(domain)
            except KeyError as e:
                if domain not in (d1, d2):
                    raise e
            if len(self.domain_state_subspecies[old_domspec]) == 0:
                depleted_domspecs.add(old_domspec)
                del self.domain_state_subspecies[old_domspec]
            else:
                changed_domspecs.add(old_domspec)
            domain.state_change_reset()
            new_domspec = domain.domain_state_fingerprint()
            if new_domspec not in self.domain_state_subspecies:
                new_domspecs.add(new_domspec)
            else:
                changed_domspecs.add(new_domspec)
            self.domain_state_subspecies[new_domspec].add(domain)

        ## Do we need to update hybridized domains (duplexes)?
        ## - As long as we don't have bending that affects k_off rates: No.
        ## But then... how to keep track of hybridized domains??
        ## Things like stacking/unstacking will affect their k_off rate constants, so the count depends on
        ## the state. Unless I do it completely different and just say that
        ## the domspec "fingerprint" for a hybridized domain does not depend on the complex state at all,
        ## but only depends on whether the domain is stacked or not.
        ## E.g. fingerprint = (domain_species_name, (end5p.stack_partner is not None, end3p.stack_partner is not None)


        ## Consideration: Do we really need to check for depleated states all the time?? Maybe just do it
        ## once every 10th run or so.

        ## Depleted domain states:
        for domspec in depleted_domspecs:
            depleted_reactions = {k: v for k, v in self.possible_hybridization_reactions.items()
                                  if domspec in k} # k is frozendict(domspec1, domspec2)
            for k, v in depleted_reactions:
                del self.possible_hybridization_reactions[k]
                del self.propensity_functions[k]
                self.depleted_hybridization_reactions[k] = v
                self.N_state_depletions += 1

        ## New domain states:
        for domspec in new_domspecs:
            repleted_reactions = {k: v for k, v in self.depleted_hybridization_reactions.items()
                                  if domspec in k and all(len(self.domain_state_subspecies[F]) > 0) for F in k}
            for k, v in repleted_reactions:
                del self.depleted_hybridization_reactions[k]
                self.possible_hybridization_reactions[k] = v
                self.N_state_repletions += 1
            d1 = next(iter(self.domain_state_subspecies[domspec]))
            if d1.partner:
                d2 = d1.partner
                d2_domspec = d2.domain_state_fingerprint()
                doms_specs = frozenset((domspec, d2_domspec))
                if doms_specs not in Rxs:
                    # tuple values are (propensity constant c_j, is_hybridizing)
                    Rxs[doms_specs] = self.calculate_c_j(d1, d2, is_hybridizing=False)
                else:
                    # Should not happen - it is supposedly a "new domain state".
                    print(("Weird: doms_specs %s for hybridized domains d1, d2 = (%s, %s) is"
                           "already in possible_reactions - that should not happen.") % (doms_specs, d1, d2))
            else:
                ## No partner -- consider hybridization reactions with all other unhybridized complementary domains.
                for d2 in (d for cname in self.domain_pairs[d1.name]
                           for d in self.unhybridized_domains_by_name[cname]):
                    d2_domspec = d2.domain_state_fingerprint()
                    doms_specs = frozenset((domspec, d2_domspec))
                    if doms_specs not in Rxs:
                        # R_j = (c_j, v_j) - propensity constant for reaction j
                        Rxs[doms_specs] = self.calculate_c_j(d1, d2, is_hybridizing=True)
                        self.hybridization_reactions_by_domspec[domspec] = doms_specs
                        self.hybridization_reactions_by_domspec[d2_domspec] = doms_specs

        ## Changed domains states (i.e. only a change in number of domains in that state)
        ## Just update propensity_functions for those domspec
        for domspec in changed_domspecs | new_domspecs:
            for doms_specs in self.hybridization_reactions_by_domspec[domspec]:
                self.recalculate_propensity_function(doms_specs)


    def recalculate_propensity_function(self, doms_specs):
        """
        Recalculate propensity function for a single doms_specs = frozenset((F1, F2)).
        Assumes that c_j, is_hybridizing is already available in
        self.possible_hybridization_reactions.
        """
        c_j, is_hybridizing = self.possible_hybridization_reactions[doms_specs]
        if is_hybridizing:
            # a_j = c_j * x₁ * x₂
            self.propensity_functions[doms_specs] = c_j * \
                np.prod([len(self.domain_state_subspecies[ds]) for ds in doms_specs])
        else:
            # a_j = c_j * x₃ ,      x₃ is number of duplexes
            self.propensity_functions[doms_specs] = c_j * len(self.domain_state_subspecies[doms_specs[0]])



    def init_all_propencity_functions(self):
        """
        Reaction propensity specifies how likely a certain reaction of one or more species is to occour.
        The product "propensity × dt" denotes the probability that a reaction will occour in infinitesimal time dt.
        This is just the reaction rate times the number of reactant species present (concentration):
            propensity = k_on [domain1] [domain2]   # if hybridizing
            propensity = k_off * [duplex]           # if dehybridizing
        However, since our simulation operates stochastically at the molecular level, we prefer to express reaction
        propencities in term of discrete molecule count, as:
            a_j = c_j * X₁ * X₂         # for hybridization reactions
            a_j = c_j * X₁              # for dehybridization reactions

        Where c_j is the "stochastic rate constant" (which is essentially k_j * stochastic_activity), and
        X₁, X₂ are the number of domains of species X₁, X₂ or duplexes of species X₁.

        Gillespie calls it "propensity function", since propensity is a function of state (concentration/species counts)
        and uses the symbol a_j to indicate the propensity of reaction j. This should not be confused with "activity".
        Gillespie uses a_0 to indicate sum(a_j for a_j in a) - I just use "a_sum".

        The propencify function sum, a_0 or a_sum, is simply: a_sum = sum(a)

        Calculation of initial propencities is done here as:
            a = [c_j * (prod(len(domain_subspecies[ds]) for ds in doms_spec)  # equivalnt to: c_j * X₁ * X₂, or c_j * X₁
                        if is_hybridizing else len(self.domain_state_subspecies[doms_specs[0]]))
                 for doms_spec, c_j, is_hybridizing in possible_reactions_lst]
        Where doms_spec is a set of Fᵢ = (domain-specie, cstate) for each domain in the reaction, used to specify
        all domains of a particular species (e.g. "domain A") in a particular state (e.g. free or in any complex state).
        The "prod(len(domain_subspecies[ds]) for ds in doms_spec)" expression is equivalent to c_j * X₁ * X₂,
        for hybridization reactions, while "len(self.domain_state_subspecies[doms_specs[0]]))" is equivalent to c_j * X₁
        for dehybridzation reactions.

        For uni-molecular reactions, c_j is simply the reaction rate constant, k_j.
        For bi-molecular reactions, c_j is k_j *times* the concentration of a single reactant molecule.
        E.g. for the bi-molecular reaction j between two species, S₁ and S₂, with population/count x₁ and x₂
        and rate constant k_j, we calculate c_j as:
            c_j = k_j / N_AVOGADRO / volume    # c_j should be in units of s⁻¹
        At a constant concentration, as we increase the simulation volume, we decrease the concentration of each
        individual molecule (instance), so c_j decreases. But since the number of specie molecules/instances
        increases by the same factor, the product c_j * x₁ is constant, and the propensity a_j = c_j * x₁ * x₂
        actually *increases*. This is intuitive: If we have a larger reaction volume at a constant concentration,
        we expect more reactions to occour within a given timespan Δt. This is even simpler to see for uni-molecular
        reactions (Sᵢ -> ...), where c_j = k_j and a_j increases linearly with xᵢ, which again is linearly proportional
        to the volume. Meaning that in a given timespan Δt, the larger the volume, the more Sᵢ will react.

        We typically save the "unit-normalized" c_j instead of k_j, because:
            (1) it prevents repeatedly doing the floating-point calculation c_j = k_j / N_AVOGADRO / volume, and
            (2) we don't need to worry about uni- vs bi- vs tri-molecular etc, we just need to multiply c_j with
                the reactant species population count, x₁ (x₂, x₃, etc) to get the propensity function a_j.
        """
        # doms_specs (plural) is: frozenset({(domain1, cstate), (domain2, cstate)})
        a = {doms_specs: c_j * (np.prod([len(self.domain_state_subspecies[ds]) for ds in doms_specs])
                                if is_hybridizing else len(self.domain_state_subspecies[doms_specs[0]]))
             for doms_specs, (c_j, is_hybridizing) in self.possible_hybridization_reactions.items()}
        self.propensity_functions = a


    def calculate_c_j(self, d1, d2, is_hybridizing):
        """
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
        * a_j = c_j x₁ x₂       # propencity function, has unit of s⁻¹

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
        if is_hybridizing:
            if d1.complex == d2.complex != None:
                # Intra-complex reaction:
                stochastic_activity = self.intracomplex_activity(d1, d2)
            else:
                # Inter-complex reaction:
                stochastic_activity = self.specific_bimolecular_activity # Same as 1/self.volume/N_AVOGADRO × M
            if self.include_steric_repulsion:
                steric_activity_factor = np.prod([1 if d.complex is None else self.steric_activity_factor(d)
                                                  for d in (d1, d2)])
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
            k_off = self.dehybridization_rate_constant(d1, d2)
            # For uni-molecular reactions, c_j = k_j
            c_j = k_off
        return c_j, is_hybridizing

    def steric_activity_factor(self, d):
        """
        """
        # TODO: Implement steric_activity_factor
        return 1

    def intracomplex_activity(self, d1, d2):
        """
        Return the activity for hybridization of two domains within a complex.
        TODO: Add steric/electrostatic/surface repulsion (individual activity contribution for d1 and d2,
              so that activity = intercomplex_activity*d1_activity*d2_activity)
        """
        assert d1.complex == d2.complex != None
        cache_key = frozenset([d.domain_state_fingerprint() for d in (d1, d2)])
        cache = self.cache['intracomplex_activity']
        if cache_key in cache:
            return cache[cache_key]
        activity = d1.complex.intracomplex_activity(d1, d2)
        cache[cache_key] = activity
        return activity


    def hybridization_rate_constant(self, d1, d2):
        #T = self.temperature
        # At T >> Tm (ΔG° ~ 0), this depends on ΔG°.
        # Also depends on:
        # * salt / ionic strength
        # * length of d1, d2 (and the full length of their strands)
        # But for now, just return default rate constant of 1e5 /s/M
        # Another question: What about intra-complex hybridization?
        #   Should this be accounted for here at the hybridization rate constant,
        #   or at the propensity constant?
        return 1e5


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

        stack_dH, stack_dS = self.stacking_energy(d1, d2)
        dH += stack_dH
        dS += stack_dS

        K = exp(dS - dH/T)
        k_on = 1e5  # Is only valid for if K > 1
        k_off = k_on/K
        return k_off

    def stacking_energy(self, d1, d2):
        """ Calculate stacking energy for d1 and d2. """
        avg_stack_dH, avg_stack_dS = -4000, -10
        # TODO: Add base-specific stacking interactions
        stack_dH, stack_dS = 0, 0
        for d in (d1, d2):
            if d.end5p.stack_partner:
                stack_dH += avg_stack_dH
                stack_dS += avg_stack_dS
            if d.end3p.stack_partner:
                stack_dH += avg_stack_dH
                stack_dS += avg_stack_dS
        return stack_dH, stack_dS

    def dHdS_from_state_cache(self, d1, d2):
        """
        Cache hybridization energy
        Returns
            dH, dS
        where dH is in units of gas constant R, and dS in units of R/K
        """
        doms_state_hash = frozenset(d.domain_state_fingerprint() for d in (d1, d2))
        if doms_state_hash not in self._statedependent_dH_dS:
            dH, dS = self.domain_dHdS[(d1.name, d2.name)]
            # TODO: make adjustments, e.g. due to:
            # * dangling ends
            # * stacking
            # * mechanical bending of the helix
            # * electrostatic repulsion
            self._statedependent_dH_dS[doms_state_hash] = (dH, dS)
        return self._statedependent_dH_dS[doms_state_hash]


    def get_relative_activity(self, d1, d2):
        # Historical domain reaction rate constants
        # [T][{(domain1-specie, complex-state), (domain2-specie, complex-state)}]
        # [T][{(domain1-specie, complex-state), (domain2-specie, complex-state)}]
        # Question: Is it really needed to save all permutations of free and complexed rate constants?
        # At least when the two domains are in separate complexes, it seems excessive.
        # But is is probably nice to have for when the two domain species occupy the same complex
        # (where a detailed calculation actually matter)
        # maybe just index by:
        #  - For intra-complex reactions:
        #    relative_intra-complex-reactivity for d1<->d2 =
        #        [({domain1-specie, domain2-specie}, complex-state-fingerprint)]
        #  - For inter-complex reactions, you extract the relative reactivity of each:
        #     d1_relative_activity = [(domain1-specie, domain1-complex-state-fingerprint)]
        #     d2_relative_activity = [(domain2-specie, domain2-complex-state)]
        #  - If the domain strand is not in a complex, relative activity is set to 1.
        #      [{domain1-specie, complex-statedomain2-specie}, complex-fingerprint, intra-complex-reactivity)]
        if d1.complex is None:
            d1_rel_act = 1
        if d2.complex is None:
            d2_rel_act = 1
        if d1.complex == d2.complex != None:
            rel_act = 1


    def n_hybridized_domains(self):
        """ Count the number of hybridized domains. """
        count = sum(1 for domain in self.domains_list if domain.partner)
        if not count % 2 == 0:
            print("Weird - n_hybridized_domains counts to %s (should be an even number)" % count)
            print("Hybridized domains:", ", ".join(str(domain) for domain in self.domains_list if domain.partner))
        return count

    def n_hybridized_strands(self):
        """ Count the number of hybridized strands. """
        return sum(1 for oligo in self.strands if oligo.is_hybridized())


    def hybridize(self, domain1, domain2):
        """
        Splitting out logic in preparation for Julia implementation.
        returns
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        Where changed_complexes includes new_complexes but not obsolete_complexes.
        """

        assert domain1 != domain2
        assert domain1.partner is None
        assert domain2.partner is None

        #dset = frozenset((domain1, domain2))
        #sset = frozenset(domain1.strand, domain2.strand)
        domain1.partner = domain2
        domain2.partner = domain1
        strand1 = domain1.strand
        strand2 = domain2.strand
        c1 = strand1.complex
        c2 = strand2.complex

        # Update system-level graphs:
        edge_kwargs = {"key": HYBRIDIZATION_INTERACTION, "interaction": HYBRIDIZATION_INTERACTION}
        key = (domain1.universal_name, domain2.universal_name, HYBRIDIZATION_INTERACTION)
        self.strand_graph.add_edge(strand1, strand2, key=key, interaction=HYBRIDIZATION_INTERACTION)
        self.domain_graph.add_edge(domain1, domain2, **edge_kwargs)
        self.ends_graph.add_edge(domain1.end5p, domain2.end3p, **edge_kwargs)
        self.ends_graph.add_edge(domain2.end5p, domain1.end3p, **edge_kwargs)

        changed_complexes, new_complexes, obsolete_complexes, free_strands = None, None, None, []

        if strand1 == strand2:
            # If forming an intra-strand connection, no need to make or merge any Complexes
            if strand1.complex:
                assert strand1.complex == strand2.complex
                strand1.complex.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION)
                changed_complexes = [c1]
            return changed_complexes, new_complexes, obsolete_complexes, free_strands


        ## Update complex:
        if c1 and c2:
            if c1 == c2:
                # Intra-complex hybridization
                assert strand1 in c2.strands and strand2 in c1.strands
                # c1.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION) # done below
                c_major = c1
            else:
                # Merge the two complexs:
                #obsolete_complexes = merge_complex()
                c_major, c_minor = (c1, c2) if (len(c1.nodes()) >= len(c2.nodes())) else (c2, c1)
                for strand in c_minor.strands:
                    strand.complex = c_major
                c_major.strands |= c_minor.strands
                c_major.N_strand_changes += c_minor.N_strand_changes
                c_major.add_edges_from(c_minor.edges())
                c_minor.strands = set()
                c_minor.adj = {}        # For networkx, the graph.adj attribute is the main data store
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
            new_complexes = [new_complex]

        c_major.domain_distances = {} # Reset distances
        c_major.add_edge(domain1, domain2, interaction=HYBRIDIZATION_INTERACTION)
        changed_complexes = [c_major]

        assert strand1.complex == strand2.complex != None

        return changed_complexes, new_complexes, obsolete_complexes, free_strands



    def dehybridize(self, domain1, domain2):
        """
        Dehybridize domain2 from domain1.
        Returns:
            changed_complexes, new_complexes, obsolete_complexes, free_strands
        """

        assert domain1 != domain2
        assert domain1.partner == domain2 != None
        assert domain2.partner == domain1 != None

        # dset = frozenset((domain1, domain2))
        #sset = frozenset(domain1.strand, domain2.strand)
        domain1.partner = None
        domain2.partner = None

        strand1 = domain1.strand
        strand2 = domain2.strand
        c = strand1.complex
        c.domains_distances = {}    # Reset distances.
        assert c == strand2.complex

        # Update system-level graphs:
        s_edge_key = (frozenset((domain1.universal_name, domain2.universal_name)),
                      HYBRIDIZATION_INTERACTION)
        self.strand_graph.remove_edge(strand1, strand2, key=s_edge_key)
        self.domain_graph.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        self.ends_graph.remove_edge(domain1.end5p, domain2.end3p, key=HYBRIDIZATION_INTERACTION)
        self.ends_graph.remove_edge(domain2.end5p, domain1.end3p, key=HYBRIDIZATION_INTERACTION)

        if strand1 == strand2 and c is None:
            # The domains are on the same strand and we don't have any complex to update
            return None, None, None, [strand1]
        assert c is not None

        c.remove_edge(domain1, domain2, key=HYBRIDIZATION_INTERACTION)
        c.strand_graph.remove_edge(strand1, strand2, s_edge_key)

        changed_complexes, new_complexes, obsolete_complexes, free_strands = [c], None, None, []

        ## Break the complex:
        # Update distance:
        c.domain_distances = {}
        #dist = distance_cached(domain1, domain2) # We have just reset, so needs to be re-calculated
        dom1_cc_oligos = nx.node_connected_component(self.strand_graph, strand1)
        # Could also use nx.connected_components_subgraphs(c)

        if strand2 in dom1_cc_oligos:
            ## The two strands are still connected: No need to do anything further
            return changed_complexes, new_complexes, obsolete_complexes, free_strands


        #### The two strands are no longer connected: ####

        ## Need to split up. Three cases:
        ## Case (a) Two smaller complexes - must create a new complex for detached domain:
        ## Case (b) One complex and one unhybridized strand - no need to do much further
        ## Case (c) Two unhybridized strands

        dom2_cc_oligos = nx.node_connected_component(self.strand_graph, strand1)
        assert strand2 not in dom1_cc_oligos
        assert strand1 not in dom2_cc_oligos
        dom1_cc_size = len(dom1_cc_oligos)  # cc = connected component
        dom2_cc_size = len(dom2_cc_oligos)
        print("dom1_cc_size=%s, dom2_cc_size=%s, len(c.strands)=%s" %
              (dom1_cc_size, dom2_cc_size, len(c.strands)))

        assert len(c.strands) == dom1_cc_size + dom2_cc_size

        if dom2_cc_size > 1 or dom1_cc_size > 1:
            # Case (a) Two smaller complexes or case (b) one complex and one unhybridized strand
            # Determine which of the complex fragments is the major and which is the minor:

            if dom2_cc_size > 1 and dom1_cc_size > 1:
                # Case (a) Two smaller complexes - must create a new complex for detached domain:
                graph_minor, graph_major = sorted(connected_component_subgraphs(c),
                                                  key=lambda g: len(g.nodes()))
                # Remember to add the departing domain.strand to the new_complex_oligos list:
                new_complex_oligos = set(dom2_cc_oligos if domain1 in graph_major.nodes() else dom1_cc_oligos)
                # Use graph_minor to initialize; then all other is obsolete
                c_new = Complex(data=graph_minor, strands=new_complex_oligos)
                c.strands -= new_complex_oligos
                c.remove_nodes_from(graph_minor.nodes())
                new_complexes = [c_new]
                changed_complexes.append(c_new)
            else:
                # case (b) one complex and one unhybridized strand - no need to do much further
                domain_minor = domain1 if dom1_cc_size == 1 else domain2
                free_strands = [domain_minor.strand]
                c.strands.remove(domain_minor.strand)
                c.strand_graph.remove_node(domain_minor.strand)
        else:
            # Case (c) Two unhybridized strands
            free_strands = [d.strand for d in (domain1, domain2)]
            c.strands -= {strand1, strand2}
            assert c.strands == set()
            obsolete_complexes = [c]
            strand1.complex, strand2.complex = None, None

        assert domain1.partner is None
        assert domain2.partner is None
        return changed_complexes, new_complexes, obsolete_complexes, free_strands
