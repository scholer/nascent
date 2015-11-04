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



class StructureAnalyzer():

    def __init__(self, ends5p3p_graph=None, domain_graph=None, strand_graph=None):
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
        self.ends5p3p_graph = ends5p3p_graph or nx.MultiGraph()
        self.domain_graph = domain_graph or nx.MultiGraph()
        self.strand_graph = strand_graph or nx.MultiGraph()


    def domain_path_elements(self, path):
        """
        Returns a list of structural elements based on a domain-level path (list of domains).
        """
        #path_set = set(path)
        remaining_domains = path[:] # deque(path)
        elements = [] # list of structural elements on path
        while remaining_domains:
            domain = remaining_domains.pop(0) # .popleft() # use popleft if using a deque
            if not domain.partner:
                elem = SingleStrand(domain)
            else:
                elem = DsHelix(domain)
                # Determine if DsHelix is actually a helix bundle:
                # Uhm...
                # Would it be better to have a "complete" helix structure description of the complex
                # at all time?
                # Whatever... for now
            elements.append(elem)
            if not remaining_domains:
                break
            i = 0
            # while i < len(remaining_domains) and remaining_domains[i] in elem.domains:
            for i, domain in remaining_domains:
                if domain not in elem.domains:
                    break
            else:
                remaining_domains = []
            if i > 0:
                remaining_domains = remaining_domains[i:]

    def ends5p3p_path_partial_elements(self, path, length_only=False, summarize=False):
        """
        Returns a list of structural elements based on a 5p3p-level path (list of 5p3p ends).

        For this, I will experiment with primitive list-based representations rather than making full
        element objects.
        path_edges = [(1, [(length, length_sq, source, target), ...]),
                      (3, [(length, length_sq, source, target), ...]),
                      (2, [(length, length_sq, source, target)]
        Where 1 indicates a list of single-stranded edges,
        2 indicates a hybridization edge, and 3 indicates a list of stacked edges.
        Since only the interaction type and lenghts are important, maybe just
        path_edges = [(1, [length, length, ...]), (3, [length, length, ...], ...)]
        Edit: Instead of 1/2/3, use the standard INTERACTION constants values.
        """
        path_edges = []
        last_interaction = None
        interaction_group = None
        for i in range(len(path)-1):
            source, target = path[i], path[i+1]
            # Method 1: Use ends5p3p_graph edge attributes (which must include stacking edges!):
            edge = self.ends5p3p_graph[source][target]
            length = edge.get('length')     # TODO: Reach consensus on "length" vs "len" vs "weight"
            interaction = edge.get('interaction')
            if interaction is None:
                # Fallback method: Manually determine what type of interaction we have based on source, target:
                if target == source.stack_partner:
                    interaction = STACKING_INTERACTION
                    length = 1
                    length_sq = 1
                elif target == source.hyb_partner:
                    interaction = HYBRIDIZATION_INTERACTION
                    length = 1 # one nm from one helix to the next. We really don't know for sure because it turns.
                    length_sq = 1
                elif target in (source.bp_upstream, source.pb_downstream):
                    if source.hyb_partner and \
                        ((source.end == "5p" and target == source.pb_downstream) or
                         (source.end == "3p" and target == source.bp_upstream)):
                        # Above could probably be done with an XOR:
                        # (source.end == "5p") != (target == source.pb_downstream) # boolean xor
                        # (source.end == "5p") ^ (target == source.pb_downstream)  # bitwise xor
                        # We have a stacked duplex:
                        interaction = STACKING_INTERACTION
                        length = source.domain.ds_length_nm
                        length_sq = source.domain.ds_length_sq
                    else:
                        # We have a single-stranded domain:
                        interaction = PHOSPHATEBACKBONE_INTERACTION
                        length = source.domain.ss_length_nm     # mean end-to-end length; not contour length
                        length_sq = source.domain.ss_length_sq
                else:
                    raise ValueError("Could not determine interaction between %s and %s" % (source, target))

            if interaction == HYBRIDIZATION_INTERACTION:
                # We only distinguish between ds-helix vs single-strand; use STACKING_INTERACTION to indicate ds-helix:
                interaction = STACKING_INTERACTION

            if interaction != last_interaction:
                path_edges.append((last_interaction, interaction_group))
                interaction_group = []
                last_interaction = interaction
            if length_only:
                if length_only == 'sq':
                    interaction_group.append(length_sq)
                elif length_only == 'both':
                    interaction_group.append((length, length_sq))
                else:
                    interaction_group.append(length)
            else:
                interaction_group.append((length, length_sq, source, target))
        if summarize and length_only:
            if length_only == 'both':
                # Return a list of (interaction, (length, length_squared)) tuples:
                return [(interaction, (sum(lengths) for lengths in zip(*lengths_tup))) # pylint: disable=W0142
                        for interaction, lengths_tup in path_edges]
            else:
                return [(interaction, sum(lengths)) for interaction, lengths in path_edges]
        return path_edges



    def domains_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        """
        return shortest_path(self.domain_graph, domain1, domain2)

    def ends5p3p_shortest_path(self, domain1, domain2):
        """
        TODO: This should certainly be cached.
        TODO: Verify shortest path for end3p as well?
        """
        return shortest_path(self.ends5p3p_graph, domain1.end5p, domain2.end5p)


    def intracomplex_activity(self, domain1, domain2):
        r"""
        Returns
            :intracomplex_activity:
        between domain1 and domain2, so that
            c_j = k_j * intracomplex_activity
        The intracomplex activity is basically just:
            activity = 1 / (N_A * effective_volume) = N_A⁻¹ * Ω⁻¹    [unit: M = mol/L]
        where NA is Avogadro's constant, 6.022e23/mol.

        The activity has the same value as the unitless (P_loop/P_v0) just multiplied with "× M" to get unit of M.
        Thus, the activity returned by this function can be interpreted as a relative probability ratio
        denoting the probability that two domains/reactants will be in sufficient proximity to react,
        relative to the probability of two reactants confined within an Avogadro volume, v0 will react, with
            v0 = 1/(NA M) = 1/(6.022e23 mol⁻¹ mol L⁻¹) = 1.6e-27 m³ = 1.6 nm³
        Two molecules confined in v0 is equivalent to a solution of reactants with a concentration of 1 M
        (i.e. standard conditions and a standard activity of 1). The reason we want to specify a relative collision
        probability / activity (instead of an absolute) is that determining reaction probabilities from
        first principles is very hard. Calculating a relative probability/activity allows us to use empirical data
        such as k_on, and energies (ΔH, ΔS) at standard conditions.

        For a walkthrough of this argument, see Dannenberger et al, 2015.

        In addition to determining the stochastic rate constant, activity (or volume or inverse volume)
        can also be used to calculate the loop energy:
        dG  = R T ln(effective_volume/avogadro_volume)      # avogadro_volume = 1/(NA × M) = 1.6 nm³
            = R T ln((NA × M)/molecular_activity)           # molecular_activity = 1/effective_volume
            = - R T ln(loop_activity × M⁻¹)                 # activity = 1/(effective_volume × NA)

        Initially, I considered returning inverse_volume (L-1) or mean-root-square-radius (m2), but I think
        activity is the most straight-forward to use and it can always be converted to either of the other values.

        Regarding names:
        - "activity" -  should be have unit of M=mol/L...
                        although at the single molecule level it could be argued to be just 1/L.
        - volume    -   Unit of L, obviously. (Or maybe m³ or nm³.)
        - probability   Should be unit-less. Although we can just say that we return a relative probability factor
                        so that c_j = k_j × rel_prob_factor × M.

        Alternative names:
        - loop_activity
        - intracomplex_activity             [could be return value as 1/L or mol/L=M]
        - intracomplex_stochastic_activity
        - molar_collision_probability  (except it is relative)
        - spatial_overlap_factor
        - localization_cost  (except we usually associate "cost" with energy)

        ## IMPLEMENTATION: ##

        a: Flow-chart style:
        1. Are the two domains connected by a single ds helix?
        2. (...)

        b: More direct approach:
        1. Determine one (or multiple?) path(s) connecting domain 1 and 2.
        2. Determine the structural elements that make up this path:
            Flexible, single-stranded connections.
            Semi-rigid double-stranded helices.
            Rigid, multi-helix bundles.
             * Question: How about stacked/rigid interface between a ds-helix and multi-helix bundle?
                - Consider as a single, hard-to-calculate element?
        4. Determine the length of longest rigid element (LRE),
            and the length of all other elements plus half of the length of domain1/2.
            (sum remaining elements, SRE)
        5. If the SRE is as long or longer than the LRE, then the domains can immediately hybridize,
            under loop energy EX
        6. If the SRE is shorter than the LRE, but the LRE is a single ds helix, then the
            domains can hybridize under helix bending energy EB.
        7. Otherwise, if the LRE is longer than the SRE and the LRE is a rigid multi-bundle element,
            then for now we assume that the domains cannot bind.

        Regarding mixed-level optimization (APPROXIMATION!):
        1. Get the path at the strand level
        2. Get domain-level subgraph only for domains with strand in the strand-level path
        3. Get domain-level shortest path using the subgraph from (2).
        4. Get 5p3p-level subgraph with 5p3p ends whose domain is in the domain-level shortest path.
        5. Get 5p3p-level shortest path using the subgraph from (4).
        Critizism:
        * May be far from the "proper" 5p3p shortest path.
        * The strand-level graph has a hard time evaluating edge distances (weights).
            For instance, all staples are connected to the scaffold. What is the distance between two staples?

        Edit: I really feel the only way to properly represent the path is to use 5p3p representation.
        Domain-level representation is simply not sufficient.
        For instance, in the "bulge" structure below, domains C and c are certainly within reach:
                                            _ _ _C_ _ _  3'
        5'------------A---------------B----/
        3'------------a----------
                                 \ _ _ _c_ _ _ 5'
        However, the domain level representation:
            A -- B -- C
            |
            a
              \- c
        Is equivalent to the would-be circular structure:
                3'_ _ _C_ _ _
                             \----B----------------A----------5'
                                       ------------a----------
                                                              \ _ _ _c_ _ _ 3'
        Where, depending in helix Aa, domains C and c may not be within reach.
        The domain-level graph looses some detail.
        This *can* be recovered via edge attributes, e.g. edge A-B can be directed or otherwise
        inform that B is on the same side of A as c is. But that might be just as much work
        as just using the 5p3p-ends graph.

        Note: Previously I also intended this function to determine whether domain1 and domain2 can hybridize
        and what the energy penalty is, i.e. loop energy, helix-bending energy, etc.


        ## GLOBAL vs LOCAL model: Using only the minimum loop (shortest path) vs all paths ##

        Consider the following two model cases:
         ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
         |           ⁞‾‾‾‾‾‾‾‾‾‾‾               -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾
         |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃                      |           ⁞___________
         ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'
         ˏ_____A_____ˍ_____B_____ˍ_____C_____₅       ˏ_____A_____ˍ_____B_____ˍ_____C_____₅
         |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|  -->  |           ⁞‾‾‾‾‾‾‾‾‾‾‾ ‾‾‾‾‾‾‾‾‾‾‾|
         |           ˋ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃   ₃__________⌡       |           ⁞___________ ₃__________⌡
         ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'      ˋ‾‾‾‾‾D‾‾‾‾‾ˉ‾‾‾‾‾E‾‾‾‾‾ˉ‾‾‾‾‾F‾‾‾‾‾³'

        In both cases, when connecting (B3p-E5p), we would only consider the shortest path, (B3p-A5p-A3p-D5p-D3p-E5p).
        However, the second-shortest path (B3p-B5p-C3p-C5p-F3p-F5p-E3p-E5p) would clearly also have an influence
        on the PDF overlap (or, activity/effective_volume) of domain B and E.

        To make the energy calculation more precise, you could do a search for secondary and tertiary loops.
        If these are present, that would increase the activity.

        More refs:
        * https://en.wikipedia.org/wiki/Loop_entropy

        TODO: Check for secondary loops
        TODX: Saving secondary loops/paths might be be useful when determining if a dehybridization will split
              up a complex. Although that is pretty much the reverse process and can't really be used. Nevermind.

        """
        #path = self.domains_shortest_path(domain1, domain2)
        #path_elements = self.domain_path_elements(path)
        # NOTE: The path does not have to span the full length of the element!
        # The element could be really long, e.g. a long DsHelix or the full ss scaffold,
        # with the path only transversing a fraction of the element.
        # Also, for domain-level helices, it is hard to know if the path traverses
        # the domain, or just uses it for the hybridization, traversed at one end only.
        #      To here
        #         |
        # -------´ ---------
        # ------------------ ------from here
        #LRE_len, LRE = max((elem.length_nm, elem) for elem in path_elements if not isinstance(elem, SingleStrand))

        ## 5p3p-level shortest path:
        path = self.ends5p3p_shortest_path(domain1, domain2)
        path_elements = self.ends5p3p_path_partial_elements(path, length_only='both', summarize=True)

        ## TODO: Check for secondary loops!

        # list of [(interaction, total-length), ...]
        # For rigid, double-helical elements, element length, l, is N_bp * 0.34 nm.
        # For single-stranded elements, we estimate the end-to-end distance by splitting the strand into
        #   N = N_nt*0.6nm/1.8nm segments, each segment being the Kuhn length 1.8 nm of ssDNA,
        # where 0.6 nm is the contour length of ssDNA and N_nt is the number of nucleotides in the strand.
        #   E[r²] = ∑ Nᵢbᵢ² for i ≤ m = N (1.8 nm)²
        #         = round(N_nt*0.6nm/1.8nm) (1.8 nm)² = N_nt * 0.6/1.8 * 1.8*1.8 nm² = N_nt 0.6*1.8 nm²
        #         = N_nt * lˢˢ * λˢˢ = N_nt * 1.08 nm²
        # Why use interaction as first maximum criteria??
        _, LRE_len, LRE_len_sq, LRE_idx = max((interaction, elem_length, elem_len_sq, i)
                                              for i, (interaction, (elem_length, elem_len_sq))
                                              in enumerate(path_elements)
                                              if interaction in (STACKING_INTERACTION, HYBRIDIZATION_INTERACTION))
        #LRE = path_elements[LRE_idx] # .pop(LRE_idx)
        # Exclude LRE when calculating SRE length:
        SRE_lengths, SRE_sq_lengths = zip(*[(elem_length, elem_len_sq)
                                            for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
                                            for interaction, (elem_length, elem_len_sq) in sub_path])
        # SRE_len_sq = sum(elem_len_sq for sub_path in (path_elements[:LRE_idx], path_elements[LRE_idx+1:])
        #               for interaction, (elem_length, elem_len_sq) in sub_path)
        SRE_len, SRE_len_sq = sum(SRE_lengths), sum(SRE_sq_lengths)

        # Comparing mean-end-to-end-squared values vs mean end-to-end lengths vs full contour length?
        # There is a difference that sum of squares does not equal the square of sums, so
        # even if LRE_len_sq is > SRE_len_sq, LRE_len could be less than SRE_len.
        # Also, while ds duplexes has full contour length equal to mean end-to-end length,
        # this is not true for ssDNA. Indeed if calculating whether two domains *can* reach,
        # it is probably better to use domain.ds_length.

        if LRE_len > SRE_len:
            # The domains cannot reach each other.
            # Hybridization requires helical bending; Not implemented yet; just returning 0 meaning "impossible".
            # TODO: Implement hybridization via helical bending.
            #   Persistance length 50 nm (physiological salt)
            #   - Depends on ionic strength and cationic valency
            # TODO: Look at formulas for k_on and k_off rates under stress.
            #   For DNA, there is certainly a difference between axial "ripping" and perpendicular "zipping".
            #   - Zippering occours at about 10-15 pN (sequence dependent).
            #   -
            return 0
        ## There is probably some profound relation between the elements and the gamma factor.
        ## E.g. if the contour length is long enough for the domains to reach, but the
        ## SRE mean squared end-to-end distance is less than the LRE, then the SRE will rarely
        ## be sufficiently extended for the domains to hybridize. This decrease in spatial pdf
        ## can be considered equivalent to an increase in effective volume.
        ## Another, more complex case, is when the SRE has only (a) one, or (b) a few links,
        ## in which case the mean squared end-to-end distance is not a good measure of spatial pdf.
        ## In the case where SRE is a rigid 1-element chain of same length as the LRE, the pdf
        ## is essentially a sphere centered at the joint between LRE and SRE. (Similar case when
        ## the LRE is flanked by two rigid elements.)

        # Example 1: SRE_len = 10 * 4 nm = 40 nm; SRE_len_sq = 10 * (4 nm)**2 = 160 nm2.
        #            LRE_len = 20 nm,             LRE_len_sq = (20 nm)**2 = 400 nm2.
        #            LRE_len < SRE_len, but LRE_len_sq > SRE_len_sq
        #            SRE_len/LRE_len = 2 -- higher => lower gamma_corr.
        #            LRE_len_sq/SRE_len_sq = 2.5 -- higher => higher gamma_corr.
        # We could, for instance, say: gamma_corr = 1 + ln(LRE_len_sq/SRE_len_sq)
        # Hmm... probably need to do some further analysis of different examples and try to figure out
        # a proper relationship between link-elements and gamma_corr... And it might not be as simple
        # as a simple exponential correction to (P_loop/P_v0).

        # If LRE_len_sq > SRE_len_sq, then the approximation assumption "we many links of length l_i"
        # is certainly not valid (we only have 1 link of length LRE_len).
        # Instead of considering P_loop(r<rc)/P_v0(r<rc), we have to consider P_SRE(r=LRE+/-rc)/V(r=LRE)/P_v0(r<rc).
        # that is, the probability of the SRE end-end distance equaling LRE length, normalized by the shell
        # volume at r=LRE.
        # This gives us a factor that seems to be:
        # 1/(4 π LRE_len_sq) exp(-3*LRE_len_sq / (2*SRE_len_sq))
        # although the first part would give us a non-unitless factor which is not acceptable. It probably has to be
        # normalized in some way, but for now just use the exponential part.
        # Edit: Actually, the first part is probably something like LRE_len/rc
        #
        # For LRE_len_sq == SRE_len_sq, this gives us exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 1/e = 0.22.
        # For example 1, this will give us:
        # exp(-3*LRE_len_sq / (2*SRE_len_sq)) = 0.02.
        LRE_factor = math.exp(-3*LRE_len_sq / (2*SRE_len_sq))

        gamma_corr = 1
        if LRE_len_sq > SRE_len_sq:
            # Domains can reach, but requires the SRE to extend beyond the mean squared end-to-end distance.
            # Should probably be approximated with some continuous function.
            # gamma = (3/2)*gamma_corr;
            gamma_corr += 0.5
        # Mean end-to-end squared distance between the two domains, aka E_r_sq:
        # Mean end-to-end squared distance.
        # We already have the squared length, Nᵢbᵢ², so we just need to sum:

        mean_sq_ee_dist = LRE_len_sq + SRE_len_sq       # unit of nm

        ## Regarding "effective volume" vs P_loop/P_v0: ##
        # Using "effective volume" is more intuitive than P^rc_loop/P^rc_v0.
        # However, stricly speaking we ARE considering reactant proximity/localization/collision probability
        # relative to standard conditions (1 M) using probability distribution function of the loop-connected reactants.
        # This happens to reduce to something that looks like P_loop/P_v0 = v_0/v_eff,
        # where we can use the mean squared end-to-end distance (mean_sq_ee_dist aka E[r²]) to calculate v_eff.
        # However, this is just an approximation valid only for E[r²] >> rc (rc = critical interaction distance),
        # which happens to reduce to a simple expression for P_loop and an expression for P_loop/P_v0 that does not
        # depend on rc. In general, especially if rc is significant, P_loop/P_v0 could very well depend on rc!
        # There is also no general guarantee that the expression for P_rel = P_loop/P_v0 will reduce to something with
        # an easily-identifiable v_eff subexpression. In that case we can simply define v_eff as
        #   v_eff = v_0 * P_v0(rc) / P_loop(rc) = v_0 / P_rel(rc)       # P_rel < 1 (unless they are really close)
        # although in practice we would not need to calculate v_eff, we would just use
        #   activity = P_v0(rc) / P_loop(rc) × M

        effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2)
        #effective_volume = (2/3*math.pi*mean_sq_ee_dist)**gamma * 1e-24 # 1e-24 to convert nm3 to L.
        #activity = (1/N_AVOGADRO)*(1/effective_volume)
        # Using AVOGADRO_VOLUME_NM3 to avoid redundant conversions:
        activity = AVOGADRO_VOLUME_NM3/effective_volume_nm3
        if gamma_corr > 1:
            activity = activity**gamma_corr

        ## When/where to apply extra gamma? ##
        # Note: The extra gamma really should be applied to the activity, not just the effective volume:
        # That is, we would always use exponent of 3/2 for calculating effective volume from mean_sq_ee_dist,
        # and apply the extra gamma to the unitless activity (v0/effective_volume) as:
        #   activity = (v0/effective_volume) ** gamma_corr
        # where gamma_corr = 1+x for gamma_org = 3/2+x
        # This is also what is effectively done in (Dannenberger, 2015), where γ is outside the expression
        #     ΔG = - R T γ ln(C/E[r2])
        # TODO: Check that this is also how Dannenberger actually does it when calculating k₊ in the java code.
        # To keep ΔG "reasonable" for large E[r2], (Dannenberger et al, 2015) adjusts C approximately as:
        #   C = 2.2e-18 m² γ - 2.7e-18 m² = 3.34 C0 γ - 4.0 C0,  C(γ=1.5) = C0
        # Where C0 = 3/(2π) v0**(2/3) = 6.7e-19 m²
        # This corresponds to increasing Avogadro volume, v0.

        # Note: The "activity" as calculated above appears on paper to be a unitless ratio (P_loop/P_v0).
        # However, since it is relative to standard conditions (1 M), we just have to implicitly multiply
        # the unitless ratio with "× 1 M" to get a proper molar activity.

        # TODO: Currently not accounting for bending or zipping energy.
        # TODO: Account for secondary (and tertiary?) loops
        return activity
