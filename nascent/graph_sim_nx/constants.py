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


"""

Package-wide constants.

"""

from __future__ import absolute_import, print_function, division
from math import pi


PHOSPHATEBACKBONE_INTERACTION = 'b'     # 1, b, pb, p, backbone,
HYBRIDIZATION_INTERACTION = 'h'         # 2, h, hh, hyb, hybridization, hybridized
STACKING_INTERACTION = 's'              # 3, s, ss, bs, stacked, stacking, basestacking

BACKBONE_INTERACTION_ALTERNATIVES = {1, 'pb', 'p', 'backbone'}
HYBRIDIZATION_INTERACTION_ALTERNATIVES = {2, 'hh', 'h', 'hyb', 'hybridization', 'hybridized'}
STACKING_INTERACTION_ALTERNATIVES = {3, 's', 'ss', 'bs', 'stacked', 'stacking', 'basestacking'}

REACTION_NAMES = {}
# REACTION_NAMES[<is_forming>]
REACTION_NAMES[True] = {PHOSPHATEBACKBONE_INTERACTION: 'backbone-nicking',
                        HYBRIDIZATION_INTERACTION: 'hybridization',
                        STACKING_INTERACTION: 'stacking'}
REACTION_NAMES[False] = {PHOSPHATEBACKBONE_INTERACTION: 'backbone-ligation',
                         HYBRIDIZATION_INTERACTION: 'de-hybridization',
                         STACKING_INTERACTION: 'un-stacking'}

interactions = [None, PHOSPHATEBACKBONE_INTERACTION, HYBRIDIZATION_INTERACTION, STACKING_INTERACTION]
valid_interactions = set(interactions)
interactions_dict = {}
for i, interaction in enumerate(interactions):
    interactions_dict[i] = interaction
    interactions_dict[interaction] = interaction


N_AVOGADRO = 6.022e23   # /mol
R = 1.987  # universal gas constant in cal/mol/K
LITER_TO_M3 = 0.001
NM3_TO_LITER = 1e-24
AVOGADRO_VOLUME_NM3 = 1/(N_AVOGADRO*NM3_TO_LITER)


#### DNA STRUCTURAL PARAMETERS: ####

ss_kuhn_length = 1.8 # nm
ss_rise_per_nt = 0.60 # nm per nt
ds_rise_per_bp = 0.34 # nm per bp
ds_kuhn_length = 100 # nm

# Distance between the 3p end of one base and the 5p end of the other stack.
HELIX_STACKING_DIST = ds_rise_per_bp
# Distance / backbone contour length from the 5p end of one domain to the 3p end of an unstacked downstream domain.
# A distance of 0.66954 nm gives an activity of 1 (M):
# Uh, calculated how? Currently, graph_manager calculates effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2),
# which for mean_sq_ee_dist = 0.66954**2 = 0.44828 gives effective_volume_nm3 = 0.90974 and an activity of 1.825 !
# effective_volume_nm3 = (2/3*math.pi*mean_sq_ee_dist)**(3/2) <=>
# 3/(2 * math.pi) * effective_volume_nm3**(2/3) = mean_sq_ee_dist,
# For activity of 1, effective_volume_nm3 = AVOGADRO_VOLUME_NM3, so:
# mean_sq_ee_dist = 3/(2 * math.pi) * AVOGADRO_VOLUME_NM3**(2/3) = 0.669546718045122
# - and HELIX_XOVER_DIST = sqrt(mean_sq_ee_dist) = 0.818
# Unless, of course, we include both ss-backbone connections in the path, and get:
# mean_sq_ee_dist = 2 * HELIX_XOVER_DIST**2, thus:
# HELIX_XOVER_DIST = sqrt(mean_sq_ee_dist/2) = 0.57859

HELIX_XOVER_DIST = 0.66954  # 1 # nm
HELIX_WIDTH = 2      # nm.




# "Cgamma prefactor":
loop_Cgamma0 = 3/(2*pi)*(LITER_TO_M3/N_AVOGADRO)**(2/3)



def loop_prefactor_corr(gamma):
    """
    Return a corrected loop prefactor (C_gamma), using the same model as Dannenberg et al,
    where the loop energy of a 18 nt single-stranded loop is constant for variable gamma exponents.
    The condition for this is simply "gamma*ln(C/E18) = const" for all (gamma, C), yielding
    the relation:
        C = E18 * (C0/E18)**(1.5/gamma)
    where E18 is the mean squared end-to-end distance of a 18 nt single-stranded loop, defined by:
        E18 = 18 nt * 0.6 nm/nt * 1.8 nm = 19.44 nm² = 1.944e-17 m²
    Where
        0.6 nm/nt is the contour length per nt and 1.8 nm is the Kuhn length of ssDNA.
    """
    E18 = 18 * ss_rise_per_nt * ss_kuhn_length * 1e-18 # 1 nm² = 1e-18 m²
    C0 = loop_Cgamma0
    C = E18 * (C0/E18)**(1.5/gamma)
    return C



# Connection types constants:
#PHOSPHATE_BACKBONE_INT = 1
#DUPLEX_HYBRIDIZATION_INT = 2
#STACKING_INTERACTION_INT = 3
