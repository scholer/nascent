# -*- coding: utf-8 -*-
##    Copyright 2015-2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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

Module with various utility functions for reactions.

"""

from .constants import HYBRIDIZATION_INTERACTION, STACKING_INTERACTION

def get_reaction_spec_pair(elem1, elem2, reaction_attr):
    if reaction_attr.reaction_type is HYBRIDIZATION_INTERACTION:
        reaction_spec_pair = frozenset((elem1.state_fingerprint(), elem2.state_fingerprint()))
    elif reaction_attr.reaction_type is STACKING_INTERACTION:
        reaction_spec_pair = frozenset(((elem1[0].state_fingerprint(), elem2[1].state_fingerprint()),
                                        (elem2[0].state_fingerprint(), elem1[1].state_fingerprint())))
    return reaction_spec_pair

