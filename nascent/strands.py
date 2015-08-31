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

# pylint: disable=C0103,C0301


"""

At the core, this library models DNA/RNA as a sequence of vectors with rotation.


For now, each nucleobase/nucleotide has a single vector, which points in the 5'->3' direction.
The vector does NOT represent the phosphate-to-phosphate displacement.
Instead, the vector is in the middle of the base-pair (base-pairing center: bpc - the point in
the middle between a base-pair).
For B-DNA, this means that the vector aligns approximately with the double-helix' center and direction.

Question: Where is the vector's base?
* The vectors base is at the base-pairing center (bpc). - Because this seems the most obvious.
* The vectors *end* is at the bpc - because in this way the nucleotide's vector can affect
    it's position! (Which seems natural).
* In the middle (vector goes through bpc) - because this represents the physical extend of the base.

Neither of these solutions alone are obviously ideal.
* Having the vector's base (start) at the bpc means the nucleotide's vector does not affect it's
    position. (Only subsequent ones).
* Having the vector's end at the bpc makes it hard to calculate (and the location of the first is unintuitive).
* Having the middle of the vector at the bpc makes all calculations much harder.

To resolve this, especially if we need a more accurate model, we can have two
vectors per nucleotide. (One vector before and one after).

For now, I go with what seems the simplest solution and have the vector's base at the bpc.


Question: Is the nucleobase in the xy-plane or yz-plane?
* If xy-plane, then the z-axis is the axis of the helix.
* Actually two questions: Is the helical axis along the strand's z-axis or x-axis?
* And does the nucleotide's axis (perpendicular to the nucleobase plane) follow?


"""

import numpy as np
import pyrr

from .nucleotides import Nucleotide



def sequence_gen(sequence):
    """
    For now the sequence generator is super simple,
    but we probably have to make it more complex in the future to account for mixed
    nucleotides, specified with dA, rA, lA, pA (for deoxy-, ribose-, locked-, peptide-, etc) specifications.
    """
    return (letter for letter in sequence)


#def calc_end_pos_dir_rot(nucleotides, position, direction, rotation):
#    for nuc in nucleotides:
#        position = nuc.position_after(position)
#        direction = nuc.direction_after(direction)
#        rotation = nuc.rotation_after(rotation)
#    return position, direction, rotation


class Strand(object):

    def __init__(self, position, orientation, sequence):

        # Position, direction and rotation of the 5' nucleotide:
        self.Position = position if position is not None else pyrr.Vector3()
        self.Orientation = orientation if orientation is not None else pyrr.Quaternion()
        #self.Rotation = rotation
        self.Sequence = sequence

        self.Locked = {}    # dict of locked attributes
        self.Locked_temp = {} # dict of temporarily locked attributes

        # Nucleotide objects:
        self.Nucleotides = None

        if sequence:
            self.build_from_sequence(sequence)


    def build_from_sequence(self, sequence=None):
        """ Build nucleotides from string sequence. """
        if sequence is None:
            sequence = self.Sequence
        if self.Nucleotides:
            raise AttributeError("")
        assert not self.Nucleotides
        position = self.Position
        orientation = self.Orientation
        #rotation = self.Rotation

        #for spec in sequence_gen(sequence):
        #    nuc = Nuclotide(position=position, orientation=orientation, spec=spec)
        #    self.Nucleotides.append(nuc)
        #    position = nuc.position_after()
        #    direction = nuc.direction_after()
        #    rotation = nuc.rotation_after()

        # Recursively:
        spec, rest = sequence[0], sequence[1:]
        nuc5p = Nucleotide(position=position, orientation=orientation, spec=spec)
        nucs = nuc5p.append_nucs_recursive(sequence=rest)
        self.Nucleotides = [nuc5p] + nucs
        return self.Nucleotides
