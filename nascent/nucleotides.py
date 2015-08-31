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


Transformation libraries:
* http://www.lfd.uci.edu/~gohlke/code/transformations.py.html
** Unofficial Github mirror: https://github.com/malcolmreynolds/transformations
** https://pythonhosted.org/MDAnalysis/documentation_pages/core/transformations.html
* https://code.google.com/p/gameobjects/
* Sage, http://www.sagemath.org/, https://vnoel.wordpress.com/2008/05/03/bye-matlab-hello-python-thanks-sage/

* math3d:
** http://git.automatics.dyndns.dk/?p=pymath3d.git
** Somewhat similar to Pyrr, uses numpy to provides intuitive wrappers for geometric computations.

* http://code.enthought.com/projects/mayavi/, http://docs.enthought.com/mayavi/mayavi/

* https://github.com/adamlwgriffiths/Pyrr  This actually seems pretty intuitive!
** Basically just provides wrappers around numpy for geometric computations.
** https://pyrr.readthedocs.org/en/latest/api_geometry.html


* https://pypi.python.org/pypi/transforms3d
** https://github.com/matthew-brett/transforms3d
** Uses a lot of gohlke's transformations module.

* pypi.python.org/pypi/qmath (old, obsolete?)
** From 2012, but no home page or repository.
* pypi.python.org/pypi/se3
** https://github.com/ccorcos/se3
** New, from 2014. Also uses numpy, but much less mature thatn Pyrr/transforms3d/math3d
* http://cxc.harvard.edu/mta/ASPECT/tool_doc/pydocs/Quaternion.html
** Quaternion, from 2010.
* pypi.python.org/pypi/quaternion-algebra, bitbucket.org/sirex/quaternion
** From 2012.
* pypi.python.org/pypi/Pyternion
** github.com/philip-peterson/pyternion
** Uses standard python math module, no numpy.
* github.com/olsoneric/pedemath
* https://github.com/zonca/quaternionarray, pypi.python.org/pypi/quaternionarray
** From 2010,
* pypi.python.org/pypi/euclid, and pypi.python.org/pypi/euclid3 (py3k)
** Super old, from 2006.


2D only libs:
* Affine, https://pypi.python.org/pypi/affine, https://github.com/sgillies/affine  (ONLY 2D!)
* http://toblerity.org/shapely/project.html

Stackoverflow:
* http://stackoverflow.com/questions/16083258/python-implementation-of-3d-rigid-body-translation-and-rotation

Educational:
* http://nghiaho.com/?page_id=671
* http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/Transformation3D.ipynb
** Very educational notebook with examples.
* http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/MARBLE/high/pose/express.htm
* http://planning.cs.uiuc.edu/node99.html
* http://en.wikipedia.org/wiki/Transformation_matrix, http://en.wikipedia.org/wiki/Transformation_(function),
* http://demonstrations.wolfram.com/Understanding3DTranslation/


DNA python code:
* http://www.lfd.uci.edu/~gohlke/code/dnacurve.py.html


3D plotting:
* http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
** ()



Other code:
* ~/ipython-notebooks/3D_plot_test.ipynb


"""


import numpy as np
from numpy.core import pi
import pyrr
from pyrr import Vector3 #, Quaternion, Matrix44,
#from matplotlib import pyplot

# Units:
A = Angstrom = 0.1
nm = Nanometer = 1
degree = pi/180

import random


def make_operation_hash():
    """ Cheap way to make a random hash. """
    return random.randint(0, 2**32)



class Nucleotide(object):

    """
    How to keep track of world position and orientation?
    a) Using a 4x4 transformation matrix
    b) Using a 3x3 transformation matrix plus a vector for the world position
        A 3x3 transform includes rotation+shear+reflection, but no translation.
    c) Using a position vector plus orientation quaternion.
    http://en.wikibooks.org/wiki/GLSL_Programming/Applying_Matrix_Transformations

    Considering the coordinate frame, you may want to refer to the current standards, e.g.
        Olson, (...), Berman's "standard reference frame" from J.Mol.Bio 2001. (10.1006/jmbi.2001.4987)
        * This seems to be focusing on the nucleobase and the base-pair. I'm not sure this standard's
          standard frame origin aligns with the helix center.
          NB: For B-DNA, shift and slide step parameters are both very close to zero: -0.02 Å and 0.23 Å
          Local helical parameters x- and y-displacement are also small, 0.05 and 0.02.

    Also:
    * http://x3dna.org/articles/seeing-is-understanding-as-well-as-believing

    """

    # Specify what is "up". I currently use z-up (to match the standard from [Berman 2001]).
    # Previously, I used y-up.
    helix_up = Vector3([0., 0., 1.])

    # Defaults for how to go to the next nucleotide. Could also be part of the strand model...
    params = {'rise': {'simple': 0.45 * Nanometer,
                       'ss': 0.42 * Nanometer,
                       'ds': 0.36 * Nanometer},
              'twist': {'simple': 0 * degree,
                        'ss': 10 * degree,
                        'ds': 36 * degree},
             }
    #standard_ss_translation = helix_up * 0.42 * Nanometer
    #standard_ds_translation = helix_up * 0.36 * Nanometer
    #standard_ss_rotation = pyrr.Quaternion.from_z_rotation(np.pi*0/180)
    #standard_ds_rotation = pyrr.Quaternion.from_z_rotation(np.pi*36/180)
    # The transform between one nucleotide to the other.
    # B-DNA has essentially no slide or roll, so just do a 180 rotation around the x-axis:
    standard_bp_transform = pyrr.Quaternion.from_x_rotation(pi)
    # [Berman 2001, Table 1, average of ATGC]
    # http://www.sciencedirect.com/science/article/pii/S0022283601949873
    atoms_coords = {'C1p': [-2.478, 5.375, 0.]}
    standard_bpc_to_C1p_vector = Vector3([-2.478, 5.375, 0.]) * Angstrom

    def standard_step_param(self, param, conformation=None, direction='3p'):
        """ Return standard step parameter. """
        if conformation is None:
            conformation = 'simple'
        if param in self.params:
            return self.params[param][conformation]
        if param == 'translation':
            return (-1 if '5' in direction else 1) * self.params['rise'][conformation] * self.helix_up
        if param == 'orientation':
            # For now we just assume the rotation is a simple twist around the z-axis:
            twist = (-1 if '5' in direction else 1) * self.params['twist'][conformation]
            return pyrr.Quaternion.from_z_rotation(twist)



    def __init__(self, position, orientation, strand=None, spec='N', conformation='simple'):

        self.Strand = strand
        self.Spec = spec
        # Base position, direction and rotation of the nucleotide:
        # Note: These could be calculated on-the-fly from local values.
        self.Position = position    # The position of the base (currently defined as the base-pair center, bpc)
        # Currently, the base's direction is the same as it's vector.
        #self.World_direction = direction  # Vector of translation to the next bpc.
        self.Orientation = orientation    # Orientation quaternion.
        # (x, y, z) Vector translation from the base of this nucleotide to the base of the next.
        # Helix axis is z, unit is nm:

        #if vector_trans is None:
        #    vector_trans = (0, 0, 0.36) if basepaired else (0, 0, 0.40)
        #if vector_rot is None:
        #    vector_rot = 36 if basepaired else 25
        #self._next_nuc_transform = next_nuc_transform

        self.Conformation = conformation
        self.Basepaired_with = None
        # When moving parts (strands or nucleotides), move whichever
        # has the lowest fixation level.
        self.fixation_level = 0     # 0: Free to move. (Generally affects translation).
        self.fixed_orientation = 0     # 0: Free to move

        self.Locked = {}    # dict of locked attributes
        self.Locked_temp = {} # dict of temporarily locked attributes

        #self.Pivot_local = None    # Local pivot is always the bpc
        #self.Pivot_world = None    # Global pivot should be avoided.

        self.Nuc3p = None
        self.Nuc5p = None

        self.Operation_hashes = set()

    @property
    def Is_basepaired(self):
        """ Returns True if nucleotide has a (watson-crick) base-pair partner. """
        return self.Basepaired_with is not None

    @property
    def World_direction(self):
        """ Just the unit vector of vector_trans. """
        return self.Orientation * self.helix_up

    def next_nuc_translation_default(self, conformation=None, direction='3p'):
        """ Translation from this nucleotide to the next. """
        #if self._next_nuc_transform:
        #    # This should be the translation part of the matrix,
        #    # as long as it is a 4x4 Matrix.
        #    # If we have a 3x3 matrix, it only includes the rotation part.
        #    return self._next_nuc_transform[3,0:3]
        return self.standard_step_param('translation', conformation=conformation, direction=direction)

    def next_nuc_position_default(self, conformation=None, direction='3p'):
        """ Translation from this nucleotide to the next. """
        return self.Position + self.next_nuc_translation_default(conformation=conformation, direction=direction)

    def next_nuc_orientation_default(self, conformation=None, direction='3p'):
        """ Return default orientation for the next nucleotide. """
        step_rotation = self.standard_step_param('orientation', conformation=conformation, direction=direction)
        return step_rotation * self.Orientation

    #
    #def next_nuc_transform(self):
    #    """
    #    The full 4x4 transformation matrix for the next nucleotide.
    #    Q: Is this the transformation matrix between this and the next,
    #    or the 'world' tranform?
    #    World is generally 'position' + 'orientation'.
    #    """
    #    #if self._next_nuc_transform:
    #    #    return self._next_nuc_transform
    #    # We have to generate it from our own parameters:
    #    pass
    #
    #def next_nuc_position(self):
    #    """
    #    World position after applying vector.
    #    "Next" is the 3p nucleotide.
    #    """
    #    pass

    def get_neighbour(self, direction='3p'):
        """ Return neighbourh in direction. """
        if '5' in direction:
            return self.Nuc5p
        else:
            return self.Nuc3p


    def walk(self, direction, operation_hash=None, extend=None, follow_bp=True,
             strand_blacklist=None, visited_set=None, include_self=True, bp_strand_blacklist=None):
        """
        In theory, this should just be a generator over all nucleotides.
        An alternative to having a operation_hash would be to add
        all visited nucleotides to a temporary set.

        Note: The order of returned nucleotides cannot currently be guaranteed.
        Also, this might produce funny comparisons for pseudo-knot structures.
        """
        # Break conditions:
        if extend is not None:
            if extend <= 0:
                return
            extend -= 1
        if visited_set is not None:
            if self in visited_set:
                return
            visited_set.add(self)
        if strand_blacklist is not None and self.Strand in strand_blacklist:
            return
        if operation_hash is not None:
            if operation_hash in self.Operation_hashes:
                return
            self.Operation_hashes.add(operation_hash)

        # Yields
        if include_self:
            yield self
        if direction in ('both', '3p'):
            yield from self.Nuc3p.walk('3p', operation_hash=operation_hash, extend=extend,
                                       follow_bp=follow_bp, strand_blacklist=strand_blacklist,
                                       bp_strand_blacklist=bp_strand_blacklist,
                                       visited_set=visited_set, include_self=True)
        if direction in ('both', '5p'):
            yield from self.Nuc5p.walk('5p', operation_hash=operation_hash, extend=extend,
                                       follow_bp=follow_bp, strand_blacklist=strand_blacklist,
                                       bp_strand_blacklist=bp_strand_blacklist,
                                       visited_set=visited_set, include_self=True)
        if follow_bp:
            if follow_bp is True:
                bp_walk_dir = 'both'
            elif follow_bp == 'continue':
                bp_walk_dir = '3p' if '3' in direction else '5p'
            else:
                bp_walk_dir = follow_bp
            nuc = self.Basepaired_with
            if nuc.Strand in strand_blacklist or nuc.Strand in bp_strand_blacklist:
                return
            yield from nuc.walk(bp_walk_dir, operation_hash=operation_hash, extend=extend,
                                follow_bp=follow_bp, strand_blacklist=strand_blacklist,
                                bp_strand_blacklist=bp_strand_blacklist,
                                visited_set=visited_set, include_self=True)


    def gen_fixation_level(self, direction, **kwargs):
        """ Generator of fixation level using walk. """
        nucs = self.walk(direction, **kwargs)
        return (nuc.fixation_level if nuc else 0 for nuc in nucs)

    def max_fixation_level(self, direction, **kwargs):
        """
        kwargs are passed on to self.walk(...)
        """
        #nucs = self.walk(direction, **kwargs)
        #max_level = 0
        #for nuc in nucs:
        #    if nuc and nuc.fixation_level > max_level:
        #        max_level = nuc.fixation_level
        # Alternatively:
        #max_level = max(map(lambda x: x.fixation_level if x else 0, nucs))
        max_level = max(self.gen_fixation_level(direction, **kwargs))
        return max_level

    def max_fixated_nuc(self, direction, **kwargs):
        """ Return the most fixated nucleotide, using walk. """
        nucs = self.walk(direction, **kwargs)
        keyfunc = lambda x: x.fixation_level if x else 0
        nuc = max(nucs, key=keyfunc)
        return nuc

    def cumulative_fixation_level(self, direction, **kwargs):
        """
        kwargs are passed on to self.walk(...)
        """
        #nucs = self.walk(direction, **kwargs)
        #cum_level = np.cumsum(map(lambda x: x.fixation_level if x else 0, nucs))
        cum_level = np.cumsum(self.gen_fixation_level(direction, **kwargs))     # pylint: disable=E1101
        return cum_level

    #def max_fixation_level(self, direction, operation_hash=None, extend=None, follow_bp=True,
    #                       strand_blacklist=None):
    #    if operation_hash in self.Operation_hashes:
    #        if operation_hash is None:
    #            print("WARNING: %s.Operation_hashes contains None!!" % self)
    #        return 0
    #    if extend is not None and extend <= 0:
    #        return self.fixation_level
    #    if operation_hash is None:
    #        operation_hash = make_operation_hash()
    #        self.Operation_hashes.add(operation_hash)
    #    nuc = self.get_neighbour(direction)
    #    this_strand_fixation_level = max((self.fixation_level, nuc.max_fixation_level(direction)))
    #    if self.Is_basepaired and follow_bp:
    #        other_dir = '3p' if '5' in direction else '5p'
    #        other_strand_fixation_level = self.Basepaired_with.max_fixation_level(other_dir, operation_hash, extend)
    #    return max((this_strand_fixation_level, other_strand_fixation_level))


    #def cumulative_fixation_level(self, direction, operation_hash=None, extend=None):
    #    if extend is not None and extend <= 0:
    #        return self.fixation_level
    #    if operation_hash in self.Operation_hashes:
    #        if operation_hash is None:
    #            print("WARNING: %s.Operation_hashes contains None!!" % self)
    #        return
    #    if operation_hash is None:
    #        operation_hash = make_operation_hash()
    #        self.Operation_hashes.add(operation_hash)
    #
    #    cum_level = self.fixation_level
    #    #return self.fixation_level + nuc.cumulative_fixation_level(direction)
    #    #nuc = self.get_neighbour(direction)
    #    if '3' in direction or direction == 'both':
    #        cum_level = self.Nuc3p.cumulative_fixation_level(direction='3p',
    #                                                         operation_hash=operation_hash,
    #                                                         extend=extend)
    #    if '5' in direction or direction == 'both':
    #        cum_level = self.Nuc5p.cumulative_fixation_level(direction='5p',
    #                                                         operation_hash=operation_hash,
    #                                                         extend=extend)



    def strand_fixation_level(self):
        """ Return the fixation level of this strand. """
        return self.cumulative_fixation_level(direction='both', include_self=True)

    def basepair_with(self, nuc):
        """
        To make a basepair:
        1) Align up bpc (moving the nuc with least fixation, or both if equal fixation).
        2) If either 5p or 3p is fixed, adjust helical rotation (twist) and translation (rise)
            to dsDNA conformation.
        3) Check that the base-pair rotation is correct. Move if required.

        For each step, when moving this or bp partner nucleotide, you also have
        to move the remaing nucleotides on the strand (on the least-fixed side).
        For translation, this is simple and should match up.
        For rotation, you probably have to rotate the remaining nucleotides
        around this nucleotide's bpc.
        Can this be done simply by modifying the orientation with the same rotation,
        and then rotating the position (separately)?

        Note: What do you do if this is already connected, e.g. you are making a pseudo-knot
        (you ALWAYS make lots of pseudo-knots when you make DNA nanostructures...)
        EFFECTS:
         1) You cannot perform simple translations. You have to move by rotating the different parts.
        Transformation by rotation:
         0) Objective is to move self and nuc together. Two ways to solve this:
             a) Start from current point and find rotations to get the two nucleotides together.
             b) Move the two nucleotides together and relax the position from there.
         1) Treat all helices as rigid rods.
         2) Use regions with 1 or more un-paired nucleotides as hinges.
            A good starting point might be to distribute the rotations equally on all un-paired nucleotides,
            (perhaps straightening ss-stretches at the ends).
        """
        this_strand_fixation = self.strand_fixation_level()
        other_strand_fixation = nuc.strand_fixation_level()
        # 1: Translate to aligh the bpc of this and the other nuc:
        translation = nuc.Position - self.Position
        # If we are already aligned, don't go through the process of translating the bpc:
        if any(translation):
            isequal = this_strand_fixation == other_strand_fixation
            if this_strand_fixation >= 2**32 and other_strand_fixation >= 2**32:
                # The strands are both fixed, either completely or at one point. Special case!
                raise NotImplementedError("Base-pairing two completely fixed strands is not yet supported.")
            elif isequal:
                # Move both strands half way: (You could scale by the total mass or inertia)
                self.translate(0.5 * translation)
                nuc.translate(-0.5 * translation)
            else:
                # Move the strand with least fixation:
                if self.strand_fixation_level() < nuc.strand_fixation_level():
                    # translate_nuc = self
                    self.translate(translation)
                else:
                    # translate_nuc = nuc
                    nuc.translate(-translation)

        # 2: Adjust helical parameters, transforming from ss to ds:
        # http://www.nature.com/ncomms/journal/v3/n6/fig_tab/ncomms1903_T1.html
        # Cases:
        # - This nucleotide xor the other is fixed
        # - None of the strands are fixed (yet).
        #if not this_strand_fixation and not other_strand_fixation:
        if not this_strand_fixation > self.fixation_level \
        and not other_strand_fixation > self.fixation_level:
            # None of the strands are fixed (yet).
            # This bp forms the basis of the helix.
            # Let's just say we transform the other nuc:
            #new_orientation = self.standard_step_param['orientation']['ds'] * self.Orientation
            if self.fixation_level > nuc.fixation_level:
                new_orientation = self.standard_bp_transform * self.Orientation
                nuc.rotate_to(new_orientation)
            else:
                new_orientation = nuc.standard_bp_transform * nuc.Orientation
                self.rotate_to(new_orientation)
        elif this_strand_fixation > other_strand_fixation:
            # This strand is fixed.
            fixation_3p = self.cumulative_fixation_level('3p')
            fixation_5p = self.cumulative_fixation_level('5p')
            if fixation_5p >= fixation_3p:
                # Adjust this nucleotide according to the 5' neighbour.
                new_position = self.Nuc5p.next_nuc_position_default(conformation='ds', direction='3p')
                new_orientation = self.Nuc5p.next_nuc_orientation_default(conformation='ds', direction='3p')
                self.translate_to(new_position, neighbours='3p')
                self.rotate_to(new_orientation, neighbours='3p')
            else:
                # Adjust this nucleotide according to the 3' neighbour.
                new_position = self.Nuc3p.next_nuc_position_default(conformation='ds', direction='5p')
                new_orientation = self.Nuc3p.next_nuc_orientation_default(conformation='ds', direction='5p')
                self.translate_to(new_position, neighbours='5p')
                self.rotate_to(new_orientation, neighbours='5p')
            # Adjust partner nucleotide:
            new_orientation = self.standard_bp_transform * self.Orientation
            nuc.rotate_to(new_orientation)
        else:
            # Other strand is fixed.
            fixation_3p = nuc.cumulative_fixation_level('3p')
            fixation_5p = nuc.cumulative_fixation_level('5p')
            if fixation_5p >= fixation_3p:
                # Adjust partner nucleotide according to its 5' neighbour.
                new_position = self.Nuc5p.next_nuc_position_default(conformation='ds', direction='3p')
                new_orientation = self.Nuc5p.next_nuc_orientation_default(conformation='ds', direction='3p')
                nuc.translate_to(new_position, neighbours='3p')
                nuc.rotate_to(new_orientation, neighbours='3p')
            else:
                # Adjust partner nucleotide according to its 3' neighbour.
                new_position = nuc.Nuc3p.next_nuc_position_default(conformation='ds', direction='5p')
                new_orientation = nuc.Nuc3p.next_nuc_orientation_default(conformation='ds', direction='5p')
                nuc.translate_to(new_position, neighbours='5p')
                nuc.rotate_to(new_orientation, neighbours='5p')
            # Adjust this nucleotide according to partner:
            new_orientation = nuc.standard_bp_transform * nuc.Orientation
            self.rotate_to(new_orientation)
        self.fixation_level = self.fixation_level or 1
        nuc.fixation_level = nuc.fixation_level or 1



    def translate(self, translation, neighbours='both',
                  move_bp_partner=True, operation_hash=None):
        """
        Translate this and other strand nucleotide strands to accomodate.
        Use an operation hash to ensure that the same operation is not applied multiple times
        because they are connected by multiple connections.
        The simplest case of this is moving a hair-pin: The base-pair partner
        should be translated, but it will also eventually be translated because of the backbone
        connection.

        Returns the number of nucleotides translated by the operation.
        """
        # Special case: If any nucleotides are completely fixed, we can attempt to translate by rotating
        # around the fixed nucleotide.
        if operation_hash in self.Operation_hashes:
            if operation_hash is None:
                print("WARNING: %s.Operation_hashes contains None!!" % self)
            return
        if operation_hash is None:
            operation_hash = make_operation_hash()
            self.Operation_hashes.add(operation_hash)

        self.Position += translation
        nmoved = 1

        if neighbours and ('both' in neighbours or '3p' in neighbours):
            if self.Nuc3p:
                nmoved += self.Nuc3p.translate(translation, neighbours='3p', operation_hash=operation_hash)
        if neighbours and ('both' in neighbours or '5p' in neighbours):
            if self.Nuc5p:
                nmoved += self.Nuc5p.translate(translation, neighbours='5p', operation_hash=operation_hash)
        if move_bp_partner and self.Is_basepaired:
            nmoved += self.Basepaired_with.translate(translation, neighbours='both', operation_hash=operation_hash)

        return nmoved


    def translate_to(self, position, neighbours='both'):
        """
        Translate this nucleotide to this position and make equivalent
        translation to remaining nucleotide's in strand.
        """
        translation = position - self.Position
        return self.translate(translation, neighbours)


    def rotate(self, rotation, pivot=None, neighbours='both',
               move_bp_partner=True, operation_hash=None):
        """
        Rotate this nucleotide by transform, and perform equivalent rotation
        to remaining nucleotides in the strand.
        Rotation is applied to orientation and to self.Position by rotating around pivot.
        If no pivot is given, bpc is used.
        Note: rotation must be a quaternion, which has axis and rotation (radians).
        """
        if operation_hash in self.Operation_hashes:
            if operation_hash is None:
                print("WARNING: %s.Operation_hashes contains None!!" % self)
            return
        if operation_hash is None:
            operation_hash = make_operation_hash()
            self.Operation_hashes.add(operation_hash)

        self.Orientation += rotation * self.Orientation
        # To rotate position around a pivot point: subtract the pivot point from position,
        # perform rotation, and add pivot point again:
        # http://www.euclideanspace.com/maths/geometry/affine/aroundPoint/
        if pivot:
            self.Position = rotation*(self.Position-pivot) + pivot
        else:
            self.Position = rotation * self.Position
            pivot = self.Position   # Later rotation use this nuc's bpc as pivot
        nmoved = 1

        if neighbours and ('both' in neighbours or '3p' in neighbours):
            if self.Nuc3p:
                nmoved += self.Nuc3p.rotate(rotation, pivot=pivot, neighbours='3p', operation_hash=operation_hash)
        if neighbours and ('both' in neighbours or '5p' in neighbours):
            if self.Nuc5p:
                nmoved += self.Nuc5p.rotate(rotation, pivot=pivot, neighbours='5p', operation_hash=operation_hash)
        if move_bp_partner and self.Is_basepaired:
            nmoved += self.Basepaired_with.rotate(rotation, pivot=pivot, neighbours='both', operation_hash=operation_hash)

        return nmoved


    def rotate_to(self, orientation, neighbours='both'):
        """
        Rotate this nucleotide to given orientation, and perform equivalent rotation
        to remaining nucleotides in the strand.
        """
        # new_orientation = rotation * orientation
        # <=>    rotation = new_orientation * orientation.inverse
        rotation = orientation * self.Orientation.inverse
        return self.rotate(rotation, neighbours)


    def append_nuc(self, direction='3p', nuc=None, spec='N'):
        """
        Append nucleotide to self.
        """
        assert (self.Nuc5p if '5' in direction else self.Nuc3p) is None
        position = self.standard_step_param('translation', direction=direction) + self.Position
        orientation = self.standard_step_param('orientation', direction=direction) * self.Orientation
        nuc = Nucleotide(position=position, orientation=orientation, spec=spec)
        if '5' in direction:
            self.Nuc5p = nuc
            nuc.Nuc3p = self
        else:
            self.Nuc3p = nuc
            nuc.Nuc5p = self
        return nuc

    def append_nucs_recursive(self, sequence, direction='3p'):
        """
        Note: The order of the returned nucs may be [3'-5'] if direction is '5p'.
        """
        #assert sequence
        if not sequence:
            return []
        spec, rest = sequence[0], sequence[1:]
        nuc = self.append_nuc(direction, spec=spec)
        if rest:
            return [nuc] + nuc.append_nucs_recursive(rest, direction)
        else:
            return [nuc]


    def __repr__(self):
        return "Nucleotide({}, {})".format(self.Position, self.Orientation)
