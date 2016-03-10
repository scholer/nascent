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

# pylint: disable=C0103,E1101

"""

As soon as this module is imported, vispy.app is imported, which may cause side-effects.

Only import this module if you intend to run the graphviewer application!

All other abstractions are in the other file.


Credits:
* Vispy Development Team
* github.com/campagnola, github.com/Eric89GXL, github.com/bollu,

Examples used:
* https://github.com/vispy/vispy/blob/master/examples/demo/gloo/galaxy.py (and galaxy/galaxy.py)
* https://github.com/vispy/vispy/blob/master/examples/demo/gloo/cloud.py
* https://github.com/vispy/vispy/blob/master/examples/demo/gloo/molecular_viewer.py
* https://github.com/vispy/vispy/blob/master/examples/demo/scene/force_directed_graph.py (2D)
** Using the "True Impostors" technique? http://cyrille.rossant.net/back-from-our-first-vispy-code-camp-at-esrf/
** Original: https://github.com/vispy/codecamp-esrf-2014/blob/master/gg/molecular_viewer.py

Alternatives:
* Mayavi, e.g. PDB protein rendering: http://docs.enthought.com/mayavi/mayavi/auto/example_protein.html

What kind of representations/visuals would be suitable for nodes?

 * Spheres (obviously)
 * Cubes/Boxes
 * PointVisual
 * Markers (vispy.scene.visuals.Markers - typically used for 3D scatter plots)
 * Some variation of true impostor?

What camera? Turntable or Arcball (trackball) ?
* To see the difference, go to https://www.jasondavies.com/maps/rotate/, go down to "Naive vs improved",
    view artica dead on and then turn side-ways.
* Trackball has more freedom regarding the "world-up" axis, but also makes it easy to "get lost".
* Since world-up is important to us, we use turntable camera.


"""

import sys
import os
import time
import numpy as np
sys.path.insert(1, '/Users/rasmus/Dev/repos-others/vispy')

import vispy
# from vispy.util.transforms import perspective
# from vispy.util import transforms
# from vispy import gloo
from vispy import app
from vispy import scene
from vispy.scene.canvas import SceneCanvas
from vispy import io
from vispy.color import Color
# from vispy.util.transforms import translate # Returns a translation-transformation matrix.
from vispy.visuals.transforms import STTransform  # STTransform(scale, translation)
# The vispy.visuals.transforms
# STTransform have .move() and .zoom() methods, while .transform and .scale are properties;
# MatrixTransform have generic .translate(), .rotate(), .scale() methods, applied to a generic .matrix property.
# vispy/visuals/transforms/linear.py


def increase_adjmatrix(adj_matrix, inc=1):
    """ Not sure if this is the optimal way to do it, but it works... """
    # adj_matrix = np.vstack((adj_matrix, np.empty((1, adj_matrix.shape[1]), dtype=bool)))
    # adj_matrix = np.hstack((adj_matrix, np.empty((adj_matrix.shape[0], 1), dtype=bool)))
    # The above works, but simpler to just use np.pad: (Pads all edges)
    adj_matrix = np.lib.pad(adj_matrix, (0, inc), 'constant')  # constant value, default to 0
    return adj_matrix


class GraphViewer():
    # Consider using class vispy.scene.SceneCanvas(app.Canvas, Frozen)

    def __init__(self):
        # setup initial width, height
        # app.Canvas.__init__(self, title='Nascent 3D Graph Viewer', keys='interactive', size=(1200, 800))
        self.canvas = SceneCanvas(title='Nascent 3D Graph Viewer', keys='interactive',
                                  size=(1200, 1200), show=True)
        self.main_view = self.canvas.central_widget.add_view()
        self.scene = self.main_view.scene
        self.main_view.bgcolor = '#efefef'
        self.main_view.camera = 'turntable'
        self.main_view.camera.up = '+z'  # scene.TurntableCamera(elevation=30, azimuth=30, up='+z')
        self.main_view.camera.fov = 60   # default=0 (for turntable) gives orthogonal projection
        self.main_view.padding = 100
        # Index by nodeid or sequential as list/array?
        self.nodeidx_by_id = {} # [nodeid] = idx ?
        # For now just using a dict:
        self.node_pos = np.zeros((0,3))    # [idx] = (x, y, z)  -- currently each node/Cube object has its own transform
        # self.node_obj = {}          # node object by nodeid
        self.node_obj_list = []     # [idx] = Cube
        # self.node_attrs = {}        # [nodeid] = node_attrs
        self.node_attrs_list = []   # [idx] = node_attrs
        # Edges:
        # (edge pos is given by source/target nodes)
        self.edge_idx = {}      # [source][target] = edge_idx
        self.edge_attrs = {}    # [source][target] = edge_attrs like networkx.DiGraph?
        self.edge_obj = {}     # [source][target] = Vispy arrow visuals objects used to mark edges (directional)
        self.edge_obj_list = [] # list of objects in same order as edge_list
        self.edge_attrs_list = []
        self.edge_list = []
        # edge_array: N x 2 numpy array, edge_array[N] = [i, j] connecting edge i with edge j.
        self.edge_array = np.empty((0, 2), dtype=np.uint32)
        # self.edges_pos = np.zeros((0,))  # Nope, use node_pos as edges pos.
        self.edges_line = scene.Line(pos=self.node_pos, connect=self.edge_array, color='blue', parent=self.scene)
        # Also: color, width, method='gl'/'agg', antialias=False/True
        # Regarding vispy LineVisual:
        # Does not *have* to be a contiguous line, can be a collection of discrete lines. Use 'connect' to control:
        #   connect='strip' will produce a contiguous line from the "pos" vertices;
        #   connect='segments' will interpret pos as pairs of start, end coordinates for each line.
        #   connect=np array can be used to specify exactly what vertices are connected. This is what you want.
        # Use line.set_data(pos=pos, color=color) to update

        # np.empty((0, 2), dtype='uint32')  # 'edges' in Luke's example above, pair of node idxs.
        # Well, having an adjacency matrix does make for some nice and easy numeric calculations, e.g.
        # https://github.com/vispy/vispy/blob/master/examples/demo/scene/force_directed_graph.py
        # edge_list is just a list of (idx1, idx2) tuples (more or less), not an actual adjacency matrix,
        # but is it being used to built `mask`, which is essentially an adjacency matrix.
        self.adj_matrix = np.empty((0, 0), dtype=bool)

        # self._timer = app.Timer('auto', connect=self.update, start=True)
        self._timer = app.Timer('auto', connect=self.visualize_layout, start=True)



    def add_node(self, nodeid, node_attrs):
        """ Add new node to the graph. """
        print("Adding new node: %s - %s" % (nodeid, node_attrs))
        assert nodeid not in self.nodeidx_by_id
        assert nodeid not in self.edge_attrs
        assert nodeid not in self.edge_obj

        try:
            color = Color(node_attrs['color'])
        except Exception as e:
            print("Could not extract node color from node attr:", e)
            color = Color("#3f51b5")


        # scene.visuals: Cube, Sphere, Plane, Tube, Arrow, Volume,
        # Scene visuals use objects loaded from vispy.visuals
        # Has both 'Box' and Cube.
        # Both inherits from vispy.visuals.visual.CompoundVisual
        # Box has height/depth/width segments and planes  (perhaps the only difference is Box has sub-divisions)
        # * examples/basics/visuals/box.py  vs  examples/basics/visuals/cube.py and examples/basics/scene/cube.py
        # For Cube, size=float gives a cube, size=[w, d, h] (x, y, z) gives a cuboid.
        cube = scene.visuals.Cube(size=(0.3, 0.2, 0.1), color=color, edge_color="black", parent=self.scene)


        if 'pos' in node_attrs:
            pos = np.array(node_attrs['pos'])
        else:
            pos = np.zeros((1,3))

        print(self.node_pos)
        print(pos)
        self.node_pos = np.vstack((self.node_pos, pos))  # Convenience function for concatenate. Not in-place.
        cube.transform = STTransform(scale=None, translate=pos)
        # STTransform._translate and ._scale are both vectors of length 4.
        self.nodeidx_by_id[nodeid] = len(self.node_obj_list)
        self.node_obj_list.append(cube)
        self.node_attrs_list.append(node_attrs)

        ## Prepare edge list:
        self.adj_matrix = increase_adjmatrix(self.adj_matrix)
        self.edge_idx[nodeid] = {}     # [source][target] = edge_idx
        self.edges_line.set_data(pos=self.node_pos)
        # self.edge_attrs[nodeid] = {}    # [source][target] = edge_attrs
        # self.edge_obj[nodeid] = {}     # [source][target] = edge_obj





        # cube pos can be changed with:
        # cube._mesh.vertices['position'] ?
        # visuals (and scene.visuals) can be transformed through their .transform property
        # Visuals actually have multiple transforms:
        # '.transform' is just the primary transforms.visual_transform.transforms[0]
        # self.transforms = TransformSystem()
        # TransformSystem._visual_transform = ChainTransform([NullTransform()])
        # We use a simple Scale+Translate transformation, which has a .move method that simply does:
        # self.translate = self.translate + as_vec4(move, default=(0, 0, 0, 0))

        # Can be altered by multiplying by a transformation matrix.
        # For rotation, use transformation matrix or quaternion.
        # vispy.util.quaternion/transforms for helper functions.




    def add_edge(self, source, target, edge_key, edge_attrs, directed=False):
        """ Add a new edge to the graph scene/canvas. """
        s_idx = self.nodeidx_by_id[source]
        t_idx = self.nodeidx_by_id[target]

        print("Adding edge from %s (%s) to %s (%s):" % (source, s_idx, target, t_idx))
        print(self.node_pos)
        points = self.node_pos[[s_idx, t_idx], :]
        print(points)
        color = edge_attrs.get('color', 'purple')
        # arrow = scene.visuals.Arrow(size=(0.3, 0.2, 0.1), color=color, edge_color="black", parent=self.main_view.scene)
        # arrow = scene.visuals.Tube(points=points, radius=0.1, color=color, #edge_color="black",
        #                            tube_points=6, parent=self.main_view.scene)
        # arrow = scene.visuals.Line(size=(0.3, 0.2, 0.1), color=color, edge_color="black", parent=self.main_view.scene)
        # "scene-node-aware LineVisual, is actually a collection of lines
        # We are using node_pos as pos for edges vertices, and edge_array to define line-connected vertices,
        # so it shouldn't be neccesary to "add" any new lines... hopefully...
        line = scene.visuals.Line()

        # by source/target ids:
        # self.edge_obj[source][target] = arrow
        # list (by idx):
        # self.edge_obj_list.append(arrow)

        self.adj_matrix[s_idx, t_idx] = True
        if not directed:
            self.adj_matrix[t_idx, s_idx] = True  # Directed vs undirected

        edge_idx = len(self.edge_list)
        n_edges_before = self.edge_array.shape[0]
        assert n_edges_before == edge_idx
        self.edge_list.append((s_idx, t_idx))
        # Try to update the edge_array in-place; otherwise we have to push a new edge_array copy to the GPU
        # self.edge_array.resize((n_edges_before+1,2), refcheck=False)
        # Argh, ValueError: cannot resize an array that references or is referenced
        # With refcheck=False... SegFault!
        # Bah, I say...
        # self.edge_array[edge_idx] = [s_idx, t_idx]
        # I guess I have to make a copy and update edges_line manually...
        self.edge_array = np.vstack((self.edge_array, [s_idx, t_idx]))
        print("edge_array after adding edge %s-%s" % (s_idx, t_idx))
        print(self.edge_array)

        # self.edges_line.set_data(pos=self.node_pos)
        # Make sure that edge vertices pos have been updated, otherwise the code will segfault.
        self.edges_line.set_data(connect=self.edge_array)

        # self.edges_line._changed['connect'] = True # To update connections..
        self.edge_idx[source][target] = edge_idx
        self.edge_attrs_list.append(edge_attrs)
        # self.edge_attrs[source][target] = edge_attrs


    def update_node_objs_pos(self):
        """ Update pos for all node/Cube objects to match the values in self.node_pos. """
        for pos, cube in zip(self.node_pos, self.node_obj_list):
            cube.transform.translate = pos

        # Also update edge objects
        # If using line visuals, this is as simple as:
        self.edges_line._changed['pos'] = True  # To update edge vertices positions
        # Alternatively, more correctly:
        # self.edges_line.set_data(pos=self.node_pos)
        # for eidx, (s_idx, t_idx) in enumerate(self.edge_list): #, self.edge_obj):
        #     # srcid, tgtid = self.nodeidx_by_id[s_idx], self.nodeidx_by_id[t_idx]
        #     edge_attrs = self.edge_attrs_list[eidx]
        #     # Maybe it is easier to create new tube rather than move the vertices?
        #     points = self.node_pos[[s_idx, t_idx], :]
        #     color = edge_attrs.get('color', 'purple')
        #     # new_arrow = scene.visuals.Tube(points=points, radius=0.1, color=color, #edge_color="black",
        #     #                                tube_points=6, parent=self.main_view.scene)
        #     # self.edge_obj_list[eidx] = new_arrow
        #     # Perhaps it is better to use lines instead?



    def force_directed_layout_2d(self, iterations=100):
        """ Apply force-directed layout in two dimensions (xy plane). """
        # dx = self.node_pos[:2, :] -
        pos = self.node_pos
        npts = pos.shape[0]
        # shape = (N, 3), i.e. N nodes with 3-elem position/coordinate vectors.
        # pos[:, np.newaxis, :] N x N matrix of position-vectors, each column with the same "array of pos vectors"
        # pos[np.newaxis, :, :] N x N matrix of position-vectors, each row with the same "array of pos vectors"
        while iterations > 0:
            iterations -= 1
            # Coordinate differences, only the first two (x, y) pos
            r = pos[:, np.newaxis, :2] - pos[np.newaxis, :, :2]
            # Calculate scalar distances and normalized node-to-node directionality unit vector
            dist = (r**2).sum(axis=2)**0.5     # pair-wise distance between all nodes
            dist[dist == 0] = 1.                # Prevent division-by-zero errors
            er = r / dist[..., np.newaxis]    # normalized node-node directionality vector

            # Repulsive force:  ke*(q1q2)/r^2
            # all points push away from each other
            # F_r = 0.02 / dist[..., np.newaxis]**4 * er  # Direct calculation+assignment
            F_r = 0.1 / dist[..., np.newaxis]**2 * er  # Direct calculation+assignment

            # Attractive spring force:
            # connected points pull toward each other with hookean force: F = ks*x,  ks = spring force constant
            # (Note: pulsed force may help the graph to settle faster)
            ks = 0.05
            #ks = 0.05 * 5 ** (np.sin(i/20.) / (i/100.))
            #ks = 0.05 + 1 * 0.99 ** i  # Exponentially decaying force to a constant value
            # ks = 0.5 + 1 * 0.999 ** iterations  # Exponentially decaying force to a constant value
            F_spring = ks * dist[..., np.newaxis] * self.adj_matrix[:, :, np.newaxis] * er
            # F_spring = ks * dist[..., np.newaxis]**2 * self.adj_matrix[:, :, np.newaxis] * er

            # Gravity force: F = G * m1m2/r^2 * dre,
            # where dre is normalized node-node directionality vector and G is gravitational constant
            Gm1m2 = 0.05
            F_g = Gm1m2 * er / dist[..., np.newaxis] # dist[..., np.newaxis]**2

            # All forces from all particles:
            F_all = F_spring - F_r # + F_g
            # Make sure points do not exert force on themselves:
            F_all[np.arange(npts), np.arange(npts)] = 0

            # dv = F/m*dt,  v1 = v0*dv,  dx = v*dt
            # But we don't have persistent velocities (i.e. start and end with v=0 for every step).
            # So we just use F*dt^2/m/4 for all particles, and say that dt^2/m/4 is some constant
            c_dt2m = 0.5 # 0.09   # In units of length/force
            F_sum = F_all.sum(axis=0)   # Sum all contributions from all other particles
            # print(F_sum)

            # Make sure no force is "too high":
            dx = np.clip(F_sum, -3, 3) * c_dt2m
            pos[:, :2] += dx


    def visualize_layout(self, event):
        """ Apply a single round of force-directed layout and update visuals. """
        self.force_directed_layout_2d(1)
        self.update_node_objs_pos()
        # time.sleep(0.1)


    ## OTHER THINGS:
    # visibility can be changed with .visibility property




    # From the "galaxy" and "molecular_viewer" examples:
    # create a new shader program
    # self.program = gloo.Program(VERT_SHADER, FRAG_SHADER, count=len(galaxy))

    # load the star texture
    # self.texture = gloo.Texture2D(load_galaxy_star_image(), interpolation='linear')
    # self.program['u_texture'] = self.texture

    # construct the model, view and projection matrices
    # self.view = transforms.translate((0, 0, -5))
    # self.program['u_view'] = self.view
    #
    # self.model = np.eye(4, dtype=np.float32)
    # self.program['u_model'] = self.model
    #
    # self.program['u_colormap'] = colors

    # w, h = self.size
    # self.projection = perspective(45.0, w / float(h), 1.0, 1000.0)
    # self.program['u_projection'] = self.projection
    #
    # # start the galaxy to some decent point in the future
    # galaxy.update(100000)
    # data = self.__create_galaxy_vertex_data()
    #
    # # setup the VBO once the galaxy vertex data has been setup
    # # bind the VBO for the first time
    # self.data_vbo = gloo.VertexBuffer(data)
    # self.program.bind(self.data_vbo)

    # setup blending
    # gloo.set_state(clear_color=(0.0, 0.0, 0.03, 1.0), depth_test=False, blend=True,
    #                blend_func=('src_alpha', 'one'))




    # def on_paint(self):
    #     pass



if __name__ == '__main__':
    viewer = GraphViewer()

    viewer.add_node(1, {'color': '#ff660033', 'pos': [3, 2, 0.5]})
    viewer.add_node(4, {'color': '#0000ff88', 'pos': [1, 2, 0]})
    viewer.add_node(5, {'color': '#00ff44cc'})

    viewer.add_edge(1, 4, 14, {'color': '#00ff0055'})
    viewer.add_edge(5, 4, 14, {'color': '#0000ff55'})

    # if the 'interactive' flag is set, you can interact with the program via an interactive terminal (e.g ipython)
    if sys.flags.interactive == 0:
        app.run()
