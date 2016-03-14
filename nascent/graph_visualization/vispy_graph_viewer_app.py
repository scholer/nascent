# -*- coding: utf-8 -*-
##    Copyright 2016 Rasmus Scholer Sorensen, rasmusscholer@gmail.com
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
import json
import numpy as np
sys.path.insert(1, '/Users/rasmus/Dev/repos-others/vispy')
if sys.platform.startswith('win'):
    from time import clock as system_time
    cstart = system_time.clock()
else:
    from time import time as system_time


import vispy
# from vispy.util.transforms import perspective
# from vispy.util import transforms
# from vispy import gloo
from vispy import app
from vispy import scene
from vispy.scene.canvas import SceneCanvas
# from vispy import io   # use vispy.io to avoid confusion with std lib of same name
from vispy.color import Color
# from vispy.util.transforms import translate # Returns a translation-transformation matrix.
from vispy.visuals.transforms import STTransform  # STTransform(scale, translation)
# The vispy.visuals.transforms
# STTransform have .move() and .zoom() methods, while .transform and .scale are properties;
# MatrixTransform have generic .translate(), .rotate(), .scale() methods, applied to a generic .matrix property.
# vispy/visuals/transforms/linear.py

from nascent.graph_visualization.graph_layout import force_directed_layout, force_directed_layout_2d


def increase_adjmatrix(adj_matrix, inc=1):
    """ Not sure if this is the optimal way to do it, but it works... """
    # adj_matrix = np.vstack((adj_matrix, np.empty((1, adj_matrix.shape[1]), dtype=bool)))
    # adj_matrix = np.hstack((adj_matrix, np.empty((adj_matrix.shape[0], 1), dtype=bool)))
    # The above works, but simpler to just use np.pad: (Pads all edges)
    adj_matrix = np.lib.pad(adj_matrix, (0, inc), 'constant')  # constant value, default to 0
    # Alternatively, you could use ndarray.resize((n_edges_before+1,2) to resize in-place
    return adj_matrix


class GraphViewer():
    # Consider using class vispy.scene.SceneCanvas(app.Canvas, Frozen)

    def __init__(self, config=None, **kwargs):
        if config is None:
            config = kwargs
        else:
            config.update(kwargs)
        self.config = config
        # setup initial width, height
        # app.Canvas.__init__(self, title='Nascent 3D Graph Viewer', keys='interactive', size=(1200, 800))
        self.canvas = SceneCanvas(title='Nascent 3D Graph Viewer', keys='interactive',
                                  size=(1200, 1200), show=True)
        self.main_view = self.canvas.central_widget.add_view()
        self.scene = self.main_view.scene
        self.main_view.bgcolor = config.get('bgcolor', '#efefef')
        self.main_view.camera = config.get('camera', 'turntable')
        self.main_view.camera.up = config.get('camera.up', '+z')  # scene.TurntableCamera(elevation=30, azimuth=30, up='+z')
        self.main_view.camera.fov = config.get('camera.fov', 60)   # default=0 (for turntable) gives orthogonal projection
        self.main_view.padding = config.get('view.padding', 100)
        # Index by nodeid or sequential as list/array?
        self.nodeidx_by_id = {} # [nodeid] = idx ?
        # For now just using a dict:
        self.node_pos = np.zeros((0, 3), dtype=float) # [idx] = (x, y, z) - each node/Cube object has its own transform
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
        # self.adj_matrix = np.empty((0, 0), dtype=bool)
        self.adj_matrix = np.empty((0, 0), dtype=float)
        # Do we want adj_matrix to also be weights? (an edge with 0 weight is effectively the same as no edge)
        # Or do we want a separate self.edges_weights data structure?
        # self.edges_weights = np.empty_like(self.adj_matrix)

        # self.layout_method = force_directed_layout
        self.node_grid_enabled = config.get("node_grid_enabled", True)   # Set to False if you are not using node_grid.
        self.node_grid_force = config.get("node_grid_force", np.array([1.], dtype=float))  # Disable node_grid_forces.
        # self.node_grid_force = np.array([1.])  # Use the same node_grid force scalar constant in all directions
        # self.node_grid_force = np.array([0.1, 1.0]) # Apply different forces in x and y:
        # Use node_grid_force = np.array([0.2, 0.]) to only apply node_grid forces on x-axis.
        self.node_grid_force_per_node = config.get("node_grid_force_per_node")
        # set node_grid_force_per_node=True if you want to have per-node force constants.
        # (the primary effect is that node_grid_force will be re-sized automatically when new nodes are added).
        # If re-assigning self.node_grid_force, make sure to update self.force_params!
        # TODO: Consider making force_params a dict to make it easier to update, and pass force_params.values()

        self.force_params = [
            # scale, exponent, use_adj_matrix, use_node_grid,  F = scale * dist**exponent
            [-0.1, -2, False, False],                   # Repulsion force
            [0.05, 1, True, False],       # scalar edge spring force, use floating point adj_matrix for edge weight
            # [self.edges_weights, 1, True, False],       # Per-edge spring force (if boolean adj_matrix)
            # [self.node_grid_force, 1, False, True],     # node_grid force
        ]
        print("force_params:", self.force_params, sep='\n')

        # Node grid: Direct certain nodes towards certain positions.
        self.node_grid = np.empty_like(self.node_pos)
        # node_grid_force can be either a scalar, or a per-node N x 2 array (for N nodes, doing )

        self.graph_stream = None
        self.graph_stream_stepping_gen = None
        self.sleep_until = None
        self.render_filename_fmt = config.get('render_filename_fmt', "graphviewer_{i}_{time:0.02f}.png")
        self.render_file_enumeration = 0


        ### Timers: ###
        self._app_start_time = system_time()
        # self._timer = app.Timer('auto', connect=self.update, start=True)
        # self._timer = app.Timer('auto', connect=self.visualize_layout, start=True)
        self.refresh_rate = 0.05 # default='auto'==1/60
        self._timer = app.Timer(self.refresh_rate, connect=self.read_stream_and_update, start=True)
        print("Timer type:", type(self._timer))
        # print("Timer.events.ignore_callback_error:", self._timer.events.ignore_callback_errors)  # Doesn't work, bug.
        print("Timer.events.ignore_callback_error:", self._timer.events._ignore_callback_errors)
        # Render timer: Renders the canvas and saves to file in 2 seconds interval:
        self._render_timer = app.Timer(2.0, connect=self.render_and_save, start=False, iterations=10)
        # app.Timer.events is an EmitterGroup(EventEmitter) with three event emitters: start, stop and timeout.
        for timer in (self._timer, self._render_timer):
            timer.events.print_callback_errors = "always"
            timer.events.ignore_callback_errors = False
        # iterations=20, )
        # Note: Setting iterations will not stop the app, it will just stop the timer event from triggering.



    def render(self):
        """ Return data from canvas.render() """
        return self.canvas.render()


    def render_and_save(self, event, fnfmt=None):
        """ Render canvas and save to file. """
        if fnfmt is None:
            fnfmt = self.render_filename_fmt
        self.render_file_enumeration += 1
        filename = fnfmt.format(
            i=self.render_file_enumeration,
            time=system_time()-self._app_start_time)
        image = self.render()
        vispy.io.write_png(filename, image)
        print("Rendered image saved to file:", os.path.abspath(filename))


    def get_node_visual_attrs(self, node_attrs):
        try:
            color = Color(node_attrs['color'])
        except Exception as e:
            print("Could not extract node color from node attr:", e)
            color = Color("#3f51b5")
        # if 'scale' in node_attrs:
        #     scale = node_attrs['shape']
        # elif 'size' in node_attrs:
        #     scale = node_attrs['size']**(1/3)
        scale = node_attrs.get('shape', node_attrs.get('size', 1)**(1/3))
        return color, scale


    def add_node(self, nodeid, node_attrs):
        """ Add new node to the graph. """
        print("Adding new node: %s - %s" % (nodeid, node_attrs))
        assert nodeid not in self.nodeidx_by_id
        assert nodeid not in self.edge_attrs
        # assert nodeid not in self.
        assert len(self.nodeidx_by_id) == len(self.node_attrs_list) == len(self.node_obj_list) == len(self.node_pos)
        if self.node_grid_enabled:
            assert len(self.node_pos) == len(self.node_grid)
        new_node_idx = len(self.node_obj_list)

        try:
            color = Color(node_attrs['color'])
        except ValueError as e:
            print("Could not extract node color from node attr:", e)
            color = Color("#3f51b5")
        scale = node_attrs.get('shape', node_attrs.get('size', 1)**(1/3))


        # scene.visuals: Cube, Sphere, Plane, Tube, Arrow, Volume,
        # Scene visuals use objects loaded from vispy.visuals
        # Has both 'Box' and Cube.
        # Both inherits from vispy.visuals.visual.CompoundVisual
        # Box has height/depth/width segments and planes  (perhaps the only difference is Box has sub-divisions)
        # * examples/basics/visuals/box.py  vs  examples/basics/visuals/cube.py and examples/basics/scene/cube.py
        # For Cube, size=float gives a cube, size=[w, d, h] (x, y, z) gives a cuboid.
        cube = scene.visuals.Cube(size=(0.14142, 0.1, 0.05), color=color, edge_color="black", parent=self.scene)


        if 'pos' in node_attrs:
            pos = np.array(node_attrs['pos'], dtype=float)
        else:
            pos = np.zeros((1, 3))   # Or, perhaps, np.random.rand(1, 3)?
            pos[:2] += np.random.rand(1, 2) * new_node_idx**0.5

        # print(pos)
        # self.node_pos = np.vstack((self.node_pos, pos))  # Convenience function for concatenate. Not in-place.
        # new_pos_size = (self.node_pos.shape[0]+1, self.node_pos.shape[1])
        self.increment_node_pos(increment=1)  # Increment node_pos and node_grid arrays by 1
        self.node_pos[new_node_idx, :len(pos)] = pos
        print("node_pos array and self.node_grid after adding %s" % pos)
        print(self.node_pos)
        print(self.node_grid)
        # Use self.node_pos.resize(new_pos_size) to resize in-place..
        # Use np.lib.pad(self.node_pos, padding) to get same effect as np.vstack

        cube.transform = STTransform(scale=None, translate=pos)
        # STTransform._translate and ._scale are both vectors of length 4.
        if scale != 1:
            cube.scale = scale
        self.nodeidx_by_id[nodeid] = len(self.node_obj_list)
        self.node_obj_list.append(cube)
        self.node_attrs_list.append(node_attrs)

        ## Prepare edge list:
        # self.adj_matrix = increase_adjmatrix(self.adj_matrix) # Done by self.increment_node_pos
        self.edge_idx[nodeid] = {}     # [source][target] = edge_idx
        self.edges_line.set_data(pos=self.node_pos)
        # self.edge_attrs[nodeid] = {}    # [source][target] = edge_attrs
        # self.edge_obj[nodeid] = {}     # [source][target] = edge_obj


        if 'grid_pos' in node_attrs:
            grid_pos = node_attrs['grid_pos']
            self.node_grid[new_node_idx, :len(grid_pos)] = grid_pos



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

    def increment_node_pos(self, increment=1):
        """
        Increment all numpy data structures which are depending on the number of nodes.
        It is an advantage that we are using C-style memory layout for node_pos, i.e.
            (node1_x, node1_y, node1_z, node2_x, node2_y, node2_z, ...)
        because then we can easily resize/extend node_pos and node_grid in-place.
        """
        # new_n_nodes = len(self.node_pos) + increment
        new_pos_size = (self.node_pos.shape[0]+increment, self.node_pos.shape[1])
        try:
            self.node_pos.resize(new_pos_size)
        except ValueError:
            print("Failed in-place node_pos resize, creating new...")
            # Pad one row at the end of the first axis (0, 1), pad zero values on the second axis (0, 0):
            self.node_pos = np.lib.pad(self.node_pos, ((0, increment), (0, 0)), mode='constant')
            # np.vstack can be used both to add a single row to the end of an array,
            # or concatenate two arrays of the same shape:
            # self.node_pos = np.vstack((self.node_pos, np.zeros(self.node_pos.shape[1])))  # Adds a single row
            # self.node_pos = np.vstack((self.node_pos, np.zeros((increment, self.node_pos.shape[1]))))  # Concatenate

        # try:
        self.node_grid.resize(new_pos_size)
        # except ValueError:
        #     # If you are really sure, you can disable refcheck:
        #     self.node_grid.resize(new_pos_size, refcheck=False)  # Fills with zeros
        #     # Otherwise, it is better to create a new ndarray:
        #     self.node_pos = np.lib.pad(self.node_pos, ((0, increment), (0, 0)), mode='constant')

        if self.node_grid_force_per_node:
            print("Must resize node_grid_force in-place, OR update self.force_params")
            if len(self.node_grid_force.shape) == 1:
                # A single scaler for each node, applied evenly in all directions
                pass
            elif len(self.node_grid_force.shape) == 2:
                self.node_grid_force.resize(new_pos_size)
            else:
                raise NotImplementedError("node_grid_force_per_node is not implemented for node_grid_force of shape %s."
                                          % self.node_grid_force.shape)

        # Don't try to resize adj_matrix in-place; you would have to move a lot of data to keep the edges
        # between the right nodes.
        # self.adj_matrix.resize((new_n_nodes, new_n_nodes))
        self.adj_matrix = np.lib.pad(self.adj_matrix, (0, increment), 'constant')  # constant value, default to 0




    def add_edge(self, source, target, edge_key, edge_attrs, directed=False):
        """ Add a new edge to the graph scene/canvas. """
        if source is None:
            source = edge_attrs.pop("source")
        if target is None:
            target = edge_attrs.pop("target")
        s_idx = self.nodeidx_by_id[source]
        t_idx = self.nodeidx_by_id[target]

        print("Adding edge from %s (%s) to %s (%s):" % (source, s_idx, target, t_idx))
        print("Node pos:")
        print(self.node_pos)
        points = self.node_pos[[s_idx, t_idx], :]
        print(points)
        # color = edge_attrs.get('color', 'purple')
        # arrow = scene.visuals.Arrow(size=(0.3, 0.2, 0.1), color=color, edge_color="black",
        #                             parent=self.main_view.scene)
        # arrow = scene.visuals.Tube(points=points, radius=0.1, color=color, #edge_color="black",
        #                            tube_points=6, parent=self.main_view.scene)
        # arrow = scene.visuals.Line(size=(0.3, 0.2, 0.1), color=color, edge_color="black", parent=self.main_view.scene)
        # "scene-node-aware LineVisual, is actually a collection of lines
        # We are using node_pos as pos for edges vertices, and edge_array to define line-connected vertices,
        # so it shouldn't be neccesary to "add" any new lines... hopefully...
        # line = scene.visuals.Line()
        # TODO: Figure out if we can have per-edge edge color and size without creating separate line collections...
        # TODO: You could have one set of edges for "forming" reactions and another set for "breaking" reactions?

        # by source/target ids:
        # self.edge_obj[source][target] = arrow
        # list (by idx):
        # self.edge_obj_list.append(arrow)
        weight = edge_attrs.get('weight', 1)   # can also be spring force constant...

        # adj_matrix_value = True if self.adj_matrix.dtype == bool else weight # Not needed, is converted.
        self.adj_matrix[s_idx, t_idx] = weight
        if not directed:
            self.adj_matrix[t_idx, s_idx] = weight  # Directed vs undirected

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
        if not directed:
            self.edge_idx[target][source] = edge_idx
        self.edge_attrs_list.append(edge_attrs)
        # self.edge_attrs[source][target] = edge_attrs


    def update_node_objs_pos(self):
        """ Update pos for all node/Cube objects to match the values in self.node_pos. """
        for pos, cube in zip(self.node_pos, self.node_obj_list):
            cube.transform.translate = pos

        # Also update edge objects
        # If using line visuals, this is as simple as: (Yes, this works just fine)
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


    def step_event(self, event_info):
        return event_info


    def change_node(self, node, node_attrs):
        """ Change/update node_attrs for a single node. """
        nodeidx = self.nodeidx_by_id[node]
        old_node_attrs = self.node_attrs_list[nodeidx]
        cube = self.node_obj_list[nodeidx]
        changed_attrs = {}
        for k, new_value in node_attrs:
            if k not in old_node_attrs or new_value != old_node_attrs[k]:
                changed_attrs[k] = new_value
                # if k in ('color', 'scale', 'size'):
        self.node_attrs_list[nodeidx].update(node_attrs)
        if 'color' in changed_attrs:
            cube.mesh.color = changed_attrs['color']
            # Edit, probably have to change cube.mesh.color (property setter, uses set_data(color=color)
        if 'scale' in changed_attrs:
            cube.transform.scale = changed_attrs['scale']
        elif 'size' in changed_attrs:
            cube.transform.scale = changed_attrs['size']**(1/3)
        if 'marker_size' in changed_attrs:
            pass
        if 'marker_color' in changed_attrs:
            pass


    def change_edge(self, source, target, edge_key, edge_attrs, directed=None):
        """ Change/update edge attrs for a single edge. """
        if source is None:
            source = edge_attrs.pop("source")
        if target is None:
            target = edge_attrs.pop("target")
        s_idx = self.nodeidx_by_id[source]
        t_idx = self.nodeidx_by_id[target]
        edge_idx = self.edge_idx[source][target]
        old_edge_attrs = self.edge_attrs_list[edge_idx]
        # cube = self.node_obj_list[nodeidx]
        changed_attrs = {}
        for k, new_value in edge_attrs:
            if k not in old_edge_attrs or new_value != old_edge_attrs[k]:
                changed_attrs[k] = new_value

        print("Changing edge from %s (%s) to %s (%s):" % (source, s_idx, target, t_idx))
        # I would like to be able to change:
        # - edge width
        # - edge color
        # -


    def parse_graphstreaming_msg(self, msg):
        return json.loads(msg)


    def read_graph_stream(self, stream):
        """
        Returns a generator that parses a stream of graph streaming events.
        The generator will return every time a "st" graph steaming step event is encountered in the stream.
        """
        # stream is a file-like object where each next(stream) gives an event
        # <event>      ::= ( <an> | <cn> | <dn> | <ae> | <ce> | <de> | <cg> | <st> | <cl> ) ( <comment> | <EOL> )
        event_methods = {
            "an": self.add_node,
            "cn": self.change_node,
            #"dn": self.delete_node,
            "ae": self.add_edge,
            "ce": self.change_edge,
            #"de": self.delete_edge,
            # "cg": self.change_graph,
            "st": self.step_event, # set time / clock tick
        }
        for msg in stream:
            msg = msg.strip()
            if not msg:
                continue
            print("Parsing msg:", msg)
            # msg is a string in the form of:
            # {"an": {<node id>: {"weight": 1}}}, {"ae": {edge_key: {"source", "target", "weight": 1}}}
            event = self.parse_graphstreaming_msg(msg)
            print(" processing event:", event)
            for event_type, event_info in event.items():
                method = event_methods[event_type]
                # res = method(event_info)
                if event_type == "st":
                    yield event_info
                else:
                    # All other events are a one or more dicts:
                    if event_type[1] == "n":
                        # Node events:
                        for nodeid, attrs in event_info.items():
                            method(nodeid, attrs)
                    elif event_type[1] == "e":
                        for edge_key, attrs in event_info.items():
                            source = attrs.pop('source')
                            target = attrs.pop('target')
                            method(source, target, edge_key, attrs)
                    else:
                        raise ValueError("Could not understand %s event: %s", (event_type, event_info))



    def read_stream_and_update(self, event):
        """
        Read stream and update.
        """
        # This will fully exhaust the stream before passing control back to the event loop:
        # if self.graph_stream is not None:
        #     for step in self.read_graph_stream(self.graph_stream):
        #         print("(read_stream_and_update) encountered step:", step)
        #         print("(read_stream_and_update) laying out graph...")
        #         # self.force_directed_layout_2d()
        #         self.visualize_layout(None)
        #         print("(read_stream_and_update) sleeping...")
        #         time.sleep(step)
        # alternatively, read stream one step at a time and pass control back to event loop after each step:
        now = system_time()
        if self.sleep_until is not None and now < self.sleep_until:
            ## TODO: Consider having a separate timer running visualize_layout
            self.visualize_layout(None)
            return
        if self.graph_stream_stepping_gen is None and self.graph_stream is not None:
            self.graph_stream_stepping_gen = self.read_graph_stream(self.graph_stream)
        if self.graph_stream_stepping_gen:
            try:
                step = next(self.graph_stream_stepping_gen)
                print("(read_stream_and_update) sleeping for %s secs" % step)
                # time.sleep(step)
                self.sleep_until = now + step
            except StopIteration:
                pass
                # TODO: Is there a better alternative, e.g. just postpone the next call to read_stream?
        # print("(read_stream_and_update) stream input read (for now..)")
        # time.sleep(0.1)
        self.visualize_layout(None)




    def parse_graphstreaming_list(self, messages):
        """ Parse a list of graph streaming events: """
        # {event: {identifier: attrs}}
        for msg in messages:
            pass


    def force_directed_layout(self, iterations=10, dim=2):
        """ Apply force-directed layout in two dimensions (xy plane). """
        # dx = self.node_pos[:2, :] -
        # force_directed_layout(pos, adj_matrix, node_grid=None, iterations=100, dim=None,
        #                       force_params=None, check_params=False
        # print("force_params:", self.force_params, sep='\n')
        # print("node_grid:", self.node_grid, sep='\n')
        force_directed_layout(self.node_pos, self.adj_matrix,
                              node_grid=self.node_grid if self.node_grid_enabled else None,
                              iterations=iterations, dim=dim, force_params=self.force_params)


    def visualize_layout(self, event):
        """ Apply a single round of force-directed layout and update visuals. """
        self.force_directed_layout(1)
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
    # self._timer = app.Timer('auto', connect=self.read_stream_and_update, start=False)

    # viewer.add_node(1, {'color': '#ff660033', 'pos': [3, 2, 0.5]})
    # viewer.add_node(4, {'color': '#0000ff88', 'pos': [1, 2, 0]})
    # viewer.add_node(5, {'color': '#00ff44cc'})
    #
    # viewer.add_edge(1, 4, 14, {'color': '#00ff0055'})
    # viewer.add_edge(5, 4, 14, {'color': '#0000ff55'})

    import io
    # Note: JSON keys must be enclosed in double-quotes; ints or similar are not allowed as keys.
    # For numeric values, + is not an allowed prefix.
    test_stream = io.StringIO("""
{"an": {"1": {"size": 1,     "pos": [-1, 0, 0], "color": "red", "grid_pos": [0, 0]}}}
{"st": 1.0}
{"an": {"2": {"x_center": 1, "pos": [1, 0, 0], "color": "green", "grid_pos": [0.5, -0.5]}}}
{"an": {"3": {"x_center": 1, "pos": [0, -1, 0], "color": "blue", "grid_pos": [0.5, 1]}}}
{"st": 1.0}
{"an": {"4": {"x_center": 1, "pos": [0, 1, 0], "color": "yellow", "grid_pos": [1, 1.5]}}}
{"an": {"5": {"x_center": 1, "pos": [0, 0, -1], "color": "cyan", "grid_pos": [1, 0]}}}
{"an": {"6": {"x_center": 1, "pos": [0, 0, 1], "color": "magenta", "grid_pos": [1, -1]}}}
{"st": 3.0}

{"ae": {"1-2": {"source": "1", "target": "2", "weight": 1}}}
{"st": 1.0}
{"ae": {"1-4": {"source": "1", "target": "4", "weight": 1}}}
{"st": 1.0}
{"ae": {"1-5": {"source": "1", "target": "5", "weight": 1}}}
{"ae": {"2-4": {"source": "2", "target": "4", "weight": 1}}}
{"st": 1.0}
{"ae": {"5-4": {"source": "5", "target": "4", "weight": 1}}}
{"st": 1.0}
{"ae": {"3-4": {"source": "3", "target": "4", "weight": 1}}}
{"st": 1.0}
{"ae": {"5-6": {"source": "5", "target": "6", "weight": 1}}}
{"st": 1.0}
{"ae": {"5-3": {"source": "5", "target": "3", "weight": 1}}}
{"ae": {"1-6": {"source": "1", "target": "6", "weight": 1}}}
    """)
    viewer.graph_stream = test_stream
    # print("node_grid_force before resizing:", viewer.node_grid_force, sep='\n')
    # viewer.node_grid_force.resize((2,), refcheck=False)
    # print("node_grid_force after resizing:", viewer.node_grid_force, sep='\n')
    # viewer.node_grid_force[:] = [0.5, 0.8] # Only apply node_grid force in the x-direction (on x-coordinates).
    # print("node_grid_force after assignment:", viewer.node_grid_force, sep='\n')
    # if the 'interactive' flag is set, you can interact with the program via an interactive terminal (e.g ipython)
    # viewer.read_stream_and_update(None)
    if sys.flags.interactive == 0:
        app.run()
        # viewer.read_graph_stream(viewer.graph_stream)
        # pass
        # filename = 'test.png'
        # print("App run done; saving rendered image to file:", os.path.abspath(filename))
        # image = viewer.render()
        # vispy.io.write_png(filename, image)
        print("Done app.run() !")
    else:
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)

    print("Final node_pos:", viewer.node_pos, sep='\n')
    print("viewer.node_grid:", viewer.node_grid, sep='\n')
    print("node_pos - node_grid:", viewer.node_pos - viewer.node_grid, sep='\n')
