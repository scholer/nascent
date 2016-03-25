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
import numbers
import numpy as np
sys.path.insert(1, '/Users/rasmus/Dev/repos-others/vispy')
if sys.platform.startswith('win'):
    from time import clock as system_time
    cstart = system_time()
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

from nascent.graph_visualization.graph_layout import (
    force_directed_layout, force_directed_layout2, force_directed_layout_2d)


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
        self.n_events_processed = 0
        self.n_steps_processed = 0
        self.n_layout_iterations_count = 0
        # continuous_iterations: True=Keep track of iterations between calls to layout; False=Reset it for every call.
        self.continuous_iterations = True
        self.config = config
        self.auto_center_view_factor = 1.0  # Set to False or None to disable view centering..
        # setup initial width, height
        # app.Canvas.__init__(self, title='Nascent 3D Graph Viewer', keys='interactive', size=(1200, 800))
        self.canvas = SceneCanvas(title='Nascent 3D Graph Viewer', keys='interactive',
                                  size=(1400, 1000), show=True)
        self.main_view = self.canvas.central_widget.add_view()
        self.scene = self.main_view.scene
        self.main_view.bgcolor = config.get('bgcolor', '#efefef')
        camera_type = config.get('camera', 'turntable')
        self.main_view.camera = camera_type
        # 'arcball, 'turntable', 'panzoom', 'fly',
        # 'panzoom' is great for 2D projections, but not 3D. 'turntable' is easy to use but you cannot re-center.
        # 'arcball' is advanced and the most flexible, but hard to zoom and a bit confusing.
        # 'fly' is ...
        # scene.TurntableCamera(elevation=30, azimuth=30, up='+z')
        self.main_view.camera.up = config.get('camera.up', '+z')
        self.main_view.camera.fov = config.get('camera.fov', 30) # default=0 (for turntable) = orthogonal projection
        self.main_view.padding = config.get('view.padding', 0)
        if camera_type == 'turntable':
            # if 'camera.azimuth' in config:
            self.main_view.camera.azimuth = config.get('camera.azimuth', 90)
            # if 'camera.elevation' in config:
            self.main_view.camera.elevation = config.get('camera.elevation', 45)
            # if 'camera.distance' in config or True:
            self.main_view.camera.distance = config.get('camera.distance', 10)


        # Index by nodeid or sequential as list/array?
        self.nodeidx_by_id = {} # [nodeid] = idx
        self.nodeid_list = []
        # nodeid_list == [nodeid for idx, nodeid in sorted([v, k for k, v in self.nodeidx_by_id.items()])]
        # nodeidx_by_id == {nodeid: idx for idx, nodeid in enumerate(self.nodeid_list)}

        self.node_pos = np.zeros((0, 3), dtype=float) # [idx] = (x, y, z) - each node/Cube object has its own transform
        # self.node_obj = {}          # node object by nodeid
        self.node_obj_list = []     # [idx] = Cube
        # self.node_attrs = {}        # [nodeid] = node_attrs
        self.node_attrs_list = []   # [idx] = node_attrs
        self.node_labels_list = []  # [idx] = label
        self.enable_node_labels = config.get('enable_node_labels', True)

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
        ## TODO: BUG: Setting width also affects cubes (and ground): !!
        self.edges_line = scene.Line(pos=self.node_pos, connect=self.edge_array, parent=self.scene,
                                     color='blue', width=1, method='gl')  # method='agg' or 'gl' (default)
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


        ### LAYOUT parameters and forces: ###

        # If you do not like the original node positions, you can translate them when adding nodes:
        # node_pos_transform_input = (scale, translate) # e.g. ([0.5, 2, 0.1], [-4, 0, 0])):
        self.node_pos_transform_input = config.get("node_pos_transform_input", ([0.5, 2, 0.2], [-2, -1, 0]))

        self.force_directed_dimensions = 3   # 2 to only apply forces in xy plane

        # number of force iterations:
        self.force_directed_iterations_per_call = 3 # default
        self.force_directed_iterations_per_callback = 3 # when called as part of the regular event loop
        self.force_directed_iterations_per_step_cont = 0  # a step was processed, but we did not break
        self.force_directed_iterations_per_step_break = 0 # Every time we break out of the stream
        self.force_directed_iterations_per_msg = 0  # 1 msg can have multiple event types (but usually just one)
        self.force_directed_iterations_per_event_type = 0  # 1 msg can have multiple event types (but usually just one)
        # event_type includes steps, node/edge changes, which are typically not relevant to layout:
        self.force_directed_iterations_per_added_node = 10  # 1 msg can have multiple event types (but usually just one)
        self.force_directed_iterations_per_added_edge = 5  # 1 msg can have multiple event types (but usually just one)

        self.do_grid_calculation = True  # Set to None to detect automatically.
        self.force_dx_factor = 0.02 # How much to move each node for a given force.
        # Typical motion simulations uses dv = F/m*dt,  v1 = v0*dv,  dx = v*dt
        # But we don't have persistent velocities (i.e. start and end with v=0 for every step).
        # So we just use F*dt^2/m/4 for all particles (where force_dx_factor = dt^2/m/4)


        # Avoid re-creating a lot of numpy arrays all the time;
        # instead, reuse old ones and clear cache when the number of nodes changes:
        self.cache = {}


        # Node grid: Direct certain nodes towards certain positions.
        self.node_grid = np.empty_like(self.node_pos)
        self.pos_is_default_grid_pos = True # True = If no 'grid_pos' in node_attrs, use pos as grid_pos.


        # self.layout_method = force_directed_layout
        self.node_grid_enabled = config.get("node_grid_enabled", True)   # Set to False if you are not using node_grid.
        # node_grid_force_k can be either a scalar, or a per-node N x 2 array (for N nodes, doing ):
        # Make sure grid_k matches self.force_directed_dimensions!
        grid_k = config.get("node_grid_force_k")  # Disable node_grid_forces.
        if grid_k is None:
            grid_k = np.array([1.], dtype=float) # isotropic scaling
            grid_k = np.array([0.02, 0, 2.5]) # high in z, little in x, none in y.
        if len(grid_k) > self.force_directed_dimensions:
            grid_k = grid_k[:self.force_directed_dimensions]
        self.node_grid_force_k = grid_k

        # self.node_grid_force_k = np.array([1.])  # Use the same node_grid force scalar constant in all directions
        # self.node_grid_force_k = np.array([0.1, 1.0]) # Apply different forces in x and y:
        # Use node_grid_force_k = np.array([0.2, 0.]) to only apply node_grid forces on x-axis.
        self.node_grid_force_per_node = config.get("node_grid_force_per_node")
        # set node_grid_force_per_node=True if you want to have per-node force constants.
        # (the primary effect is that node_grid_force_k will be re-sized automatically when new nodes are added).
        # If re-assigning self.node_grid_force_k, make sure to update self.force_params!
        # TODO: Consider making force_params a dict to make it easier to update, and pass force_params.values()

        ## OLD FORCE PARAMS: (for the record)
        # self.force_params = [
        #     # scale, exponent, use_adj_matrix, use_node_grid,  F = scale * dist**exponent
        #     ## Repulsion forces:
        #     # [-0.1, -2, False, False],     # Repulsion force
        #     # [-0.6, -3, False, False],     # Repulsion force
        #     [-0.2, -2, False, False],       # Repulsion force
        #     ## Spring forces (edges)
        #     [0.5, 1, True, False],       # scalar edge spring force, use floating point adj_matrix for edge weight
        #     # [0.1, 1, True, False],       # scalar edge spring force, use floating point adj_matrix for edge weight
        #     # [self.edges_weights, 1, True, False],       # Per-edge spring force (if boolean adj_matrix)
        #     ## Node grid pos forces:
        #     # [self.node_grid_force_k, 1, False, True],     # node_grid force
        #     # [np.array([1.]), 1, False, True],     # node_grid force
        #     # [np.array([1., 0, 0]), 1, False, True],     # node_grid force in x-direction only..
        #     # [np.array([0.1, 0, 0.4]), 1, False, True],     # node_grid force in x and z directions..
        #     [grid_k, 1, False, True],     # node_grid force in x and z directions..
        # ]
        # New-style force_params with support for fluctuations
        # dist_enum: 1=node-node, 2=node-grid, 0=do-not-use-distances
        # mask_enum: 0=no-mask, 1=adj_matrix
        # fluctuations: (additive, scaling) functions or None to not add fluctuations.
        # As always, negative k = repulsive forces: F = k * d^a
        # Nice fluctuation distributions:
        additive = lambda shape, dim, i: 0.3*np.random.normal(size=(shape[0], dim))
        scale1 = lambda shape, dim, i: np.random.exponential(size=(shape[0], dim))
        # decaying scaling fluctuations: (towards 1)
        sdecay = lambda shape, dim, i: 1 + (0.9999**i * np.random.exponential(size=(shape[0], dim)))
        scale5 = lambda shape, dim, i: 5*np.random.exponential(size=(shape[0], dim))
        randunit1 = lambda shape, dim, i: 2*np.random.uniform(size=(shape[0], dim))
        randunit0 = lambda shape, dim, i: np.random.uniform(low=-0.5, high=0.5, size=(shape[0], dim))

        self.force_params = [
            # k,  a, dist, mask, fluctuations
            # Repulsive forces, e.g. electrostatic F = k / dist**2 == k * dist**(-2)
            [-0.20, -2, 1, 0, None],
            # [-0.20, -2, 1, 0, (None, sdecay)], # Scaled by randomly-fluctuating values
            # Attractive edge forces, e.g. spring force F = k * dist
            [0.500, +1, 1, 1, None],
            # [0.010, -1, 1, 0, None], # "Gravity" (but longer range, F=Gmm/r, not F=Gmm/r^2)
            # Grid forces:
            # [grid_k, 1, 2, 0, None], # node_grid force in x and z directions..
            # [grid_k, 1, 2, 0, (None, scale5)], # node_grid force in x and z directions..
            [grid_k, 1, 2, 0, (None, randunit1)], # node_grid force in x and z directions..
            # [grid_k, 1, 2, 0, (additive, None)], # node_grid force in x and z directions..
            # [grid_k, 1, 2, 0, (additive, scale)], # node_grid force in x and z directions..
            # Adding a randomly-fluctuating scaling factor to the grid forces produces a
            # really nice effect where nodees "look unhappy/unsettled" if they are far from their
            # preferred grid location.
        ]
        print("force_params:", self.force_params, sep='\n')



        ### VISUAL representation of objects: ###

        # Node things:
        self.state_time_max = 0 # Is used to normalize observed state partition functions (cummulated state time)
        # self.state_times = []  # edit: get dynamically using [attr.get(tau_cum) for attr in self.node_attrs_list]

        ## Edge visuals:
        self.edges_line_groups = {
            # group_name: {key_value1: [list-of-edge-lines], key_value2: [...], ...}
            'is_forming': {True: [], False: []}, # group by is_forming; should be either True or False
        }
        ## TODO: Implement different line groups, e.g. "forming" vs "breaking" reactions.
        # Key-funcs, returning a value to group by key func.
        self.edges_line_group_key_funcs = {
            # group_name: group-by-key-func
            "is_forming": lambda source, target, edge_attr: edge_attr.get('is_forming')
        }
        # If key_func returns None, do not add a line_group line for this edge.
        # Otherwise, add line to the line group.
        # By example, in add_edge:
        # for group_name, keyfunc in self.edges_line_group_key_funcs.items():
        #     group_by_value = keyfunc(source, target, edge_attr)
        #     if group_by_value is None: continue
        #     edge_line_group = self.edges_line_groups[group_name].setdefault(group_by_value, {})

        # Per-edge lines:
        self.per_edge_lines = []


        ### VISUAL ENVIRONMENT: ###
        self.scene_draw_ground = True
        self.scene_ground_obj = None
        if self.scene_draw_ground:
            self.draw_ground()


        ### Streaming: ###
        self.graph_stream = None
        self.graph_stream_stepping_gen = None
        self.sleep_until = None
        self.step_sleep_scalefactor = config.get('step_sleep_scalefactor', 0.01)
        self.step_sleep_minimum_cutoff = 1e-1 # If step is < 1e-3 secs, then do not sleep...
        self.stream_collect_steps_until = 10  # Avoid repeated steps, only continue when step sum > N seconds.
        self.stream_collected_step_sum = 0

        ### Rendering: ###
        self.render_filename_fmt = config.get('render_filename_fmt', "graphviewer_{i}_{time:0.02f}.png")
        self.render_file_enumeration = 0

        ### Timers: ###
        self._app_start_time = system_time()
        # self._timer = app.Timer('auto', connect=self.update, start=True)
        # self._timer = app.Timer('auto', connect=self.visualize_layout, start=True)
        self.refresh_rate = 1/1000 # 0.05 # default='auto'==1/60
        self._timer = app.Timer(self.refresh_rate, connect=self.read_stream_and_update, start=True)
        # print("Timer type:", type(self._timer))
        # print("Timer.events.ignore_callback_error:", self._timer.events.ignore_callback_errors)  # Doesn't work, bug.
        # print("Timer.events.ignore_callback_error:", self._timer.events._ignore_callback_errors)
        # Render timer: Renders the canvas and saves to file in 2 seconds interval:
        self._render_timer = app.Timer(2.0, connect=self.render_and_save, start=False, iterations=10)
        # app.Timer.events is an EmitterGroup(EventEmitter) with three event emitters: start, stop and timeout.

        if self.auto_center_view_factor:
            self._auto_center_view_timer = app.Timer(0.2, connect=self.center_view, start=True)

        for timer in (self._timer, self._render_timer, self._auto_center_view_timer):
            timer.events.print_callback_errors = "always"
            timer.events.ignore_callback_errors = False
        # iterations=20, )
        # Note: Setting iterations will not stop the app, it will just stop the timer event from triggering.



    def center_view(self, event, center_on=None, animate=True):
        """
        Center view on the specified position. If no position is given,
        center on the mean of all nodes.
        """
        # Note: translate is a 4-vector.
        if center_on is None:
            center_on = self.node_pos.mean(axis=0) # x,y,z mean for all nodes
        # if len(center_on) != 4:
        #     center_on.resize((4,))
        if animate and self.auto_center_view_factor:
            r = center_on - self.main_view.camera.center # self.main_view.transform.translate
            print("\n\n\ncenter_view: Currently %s from center_on %s" % (r, center_on))
            if True or r.sum() > 0.1:
                dx = self.auto_center_view_factor*r
                # self.main_view.transform.translate += dx
                self.main_view.camera.center += dx
                print("\ncenter_view: Moving main_view by %s to %s" % (
                    center_on, self.main_view.camera.center))
            # self.main_view.transform.translate += (
            #     self.auto_center_view_factor * (center_on - self.main_view.transform.translate))
        else:
            print("\ncenter_view: Centering directly on %s" % center_on)
            self.main_view.camera.center = center_on


    def draw_ground(self):
        """ Draw a ground underneath the graph. """
        size = self.config.get('scene_ground_size', (10, 14, 0.01))
        color = self.config.get('scene_ground_color', 'grey')
        scale = self.config.get('scene_ground_scale')
        translate = self.config.get('scene_ground_translate', [0., 0., -4])
        cube = scene.visuals.Cube(size=size, color=color, edge_color="black", parent=self.scene)
        cube.transform = STTransform(scale=scale, translate=translate)
        self.scene_ground_obj = cube


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
        self.nodeidx_by_id[nodeid] = new_node_idx
        self.nodeid_list.append(nodeid)
        self.node_attrs_list.append(node_attrs)

        try:
            color = Color(node_attrs['color'])
        except (KeyError, ValueError) as e:
            print("Could not extract node color from node attr:", e)
            color = Color("#3f51b5")


        if "pos" in node_attrs:
            pos = np.array(node_attrs["pos"], dtype=float)
            if self.node_pos_transform_input:
                pos_scale, pos_translate = self.node_pos_transform_input
                if pos_scale and pos_scale != 1:
                    pos *= pos_scale
                if pos_translate and pos_translate != 0:
                    pos += pos_translate
        else:
            # Center the new node around new_node_idx (in x/y plane):
            pos = np.zeros((3,))
            pos[:2] += np.random.normal(new_node_idx**0.5, scale=(1+new_node_idx)**0.5, size=(2,))
        if "offset_relative_to" in node_attrs:
            # [(0, None), (random.random(), reaction_spec_source_states_list), (0, None)]
            offset = [offset +
                      (self.node_pos[[self.nodeidx_by_id[str(node)] for node in node_list], i].mean()
                       if node_list is not None else 0)
                      for i, (offset, node_list) in enumerate(node_attrs["offset_relative_to"])]
            print("Offsetting new node %s pos %s by %s" % (nodeid, pos, offset))
            pos[:len(offset)] += offset
            print(" - New node pos after offset: %s" % pos)

        ## Expand data structures and add node pos: ##
        # self.node_pos = np.vstack((self.node_pos, pos))  # Convenience function for concatenate. Not in-place.
        # new_pos_size = (self.node_pos.shape[0]+1, self.node_pos.shape[1])
        self.increment_node_pos(increment=1)  # Increment node_pos and node_grid arrays by 1
        self.node_pos[new_node_idx, :len(pos)] = pos
        print("node_pos array and self.node_grid after adding %s" % pos)
        print(self.node_pos)
        print(self.node_grid)
        # Use self.node_pos.resize(new_pos_size) to resize in-place..
        # Use np.lib.pad(self.node_pos, padding) to get same effect as np.vstack

        ## Translate node cube visual object: ##
        if 'scale' in node_attrs:
            scale = node_attrs['scale']
        elif 'size' in node_attrs:
            print("Note: Using size for scale as (%s)^0.3" % node_attrs['size'])
            scale = node_attrs['size']**(0.33)
        else:
            scale = None
        if scale is not None:
            # print("Scaling node cube visual by %s" % scale)
            # time.sleep(4)
            if isinstance(scale, numbers.Number):
                scale = np.array([scale, scale, scale], dtype='float32')


        ### Add node visual object representation: ###
        # scene.visuals: Cube, Sphere, Plane, Tube, Arrow, Volume,
        # Scene visuals use objects loaded from vispy.visuals
        # Has both 'Box' and Cube.
        # Both inherits from vispy.visuals.visual.CompoundVisual
        # Box has height/depth/width segments and planes  (perhaps the only difference is Box has sub-divisions)
        # * examples/basics/visuals/box.py  vs  examples/basics/visuals/cube.py and examples/basics/scene/cube.py
        # For Cube, size=float gives a cube, size=[w, d, h] (x, y, z) gives a cuboid.
        cube_base_size = (0.14142, 0.1, 0.05)
        # cube_base_size = (0.14142/2, 0.1/2, 0.05/2)
        # cube_base_size = (0.05, 0.04, 0.01)  # too small..
        cube = scene.visuals.Cube(size=cube_base_size, color=color, edge_color="black", parent=self.scene)
        cube.transform = STTransform(scale=scale, translate=pos)
        # STTransform._translate and ._scale are both vectors of length 4.

        self.node_obj_list.append(cube)
        if True or self.enable_node_labels:
            # node_label = scene.Label(nodeid, rotation=0, color='black')
            # node_label.stretch = (0.1, 1)
            # scene.Label is intended for 2D display, so maybe use a TextVisual? (from vispy.visuals.text)
            # from vispy.visuals.text.text import TextVisual
            # from vispy.visuals import TextVisual  # Also here
            # node_label = TextVisual()
            # or maybe use the scene-enabled version?
            node_label = scene.Text(nodeid, pos=pos, parent=self.scene, font_size=72, anchor_y='bottom')
            self.node_labels_list.append(node_label)
            # defaults: text=None, color='black', bold=False, italic=False, face='OpenSans', font_size=12,
            #           pos=[0, 0, 0], rotation=0., anchor_x='center', anchor_y='center', font_manager=None
            # move around using text.pos property

        ## Prepare edge list:
        # self.adj_matrix = increase_adjmatrix(self.adj_matrix) # Done by self.increment_node_pos
        self.edge_idx[nodeid] = {}     # [source][target] = edge_idx
        self.edges_line.set_data(pos=self.node_pos)
        # self.edge_attrs[nodeid] = {}    # [source][target] = edge_attrs
        # self.edge_obj[nodeid] = {}     # [source][target] = edge_obj

        grid_pos = node_attrs.get('grid_pos')
        if grid_pos is None:
            if self.pos_is_default_grid_pos:
                grid_pos = pos  # pos has already been transformed by node_pos_transform_input if needed
        elif self.node_pos_transform_input:
            pos_scale, pos_translate = self.node_pos_transform_input
            if pos_scale and pos_scale != 1:
                grid_pos *= pos_scale
            if pos_translate and pos_translate != 0:
                grid_pos += pos_translate
        if grid_pos is not None:
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
        self.cache.clear()
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
            print("Must resize node_grid_force_k in-place, OR update self.force_params")
            if len(self.node_grid_force_k.shape) == 1:
                # A single scaler for each node, applied evenly in all directions
                pass
            elif len(self.node_grid_force_k.shape) == 2:
                self.node_grid_force_k.resize(new_pos_size)
            else:
                raise NotImplementedError("node_grid_force_per_node is not implemented for node_grid_force_k of shape %s."
                                          % self.node_grid_force_k.shape)

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
        for pos, cube, text in zip(self.node_pos, self.node_obj_list, self.node_labels_list):
            ## TODO: See if you can do something similar to the edge lines, where they all
            ## just share the same positions = self.node_pos.
            cube.transform.translate = pos
            text.pos = pos
        # if self.enable_node_labels:
        #     for text in self.node_labels_list:


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
        """ Parse/process a single step event. """
        if isinstance(event_info, numbers.Number):
            return event_info
        else:
            return float(next(iter(event_info.keys())))


    def change_node(self, node, node_attrs):
        """ Change/update node_attrs for a single node. """
        nodeidx = self.nodeidx_by_id[node]
        old_node_attrs = self.node_attrs_list[nodeidx]
        cube = self.node_obj_list[nodeidx]
        changed_attrs = {}
        for k, new_value in node_attrs.items():
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
        if 'tau_cum' in changed_attrs:
            if changed_attrs['tau_cum'] > self.state_time_max:
                self.state_time_max = changed_attrs['tau_cum']


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
        """ Parse a single graphstreaming message (i.e. a JSON-formaated str). """
        return json.loads(msg)


    def read_graph_stream(self, stream):
        """
        Returns a generator that parses a stream of graph streaming events.
        The generator will return every time a "st" graph steaming step event is encountered in the stream.
        """
        # TODO: Move to separate graph-stream interpreter class/function.
        # stream is a file-like object where each next(stream) gives an event
        # <event>      ::= ( <an> | <cn> | <dn> | <ae> | <ce> | <de> | <cg> | <st> | <cl> ) ( <comment> | <EOL> )
        # Argh, JSON sucks: You cannot use integers as keys (but they can be values) - make sure to
        # convert all nodes to str before adding them!
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
            if not msg or msg[0] == "#":
                continue
            # print("Parsing msg:", msg)
            # msg is a string in the form of:
            # {"an": {<node id>: {"weight": 1}}}, {"ae": {edge_key: {"source", "target", "weight": 1}}}
            event = self.parse_graphstreaming_msg(msg)
            # print(" processing event:", event)
            for event_type, event_info in event.items():
                self.n_events_processed += 1
                method = event_methods[event_type]
                # res = method(event_info)
                if event_type == "st":
                    # Typically, event info will be a float, but it may also be a single elem dict: {time: {}}
                    step_time = self.step_event(event_info)
                    yield step_time
                else:
                    # All other events are a one or more dicts:
                    if event_type[1] == "n":
                        # Node events:
                        for nodeid, attrs in event_info.items():
                            method(str(nodeid), attrs)  # nodeid will always be a str (for JSON)
                            if self.force_directed_iterations_per_added_node:
                                self.force_directed_layout(self.force_directed_iterations_per_added_node)
                    elif event_type[1] == "e":
                        for edge_key, attrs in event_info.items():
                            source = str(attrs.pop('source'))  # so make sure source, target are also strs
                            target = str(attrs.pop('target'))
                            method(source, target, edge_key, attrs)
                            if self.force_directed_iterations_per_added_edge:
                                self.force_directed_layout(self.force_directed_iterations_per_added_edge)
                    else:
                        raise ValueError("Could not understand %s event: %s", (event_type, event_info))
            if self.force_directed_iterations_per_event_type:
                self.force_directed_layout(self.force_directed_iterations_per_event_type)
        if self.force_directed_iterations_per_msg:
            self.force_directed_layout(self.force_directed_iterations_per_msg)

        print("%s stream exhausted..." % stream)
        np.set_printoptions(precision=3) # linewidth=160, suppress=True,
        print("   Final node_pos (transposed):", self.node_pos.transpose(), sep='\n')
        print("   viewer.node_grid (transposed):", self.node_grid.transpose(), sep='\n')
        print("   node_pos - node_grid (transposed):", (self.node_pos - self.node_grid).transpose(), sep='\n')
        print("   nodeid_list for reference:", self.nodeid_list, sep='\n')


    def read_stream_and_update(self, event):
        """
        Read stream and update, intended to be called as a callback by an event timer.
        """
        # Read stream one step at a time and pass control back to event loop after each step:
        now = system_time()
        events_read = False
        if self.sleep_until is not None and now < self.sleep_until:
            ## TODO: Consider having a separate timer running visualize_layout
            self.visualize_layout(None)
            return None
        if self.graph_stream_stepping_gen is None and self.graph_stream is not None:
            self.graph_stream_stepping_gen = self.read_graph_stream(self.graph_stream)
        if self.graph_stream_stepping_gen:
            # self.stream_collect_steps_until = 10  # Avoid repeated steps, only continue when step sum > N seconds.
            # self.stream_collected_step_sum = 0
            # try:
            #     step = next(self.graph_stream_stepping_gen)
            for step in self.graph_stream_stepping_gen:
                self.n_steps_processed += 1
                events_read = True
                # time.sleep(step)
                self.stream_collected_step_sum += step
                if (not self.stream_collect_steps_until or
                    self.stream_collected_step_sum > self.stream_collect_steps_until):
                    # Break out of for loop..
                    self.stream_collected_step_sum = 0
                    if self.step_sleep_scalefactor and step > self.step_sleep_minimum_cutoff:
                        print("%s events in %s steps - step=%0.03f s, sleeping for %0.03g s" %
                              (self.n_events_processed, self.n_steps_processed,
                               step, step*self.step_sleep_scalefactor),
                              end='\r')
                        self.sleep_until = now + step*self.step_sleep_scalefactor
                    if self.force_directed_iterations_per_step_break:
                        self.force_directed_layout(self.force_directed_iterations_per_step_break)
                    break
                # continue to next event:
                if self.force_directed_iterations_per_step_cont:
                    # Make sure we do a bit of graph layout before reading the next step...
                    # (Or maybe do this unconditionally in self.read_graph_stream's generator?)
                    self.force_directed_layout(self.force_directed_iterations_per_step_cont)
            # except StopIteration:
            #     pass
            else:
                # Did not break out of for loop
                if events_read:
                    print("for step in self.graph_stream_stepping_gen completed (with events read)!")
                # pass
        # print("(read_stream_and_update) stream input read (for now..)")
        # time.sleep(0.1)
        self.visualize_layout(None)
        return events_read



    # def parse_graphstreaming_list(self, messages):
    #     """ Parse a list of graph streaming events: """
    #     # {event: {identifier: attrs}}
    #     for msg in messages:
    #         pass



    def visualize_layout(self, event):
        """ Apply a single round of force-directed layout and update visuals. """
        if self.force_params and self.force_directed_iterations_per_callback:
            self.force_directed_layout(self.force_directed_iterations_per_callback)
            # self.force_directed_layout_external()
        self.update_node_objs_pos()
        # time.sleep(0.1)

    #
    # def force_directed_layout(self, iterations=10, dim=2):
    #     """ Apply force-directed layout in two dimensions (xy plane). """


    def force_directed_layout_external(self, iterations=20, dim=3):
        """
        Use external force_directed_layout algorithm.
        """
        # dx = self.node_pos[:2, :] -
        # force_directed_layout(pos, adj_matrix, node_grid=None, iterations=100, dim=None,
        #                       force_params=None, check_params=False
        # print("force_params:", self.force_params, sep='\n')
        # print("node_grid:", self.node_grid, sep='\n')
        ## New-style force-params, has dist_enum, mask_enum, fluctuations_enum
        ## instead of
        additive = lambda shape, dim, i: 0.2*np.random.normal(size=(shape[0], dim))
        scale = lambda shape, dim, i: np.random.exponential(size=(shape[0], dim))
        force_params0 = [
            #  k,  a, dist, mask, fluctuations
            (-0.1, -2, 1, 0, None), # Repulsion, F = k / dist**2 == k * dist**(-2)
            (0.10, +1, 1, 1, None), # Spring force, F = k * dist
            # (0.10, +1, 2, 0, None), # grid force (2=node-grid-dist, 0=no-mask, None=no-fluctuation)
            ([.1, 0., 1.], +1, 2, 0, None), # grid force (2=node-grid-dist, 0=no-mask, None=no-fluctuation)
        ]
        # with random scaled fluctuations: (additive, scale)
        force_params1 = [
            #  k,  a, dist, mask, fluctuations
            # Repulsion, F = k / dist**2 == k * dist**(-2)
            (-0.1, -2, 1, 0, (None, scale)),
            # Spring force, F = k * dist
            (0.10, +1, 1, 1, (None, scale)),
            # grid force (2:node-grid-dist, 0:no-mask, None:no-fluctuation)
            # (0.10, +1, 2, 0, None),
            ([.1, 0., 1.], +1, 2, 0, (None, scale)), # grid force (2=node-grid-dist, 0=no-mask, None=no-fluctuation)
        ]
        # with randomly added fluctuations: (additive, scale)
        force_params2 = [
            #  k,  a, dist, mask, fluctuations
            # Repulsion, F = k / dist**2 == k * dist**(-2)
            (-0.1, -2, 1, 0, (additive, None)),
            # Spring force, F = k * dist
            (0.10, +1, 1, 1, (additive, None)),
            # grid force (2:node-grid-dist, 0:no-mask, None:no-fluctuation)
            # (0.10, +1, 2, 0, None),
            ([.1, 0., 1.], +1, 2, 0, (additive, None)), # grid force (2=node-grid-dist, 0=no-mask, None=no-fluctuation)
        ]
        # with randomly added fluctuations: (additive, scale)
        force_params3 = [
            #  k,  a, dist, mask, fluctuations
            # Repulsion, F = k / dist**2 == k * dist**(-2)
            (-0.1, -2, 1, 0, (additive, scale)),
            # Spring force, F = k * dist
            (0.10, +1, 1, 1, (additive, scale)),
            # grid force (2:node-grid-dist, 0:no-mask, None:no-fluctuation)
            # (0.10, +1, 2, 0, None),
            ([.1, 0., 1.], +1, 2, 0, (additive, scale)), # grid force (2=node-grid-dist, 0=no-mask, None=no-fluctuation)
        ]
        force_params = force_params3
        force_directed_layout2(self.node_pos, self.adj_matrix,
                               node_grid=self.node_grid if self.node_grid_enabled else None,
                               iterations=iterations, dim=dim, force_params=force_params)



    def force_directed_layout(self, iterations=None, check_params=False):
        """
        Apply two-dimensional force-directed layout to the nodes in pos.
            :pos:   Node positions, as N x 3 array for N nodes. *Will be updated in-place*.
            :adj_matrix: is the weighted adjacency matrix specifying edge-connections between nodes;
                    Can be directed (although that may produce unstable layouts): adj_matrix[src, tgt] = weight
            :node_grid: per-node "preferred" positions. Can be used to have certain nodes attract to one area
                    and other nodes attract to another area.
            :iterations: The number of iterations to do before returning.
            :dim:   The dimensions of pos to use. E.g. if you only want to adjust
                    xy-coordinates of nodes positions, use dim=slice(0, 2)
            :force_params: A list of "generalized force parameters" descring the forces in play.
            :check_params: Will check the input, and e.g. convert a dictionary-like force_params to the required tuple.
            :do_grid_calculation: Can be set to False to force-disable node_grid calculations.
                    default=None will enable grid-calculation if specified by any force_params.

        The algorithm uses generalized force parameters, i.e. a list of forces (k, a, use_adj_matrix, use_grid).
        Each entry in :force_params: applied as:
            if use_grid:
                F = k * grid_dist**a * erg
            else:
                F = k * dist**a * (adj_matrix if use_adj_matrix) * er

        Where <er> is a normalized node-to-node (unit) vector (i.e. of length 1).
        and  <erg> is a normalized node-to-grid vector pointing towards the grid position given by node_grid.

        Example: To get an "electrostatic repulsion", you could use
            (k=0.1, )
        You can, of course, use per-node forcce scaaling ("k" constants).
        E.g. if you want the per-node "grid matrix" to only apply force in the x-direction, you can use:
            k_grid = np.zeros_like(pos)
            k_grid[:,0] = 0.3

        """
        if self.force_params is None:
            return
        pos = self.node_pos
        adj_matrix = self.adj_matrix
        node_grid = self.node_grid_enabled and self.node_grid
        force_params = self.force_params
        dim = self.force_directed_dimensions
        if iterations is None:
            iterations = self.force_directed_iterations_per_call

        do_grid_calculation = self.do_grid_calculation

        if check_params:
            if any(isinstance(p, dict) for p in force_params):
                force_params = [(p['k'], p['a'], p['use_adj_matrix'], p['use_grid']) if isinstance(p, dict) else p
                                for p in force_params]
        if dim is None:
            dim = pos.shape[1]
        elif isinstance(dim, int):
            if dim < pos.shape[1]:
                pos = pos[:, :dim]
            if node_grid is not None:
                node_grid = node_grid[:, :dim]
        elif isinstance(dim, tuple) and len(dim) == 2:
            pos = pos[:, dim[0]:dim[1]]
            if node_grid is not None:
                node_grid = node_grid[:, dim[0]:dim[1]]
        else:
            pos = pos[:, dim]
            if node_grid is not None:
                node_grid = node_grid[:, dim]
        if node_grid is None:
            do_grid_calculation = False
        elif do_grid_calculation is None:
            do_grid_calculation = any(p[3] for p in force_params)
            self.do_grid_calculation = do_grid_calculation

        npts = pos.shape[0]  # number of nodes/points
        if 'r' in self.cache:
            r = self.cache['r']
        else:
            r = np.empty((npts, npts, pos.shape[1]), dtype='float32')
        if 'er' in self.cache:
            er = self.cache['er']
        else:
            er = np.empty_like(r)
        if 'dist' in self.cache:
            dist = self.cache['dist']
        else:
            dist = np.empty((npts, npts), dtype='float32')
        if 'F' in self.cache:
            F = self.cache['F']
        else:
            F = np.empty((len(force_params), npts, pos.shape[1]))
        if 'F_sum' in self.cache:
            F_sum = self.cache['F_sum']
        else:
            F_sum = np.empty_like(pos)

        if do_grid_calculation:
            # Initialize data structures for grid-force calculations:
            rg = np.empty_like(pos)
            erg = np.empty_like(rg)
            grid_dist = np.empty_like(rg.shape[0])

        if self.continuous_iterations:
            it = self.n_layout_iterations_count
            it_end_at = iterations + self.n_layout_iterations_count
        else:
            it = 0
            it_end_at = iterations
        while it < it_end_at:
            it += 1
            # Coordinate differences:
            r[:] = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
            # Calculate node-node distances and normalized node-to-node directionality unit vector
            dist = (r**2).sum(axis=2)**0.5     # pair-wise distance between all nodes
            dist[dist == 0] = 1.                # Prevent division-by-zero errors
            er[:] = r / dist[..., np.newaxis]    # normalized node-node directionality vector
            if do_grid_calculation:
                rg[:] = node_grid - pos
                # Calculate scalar distances and normalized node-to-node directionality unit vector
                grid_dist = (rg**2).sum(axis=1)**0.5     # pair-wise distance between all nodes
                grid_dist[grid_dist == 0] = 1.                # Prevent division-by-zero errors
                erg[:] = rg / grid_dist[..., np.newaxis]    # normalized node-node directionality vector

            # New-style force_params with fluctuations support:
            for Fidx, (k, a, dist_enum, mask_enum, fluctuations) in enumerate(force_params):
                # F[<type>, <node>, <force vector>]
                if dist_enum == 0:
                    # Do not use any distance dependence, only e.g. thermal fluctuations:
                    F[Fidx, :, :] = k
                elif dist_enum == 2:
                    # Use grid: grid_dist and erg unit vector
                    F[Fidx, :, :] = k * grid_dist[..., np.newaxis]**a * erg
                else:
                    # Use node-node distances (default): dist and er unit vector
                    if mask_enum == 1:
                        # Use adj_matrix mask:
                        F[Fidx, :, :] = k * (dist[..., np.newaxis]**a * adj_matrix[:, :, np.newaxis] * er).sum(axis=0)
                    else:
                        # No mask (default):
                        F[Fidx, :, :] = k * (dist[..., np.newaxis]**a * er).sum(axis=0)
                if fluctuations:
                    # F = fluc_scale * (k dist^a + fluc_additive)
                    additive, scaling = fluctuations
                    if additive:
                        F[Fidx, :, :] += additive(pos.shape, dim, it)
                    if scaling:
                        F[Fidx, :, :] *= scaling(pos.shape, dim, it)


            # Calculate displacement, making sure forces are not excessively high:
            dx = np.clip(F.sum(axis=0), -3, 3) * self.force_dx_factor
            pos[:, :] += dx
        self.n_layout_iterations_count += self.force_directed_iterations_per_call


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



def run_tests():
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
    # graph_stream = test_stream

    run = "2016-03-14 161538"
    run = "2016-03-14 182118"
    run = "2016-03-14 184805"
    run = "2016-03-14 192212"
    run = "2016-03-14 193519"
    graph_stream = open(("/Users/rasmus/Dev/nascent/examples/single_duplex/simdata/fourway_junction_1/"
                         "%s/complexes/reaction_graph_eventstream.json" % run))
    ## How long does it take to read the whole file and decode it, line by line?
    # t_start = system_time()
    # print("Test-reading and parsing whole file stream, line by line....")
    # n_read = 0
    # for i, line in enumerate(graph_stream):
    #     msg = line.strip()
    #     if msg and msg[0] != "#":
    #         # viewer.parse_graphstreaming_msg(line)
    #         data = json.loads(msg)
    #         n_read += 1
    #         print(data)
    #     print("%s lines read and parsed..." % (i+1), end='\r')
    # t_end = system_time()
    # print("\n %s lines read and parsed line by line in %s seconds..." % (n_read, t_end-t_start))
    # ## Takes only 0.06 seconds to read 6000 lines..
    # ## if we did that at 60 reads per second it would be 100 seconds..
    # answer = input("Press any key to continue...")
    # if answer in ('q', 'n'):
    #     sys.exit()
    # graph_stream.seek(0)


    # print("node_grid_force_k before resizing:", viewer.node_grid_force_k, sep='\n')
    # viewer.node_grid_force_k.resize((2,), refcheck=False)
    # print("node_grid_force_k after resizing:", viewer.node_grid_force_k, sep='\n')
    # viewer.node_grid_force_k[:] = [0.5, 0.8] # Only apply node_grid force in the x-direction (on x-coordinates).
    # print("node_grid_force_k after assignment:", viewer.node_grid_force_k, sep='\n')
    # if the 'interactive' flag is set, you can interact with the program via an interactive terminal (e.g ipython)
    # viewer.read_stream_and_update(None)


    viewer = GraphViewer()
    viewer.graph_stream = graph_stream

    if sys.flags.interactive == 0: #and False:
        app.run()
        # viewer.read_graph_stream(viewer.graph_stream)
        # pass
        # filename = 'test.png'
        # print("App run done; saving rendered image to file:", os.path.abspath(filename))
        # image = viewer.render()
        # vispy.io.write_png(filename, image)
        print("Done app.run() !")
    else:
        test_events_read = True # and False
        n_steps = 0
        while test_events_read is not False:
            test_events_read = viewer.read_stream_and_update(None)
            n_steps += 1
        else:
            print("\n\nStream manually exhausted, %s steps encountered.\n" % n_steps)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)
        viewer.read_stream_and_update(None)

    print("Final node_pos:", viewer.node_pos, sep='\n')
    print("viewer.node_grid:", viewer.node_grid, sep='\n')
    print("node_pos - node_grid:", viewer.node_pos - viewer.node_grid, sep='\n')
    print("nodeid_list for reference:", viewer.nodeid_list, sep='\n')


if __name__ == '__main__':
    run_tests()

