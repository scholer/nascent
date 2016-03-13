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

Module with graph layout algorithms.


"""

import numpy as np



def force_directed_layout(pos, adj_matrix, node_grid=None, iterations=100, dim=None,
                          force_params=None, check_params=False, do_grid_calculation=None):
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

    if force_params is None:
        force_params = [
            (-0.1, -2, False, False),    # Repulsion, F = k / dist**2 == k * dist**(-2)
            (0.05,  1, True, False),     # Spring force, F = k * dist
            (0.05, -1, False, False),   # "Gravity" (but long range, F=Gmm/r, not F=Gmm/r^2)
        ]
    if check_params:
        if any(isinstance(p, dict) for p in force_params):
            force_params = [(p['k'], p['a'], p['use_adj_matrix'], p['use_grid']) if isinstance(p, dict) else p
                            for p in force_params]
    if dim is None:
        dim = pos.shape[1]
    elif isinstance(dim, int):
        pos = pos[:, :dim]
        if node_grid is not None:
            node_grid = node_grid[:, :dim]
    elif isinstance(dim, tuple) and len(dim) == 2:
        pos = pos[:, dim[0]:dim[1]]
        if node_grid is not None:
            node_grid = node_grid[:, dim[0]:dim[1]]
    else:
        pos = pos[:, dim[0]:dim[1]]
        if node_grid is not None:
            node_grid = node_grid[:, dim[0]:dim[1]]
    if do_grid_calculation is None:
        do_grid_calculation = any(p[3] for p in force_params)
    if do_grid_calculation:
        assert node_grid is not None

    npts = pos.shape[0]  # number of nodes/points
    r = np.empty((npts, npts, pos.shape[1]), dtype='float32')
    er = np.empty_like(r)
    dist = np.empty((npts, npts), dtype='float32')
    F = np.empty((len(force_params), npts, pos.shape[1]))
    F_sum = np.empty_like(pos)

    # Grid stuff:
    # TODO: Make a short-circuit while loop that doesn't check "if use_grid" for every iteration!
    if do_grid_calculation:
        rg = np.empty_like(pos)
        erg = np.empty_like(rg)
        grid_dist = np.empty_like(rg.shape[0])


    while iterations > 0:
        iterations -= 1
        # Coordinate differences, only the first two (x, y) coordinates:
        # TODO: Use a persistent "r", "er", "dist", "F_r", "F_s", etc datastructures and update in-place.
        r[:] = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
        # Calculate scalar distances and normalized node-to-node directionality unit vector
        dist = (r**2).sum(axis=2)**0.5     # pair-wise distance between all nodes
        dist[dist == 0] = 1.                # Prevent division-by-zero errors
        er[:] = r / dist[..., np.newaxis]    # normalized node-node directionality vector
        if do_grid_calculation:
            rg[:] = node_grid - pos
            # Calculate scalar distances and normalized node-to-node directionality unit vector
            grid_dist = (rg**2).sum(axis=1)**0.5     # pair-wise distance between all nodes
            grid_dist[grid_dist == 0] = 1.                # Prevent division-by-zero errors
            erg[:] = rg / grid_dist[..., np.newaxis]    # normalized node-node directionality vector


        # # Repulsive force:  ke*(q1q2)/r^2
        # # all points push away from each other
        # # F_r = 0.02 / dist[..., np.newaxis]**4 * er  # Direct calculation+assignment
        # k_r = -0.1
        # F_r = k_r / dist[..., np.newaxis]**2 * er  # Direct calculation+assignment
        #
        # # Attractive spring force:
        # # connected points pull toward each other with hookean force: F = ks*x,  ks = spring force constant
        # # (Note: pulsed force may help the graph to settle faster)
        # ks = 0.05
        # #ks = 0.05 * 5 ** (np.sin(i/20.) / (i/100.))
        # #ks = 0.05 + 1 * 0.99 ** i  # Exponentially decaying force to a constant value
        # # ks = 0.5 + 1 * 0.999 ** iterations  # Exponentially decaying force to a constant value
        # F_spring = ks * dist[..., np.newaxis] * adj_matrix[:, :, np.newaxis] * er
        # # F_spring = ks * dist[..., np.newaxis]**2 * self.adj_matrix[:, :, np.newaxis] * er
        #
        # # Gravity force: F = G * m1m2/r^2 * dre,
        # # where dre is normalized node-node directionality vector and G is gravitational constant
        # # Edit: Make sure gravity and repulsion have different distance dependencies.
        # Gm1m2 = 0.05
        # F_g = Gm1m2 * er / dist[..., np.newaxis] # dist[..., np.newaxis]**2
        #
        for i, (k, a, use_edges, use_grid) in enumerate(force_params):
            # Should we sum here or later? I.e. should we have F[i,:,:] be a the "total force" from all
            # nodes on each node, or should we keep the details and have F[i,:,:,:] ?
            if use_grid:
                F[i, :, :] = k * grid_dist[..., np.newaxis]**a * erg
            elif use_edges:
                F[i, :, :] = k * (dist[..., np.newaxis]**a * adj_matrix[:, :, np.newaxis] * er).sum(axis=0)
                # if a == 1 and k == ks:
                #     assert np.isclose(F[i, :, :], F_spring.sum(axis=0), rtol=1e-4).all()
                # else:
                #     print("a == %s and k == %s" % (a, k))
            else:
                F[i, :, :] = k * (dist[..., np.newaxis]**a * er).sum(axis=0)
                # if a == -2 and k == k_r:
                #     assert np.isclose(F[i, :, :], F_r.sum(axis=0), rtol=1e-4).all()
                # else:
                #     print("a == %s and k == %s" % (a, k))




        # # All forces from all particles:
        # F_all = F_spring + F_r + F_g
        # # Make sure points do not exert force on themselves:
        # F_all[np.arange(npts), np.arange(npts)] = 0

        # dv = F/m*dt,  v1 = v0*dv,  dx = v*dt
        # But we don't have persistent velocities (i.e. start and end with v=0 for every step).
        # So we just use F*dt^2/m/4 for all particles, and say that dt^2/m/4 is some constant
        c_dt2m = 0.01 # 0.09   # In units of length/force
        # F_all_sum = F_all.sum(axis=0)   # Sum all contributions from all other particles
        # print(F_sum)
        F_sum = F.sum(axis=0)
        # assert F_all_sum.shape == F_sum.shape
        # if force_params == [
        #     (-0.1, -2, False, False),    # Repulsion, F = k / dist**2 == k * dist**(-2)
        #     ( 0.05, 1, True, False),     # Spring force, F = k * dist
        # ]:
        #     assert np.isclose(F_sum, F.sum(axis=0), rtol=1e-4).all()

        # Make sure no force is "too high":
        dx = np.clip(F_sum, -3, 3) * c_dt2m
        pos[:, :] += dx



def force_directed_layout_2d(pos, adj_matrix, iterations=100):
    """ Apply two-dimensional force-directed layout to the nodes in pos. """
    npts = pos.shape[0]
    # shape = (N, 3), i.e. N nodes with 3-elem position/coordinate vectors.
    # pos[:, np.newaxis, :] N x N matrix of position-vectors, each column with the same "array of pos vectors"
    # pos[np.newaxis, :, :] N x N matrix of position-vectors, each row with the same "array of pos vectors"
    while iterations > 0:
        iterations -= 1
        # Coordinate differences, only the first two (x, y) pos
        # TODO: Use a persistent "r", "er", "dist", "F_r", "F_s", etc datastructures and update in-place.
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
        F_spring = ks * dist[..., np.newaxis] * adj_matrix[:, :, np.newaxis] * er
        # F_spring = ks * dist[..., np.newaxis]**2 * self.adj_matrix[:, :, np.newaxis] * er

        # Gravity force: F = G * m1m2/r^2 * dre,
        # where dre is normalized node-node directionality vector and G is gravitational constant
        # Edit: Make sure gravity and repulsion have different distance dependencies.
        Gm1m2 = 0.05
        F_g = Gm1m2 * er / dist[..., np.newaxis] # dist[..., np.newaxis]**2

        # All forces from all particles:
        F_all = F_spring - F_r + F_g
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





if __name__ == '__main__':

    test_pos_org = np.random.rand(7,3)
    print("Starting pos:")
    print(test_pos_org)

    test_edge_list = [
        (0, 1), (0, 2), (1, 2), (2, 3),
        (3, 4), (3, 5), (3, 6),
        (0, 6), (5, 6)
    ]
    test_adj_matrix = np.empty((7, 7), dtype=bool)
    for src, tgt in test_edge_list:
        test_adj_matrix[src, tgt] = True
        test_adj_matrix[tgt, src] = True
    print("Adj matrix:")
    print(test_adj_matrix)

    test_pos = test_pos_org.copy()
    # force_directed_layout_2d(test_pos, test_adj_matrix, iterations=1)
    force_directed_layout_2d(test_pos, test_adj_matrix)
    print("\nDone force_directed_layout_2d(test_pos, test_adj_matrix), pos:")
    print(test_pos)


    # def force_directed_layout(pos, adj_matrix, node_grid=None, iterations=100, dim=None,
    #                           force_params=None, check_params=False):
    # force_params: list of forces (k, a, use_adj_matrix, use_grid)

    # Test dim=2:
    test_pos = test_pos_org.copy()
    force_directed_layout(test_pos, test_adj_matrix, dim=2)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, dim=2), pos:")
    print(test_pos)


    # Test DEFAULTS:
    test_pos = test_pos_org.copy()
    force_directed_layout(test_pos, test_adj_matrix)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix), pos:")
    print(test_pos)

    # Test force_params:
    # test_pos = np.random.rand(7,3)
    test_pos = test_pos_org.copy()
    test_force_params = [
        (-0.1, -2, False, False),
        (0.06, 1, True, False),
    ]
    force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params), pos:")
    print(test_pos)


    # Test node_grid:
    # test_pos = np.random.rand(7,3)
    test_pos = test_pos_org.copy()
    test_force_params = [
        (-0.1, -2, False, False),
        (0.05, 1, True, False),
        (0.01, 1, False, True),  # grid force. Note: Can only use both adj_matrix or node_grid, not both.
    ]
    test_node_grid = np.array([
        [0, 0],
        [1, 0],
        [2, -1],
        [2, 2],
        [3, -1],
        [3, 0],
        [3, 1],
    ])
    force_directed_layout(test_pos, test_adj_matrix, node_grid=test_node_grid, dim=2,
                          force_params=test_force_params)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params), pos:")
    print(test_pos)

    # Test node_grid with per-node force constant:
    # test_pos = np.random.rand(7,3)
    test_pos = test_pos_org.copy()
    test_grid_k = np.zeros_like(test_node_grid)
    test_grid_k[0, 0] = 1  # first node, force constant on x-axis
    test_grid_k[1, 0] = 1
    test_grid_k[2, 1] = 1
    test_force_params = [
        (-0.1, -2, False, False),
        (0.05, 1, True, False),
        (test_grid_k, 1, False, True),    # Use node_grid force
    ]
    force_directed_layout(test_pos, test_adj_matrix, node_grid=test_node_grid, dim=2,
                          force_params=test_force_params)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params)")
    print("test_grid_k:")
    print(test_grid_k)
    print("pos:")
    print(test_pos)

    # Test node_grid with per-coordinate axis force constant:
    test_pos = test_pos_org.copy()
    test_grid_k = np.array([1, 0])  # scale 1 on x-axis, 0 on y-axis
    test_force_params = [
        (-0.1, -2, False, False),
        (0.05, 1, True, False),
        (test_grid_k, 1, False, True),    # Use node_grid force
    ]
    force_directed_layout(test_pos, test_adj_matrix, node_grid=test_node_grid, dim=2,
                          force_params=test_force_params)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params)")
    print("test_grid_k:")
    print(test_grid_k)
    print("pos:")
    print(test_pos)

    # Test node_grid with 1-element array as force constant:
    test_pos = test_pos_org.copy()
    test_grid_k = np.array([1])  # scale 1 on x-axis, 0 on y-axis
    test_force_params = [
        (-0.1, -2, False, False),
        (0.05, 1, True, False),
        (test_grid_k, 1, False, True),    # Use node_grid force
    ]
    force_directed_layout(test_pos, test_adj_matrix, node_grid=test_node_grid, dim=2,
                          force_params=test_force_params)
    print("\nDone force_directed_layout(test_pos, test_adj_matrix, force_params=test_force_params)")
    print("test_grid_k:")
    print(test_grid_k)
    print("pos:")
    print(test_pos)
