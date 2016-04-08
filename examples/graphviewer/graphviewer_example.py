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



"""

Basic vispy graphviewer example.

"""

import sys
import os

sys.path.insert(0, ".")
print("os.path.abspath('.'):", os.path.abspath('.'))


import vispy
# from vispy.util.transforms import perspective
# from vispy.util import transforms
# from vispy import gloo
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



from nascent.graph_visualization.vispy_graph_viewer_app import GraphViewer



def main():

    from vispy import app


    # Select which run to analyse:
    run = "2016-03-14 161538"
    run = "2016-03-14 182118"
    run = "2016-03-14 184805"
    run = "2016-03-14 192212"
    run = "2016-03-14 193519"
    run = "2016-04-06 153346"
    run = "2016-04-07 164602"

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

    config = {
        'scene_draw_ground': False, # Disable ground
        'scene_ground_color': 'grey', # default
        'scene_ground_scale': None,
        # 'scene_ground_translate': [0., 0., -4],
        # 'scene_background_color': False,
        # 'scene_bgcolor': '#efefef', # slightly off-white
        'scene_bgcolor': 'white',
    }

    viewer = GraphViewer(config=config)
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
    main()

