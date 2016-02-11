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

Thread refs:
* https://docs.python.org/3/library/queue.html

"""



from __future__ import absolute_import, print_function, division
import sys
import os
import yaml
import re
import time
import threading
import queue # was Queue in python 2.
import random
import math


class DispatcherThread(threading.Thread):
    """
    Sub-classing approach
    """

    def __init__(self, msg_queue, config):
        self.msg_queue = msg_queue
        #self.config = config.copy()

    def run(self):
        while True:
            msg = self.msg_queue.get()  # This will hang execution until an element is available in the queue.
            print(msg)
            if msg is None:
                # Reading a None value from the queue will cause the loop to break and the thread to die.
                break


def dispatch_from_queue(msg_queue, config):
    """
    Direct threading from a function
    """
    outputfn = config.get('statechange_outputfn', 'graphstatechanges.gsc')
    mode = config.get('statechange_output_mode', 'w')
    with open(outputfn, mode) as fp:
        while True:
            msg = msg_queue.get()  # This will hang execution until an element is available in the queue.
            print(msg)
            if msg is None:
                # Reading a None value from the queue will cause the loop to break and the thread to die.
                break
            # Do whatever you need to do:
            fp.write(msg)


def init_dispatcher(config, msg_queue=None, start=False):
    """
    Mostly for example usage...
    """
    if msg_queue is None:
        # Insertions into the que will block once maxsize is reached, until someone get()s from the other end.
        msg_queue = queue.Queue(maxsize=0)  # maxsize < 1 makes an infinite queue.

    # Version 1, using a Thread subclass:
    thread = DispatcherThread(msg_queue, config)
    # Version 2, instantiating directly from the Thread class using target:
    thread = threading.Thread(target=dispatch_from_queue, args=(msg_queue, config))
    # Start the thread with:
    if start:
        thread.start()

    return thread, msg_queue



def computational_result(val):
    return math.sqrt(val+1)



def calc_chunck(data):
    return sum(computational_result(val) for val in data)



def thread_tester(data, N=2):

    assert len(data) % N == 0
    slice_size = int(len(data)/N)
    data_chunks = [data[slice_size*i:slice_size*(i+1)] for i in range(N)]
    threads = [threading.Thread(target=calc_chunck, args=(data_chunks[i],)) for i in range(N)]




def test_threads():

    data = [random.random() for i in range(10000000)]
