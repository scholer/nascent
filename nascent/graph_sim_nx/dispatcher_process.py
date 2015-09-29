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


"""

Dispatcher script. Responsible for:
 1) Saving state change events to file:
 2) Forwarding the state change to real-time graph visualizer,
    taking care of proper graph style hints in the process.


Implementation notes:


Which is better: a separate subprocess, or a thread?

Subprocesses
* For inter-process communication, we typically communicate through a pipe,
    this side doing stdin.write(...) and the subprocess reading with "for line in stdin"
    Note:

Threads:
* For threads we have a few options, e.g. a Queue:
    this side does "queue.put(msg)", the sub-thread does "while True: queue.get(msg)"
    Using Thread+Queue is very easy and convenient. However, I'm not sure how much this is
    going to help me, since threads are still subject to the Global Interpreter Lock (GIL).
    "Python threads cannot actually concurrently access state in the interpreter (there's one big lock, the
    infamous Global Interpreter Lock.) What that means in practice is that threads are useful for I/O bound
    tasks (networking, writing to disk, and so on), but not at all useful for doing concurrent computation."
    -- this is actually exactly what I want to use the dispatcher for, so that should be A-OK.

Subprocesses - with the multiprocessing library:
* If you have to do concurrent CPU-bound work on multiple processors, doing this using the multiprocess
    module makes it A LOT easier. It even supports the same "Queue" model as threading, so
    Note that spawning processes is a lot costlier than spawning threads, and communication might
    be a bit slower.

Conclusion:
 1. I should really try to stick to using threads.
 2. If I really need to offload work, use the multiprocessing library!
 3. Spawning subprocesses with subprocess module is only for non-python code,
    e.g. if you need to interact with a pre-compiled executable.


More on threads and concurrency in Python:
* http://www.slideshare.net/dabeaz/an-introduction-to-python-concurrency
* http://stackoverflow.com/questions/1190206/threading-in-python
* https://www.praetorian.com/blog/multi-core-and-distributed-programming-in-python
* https://docs.python.org/3/library/multiprocessing.html
* http://stackoverflow.com/questions/2629680/deciding-between-subprocess-multiprocessing-and-thread-in-python



Other alternatives?
* Use zeromq - http://zeromq.org/
* Use a sockets (on Unix)
* Julia has parallel computing as an integral part: http://julia.readthedocs.org/en/latest/manual/parallel-computing/



"""

import sys
import os
import argparse
import yaml
import re
import time
import subprocess


class Dispatcher():

    def __init__(self, config):
        pass




def parse_args(argv=None):
    """
    Parse args from the command line.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('-o', '--outputfn')
    parser.add_argument('-y', '--overwrite', action="store_true")
    parser.add_argument('-c', '--cfgfile')
    parser.add_argument('-a', '--appendfile', action="store_true")
    parser.add_argument('--testing', action="store_true")
    parser.add_argument('--print-input', action="store_true")
    #parser.add_argument('--sep-by-eof', action="store_true", default=None)
    parser.add_argument('--binary-input', action="store_true", default=None)
    parser.add_argument('--output-encoding', action="store_true", default=None)
    parser.add_argument('--output-buffer', type=int, default=-1)

    argns = parser.parse_args(argv)
    return argns


def process_args(args):
    """
    Process command line args, load config from file if requested and make a ready-to-use args dict.
    """

    default_args = {'outputfn': 'defaultoutputfn.gsc'} # .gsc = graph-state-change
    if args.get('cfgfile'):
        with open(args['cfgfile']) as fp:
            fileconf = yaml.load(fp)
        default_args.update(fileconf)

    for k, v in default_args.items():
        if args.get(k) is None:
            args[k] = v

    # If overwriting or appending to file, don't generate unique outputfn.
    if not (args.get('overwrite') or args.get('appendfile')):
        args['outputfn'] = find_unique_fn(args['outputfn'])

    return args





def dispatch_loop(outputfn, appendfile=False, print_input=False,
                  binary_input=None, output_encoding=None, output_buffer=-1, **kwargs):
    """
    Aka read_stdin_save_to_disk_and_forward_to_graph_visualizer.

    Note: If this process is sharing a terminal with a parent process,
    any input that this program outputs will be displayed in the terminal,
    even if the parent process is closed.

    Args:
        output_buffer : 0=unbuffered, 1=line-buffered, negative=system-default, >1=n-bytes-buffered.
                        un-buffered only works for bytes, not text.
                        Note: Use type/cat to check the content of file; don't rely on subl.

    Refs on how to get this to work:
    http://blog.codedstructure.net/2011/02/concurrent-queueget-with-timeouts-eats.html (Threading)


    Options for reading on sys.stdin:
    for loop: - while True: sys.stdin.readline()
        Infinite loop, very CPU intensive. Does not explicitly wait for new input, just
        repeatedly reads from stdin even when no new input is available.
    for line in sys.stdin:
        Reads *a single line* on stdin, then waits for new input.
        This is very cheap and cpu efficient.
    for msg in sys.stdin.readlines()
        This waits, not until a new line is received, but until EOF is received.
        Essentially, this is like any other file's readlines, which is not an iterable, but first invokes
        file.read(), which reads the whole file into memory, and then splits the text by newlines.
        Note: you can send EOF without closing stdin. This is actually quite good for sending
        full messages that may include new lines.

    Issues:
    * If parent process fails, or sys.stdin is otherwise "closed by EOF but not sys.stdin.closed=True",
        then we have a fast-cycling, cpu intensive infinite loop. Bad.

    """
    # reads and writes to stdin/stdout (i.e. msg) is bytes, not string, so open in binary mode:
    mode = ('a' if appendfile else 'w')
    if binary_input is None:
        binary_input = 'b' in sys.stdin.mode  # auto-detect
    print("binary_input: %s, sys.stdin.mode: %s" % (binary_input, sys.stdin.mode))
    if binary_input:
        mode += 'b'
    #if not binary_input and output_encoding is None:
    #    output_encoding = sys.stdin.encoding or 'utf-8'
    with open(outputfn, mode, output_buffer) as fp:
        try:
            for msg in sys.stdin:
                fp.write(msg)
                #fp.flush()  # flush()ing to disk may degrade performance, but is nice for debugging.
                if print_input:
                    print("msg: %s" % msg)
                # Do whatever else you need to do with msg:
                time.sleep(1) # Test sleeping here
                # The dispatch_loop can cause parent process to halt if the buffer is full
                # when the parent process attempts to write to the buffer.
        except KeyboardInterrupt:
            print("Process interrupted by KeybordInterrupt.")


def test_self(subargv, args):
    """
    This is just as much to give a use example for this script:

    Note: Be careful with pipes:
    * http://blogs.msdn.com/b/oldnewthing/archive/2011/07/07/10183884.aspx
    """
    exe = 'python' # or sys.executable
    if exe not in subargv[0]:
        subargv = [exe] + subargv
    print("Spawning subprocess with:", subargv)
    proc = subprocess.Popen(subargv,
                            shell=True,
                            universal_newlines=True, # False: open stdin/out/err as binary streams (default)
                            # 0=unbuffered, 1=line-buffered, -1=system default
                            # bufsize=0 only for binary streams (universal_newlines=False)
                            bufsize=-1, # doesn't matter, we flush() it after each message.
                            stdin=subprocess.PIPE)
    sleep = 0.5
    try:
        for i in range(int(10/sleep)):
            print("test cycle:", i)
            # write binary data if universal_newlines=False, else write strs
            # Probably wrap in try... except to catch IO errors if subprocess is closed.
            # e.g. if e.errno == errno.EPIPE or e.errno == errno.EINVAL
            proc.stdin.write("hej der\n\n")
            #proc.stdin.write("this line doesn't end here: -->.")
            #proc.stdin.write("this line is really, really long."*30)
            # A note on buffering: Generally, you probably
            proc.stdin.flush()  # Flush is needed on this side if you want "instant messaging"
            #time.sleep(sleep)
            # Note: proc.stdin.write will write to stdin's buffer. Once the buffer is full, it will
            # THEN WAIT here until someone has read stdin from the other end!
            # invoking stdin.flush() basically tells the buffer to start unloading to who-ever
            # is "buying" on the other end.
            # Thus, using subprocesses may not be that
            # http://blogs.msdn.com/b/oldnewthing/archive/2011/07/07/10183884.aspx
            # Perhaps
            time.sleep(0.2)
    except KeyboardInterrupt:
        print("test_self: KeyboardInterrupt requested")



def main(argv=None):
    argns = parse_args(argv)
    args = process_args(argns.__dict__)

    if args.get('testing'):
        subargv = [arg for arg in (argv or sys.argv) if 'testing' not in arg]
        test_self(subargv, args)
    else:
        dispatch_loop(**args)


if __name__ == '__main__':
    main()
