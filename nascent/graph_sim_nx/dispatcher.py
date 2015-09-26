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
    parser.add_argument('--sep-by-eof', action="store_true", default=None)
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


def touch(filepath):
    with open(filepath, 'a'):
        os.utime(filepath)


def find_unique_fn(filename):
    if not os.path.exists(filename):
        return filename
    regex_prog = re.compile(r"(.+?)(\d*)")
    fnroot, ext = os.path.splitext(filename)
    match = regex_prog.fullmatch(fnroot) # python 3.4+ only. Appends ^ before and $ after the regex.
    fnbase, num = match.groups()
    num = int(num) if num else 0
    fnfmt = "%s_%s%s"
    while os.path.exists(os.path.join(filename)):
        num += 1
        filename = fnfmt % (fnbase, num, ext)
    return filename




def dispatch_loop(outputfn, appendfile=False, print_input=False, sep_by_eof=True,
                  binary_input=None, output_encoding=None, output_buffer=-1, **kwargs):
    """
    Aka read_stdin_save_to_disk_and_forward_to_graph_visualizer.

    Note: If this process is sharing a terminal with a parent process,
    any input that this program outputs will be displayed in the terminal,
    even if the parent process is closed.

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
    mode = ('a' if appendfile else 'w') + 'b'
    fout_buf = 0 # 0=unbuffered, 1=line-buffered, negative=system-default, >1=n-bytes-buffered.
    # un-buffered only works for bytes, not text. Use type/cat to check the content of file; don't rely on subl.
    if binary_input is None:
        # auto-detect:
        binary_input = 'b' in sys.stdin.mode
    if not binary_input and output_encoding is None:
        output_encoding = sys.stdin.encoding
    with open(outputfn, mode, fout_buf) as fp:
        try:
            print("sep_by_eof:", sep_by_eof)
            #while True:
                # msg = sys.stdin.readline()  # repeated calls to stdin.read in infinite loop - very cpu intensive.
            done = False
            # I don't think sys.stdin.closed will ever be True, if spawned as PIPE with subprocess.Popen.
            # unless maybe if parent process is terminated.
            while not sys.stdin.closed and not done:
                # There is an issue if using subprocess PIPEs (io.BufferedReader on this side)
                # It is never really closed.
                for msg in (sys.stdin.readlines() if sep_by_eof else sys.stdin): # pylint: disable=C0325
                    # stdin.readlines() will read from stdin, halting here until it is closed or we receive an EOF byte.
                    # Only after EOF has been read will the input be split by lines and execution resumed.
                    # EOF can be given using ctrl+D on unix and ctrl+Z on windows terminal.
                    # you could use this by encapsulating readlines() having an infinite loop outside
                    # if using "for msg in stdin", then stdin.read() is called until reading a newline char.
                    # This is implemented to be very easy on the cpu, waiting until there is actually input.
                    # messages are separated by newline chars:
                    print("msg:", msg.strip())
                    if not msg:
                        print("Received empty msg:", msg)
                        continue
                    if not binary_input: # msg is not binary, convert
                        msg = bytes(msg, output_encoding)
                    fp.write(msg)
                if not sep_by_eof:
                    # If we are reading "for msg in sys.stdin" separated by new line, then do not loop.
                    done = True
                print('done with "for msg in (sys.stdin.readlines() if %s else sys.stdin)"' %  sep_by_eof)
                print("sys.stdin.closed:", sys.stdin.closed)
                time.sleep(0.2)
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
                            bufsize=0, # 0=unbuffered, 1=line-buffered, -1=system default
                            stdin=subprocess.PIPE)
    sleep = 0.5
    try:
        for i in range(int(3/sleep)):
            print("cycle:", i)
            # write binary data if universal_newlines=False, else write strs
            # Probably wrap in try... except to catch IO errors if subprocess is closed.
            # e.g. if e.errno == errno.EPIPE or e.errno == errno.EINVAL
            proc.stdin.write("hej der\n\n")
            proc.stdin.write("this line doesn't end here: -->.")
            #proc.stdin.write("This is a long repeated line. "*10)
            if args['sep_by_eof']:
                print("Sending EOF...")
                # Send EOF by closing the pipe?
                # proc.stdin.close()
                # That will send an EOF, but it also
                # closes the pipe completely from this side, meaning I cannot write to it again.
                # (The pipe remains open on the subprocess' side.)
                # Actually, that is exactly what ctrl+D does (ctrl+Z on Windows), it signals to the pty to close the stream.
                # However, you might be able to emulate that by sending an EOF charater:
                if 'win' in sys.platform:
                    EOF = "\x1a"
                else:
                    EOF = "\x04"
                # Edit: this is probably doomed, since proc.stdin.isatty() returns False for io.BufferedWriter objects.
                proc.stdin.write(EOF)
            #else:
                # flush the pipe:
            proc.stdin.flush()
            time.sleep(sleep)
            # if using while loop, you should have a p.wait()
    except KeyboardInterrupt:
        print("test_self: KeyboardInterrupt requested")
    print("Closing proc.stdin...")
    #proc.stdin.close()
    # this doesn't work. At least, proc.stdin.closed is still False.
    # This is probably because:
    # (1) proc.stdin is just a io.BufferedWriter(io.FileIO()) (on this side)
    # (2) proc.stdin is just a io.BufferedReader o the subprocess' side.
    # closing proc.stdin from here only affect io.BufferedWriter.
    print("proc.stdin.closed:", proc.stdin.closed)
    # proc.stdin.close() will make sure that the EOF is read when invoking stdin.read()
    print("Manually setting proc.stdin.closed=True")
    # proc.stdin.closed = True  # AttributeError: attribute 'closed' of '_io.TextIOWrapper' objects is not writable
    print("proc.stdin.closed:", proc.stdin.closed)
    # proc.wait() # wait here until the proc subprocess has terminated.
    # proc.poll() # check to see if the proc subprocess has terminated.
    # proc.communicate() # send optional input to proc, read from stdout/stderr and wait for proc to terminate.
    time.sleep(10)
    # proc.terminate()  # Not very effective if we have encountered the infinite loop. Sends SIGTERM.
    proc.kill()  # on Windows, this is just an alias for terminate()
    print("proc.poll():", proc.poll())   # None
    import signal
    proc.send_signal(signal.CTRL_BREAK_EVENT)
    time.sleep(10)
    print("proc.poll():", proc.poll())   # 1 = completed, with errors.



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
