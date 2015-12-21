"""
Debugging tips and tricks


Q:  I have the id of an object, but I don't have the variable reference
    directly available. Can I get it from the id?
A:  Yes, but it is either (a) CPython specific (implementation may vary)
    or very, very slow.
    Example 1: Get the ends objects from this exception message:
    Exception: <class 'networkx.exception.NetworkXNoPath'> No path between
    <nascent.graph_sim_nx.domain.Domain5pEnd object at 0x0000000005FA7358> and
    <nascent.graph_sim_nx.domain.Domain5pEnd object at 0x0000000005FA7940>.
"""
import ctypes
end1 = ctypes.cast(0x0000000005FA7358, ctypes.py_object).value
end2 = ctypes.cast(0x0000000005FA7940, ctypes.py_object).value
"""
    Alternatively, using gc.get_objects():
"""
import gc
end1 = next(obj for obj in gc.get_objects() if id(obj) == 0x0000000005FA7358)
end2 = next(obj for obj in gc.get_objects() if id(obj) == 0x0000000005FA7940)


"""
Q: What is the best way to use the debugger.
A: Run the program with
        python -m pdb <script/program.py>
   You need to press 'c' to continue initially, but then the programme will run normally until
   an un-caught exception is raised.

More debugger refs:
* http://www.sixfeetup.com/blog/debugging-python-code-with-pdb



Q: How do I displaly unicode characters in a Windows console?
A: `chcp 65001`  will will change the code page to UTF-8. Use e.g. Lucida Console font.





Memory profiling refs:
* http://www.slideshare.net/PiotrPrzymus/pprzymus-europython-2014
* http://chase-seibert.github.io/blog/2013/08/03/diagnosing-memory-leaks-python.html
* http://www.asimihsan.com.s3.amazonaws.com/presentations/profiling_presentation/index.html


Memory profiling libs:
* Pympler - http://pythonhosted.org/Pympler/index.html, https://github.com/pympler/pympler/
* objgraph - http://mg.pov.lt/objgraph/
* Dowser - http://www.aminus.net/wiki/Dowser
* memprof - https://github.com/jmdana/memprof
* memory_profiler - https://github.com/fabianp/memory_profiler, https://pypi.python.org/pypi/memory_profiler (most mature and standard).
*



Check stacking:
"""
print("  %s -: :- %s \n  %s -: :- %s  " % (h1end3p, h1end5p, h2end5p, h2end3p))

"""


"""
