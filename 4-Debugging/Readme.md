# Debugging

[Back To Remote Version
Control](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/3b-VersionControlRemote/)
- [Forward To
Testing](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/5-Testing/)

**Presented and Designed by Anthony Scopatz**

Below are some tools that are useful for debugging. We'll need to do
some extra installation to have all of them available:

    cd ~/UofCSCBC2012/4-Debugging/
    sudo ./install.sh

## What is Debugging (Exercise)?

Before I show you the practice (art) of debugging, separate out into
groups of 2-3 people. Follow these steps:

1.  **Come up with a definition of debugging.**
2.  **Write it down on a strip of paper.**
3.  **Put paper in hat.**
4.  **???**
5.  **Profit.**

(Bonus Challenge: Make a new friend!)

Time limit: 5 min.

## pdb

Following Python's moto of "batteries included", the language itself
comes packaged with its own aptly named Python DeBugger (pdb). From any
Python code anywhere, simply make sure that pdb is imported and then
call the `set_trace()` function.

```python
import pdb
pdb.set_trace()
```

This drops you into debugging mode, where you have *exactly* the state
of the program that the trace was set at!

Once inside the debugger, 'l(ist)' will list the commands available,
'h(elp)' will give you help on those commands, and 'q(uit)' exits the
debugger.

**Example:** `python pdb_example.py`

**Other Resources:** [PDB
Docs](http://docs.python.org/library/pdb.html),
[Wingware](http://wingware.com/doc/debug/advanced),
[O'Reily](http://onlamp.com/pub/a/python/2005/09/01/debugger.html),
[Great
Blog](http://pythonconquerstheuniverse.wordpress.com/category/the-python-debugger/).

## Profiling

Various profiling tools exist for every language out there. However, the
general idea is always the same. Different parts of your code take up
different amounts of the processing time. Your (human) time is limited.
Therefore, you should focus on optimizing/fixing/etc. only the most
important parts of your code. You discover which parts are most
important by using a profiler.

Here we will be using kernprof, a line profiler by Robert Kern.

**Example:** `kernprof.py -lv profiler_example.py`

**Other Resources:**
[kernprof](http://packages.python.org/line_profiler/)

## Linting

Linting in software is the process of discovering errors in a code
(typically typos and syntax errors) before the code is ever run or
compiled. Some people use such power automatically by not allowing their
version control to check in code unless the linter passes. Or they have
the linter run each time they exit their text editor. *This is not
crazy!*

In Python there are two main libraries that help in this regard:
**pylint** and **pyflakes**. These codes work by statically analyzing
the parse tree, rather than importing and running the module. We'll be
talking about pyflakes.

**Example:** `pyflakes lint_example.py`

**Other Resources:** [basic
definition](http://en.wikipedia.org/wiki/Lint_(software)),
[pyflakes](http://pypi.python.org/pypi/pyflakes/),
[pylint](http://www.logilab.org/857),
[comparison](http://www.doughellmann.com/articles/pythonmagazine/completely-different/2008-03-linters/).

## Segfaults

### The Scourge of C{K&R, ANSI, ISO, 99, 11, Embedded, Objective}!

Segmentation faults (*segfaults*) are some of the most obscure, most
annoying, and most difficult to debug errors in existence. This is
because they are a function of the state of the computer's RAM or
virtual memory at runtime.

Segfaults occur when the program tries to access a part of memory that
it expects to be able to get to, and for whatever reason it is not
available. At this point the code cannot continue and typically just
prints out `Segmentation fault` to the screen.

As the above error message does not indicate *where* in the execution
the segfault occurred, it very could have been anywhere. However, all
hope is not lost! Even high-level languages like Python have ways of
handling segfaults made on the C/C++/Fortran level and turning them into
standard exceptions. A great module for doing this is faulthandler,
which joined the Python standard library in v3.3.

**Example::**

> python segfault\_unhandled\_example.py python
> segfault\_handled\_example.py

**Other Resources:**
[faulthandler](https://github.com/haypo/faulthandler/wiki/),
[WAD](http://www.dabeaz.com/papers/Python2001/python.html), [HOWTO Crash
Python](http://wiki.python.org/moin/CrashingPython).

Compile each program without optimization first.

For simpleTest.cc, run this line to see errors in this code.

    valgrind --track-origins=yes --leak-check=full ./simpleTest 300 300

We also have a cache test line. Run this line to see the cache errors.

    valgrind --tool=cachegrind ./a.out 0 1000 100000

There are two paths in this code. If the first input is 1, it runs a
cache-sensitive version of the loop. If it is 0, it runs a
cache-insensitive version.

FYI: on the Trieste lab machines, this is what cache looks like:

    guy ~>dmesg | grep cache
    CPU: L1 I cache: 32K, L1 D cache: 32K
    CPU: L2 cache: 6144K
    CPU: L1 I cache: 32K, L1 D cache: 32K
    CPU: L2 cache: 6144K

You can run the same command to see cache on your linux machine. Another
way to see the exact cache setup that valgrind found is the following:

    cg_annotate --auto=yes cachegrind.out.21960

Note that your cachegrind.out will have a different number. This command
is also handy because it shows which functions caused cache misses.
