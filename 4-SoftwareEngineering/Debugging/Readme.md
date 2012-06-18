# Debugging

[Back To Remote Version
Control](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/3b-VersionControlRemote/)
- [Forward To
Testing](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/5-Testing/)

**Presented and Designed by Anthony Scopatz**

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

## How I Debug

I write perfect code the first time.

## Errors, Exceptions, & Tracebacks

Bugs start where execution ends. In many modern languages, when an
invalid operation occurs an exception is *thrown* or *raised*. These
exceptions may be *handled* or *caught*. In Python, there are \~165
built-in exceptions.

```python
try:
    a = 1.0 / 0.0
except ZeroDivisionError as e:
    print "Going from zero to hero."
    a = 1.0
```

In languages that have functions and exceptions, you can typically get a
hold of what is known as a traceback. This displays the history of
function calls leading up to the error.

**Example:** `python tb_example.py`

**Other Resources:** [Exception
Handling](http://www.doughellmann.com/articles/how-tos/python-exception-handling/index.html)

## pdb

Following Python's motto of "batteries included", the language itself
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

## Coding Standards

Much like a written natural language, there are many ways to express the
same idea. The strict syntax of languages are necessarily more forgiving
than what the *correct* way of doing things (think Oxford comma). To
make the consumption of information easier, style guides exists to
enforce particularly effective ways of writing.

Coding standards fill the same role but for programming languages. They
become absolutely essential as projects become large (\>1 person). \`
Now, the wonderful thing about standards is that there are so many to
choose from!\`\_

Python is somewhat unique in that the language itself has an approved
coding standard called [PEP8](http://www.python.org/dev/peps/pep-0008/).
The overwhelming majority (80-90%) of Python code that is available on
the internet is written in a way that is PEP8-compliant. Unfortunately,
some of that 20% is in standard library...

Thus not adhering to your coding standard is often considered "against
best practices", ie a bug. Luckily there are tools to test for
compliance:

    pep8 style_example.py

## Segfaults

### The Scourge of {K&R, ANSI, ISO, 99, 11, Embedded, Objective} C!

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
which joined the Python 3.3 standard library.

**Example**:

    python segfault_unhandled_example.py
    python segfault_handled_example.py

**Other Resources:**
[faulthandler](https://github.com/haypo/faulthandler/wiki/),
[WAD](http://www.dabeaz.com/papers/Python2001/python.html), [HOWTO Crash
Python](http://wiki.python.org/moin/CrashingPython).

## Valgrind

Valgrind is a utility for compiled codes which aids in debugging,
finding memory leaks, and profiling. This is invaluable for codes
tracking down errors that only happen at runtime, such as segfaults.

As an example, first compile the following program without optimization.
For simpleTest.cc, run this line to see errors in this code:

    g++ simpleTest.cc -o simpleTest
    valgrind --track-origins=yes --leak-check=full ./simpleTest 300 300

We also have a cache test line. Run this line to see the cache errors:

    g++ cacheTest.cc
    valgrind --tool=cachegrind ./a.out 0 1000 100000

There are two paths in this code. If the first input is 1, it runs a
cache-sensitive version of the loop. If it is 0, it runs a
cache-insensitive version. The cache should look like:

    ~ $ dmesg | grep cache
    CPU: L1 I cache: 32K, L1 D cache: 32K
    CPU: L2 cache: 6144K
    CPU: L1 I cache: 32K, L1 D cache: 32K
    CPU: L2 cache: 6144K

You can run the same command to see cache on your linux machine. Another
way to see the exact cache setup that valgrind found is the following:

    cg_annotate --auto=yes cachegrind.out.21960

Note that your cachegrind.out will have a different number. This command
is also handy because it shows which functions caused cache misses.

**Other Resources:** [Valgrind](http://valgrind.org/)
