Day 2 -- many miscellaneous topics
==================================

You should all start reading `XKCD <http://xkcd.com>`__.

Python data structures
----------------------

See `Python data structures <Python2/index.html>`__

Installing Python packages; useful Python packages
--------------------------------------------------

See :doc:`python-packages`

Useful UNIX tools
-----------------

ssh and screen; BLAST.

`SSH instructions here <SSH/index.html>`__

`Installing and using command line BLAST <SSH/blast.html>`__

Testing
-------

When we say "testing" we really mean *automated testing*.
The central problems addressed by testing are *correctness* and
*reproducibility*.  (While these are linked, they are not the
same!)

There are two basic kinds of tests that I'd like to briefly
discuss.  One kind of test is a *unit test*.  The other kind
of test is a *regression test*.  (There are also many more.)

Unit tests address *small units of code*, like functions.  They
are used to isolate and nail down and prove the functionality
of potentially complicated little functions.

Regression tests address *the overall function of code*, and
they are used to make sure that your code is doing the same
thing *today* as it was *yesterday*.

I'll show you examples of both, but quickly :).

Writing tests
~~~~~~~~~~~~~

We're going to be using the nose testing framework, which is
just a framework that makes it easy to find and execute
tests.

Basically, 'nose' creates a command 'nosetests' that finds and
runs tests.  The idea is that you won't need to register new tests.

A test function looks like this::

   def test_something():
      # run some code
      # fail loudly or succeed silently

(Notebook will go here.)

More reading
~~~~~~~~~~~~

For more reading, see:

   http://software-carpentry.org/4_0/test/

and

   http://ivory.idyll.org/articles/nose-intro.html

and also

   http://ivory.idyll.org/blog/software-quality-death-spiral.html

Version Control
---------------

(An abbreviated version of: http://ged.msu.edu/angus/git-intro.html)

The purpose of version control is to serve as a method for tacking
changes to files, which enables lots of things:

1. track changes (wanted and unwanted) in your files.
2. keep track of an entire history of changes.
3. track multiple independent "branches" of work.
4. collaborate sanely.

----

A brief introduction:

1. Web interface and editing files.

2. One-person repositories.

3. One-person repositories with multiple locations.

4. One-person repositories with branches.

5. Collaborations!

Pipelines
---------

Operating on more complex files.

Doing stuff with BLAST output.
