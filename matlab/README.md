Software Carpentry techniques in MATLAB
=======================================

So far this workshop has been teaching you general principles of good
programming practice in the context of learning the Python programming
language.  However, these principles can be applied in any language,
and so in this session we'll demonstrate how they can be used to make
you more productive when using [MATLAB][].

As a setting we will imagine performing some exploratory analysis on a
few datasets, and build towards a robust, tested program for
summarising new, similar data.


The data analysis
-----------------

Notes:
* Walk through the main steps in analysis_script.m at the command line.
* Show how to save explorations to a script and re-run whole analysis
* Show how to run particular sections in cell mode too
* Show how to use functions to split things up for testing, reuse, etc.
* Brief mention of nested functions, function handles, etc.  But they should know all that!


Unit testing in MATLAB
----------------------

The latest release of MATLAB, R2013a, has [support for testing built
in][R2013a testing].  However, since many of you will still have 2012
or earlier releases, we're going to focus on a separate library
available from MATLAB Central, called [MATLAB xUnit][xunit].  This
provides similar features to nosetests, allowing you to write tests
for MATLAB functions, automatically find tests to run, compare using
floating-point tolerances, share test setup code, and so on.

Notes:
* Go through examples similar to xUnit documentation.
* Mention suggested structure of projects, with data/tests folders as in repo.


Exercise
--------

Take the data analysis script we developed in the first part of the
session, split it into functions wherever appropriate, and write unit
tests for these.

While normally you would assume that routines provided by MATLAB are
correct, you could for instance, write tests for the basic stats
methods (mean, var, etc.) and linear regression to verify these.  You
can also test the code for finding the nearest/furthest point from the
regression line, or consider what further analyses could be performed
and how to test these.

As an alternative exercise, try adding tests to some of your own
MATLAB code.


Postscript: accessing Sqlite databases from MATLAB
--------------------------------------------------

Just as there are Python libraries that allow you to work with a
database from your programs, there are also a few MATLAB interfaces.
They do not appear to be anywhere near the maturity of the Python
versions, but may still be useful enough for you, depending on your
needs.

[mksqlite] has no dependencies itself, and has a slightly better
interface.  Unfortunately it is a MEX library, and only 32-bit Windows
binaries are provided, so if you have a 64-bit MATLAB or don't use
Windows, you need to compile it from source yourself.

Alternatively [Garrett Foster has produced a simple interface][GF sqlite]
using the Database Toolkit.  If you have a licence for the toolkit this
is the simplest to get going, but is a bit clunky for building up SQL
queries programmatically.



[MATLAB]: http://www.mathworks.co.uk/products/matlab/
[GF sqlite]: https://bitbucket.org/GarrettFoster/sqlite-matlab/src
[mksqlite]: http://mksqlite.berlios.de/mksqlite_eng.html
[xunit]: http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
[R2013a testing]: http://blogs.mathworks.com/steve/2013/03/12/matlab-software-testing-tools-old-and-new-r2013a/
