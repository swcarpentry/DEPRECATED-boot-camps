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

Start by firing up MATLAB, and navigate to the `matlab` folder within
your copy of the materials for this workshop.  Six sample datasets
have been saved in the `data` subfolder.  Let's start by loading the
first of these and see what we can discover about it.  Note that the
details of the analyses are not the most sensible approach to doing
this; I've chosen this sequence for didactic purposes!

```matlab
    load('data/sample1')
    whos
    mean(x), mean(y)
    var(x), mean(y)
    corrcoef(x, y)
```

There's a reasonably high correlation between x and y, so performing a
linear regression might not seem unreasonable.

```matlab
    p = polyfit(x, y, 1)
```

Finally, we'll do the step we should really have started with, and
plot the data itself, along with the regression line.

```matlab
    line_x = 1:19;
    line_y = polyval(p, line_x);
    plot(x, y, '+b', line_x, line_y, '-k');
```

Now, we'll want to perform this analysis for multiple datasets, and
it's tedious and error prone to type the commands repeatedly at the
command line.  Instead, we want one program we can execute to run the
full analysis pipeline and produce all the desired results figures and
statistics.  This is crucial for _reproducibility_.  Victoria Stodden
has written extensively about the idea of reproducibility in
scientific software - you may want to look up [some of her
papers][Stodden] for reference.

For our purposes, we can summarize the goal of reproducibility in two
related ways, one technical and one colloquial.

In a technical sense, your goal is to __have a complete chain of custody (i.e.,
provenance) from your raw data to your finished results and figures__. That is,
you should _always_ be able to figure out precisely what data and what code
were used to generate what result - there should be no "missing links". If you
have ever had the experience of coming across a great figure that you made
months ago and having no idea how in the world you made it, then you understand
why provenance is important. Or, worse, if you've ever been unable to recreate
the results that you once showed on a poster or (gasp) published in a
paper...

In a colloquial sense, I should be able to sneak into your lab late at night,
delete everything except for your raw data and your code, and __you should be
able to run a single command to regenerate EVERYTHING, including all of your
results, tables, and figures in their final, polished form__. Think of this as
the "push button" workflow. This is your ultimate organizational goal as a
computational scientist. Importantly, note that this rules out the idea of
manual intervention at any step along the way - no tinkering with figure axes
in a pop-up window, no deleting columns from tables, no copying data from one
folder to another, etc. All of that needs to be fully automated.

As an added bonus, if you couple this with a version control system that tracks
changes over time to your raw data and your code, you will be able to instantly
recreate your results from any stage in your research (the lab presentation
version, the dissertation version, the manuscript version, the Nobel Prize
committee version, etc.). Wouldn't that be nice?


So, let's save the commands we've run so far to a script as a starting
point.  Select an entry or entries in the command history list, and
then right-click and select Create Script from the context menu.  The
Editor opens a new file that contains the statements you selected.
Save the script with a suitable name in the `matlab` folder.  You can
now re-run all the commands just by pressing F5 (or clicking the
appropriate button).

MATLAB also has a nice feature when working with scripts, whereby you
can split a script into sections and just run one at a time.  These
are known as [code sections][].  Starting a comment line with a double
percent (`%%`) starts a section.  The current section (i.e. the one
containing your cursor) gets highlighted in yellow, and can be run
using the 'Run Section' button (Ctrl+Enter).  We'll split our script
into four sections: loading the data, summary statistics, linear
regression, and plotting.


So, now we can alter the first line stating what data to load, and run
different parts of our analysis easily, possibly making variations to
the script and re-running just those parts.  This is good for
exploring and experimenting, but not good for reproducibility: will
you remember to save the version that worked?  What if it only worked
due to running the sections in a different order?  Instead, we want to
transform the script into a function that takes the _x_ and _y_ data
to process as inputs.

### Exercise

Create a new .m file defining a function to perform the whole
analysis, returning the results.

### Structured data types

One consideration that arises from this exercise is how to return the
results.  MATLAB allows you to return multiple values from functions,
so you could use this feature to return each computed statistic as a
separate variable.  However, when the number of things to return is
more than a couple, this makes life difficult for users of the
function.  A better approach here is to return a [struct][] with
separate fields for each item.  Much like a Python dictionary, these
more structured data types are great for keeping related data together
and organised.

The data file `all_samples.mat` contains each of our test datasets
within an array of structs.  We can therefore loop over these to run
our analysis function on each in turn.

```matlab
    clear all
    load('data/all_samples'); % Contains a single struct array 'data'
    whos
    for i=1:length(data)
        figure;
        results = analysis_func(data(i).x, data(i).y);
        disp(results);
    end
```

While structs are useful for keeping related data together, they don't
allow you to associate behaviour with that data, or enforce
constraints to make sure that relationships between the data are
maintained.  This kind of functionality is the topic of [Object
Oriented Programming][OOP], and it's worth exploring MATLAB's (and
Python's) capabilities here if you're going to be writing any
moderately large programs.  We don't have time to go into it now
though.


The next topic we're going to cover is _testing_ your MATLAB code.
One of the general principles taught in Software Carpentry is that
breaking your code up into small functions is good for understanding
it.  It's also good for testing, as each small piece can be tested
individually.  When dealing with some legacy code not written with
testing in mind, identifying portions that can be split into functions
and tested is thus a good starting point, and that's what we'll do
next with our analysis program.  Firstly, however, we'll find out how
to write tests in MATLAB.


Unit testing in MATLAB
----------------------

The latest release of MATLAB, R2013a, has [support for testing built
in][R2013a testing].  However, since many of you will still have 2012
or earlier releases, we're going to focus on a separate library
available from MATLAB Central, called [MATLAB xUnit][xunit].  This
provides similar features to nosetests, allowing you to write tests
for MATLAB functions, automatically find tests to run, compare using
floating-point tolerances, share test setup code, and so on.

To make the xUnit functions available to our programs, we need to add
the library folder to MATLAB's search path.  This can be done in the
GUI, by expanding the 'matlab_xunit' folder, right-clicking on the
'xunit' folder, and selecting 'Add to path'.  Alternatively this can
be done from the command line using the `addpath` function, along the
lines of:
 ```matlab
    addpath([pwd '\matlab_xunit\xunit']);
```

Several examples are shipped with xUnit, so we'll go through a few of
these to demonstrate its abilities, and then set you loose to write
your own tests.  So, firstly we need to get into the directory
containing the examples.

```matlab
    cd matlab_xunit/doc
    cd example_quick_start
```

Tests are M-file functions that return no results, and the function
name must start or end with 'test' or 'Test'.  Test cases are
considered to pass if the function runs with no errors produced.  Your
tests should thus raise an error if the output is not as expected.
The first simple example checks that `fliplr` of a vector works.

```matlab
    type testFliplrVector
```

Raising errors yourself gets rather tedious, so as with other
frameworks xUnit provides various utility functions to perform common
checks.  `assertEqual` is demonstrated in the second example.

```matlab
    type testFliplrMatrix
```

To run tests, we use xUnit's `runtests` function.  Much like
`nosetests`, this will automatically find tests matching its rules.

```matlab
>> runtests
Test suite: C:\cygwin\home\jonc\work\git_repos\swc\boot-camps\matlab\matlab_xunit\doc\example_quick_start
02-May-2013 14:18:02

Starting test run with 2 test cases.
..
PASSED in 0.022 seconds.
```

Writing many files, one for each test, quickly makes file organisation
tricky.  As in Python, it's possible to combine multiple test cases in
a single file to make up a test suite.  The file name is the same, but
this time the function returns one variable, which must be named
`test_suite`.  It first calls the xUnit function `initTestSuite`, then
defines a subfunction for each test case, the names of which start or
end with 'test'.

```matlab
    cd ../example_subfunction_tests
    type testFliplr
    runtests
```

The `runtests` command also takes arguments allowing you to select
specific tests to run, rather than running all tests found.  This is
useful while developing code to make a particular test case or test
suite pass.

```matlab
    runtests testFliplr
    runtests testFliplr:testFliplrVector
    cd ..
    runtests example_subfunction_tests
```

### Exercise

Write tests that check the output of the `sin` function for a few
well-known values, including `pi`.  You may find the following list of
utility functions provided by xUnit helpful.
- assertTrue
- assertFalse
- assertEqual
- assertFilesEqual
- assertElementsAlmostEqual
- assertVectorsAlmostEqual
- assertExceptionThrown


Longer testing exercise
-----------------------

Take the data analysis script we developed in the first part of the
session, split it into functions wherever appropriate, and write unit
tests for these.

While normally you would assume that routines provided by MATLAB are
correct, you could for instance, write tests for the basic stats
methods (mean, var, etc.) and linear regression to verify these.
Another test is to verify the expected statistical properties of the
[Anscombe datasets][Anscombe], the first four of our samples.

Next, extend your testing to _test driven development_, where the
tests for a new feature are written _before_ implementing the feature!
Suppose we want to identify which points lie respectively nearest and
furthest from the regression line.  Write tests for a new function
`distances_from_line` which will compute all the distances, then write
the function and use it within your main analysis pipeline.  See
http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line for help
in working out the calculation.


For this exercise, it's easiest to place test files alongside your
source code, so that when running the tests they can find the
functions they're testing.  In a larger project, it makes more sense
to organise tests into their own subfolder (say 'tests') alongside
'source' and 'data' folders, and ensure that your source code is on
the MATLAB search path.  You can also use [packages][] to add further
structure to your code.


As an alternative exercise, try adding tests to some of your own
MATLAB code.


Further topics
--------------

There's obviously much more to good programming in MATLAB than we've
been able to cover in this session.  If you're going to be using it in
anger, it's worth becoming familiar with topics such as:

* The distinction between [local, nested and anonymous functions](http://www.mathworks.co.uk/help/matlab/matlab_prog/types-of-functions.html), and the use of function handles to pass functions to functions
* Defensive programming, for instance [checking that function inputs are as expected](http://www.mathworks.co.uk/help/matlab/input-and-output-arguments.html), and the related topic of [dealing with errors](http://www.mathworks.co.uk/help/matlab/error-handling.html)
* Keeping track of tasks with MATLAB's [reminder system](http://www.mathworks.co.uk/help/matlab/matlab_prog/add-reminders-to-files.html), as a simple form of issue tracker
* [Object Oriented Programming][OOP]


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
[code sections]: http://www.mathworks.co.uk/help/matlab/matlab_prog/run-sections-of-programs.html
[OOP]: http://www.mathworks.co.uk/help/matlab/object-oriented-programming.html
[struct]: http://www.mathworks.co.uk/help/matlab/structures.html
[packages]: http://www.mathworks.co.uk/help/matlab/matlab_oop/scoping-classes-with-packages.html

[GF sqlite]: https://bitbucket.org/GarrettFoster/sqlite-matlab/src
[mksqlite]: http://mksqlite.berlios.de/mksqlite_eng.html

[xunit]: http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
[R2013a testing]: http://blogs.mathworks.com/steve/2013/03/12/matlab-software-testing-tools-old-and-new-r2013a/

[Anscombe]: http://en.wikipedia.org/wiki/Anscombe%27s_quartet
[Stodden]: http://www.stanford.edu/~vcs/Papers.html
