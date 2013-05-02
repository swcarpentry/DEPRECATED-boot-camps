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

Notes for things to include somewhere:
* Structs, cell arrays, objects, function handles, etc.  Some of this at least should be familiar.
* Distinction between nested & local functions?
* http://www.mathworks.co.uk/help/matlab/matlab_prog/add-reminders-to-files.html and relation to issue tracker
* http://www.mathworks.co.uk/help/matlab/error-handling.html
* Defensive programming: http://www.mathworks.co.uk/help/matlab/input-and-output-arguments.html
* TODO: Move my versions into a solutions folder

The data analysis
-----------------

Start by firing by MATLAB, and navigate to the `matlab` folder within
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
regression might not seem unreasonable.

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
full analysis pipeline and produce all the desired results figures &
statistics.  This is crucial for _reproducibility_.  Victoria Stodden
has written extensively about the idea of reproducibility in
scientific software - you may want to look up [some of her
papers][Stodden] for reference.

For our purposes, we can summarize the goal of reproducibility in two
related ways, one technical and one colloquial.

In a technical sense, your goal is to __have a complete chain of custody (ie, 
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

## Exercise

Create a new .m file defining a function to perform the whole
analysis, returning the results.

TODO: Use as an opportunity to loop over datasets in all_samples and
talk about structs.  (Also useful for returning results.)

## Exercise

Create a 'distance from line' function and incorporate this into your
analysis program.


TODO: Talk about functions for reuse, and for testing (lead into next section).


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

To run tests, we use xUnit's `runtests` function.  Much like
`nosetests`, this will automatically find tests matching its rules.

Notes:
* Go through examples similar to xUnit documentation.
* Mention suggested structure of projects, with data/tests folders as in repo.


Testing exercise
----------------

Take the data analysis script we developed in the first part of the
session, split it into functions wherever appropriate, and write unit
tests for these.

While normally you would assume that routines provided by MATLAB are
correct, you could for instance, write tests for the basic stats
methods (mean, var, etc.) and linear regression to verify these.  You
can also test the code for finding the nearest/furthest point from the
regression line, or consider what further analyses could be performed
and how to test these.  Another test is to verify the expected
statistical properties of the [Anscombe datasets][Anscombe], the first
four of our samples.

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
[code sections]: http://www.mathworks.co.uk/help/matlab/matlab_prog/run-sections-of-programs.html

[GF sqlite]: https://bitbucket.org/GarrettFoster/sqlite-matlab/src
[mksqlite]: http://mksqlite.berlios.de/mksqlite_eng.html

[xunit]: http://www.mathworks.com/matlabcentral/fileexchange/22846-matlab-xunit-test-framework
[R2013a testing]: http://blogs.mathworks.com/steve/2013/03/12/matlab-software-testing-tools-old-and-new-r2013a/

[Anscombe]: http://en.wikipedia.org/wiki/Anscombe%27s_quartet
[Stodden]: http://www.stanford.edu/~vcs/Papers.html
