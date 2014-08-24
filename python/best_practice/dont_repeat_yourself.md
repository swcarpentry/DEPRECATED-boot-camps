[Up To Schedule](../../README.md) -
Back To [Make Incremental Changes I](../../version-control/git/local/Readme.md) - 
Forward To [Make Incremental Changes II](../../version-control/git/local/Revert_and_branch.md)

- - - - 

# Don't repeat yourself (or others): Python Modules

**Based on Lecture Materials By: Milad Fatenejad, Katy Huff, and Paul Wilson**

In this lecture we shift gears slightly.  The best practices of "Write Code
for People" continue to be important, but we'll start to place more emphasis
on the best practices of "Don't Repeat Yourself", and its twin practice of
"Don't Repeat Others".

We've already seen the importance of loops as a way to avoid repeating
yourself.  We've also already discussed the use of functions.  In this lesson,
we'll explore how to use other people's modules, and how to package pieces of
your code together into reusable modules.

## Reminder: Pasting into iPython

To paste text from another application (i.e. these lecture notes) into
iPython :

1.  select text from the wiki
2.  copy with <kbd>cntl</kbd>+<kbd>c</kbd>
3.  unfortunately pasting depends on your operating system and ssh program:

### Windows

#### Putty
Click with the right mouse button over the window.

#### Bitvise SSH
Click with the right mouse button over the window and then "Paste".

### Mac OSX
Press cmd+C.

### Linux
Click with the right mouse button over the window and then "Paste".


The code should paste and execute in iPython.

If you also type `%autocall` to turn autocall OFF, you may be able to paste
with <kbd>ctrl</kbd>+<kbd>v</kbd> though this won't work with all iPython
builds.

# Python Modules

Python has a lot of useful data type and functions built into the language, some of which you have already seen. For a full list, you can type `dir(__builtins__)`. However, there are even more functions stored in modules. An example is the sine function, which is stored in the math module. In order to access mathematical functions, like `sin`, we need to `import` the math module. Let's take a look at a simple example:

```python

print sin(3) # Error! Python doesn't know what sin is...yet

import math # Import the math module
math.sin(3)

print dir(math) # See a list of everything in the math module

help(math) # Get help information for the math module
```

It is not very difficult to use modules - you just have to know the module name and import it. There are a few variations on the import statement that can be used to make your life easier. Let's take a look at an example:

```python

from math import *  # import everything from math into the global namespace (A BAD IDEA IN GENERAL)
print sin(3)        # notice that we don't need to type math.sin anymore
print tan(3)        # the tangent function was also in math, so we can use that too
```

```python
reset # Clear everything from IPython

from math import sin  # Import just sin from the math module. This is a good idea.
print sin(3)          # We can use sin because we just imported it
print tan(3)          # Error: We only imported sin - not tan
```

```python
reset                 # Clear everything
import math as m      # Same as import math, except we are renaming the module m
print m.sin(3)        # This is really handy if you have module names that are long
```

If you intend to use python in your workflow, it is a good idea to skim the standard library documentation at the main Python documentation site, [docs.python.org](http://docs.python.org) to get a general idea of the capabilities of python "out of the box".

Let's take a look at some nice docstrings:

```python
import numpy
numpy.sum??

Type:       function
String Form:<function sum at 0x0000000002D1C5F8>
File:       c:\anaconda\lib\site-packages\numpy\core\fromnumeric.py
Definition: numpy.sum(a, axis=None, dtype=None, out=None, keepdims=False)
Source:
def sum(a, axis=None, dtype=None, out=None, keepdims=False):
    """
    Sum of array elements over a given axis.

    Parameters
    ----------
    a : array_like
        Elements to sum.
    axis : None or int or tuple of ints, optional
        Axis or axes along which a sum is performed.
        The default (`axis` = `None`) is perform a sum over all
        the dimensions of the input array. `axis` may be negative, in
        which case it counts from the last to the first axis.

        .. versionadded:: 1.7.0

        If this is a tuple of ints, a sum is performed on multiple
        axes, instead of a single axis or all the axes as before.
    dtype : dtype, optional
        The type of the returned array and of the accumulator in which
        the elements are summed.  By default, the dtype of `a` is used.
        An exception is when `a` has an integer type with less precision
        than the default platform integer.  In that case, the default
        platform integer is used instead.
    out : ndarray, optional
        Array into which the output is placed.  By default, a new array is
        created.  If `out` is given, it must be of the appropriate shape
        (the shape of `a` with `axis` removed, i.e.,
        ``numpy.delete(a.shape, axis)``).  Its type is preserved. See
        `doc.ufuncs` (Section "Output arguments") for more details.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.

    Returns
    -------
    sum_along_axis : ndarray
        An array with the same shape as `a`, with the specified
        axis removed.   If `a` is a 0-d array, or if `axis` is None, a scalar
        is returned.  If an output array is specified, a reference to
        `out` is returned.

    See Also
    --------
    ndarray.sum : Equivalent method.

    cumsum : Cumulative sum of array elements.

    trapz : Integration of array values using the composite trapezoidal rule.

    mean, average

    Notes
    -----
    Arithmetic is modular when using integer types, and no error is
    raised on overflow.

    Examples
    --------
    >>> np.sum([0.5, 1.5])
    2.0
    >>> np.sum([0.5, 0.7, 0.2, 1.5], dtype=np.int32)
    1
    >>> np.sum([[0, 1], [0, 5]])
    6
    >>> np.sum([[0, 1], [0, 5]], axis=0)
    array([0, 6])
    >>> np.sum([[0, 1], [0, 5]], axis=1)
    array([1, 5])

    If the accumulator is too small, overflow occurs:

    >>> np.ones(128, dtype=np.int8).sum(dtype=np.int8)
    -128

    """

```

We see a nice docstring seperating into several sections. A short description of the function is given, then all the input parameters are listed, then the outputs, there are some notes and examples. 
Please note the docstring is longer than the code. And there are few comments in the actual code.

* * * *
![Exercise](pics/exercise.jpg) **Short exercise: Learn about sys**

In the script we wrote to convert text files to CSV, we used the `sys` module.  Use the iPython interpreter to learn more about the `sys` module and what it does.  What is `sys.argv` and why did we only use the last n-1 elements? What is `sys.stdout`?

* * * *

## The python CSV Module

In the last lecture, we wrote code to write out some structured information in
a CSV formatted file.  In that case, we enclosed every element of the table in
quotes.  For most tools that read CSV files, this causes those elements to
always be interpretted a strings.  However, we had many numerical entries and
may want those to be recorded as numbers.

There are three ways to do this:
1. Rewrite our functions to only enclose the string elements in quotation
   marks for this particular data set.  Straightforward, but only useful in
   this particular project.
2. Rewrite our functions to detect which elements are really strings and which
   are not.  Probably quite difficult and error prone.
3. Find a module that does this for us already.  Clearly the best choice!!!

In fact, there is a CSV module already available.

![Exercise](pics/exercise.jpg) **Short Exercise** Use google to search for a `python csv module`.

Instead of having you learn about the CSV module from the documentation and
examples, we'll point you to the most important things:

```
In [25]: import csv
In [26]: dir(csv)
In [27]: help(csv)
In [28]: dir(csv.DictWriter)
In [29]: csv.DictWriter?
```

The truth is that the documentation provided here is probably not enough to
learn it by yourself.  The online documentation is better, and even better are
some examples that you can find online.

We'll start by adding the following to our file:

```python
import csv

csv_writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=column_labels)
```

This provides a way to write CSV files from python dictionaries.  More specifically:
* we will write this to the screen using `sys.stdout`
* we will use a comman (',') as the delimiter (Note: we can remove any reference to `csv_separator` now)
* we will use `column_labels` to provide the order in which to write the fields from our data

Now we can replace the two functions we added previously:
* instead of `writeCSVHeader(column_labels,csv_separator)` we have `csv_writer.writeheader()`
* instead of `writeCSVRow(column_labels,data_record,csv_separator)` we have `csv_writer.writerow(data_record)`

So that our script now ends with the following lines:

```python
import csv
csv_writer = csv.DictWriter(sys.stdout, delimiter=',', fieldnames=column_labels)

csv_writer.writeheader()
csv_writer.writerows(all_data)
```

![Exercise](pics/exercise.jpg) **Try it now!**

![Exercise](pics/exercise.jpg) **Short Exercise:** Figure out how to change the quoting behavior of the CSV writer and explore different options.  Do any of them put the strings in quotations but not the numbers?

![Exercise](pics/exercise.jpg) **Bonus Exercise:** The CSV module even has a method to write multiple rows: `writerows()`.  Try using it in stead of the loop over `all_data`.

We still need to ensure that the numerical results in our original files are
being treated as numbers.  To do so, immediately before returning our value
from `extractData()` we'll:
* define a list of columns that should be treated as numeric
* if we find one, convert it to a float

```python
    numeric_columns = ("CI type","Volume","Range","Discrimination")
    if key in numeric_columns and value != "" :
        value = float(value)
```

# Writing Your Own Module

For our set of cochlear implant data, we may now be interested in performing
some simple statistical analysis of the results.  For example, we may want to
know the mean value of the `Volume` data from all the subjects.

Since we know that taking the mean value of many numbers (and other
statistics) is something we may want to do in many different projects in the
future, let's make a module that contains those functions.  (Note: Most of
these functions are already available in existing modules and it would be best
to use those, but this is a convenient example that everyone probably can
understand equally well.)

Making a new module is as simple as defining functions in a new python file.
Let's call our file `stats.py` and start by adding a function to calculate the
mean.

```python
def mean(vals):
    """Calculate the arithmetic mean of a list of numbers in vals"""
    total = sum(vals)
    length = len(vals)
    return total/length
```

We can now use this module in our original script.

![Exercise](pics/exercise.jpg) **Exercise:**

1. Add some lines to the original script to get a list with only the `Volume` data.
2. Use this `mean()` function to calculate the mean of those numbers.

* * * *


##The General Problem##

![xkcd](http://imgs.xkcd.com/comics/the_general_problem.png "I find that when
 someone's taking time to do something right in the present, they're a
 perfectionist with no ability to prioritize, whereas when someone took time
 to do something right in the past, they're a master artisan of great
 foresight.")

From [xkcd](http://www.xkcd.com)
 
Now that you can write your own functions, you too will experience the dilemma
of deciding whether to spend the extra time to make your code more general,
and therefore more easily reused in the future.


- - - -

[Up To Schedule](../../README.md) -
Back To [Make Incremental Changes I](../../version-control/git/local/Readme.md) - 
Forward To [Make Incremental Changes II](../../version-control/git/local/Revert_and_branch.md)
