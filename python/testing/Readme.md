[Up To Schedule](../../README.md) - Back To [Make Incremental Changes](../../version-control/git/local/Readme.md) -
Forward To [Collaborate](../../version-control/git/remote/Readme.md)

# Plan for Mistakes (or: Testing)
----

**Based on materials by Katy Huff, Rachel Slaybaugh, and Anthony
Scopatz**

# What is testing?

Software testing is a process by which one or more expected behaviors
and results from a piece of software are exercised and confirmed. Well
chosen tests will confirm expected code behavior for the extreme
boundaries of the input domains, output ranges, parametric combinations,
and other behavioral **edge cases**.

# Why test software?

Unless you write flawless, bug-free, perfectly accurate, fully precise,
and predictable code **every time**, you must test your code in order to
trust it enough to answer in the affirmative to at least a few of the
following questions:

-   Does your code work?
-   **Always?**
-   Does it do what you think it does? ([Patriot Missile Failure](http://www.ima.umn.edu/~arnold/disasters/patriot.html))
-   Does it continue to work after changes are made?
-   Does it continue to work after system configurations or libraries
    are upgraded?
-   Does it respond properly for a full range of input parameters?
-   What about **edge and corner cases**?
-   What's the limit on that input parameter?
-   How will it affect your
    [publications](http://www.nature.com/news/2010/101013/full/467775a.html)?

## Verification

*Verification* is the process of asking, "Have we built the software
correctly?" That is, is the code bug free, precise, accurate, and
repeatable?

## Validation

*Validation* is the process of asking, "Have we built the right
software?" That is, is the code designed in such a way as to produce the
answers we are interested in, data we want, etc.

## Uncertainty Quantification

*Uncertainty Quantification* is the process of asking, "Given that our
algorithm may not be deterministic, was our execution within acceptable
error bounds?" This is particularly important for anything which uses
random numbers, eg Monte Carlo methods.

# Where are tests?

Let's return the the simplestats module that we began earlier:

    cd ~/simplestats

It has a file `stats.py` with the following averaging function:

```python
def mean(numlist):
    """Calculate the arithmetic mean of a list of numbers in numlist"""
    total = sum(numlist)
    length = len(numlist)
    return total/length
```

## Practice using git

Make a branch in this repository for adding tests:

    git checkout -b add_tests

or

    git branch add_tests
    git checkout add_tests

The simplest way to add a test is to add a function that calls this function
with arguments for which we already know the answer.

```python
def mean(numlist):
    """Calculate the arithmetic mean of a list of numbers in numlist"""
    total = sum(numlist)
    length = len(numlist)
    return total/length
    
def test_mean():
    """Test some standard behavior of the mean() function."""
    assert(mean([2, 4]) == 3)
```

The `assert` command will make your program stop if the condition is not true
and is common when writing tests.

You can try this test in iPython:

```
In [1]: import stats as s
In [2]: s.test_mean()
```

**Practice using git:** Commit this addition to the repository

    git add stats.py
    git commit -m "Added the first test"

Let's add one more test:

```python
def test_float_mean():
    """Test some standard behavior when the result is not an integer."""
    assert(mean([1, 2]) == 1.5)
```

and try it with:

```
In [4]: reload(s)
In [5]: s.test_float_mean()
```

The newest test fails, but we don't get much explanation of why it fails.

**Practice using git:** Commit this change to the repository

    git add stats.py
    git commit -m "Added a floating point test, but if fails"

# Separating Tests

It is more common to place tests in a different file so that they don't
clutter the module that does the real work.  Let's move our tests to a new
file called `test_stats.py`.  Now, our `stats.py` file contains only:

```python
def mean(numlist):
    """Calculate the arithmetic mean of a list of numbers in numlist"""
    total = sum(numlist)
    length = len(numlist)
    return total/length
```

and our `test_stats.py` file contains:

```python
from stats import mean

def test_mean():
    """Test some standard behavior of the mean() function."""
    assert(mean([2, 4]) == 3)

def test_float_mean():
    """Test some standard behavior when the result is not an integer."""
    assert(mean([1, 2]) == 1.5)
```

To make it even easier to test, we can add some lines at the bottom of
`test_stats.py` to run each of our tests:

```python
test_mean()
test_float_mean()
```

and then run this from the command-line:

    python test_stats.py

The same tests pass and fail, but still not much explanation.

**Practice using git:** Commit this change to the repository

    git add stats.py test_stats.py
    git commit -m "Moved tests to a separate file to declutter module."

We could start adding some lines to give us more information about each test
and why it might fail, but that could get tedious as we write basically the
same things over and over for each test.  Since we don't want to repeat
ourselves, we might write some functions to keep track of the expected result,
and report when it doesn't match the observed results.  However, that seems
like something that many people need, so maybe someone else did that already,
and we don't want to repeat others, either.

# Nose: A Python Testing Framework

The testing framework we'll discuss today is called `nose`. However, there are
several other testing frameworks available in most language. Most notably there
is [JUnit](http://www.junit.org/) in Java which can arguably attributed to
inventing the testing framework. Google also provides a [test
framework](http://code.google.com/p/googletest/) for C++ applications (note, there's
also [CTest](http://cmake.org/Wiki/CMake/Testing_With_CTest)).  There
is at least one testing framework for R:
[testthat](http://cran.r-project.org/web/packages/testthat/index.html).

## Where do nose tests live?

Nose tests are files that begin with `Test-`, `Test_`, `test-`, or
`test_`. Specifically, these satisfy the testMatch regular expression
`[Tt]est[-_]`. (You can also teach `nose` to find tests by declaring them
in the unittest.TestCase subclasses that you create in your code. You
can also create test functions which are not unittest.TestCase
subclasses if they are named with the configured testMatch regular
expression.)

## Nose Test Syntax

To write a nose test, we use the same assertions as we've used so far.

```python
assert should_be_true()
assert not should_not_be_true()
```

Since the `nose` package finds the tests and runs them automatically, we don't
need to include lines that run the tests.  We can remove the lines that call
the tests and just use our existing `test_stat.py` file like this:

    nosetests test_stat.py

We get a little more information, but still not that helpful.

**Practice using git:** Commit this change to the repository

    git add test_stats.py
    git commit -m "Introduced nose testing"

However, nose itself defines number of convenient assert functions which can
be used to test more specific aspects of the code base.

```python
from nose.tools import assert_equal, assert_almost_equal, assert_true, assert_false, assert_raises, assert_is_instance

assert_equal(a, b)
assert_almost_equal(a, b)
assert_true(a)
assert_false(a)
assert_raises(exception, func, *args, **kwargs)
assert_is_instance(a, b)
# and many more!
```

## Short Exercise

1. Change the two tests to use `assert_equal` and run with nosetests.
2. Fix the `mean()` function to resolve the test

Notice how much useful information you get from `nose` tests:
* some .... to indicate progress
* details about the failed test including the values that were not equal
* the total number of tests that were completed
* the time it took to run those tests

**Practice using git:** Commit this change to the repository

    git add test_stats.py
    git commit -m "Using nose tools to get even better output"


## What should I test?

One of the challenges of testing is to determine what the *edge cases* might
be.  Here are some cases that we might try for the mean function

* different lengths of lists:
    * [2, 4, 6]
    * [2, 4, 6, 8]
* negative numbers:
    * [-4, -2]
    * [-2, 2, 4]
* floating point numbers:
    * [2.0, 4.0, 6.0]
    * [2.5, 4.5, 6.0]

## Short Exercise

* Add tests for each of by adding either new functions or new lines to the
  existing function.  Note: some of them may fail!
* What is necessary to fix the failing tests?  Try using `assert_almost_equal`
  and make it pass the test.

**Practice using git:** Commit this change to the repository

    git add test_stats.py
    git commit -m "Added many more tests"

# Planning for bigger mistakes

What happens if someone tries to use this function with strings?

```python
mean(['hello','world'])
```

Some mistakes don't just give a wrong answer, but fail to even finish.  Python
provides a mechanism to deal with this called **exceptions**.

In this example, python automatically *raises* a TypeError exception when we
try to take the sum of a string.

In this case, we can add a test for the expected behavior: raising a TypeError
exception, by using the nose tool `assert_raises`:

```python
def test_string_mean():
    assert_raises(TypeError, mean, ['hello','world'])
```

**Practice using git:** Commit this change to the repository

    git add test_stats.py
    git commit -m "Added a test for exceptions when passing in non-numeric results."

We can provide some extra information to the user by catching the TypeError exception:

```python
def mean(numlist):
    """Calculate the arithmetic mean of a list of numbers in numlist"""
    try:
        total = sum(numlist)
        length = len(numlist)
    except TypeError:
        raise TypeError("The list contains non-numeric elements")
    return total/length
```
**Practice using git:** Commit this change to the repository

    git add test_stats.py
    git commit -m "Added extra error message for TypeError."


# When should we test?

The three right answers are:

-   **ALWAYS!**
-   **EARLY!**
-   **OFTEN!**

The longer answer is that testing either before or after your software
is written will improve your code, but testing after your program is
used for something important is too late.

If we have a robust set of tests, we can run them before adding
something new and after adding something new. If the tests give the same
results (as appropriate), we can have some assurance that we didn't
wreck anything. The same idea applies to making changes in your system
configuration, updating support codes, etc.

Another important feature of testing is that it helps you remember what
all the parts of your code do. If you are working on a large project
over three years and you end up with 200 classes, it may be hard to
remember what the widget class does in detail. If you have a test that
checks all of the widget's functionality, you can look at the test to
remember what it's supposed to do.

# Who should test?

In a collaborative coding environment, where many developers contribute
to the same code base, developers should be responsible individually for
testing the functions they create and collectively for testing the code
as a whole.

Professionals often test their code, and take pride in test coverage,
the percent of their functions that they feel confident are
comprehensively tested.

* * * * *

# Test Driven Development

Test driven development (TDD) is a philosophy whereby the developer
creates code by **writing the tests first**. That is to say you write the
tests *before* writing the associated code!

This is an iterative process whereby you write a test then write the
minimum amount code to make the test pass. If a new feature is needed,
another test is written and the code is expanded to meet this new use
case. This continues until the code does what is needed.

TDD operates on the YAGNI principle (You Ain't Gonna Need It). People
who diligently follow TDD swear by its effectiveness. This development
style was put forth most strongly by [Kent Beck in
2002](http://www.amazon.com/Test-Driven-Development-By-Example/dp/0321146530).

## A TDD Example

Say you want to write a std() function which computes the
[Standard Deviation](http://en.wikipedia.org/wiki/Standard_deviation). You
would - of course - start by writing the test, possibly testing a single set
of numbers, by adding this to `test_stats.py`:

```python
from nose.tools import assert_equal, assert_almost_equal

def test_std1():
    obs = std([0.0, 2.0])
    exp = 1.0
    assert_equal(obs, exp)
```

You would *then* go ahead and write the actual function:

```python
def std(vals):
    # you snarky so-and-so
    return 1.0
```

**Practice using git:** Since this works, commit the change to the repository

    git add stats.py test_stats.py
    git commit -m "Added a test for std() and then a function that passes the test."

And that is it, right?! Well, not quite. This implementation fails for
most other values. Adding tests we see that:

```python
def test_std1():
    obs = std([0.0, 2.0])
    exp = 1.0
    assert_equal(obs, exp)

def test_std2():
    obs = std([])
    exp = 0.0
    assert_equal(obs, exp)

def test_std3():
    obs = std([0.0, 4.0])
    exp = 2.0
    assert_equal(obs, exp)
```

These extra tests now require that we bother to implement at least a slightly 
more reasonable function:

```python
def std(vals):
    # a little better
    if len(vals) == 0:
        return 0.0
    return vals[-1] / 2.0
```

**Practice using git:** Since this works again, commit the change to the repository

    git add stats.py test_stats.py
    git commit -m "Added moore tests for std() and updated function so that is passes all tests."

However, this function still fails whenever vals has more than two elements or
the first element is not zero. Time for more tests!

```python
def test_std1():
    obs = std([0.0, 2.0])
    exp = 1.0
    assert_equal(obs, exp)

def test_std2():
    obs = std([])
    exp = 0.0
    assert_equal(obs, exp)

def test_std3():
    obs = std([0.0, 4.0])
    exp = 2.0
    assert_equal(obs, exp)

def test_std4():
    obs = std([1.0, 3.0])
    exp = 1.0
    assert_equal(obs, exp)

def test_std5():
    obs = std([1.0, 1.0, 1.0])
    exp = 0.0
    assert_equal(obs, exp)
```

At this point, we had better go ahead and try do the right thing...

```python
def std(vals):
    # finally, some math
    n = len(vals)
    if n == 0:
        return 0.0
    mu = sum(vals) / n
    var = 0.0
    for val in vals:
        var = var + (val - mu)**2
    return (var / n)**0.5
```

**Practice using git:** Since this works again, commit the change to the repository

    git add stats.py test_stats.py
    git commit -m "Added more tests for std() and updated function so that is passes all tests."

Here it becomes very tempting to take an extended coffee break or
possibly a power lunch. But then you remember those pesky infinite values!
Perhaps the right thing to do here is to just be undefined.  Infinity in 
Python may be represented by any literal float greater than or equal to 1e309.

```python
def test_std1():
    obs = std([0.0, 2.0])
    exp = 1.0
    assert_equal(obs, exp)

def test_std2():
    obs = std([])
    exp = 0.0
    assert_equal(obs, exp)

def test_std3():
    obs = std([0.0, 4.0])
    exp = 2.0
    assert_equal(obs, exp)

def test_std4():
    obs = std([1.0, 3.0])
    exp = 1.0
    assert_equal(obs, exp)

def test_std5():
    obs = std([1.0, 1.0, 1.0])
    exp = 0.0
    assert_equal(obs, exp)

def test_std6():
    obs = std([1e500])
    exp = NotImplemented
    assert_equal(obs, exp)

def test_std7():
    obs = std([0.0, 1e4242])
    exp = NotImplemented
    assert_equal(obs, exp)
```

This means that it is time to add the appropriate case to the function
itself:

```python
def std(vals):
    # sequence and you shall find
    n = len(vals)
    if n == 0:
        return 0.0
    mu = sum(vals) / n
    if mu == 1e500:
        return NotImplemented
    var = 0.0
    for val in vals:
        var = var + (val - mu)**2
    return (var / n)**0.5
```

**Practice using git:** Since this works again, commit the change to the repository

    git add stats.py test_stats.py
    git commit -m "Added tests for infinity in std() and updated function so that is passes all tests."


## Quality Assurance Exercise

Can you think of other tests to make for the ``std()`` function? I promise there
are at least two.
<!---
	1. How about std(string) or std(array)?
	2. How about std(None)?
--->

Implement one new test in test_stats.py, run nosetests, and if it fails, implement
a more robust function for that case.

And thus - finally - we have a robust function together with working
tests!

## Exercise: A different function

Try your new test-driven development chops by implementing the ``var()``
function, noting that the variance is the square of the standard devation.

# How are tests written?

The type of tests that are written is determined by the testing
framework you adopt. Don't worry, there are a lot of choices.

## Types of Tests

**Exceptions:** Exceptions can be thought of as type of runtime test.
They alert the user to exceptional behavior in the code. Often,
exceptions are related to functions that depend on input that is unknown
at compile time. Checks that occur within the code to handle exceptional
behavior that results from this type of input are called Exceptions.

**Unit Tests:** Unit tests are a type of test which test the fundamental
units of a program's functionality. Often, this is on the class or
function level of detail. However what defines a *code unit* is not
formally defined.

To test functions and classes, the interfaces (API) - rather than the
implementation - should be tested. Treating the implementation as a
black box, we can probe the expected behavior with boundary cases for
the inputs.

**System Tests:** System level tests are intended to test the code as a
whole. As opposed to unit tests, system tests ask for the behavior as a
whole. This sort of testing involves comparison with other validated
codes, analytical solutions, etc.

**Regression Tests:** A regression test ensures that new code does
change anything. If you change the default answer, for example, or add a
new question, you'll need to make sure that missing entries are still
found and fixed.

**Integration Tests:** Integration tests query the ability of the code
to integrate well with the system configuration and third party
libraries and modules. This type of test is essential for codes that
depend on libraries which might be updated independently of your code or
when your code might be used by a number of users who may have various
versions of libraries.

**Test Suites:** Putting a series of unit tests into a collection of
modules creates a test suite. Typically the suite as a whole is
executed (rather than each test individually) when verifying that the
code base still functions after changes have been made.

# Elements of a Test

**Behavior:** The behavior you want to test. For example, you might want
to test the fun() function.

**Expected Result:** This might be a single number, a range of numbers,
a new fully defined object, a system state, an exception, etc. When we
run the fun() function, we expect to generate some fun. If we don't
generate any fun, the fun() function should fail its test.
Alternatively, if it does create some fun, the fun() function should
pass this test. The expected result should known *a priori*. For
numerical functions, this is result is ideally analytically determined
even if the function being tested isn't.

**Assertions:** Require that some conditional be true. If the
conditional is false, the test fails.

**Fixtures:** Sometimes you have to do some legwork to create the
objects that are necessary to run one or many tests. These objects are
called fixtures as they are not really part of the test themselves but
rather involve getting the computer into the appropriate state.

For example, since fun varies a lot between people, the fun() function
is a method of the Person class. In order to check the fun function,
then, we need to create an appropriate Person object on which to run
fun().

**Setup and teardown:** Creating fixtures is often done in a call to a
setup function. Deleting them and other cleanup is done in a teardown
function.

**The Big Picture:** Putting all this together, the testing algorithm is
often:

```python
setup()
test()
teardown()
```

But, sometimes it's the case that your tests change the fixtures. If so,
it's better for the setup() and teardown() functions to occur on either
side of each test. In that case, the testing algorithm should be:

```python
setup()
test1()
teardown()

setup()
test2()
teardown()

setup()
test3()
teardown()
```
