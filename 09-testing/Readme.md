# Testing

* * * * *

**Based on materials by Katy Huff, Rachel Slaybaugh, and Anthony
Scopatz**

![image](https://github.com/thehackerwithin/UofCSCBC2012/raw/scopz/5-Testing/test_prod.jpg)
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
-   What about **edge or corner cases**?
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

Say we have an averaging function:

```python
def mean(numlist):
    total = sum(numlist)
    length = len(numlist)
    return total/length
```

Tests could be implemented as runtime **exceptions in the function**:

```python
def mean(numlist):
    try:
        total = sum(numlist)
        length = len(numlist)
    except TypeError:
        raise TypeError("The number list was not a list of numbers.")
    except:
        print "There was a problem evaluating the number list."
    return total/length
```

Sometimes tests they are functions alongside the function definitions
they are testing.

```python
def mean(numlist):
    try:
        total = sum(numlist)
        length = len(numlist)
    except TypeError:
        raise TypeError("The number list was not a list of numbers.")
    except:
        print "There was a problem evaluating the number list."
    return total/length


def test_mean():
    assert mean([0, 0, 0, 0]) == 0
    assert mean([0, 200]) == 100
    assert mean([0, -200]) == -100
    assert mean([0]) == 0


def test_floating_mean():
    assert mean([1, 2]) == 1.5
```

Sometimes they are in an executable independent of the main executable.

```python
def mean(numlist):
    try:
        total = sum(numlist)
        length = len(numlist)
    except TypeError:
        raise TypeError("The number list was not a list of numbers.")
    except:
        print "There was a problem evaluating the number list."
    return total/length
```

Where, in a different file exists a test module:

```python
import mean

def test_mean():
    assert mean([0, 0, 0, 0]) == 0
    assert mean([0, 200]) == 100
    assert mean([0, -200]) == -100
    assert mean([0]) == 0


def test_floating_mean():
    assert mean([1, 2]) == 1.5
```

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
wreak anything. The same idea applies to making changes in your system
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

**Regression Tests:** A regression test ensures that new code does not
change anything unexpected.

**Integration Tests:** Integration tests query the ability of the code
to integrate well with the system configuration and third party
libraries and modules. This type of test is essential for codes that
depend on libraries which might be updated independently of your code or
when your code might be used by a number of users who may have various
versions of libraries.

**Test Suites:** Putting a series of unit tests into a collection of
modules creates, a test suite. Typically the suite as a whole is
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
pass this test. The the expected result should known *a priori*. For
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

* * * * *

# Nose: A Python Testing Framework

The testing framework we'll discuss today is called
[nose](https://nose.readthedocs.org/en/latest/). However, there are
several other testing frameworks available in most languages. Most
notably there is [JUnit](http://www.junit.org/) in Java which can
arguably attributed to inventing the testing framework.
## Where do nose tests live?

Nose tests are files that begin with `Test-`, `Test_`, `test-`, or
`test_`. Specifically, these satisfy the testMatch regular expression
`[Tt]est[-_]`. (You can also teach nose to find tests by declaring them
in the unittest.TestCase subclasses that you create in your code. You
can also create test functions which are not unittest.TestCase
subclasses if they are named with the configured testMatch regular
expression.)

## Nose Test Syntax

To write a nose test, we make assertions.

```python
assert should_be_true()
assert not should_not_be_true()
```

Additionally, nose itself defines number of assert functions which can
be used to test more specific aspects of the code base.

```python
from nose.tools import *

assert_equal(a, b)
assert_almost_equal(a, b)
assert_true(a)
assert_false(a)
assert_raises(exception, func, *args, **kwargs)
assert_is_instance(a, b)
# and many more!
```

Moreover, numpy offers similar testing functions for arrays:

```python
from numpy.testing import *

assert_array_equal(a, b)
assert_array_almost_equal(a, b)
# etc.
```

## Exercise: Writing tests for mean()

There are a few tests for the mean() function that we listed in this
lesson. What are some tests that should fail? Add at least three test
cases to this set. Edit the `test_mean.py` file which tests the mean()
function in `mean.py`.

*Hint:* Think about what form your input could take and what you should
do to handle it. Also, think about the type of the elements in the list.
What should be done if you pass a list of integers? What if you pass a
list of strings?

**Example**:

    nosetests -v test_mean.py

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

To illustrate TDD, let's return to the function you wrote yesterday,
`calculate_gc`. We'll start from scratch and develop the function
by meeting test specifications. 

The beginning of the function is contained in the file `calculate_gc.py`
in this directory. It currently takes one argument as input, but does
nothing.

```python
def calculate_gc(x):
    '''
    Calculates the GC content of DNA sequence x.
    '''
    pass
```

The tests that we must pass are contained in the file
`test_calculate_gc.py`. We can run the tests using nosetests.

    nosetests -v test_calculate_gc.py

As expected, we fail all the tests! What is the bare minimum 
functionality we must add to pass the first test below?

```python
def test_only_G_and_C():
    '''
    Sequence of only G's and C's has fraction 1.0
    '''
    fixture = 'GGCGCCGGC'
    result = calculate_gc(fixture)
    assert_equal(result, 1.0)
```

And the second test?

```python
def test_half():
    '''
    Sequence with half G and C has fraction 0.5
    '''
    fixture = 'ATGC'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)
```

Test number three?

```python
def test_lower_case():
    '''
    Sequence with lower case letters
    '''
    fixture = 'atgc'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)
```

Test number four?

```python
def test_not_DNA():
    '''
    Raise TypeError if not DNA
    '''
    fixture = 'qwerty'
    assert_raises(TypeError, calculate_gc, fixture)
```

Through this cycle of writing tests and modifying the function to pass 
the tests, we have developed a function that behaves exactly as we 
expect and nothing more. And the tests not only serve as documentation 
of what the function does, but can also be easily ran again if we made 
further modifications (regression tests). What would be the next test 
you would write for our function?

