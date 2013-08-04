## Testing in practice

The example we've looked at is based on one function. Suppose we have a complex legacy code of 10000s of lines and which takes many input files and produces many output files. Exactly the same approach can be used as above - we run our code on a set of input files and check whether the output files match what you'd expect. For example, we could,

* Run the code on a set of inputs.
* Save the outputs.
* Refactor the code e.g. to optimise it or parallelise it.
* Run the code on the inputs.
* Check that the outputs match the saved outputs. 

This was the approach taken by EPCC and the Colon Cancer Genetics Group (CCGG) of the MRC Human Genetics Unit at the Western General as part of an [Oncology](http://www.edikt.org/edikt2/OncologyActivity) project to optimise and parallelise a FORTRAN genetics code.

The [Muon Ion Cooling Experiment](http://www.mice.iit.edu/) (MICE) have a large number of tests written in Python. They use [Jenkins](), a *continuous integration server* to build their code and trigger the running of the tests which are then [published online](https://micewww.pp.rl.ac.uk/tab/show/maus).

## When 1 + 1 = 2.0000001

Computers don't do floating point arithmetic too well. This can make simple tests for the equality of two floating point values problematic due to imprecision in the values being compared. 

    $ python
    >>> expected = 1 + 1 
    >>> actual = 2.0000001
    >>> assert expected == actual

We can get round this by comparing to within a given threshold, or delta, for example we may consider *expected* and *actual* to be equal if *expected - actual < 0.000000000001*.

Test frameworks such as `nose`, often provide functions to handle this for us. For example, to test that 2 numbers are equal when rounded to a given number of decimal places,

    $ python
    >>> from nose.tools import assert_almost_equal
    >>> assert_almost_equal(expected, actual, 0)
    >>> assert_almost_equal(expected, actual, 1)
    >>> assert_almost_equal(expected, actual, 3)
    >>> assert_almost_equal(expected, actual, 6)
    >>> assert_almost_equal(expected, actual, 7)
    ...
    AssertionError: 2 != 2.0000000999999998 within 7 places

What do we consider to be a suitable threshold for equality? That is application-specific - for some domains we might be happy to round to the nearest whole number, for others we may want to be far, far more accurate.

## When should we test?

We should test,

* Always!
* Early, and not wait till after we've used it to generate data for our important paper, or given it to someone else to use.
* Often, so that we know that any changes we've made to our code, or to things that our code needs (e.g. libraries, configuration files etc.) haven't introduced any bugs.

But, when should we finish writing tests? How much is enough? 

> **What we know about software development - we can't test everything**

> "It is nearly impossible to test software at the level of 100 percent of its logic paths", fact 32 in R. L. Glass (2002) [Facts and Fallacies of Software Engineering](http://www.amazon.com/Facts-Fallacies-Software-Engineering-Robert/dp/0321117425) ([PDF](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.94.2037&rep=rep1&type=pdf)).

We can't test everything but that's no excuse for testing nothing! How much to test is something to be learned by experience, so think of it as analogous to when you finish proof reading a paper, over and over, before sending it to a conference. If you find bugs when you use your code, you did too little, so consider what you might have done and how to address this next time.

Tests, like code, should ideally be reviewed by a colleague which helps avoid tests that,

* Pass when they should fail, false positives.
* Fail when they should pass, false negatives.
* Don't test anything. 

For example,

    def test_critical_correctness():
        # TODO - will complete this tomorrow!
        pass

Yes, tests like this *do* occur on projects!

Previous: [Let's start writing some tests](Writing.md) Next: [Test-driven development](TDD.md)
