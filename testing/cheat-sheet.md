Python Testing Cheat Sheet
==========================

Defensive programming
---------------------

Using an assertion to ensure input is acceptable:

    def some_function(x):
        assert x >= 0
        # ... continue safe in knowledge that x > 0

Adding an explanatory message to the assertion:

        assert x >= 0, "Function not defined for negative x."

Alternatively, raise an exception to indicate what the problem is:

    def some_function(x):
        if x < 0:
            raise TypeError, "Function not defined for negative x."
        return 0


Unittest
--------

TODO: write this section

* extending TestCase
* assertions, e.g. self.assertEquals


Nose
----

To run tests, at the shell prompt, type

    nosetests

By default, Nose will

* look for test functions that have a name starting with `test`
* look for them in files with names starting with `test`
* look for such files in the current working directory, and in subdirectories with names starting with `test`

There are some additional rules, and you can configure your own, but this should be enough to get started.

### A simple test

    from nose.tools import assert_equal

    from mystatscode import mean

    def test_single_value():
        observed = mean([1.0])
        expected = 1.0
        assert_equal(observed, expected)

### Other assertions

Nose provides a range of assertions that can be used when a test is not just checking a simple equality, e.g.

    from nose.tools import assert_items_equal

    from mycode import find_factors

    def test_6():
        observed = find_factors(6)
        expected = [2, 3]
        assert_items_equal(observed, expected) # order of factors is not guaranteed

### Floating point tests

When comparing floating-point numbers for equality, allow some tolerance for small differences due to
the way values are represented and rounded.

    from nose.tools import assert_almost_equal

    from mycode import hypotenuse

    def test_hypotenuse_345():
        observed = hypotenuse(3.0, 4.0)
        expected = 5.0
        assert_almost_equal(observed, expected)

### Testing exceptions

TODO:

* using `@raises` decorator

### Fixtures

TODO:

* setup...
* per-test fixtures with @with_setup decorator
