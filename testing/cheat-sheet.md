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

TODO: finish this

### Floating point tests

* assert_almost_equal...

### Fixtures

TODO:

* setup...
* per-test fixtures with @with_setup decorator
