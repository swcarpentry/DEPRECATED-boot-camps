Python Testing Cheat Sheet
==========================

Why testing?
------------

1. Helps you to think about expected behavior, especially boundary cases,
2. documents expected behavior,
3. confidence recent changes didn't break anything that worked before,
4. confidence code is correct.


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

* extending TestCase
* assertions, e.g. self.assertEquals


Nose
----

To run tests, at the shell prompt, type

    nosetests

TODO: finish this (what do you want to add here?)

By default, Nose will find tests in files named `test_*`.


### A simple test

    from nose.tools import assert_equal

    from mystatscode import mean

    def test_single_value():
        observed = mean([1.0])
        expected = 1.0
        assert_equal(observed, expected)

### Other assertions

TODO: finish this
* assertTrue, assertFalse
* assertIn, assertNotIn
* assertIs, assertIsNot
* assertRaises
* (what else?)


### Floating point tests

* assert_almost_equal...
* assertGreater, assertLess

### Fixtures

(todo)


Test-driven deveopment
----------------------

***Red.*** Write test function that checks one new functionality you want to add to your code. -- tests have to fail.

***Green.*** Write minimal code that implements desired features until all tests pass.

***Refactor.*** Improve code wrt. readability and speed. Constantly check that tests still pass.

***Commit.*** Commit working code to version control.

Repeat.


General advice
--------------

* Perfect test-case coverage is impossible.
* Try to test distinct functionalities.
* If you find a bug yet undiscovered by previous test, make it a new test case.


TODO:

* setup...
* per-test fixtures with @with_setup decorator

