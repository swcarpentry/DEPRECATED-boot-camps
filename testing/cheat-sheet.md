Testing Cheat Sheet
===================

Terminology
-----------

* A *unit test* acts on an isolated component within an application.
* A *test fixture* is the input data that a test acts on.
* The *interface* of a function is its public face, defined by its
  input and output.
* The *implementation* of a function is how it gets from the input to
  the output.
* A *stub* is a very simple implementation of one function that is
  used in a test of a different function.

Unit Testing
------------

* A *normal case* is a test case that reflects what is expected to be
  typical usage of a function.
* A *boundary case* is a test case that reflects a less typical but
  potentially troublesome type of usage.

Exceptions
----------

* raise
* catch
* define

Assertions
----------

* syntax
* stops execution

Unittest
--------

* extending TestCase
* assertions, e.g. self.assertEquals

Nose
----

* invocation: nosetests
* naming conventions: test_*
* fixtures: setup
* per-test fixtures with @with_setup decorator
* assertions, e.g. assert_equal, assert_almost_equal...
