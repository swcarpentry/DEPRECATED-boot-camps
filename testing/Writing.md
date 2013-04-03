## Let's start writing some tests

In the file [dna.py](python/dna/dna.py) we have a Python dictionary that stores the molecular weights of the 4 standard DNA nucleotides, A, T, C and G, 

    NUCLEOTIDES = {'A':131.2, 'T':304.2, 'C':289.2, 'G':329.2}

and a Python function that takes a DNA sequence as input and returns its molecular weight, which is the sum of the weights for each nucelotide in the sequence,
 
    def calculate_weight(sequence):
        """
        Calculate the molecular weight of a DNA sequence.
        @param sequence: DNA sequence expressed as an upper-case string. 
        @return molecular weight.
        """
        weight = 0.0
        for ch in sequence:
            weight += NUCLEOTIDES[ch]
        return weight

We can calculate the molecular weight of a sequence by,
 
    weight = calculate_weight('GATGCTGTGGATAA')
    print weight

We can add a test to our code as follows,

    def calculate_weight(sequence):
        """
        Calculate the molecular weight of a DNA sequence.

        @param sequence: DNA sequence expressed as an upper-case string.
        @return molecular weight.
        """
        weight = 0.0
        try:
            for ch in sequence:
                weight += NUCLEOTIDES[ch]
            return weight
        except TypeError:
            print 'The input is not a sequence e.g. a string or list'

If the input is not a string, or a list of characters then the `for...in` statement will *raise an exception* which is *caught* by the `except` block. For example,

    print calculate_weight(123)

This is a *runtime test*. It alerts the user to exceptional behavior in the code. Often, exceptions are related to functions that depend on input that is unknown at compile time. Such tests make our code robust and allows our code to behave gracefully - they anticipate problematic values and handle them.

Often, we want to pass such errors to other points in our program rather than just print a message and continue. So, for example we could do,

    except TypeError:
        raise ValueError('The input is not a sequence e.g. a string or list')

which raises a new exception, with a more meaningful message. If writing a complex application, our user interface could then present this to the user e.g. as a dialog box.

Runtime tests don't test our functions behaviour or whether it's implemented correctly. So, we can add some tests,

    print calculate_weight('A')
    print calculate_weight('G')
    print calculate_weight('GA')

But we'd have to visually inspect the results to see they are as expected. So, let's have the computer do that for us and make our lives easier, and save us time in checking,

    assert calculate_weight('A') == 131.2
    assert calculate_weight('G') == 329.2
    assert calculate_weight('GA') == 460.4

`assert` checks whether a condition is true and, if not, raises an exception.

We explicitly list the expected weights in each statement. But, by doing this there is a risk that we mistype one. A good design principle is to define constant values in one place only. As we already have defined them in `nucleotides` we can just refer to that,

    assert calculate_weight('A') == NUCLEOTIDES['A']
    assert calculate_weight('G') == NUCLEOTIDES['G']
    assert calculate_weight('GA') == NUCLEOTIDES['G'] + NUCLEOTIDES['A']

But this isn't very modular, and modularity is a good design principle, so let's define some test functions,

    def test_a():
        assert calculate_weight('A') == NUCLEOTIDES['A']
    def test_g():
        assert calculate_weight('G') == NUCLEOTIDES['G']
    def test_ga():
        assert calculate_weight('GA') == NUCLEOTIDES['GA'] + NUCLEOTIDES['A']

    test_a()
    test_g()
    test_ga()

And, rather than have our tests and code in the same file, let's separate them out. So, let's create

    $ nano test_dna.py

Now, our function and nucleotides data are in `dna.py` and we want to refer to them in `test_dna.py` file, we need to *import* them. We can do this as,

    from dna import calculate_weight
    from dna import NUCLEOTIDES

Then we can add all our test functions and function calls to this file. And run the tests,

    $ python test_dna.py

## `nose` - a Python test framework

`nose` is a test framework for Python that will automatically find, run and report on tests written in Python. It is an example of what has been termed an *[xUnit test framework](http://en.wikipedia.org/wiki/XUnit)*, perhaps the most famous being JUnit for Java.

To use `nose`, we write test functions, as we've been doing, with the prefix `test_` and put these in files, likewise prefixed by `test_`. The prefixes `Test-`, `Test_` and `test-` can also be used.

Typically, a test function,

* Sets up some inputs and the associated expected outputs. The expected outputs might be a single number, a range of numbers, some text, a file, a set of files, or whatever.
* Runs the function or component being tested on the inputs to get some actual outputs.
* Checks that the actual outputs match the expected outputs. We use assertions as part of this checking. We can check both that conditions hold and that conditions do not hold.

So, we could rewrite `test_a`, as the more, verbose, but equivalent,

    def test_a():
        expected = NUCLEOTIDES['A']
        actual = calculate_weight('A')                     
        assert expected == actual

Python `assert` allows us to check,

    assert should_be_true()
    assert not should_not_be_true()

`nose` defines additional functions which can be used to check for a rich range of conditions e.g..

    from nose.tools import *

    assert_equal(a, b)
    assert_almost_equal(a, b, 3)
    assert_true(a)
    assert_false(a)
    assert_raises(exception, func, *args, **kwargs)
    ...

`assert_raises` is used for where we want to test that an exception is raised if, for example, we give a function a bad input.

To run `nose` for our tests, we can do,

    $ nosetests test_dna.py

Each `.` corresponds to a successful test. And to prove `nose` is finding our tests, let's remove the function calls from `test_dna.py` and try again,

    $ nosetests test_dna.py

nosetests can output an "xUnit" test report,

    $ nosetests --with-xunit test_dna.py
    $ cat nosetests.xml

This is a standard format that that is supported by a number of xUnit frameworks which can then be converted to HTML and presented online. 

## Write some more tests

Let's spend a few minutes coming up with some more tests for `calculate_weight`. Consider,

* What haven't we tested for so far? 
* Have we covered all the nucleotides? 
* Have we covered all the types of string we can expect? 
* In addition to test functions, other types of runtime test could we add to `calculate_weight`?

Examples of tests we could add include,

* `calculate_weight('T')`
* `calculate_weight('C')`
* `calculate_weight('TC')`
* `calculate_weight(123)` 

The latter requires us to check whether an exception was raised which we can do as follows:

    try:
        calculate_weight(123) 
        assert False
    except ValueError:
        assert True

This is like catching a runtime error. If an exception is raised then our test passes (`assert True`), else if no exception is raised, it fails.

One other test we could do is `calculate_weight('GATCX')` for which we can add another runtime test,

        ...
    except KeyError:
        raise ValueError('The input is not a sequence of G,T,C,A')

Previous: [Testing](README.md) Next: [Testing in practice](RealWorld.md)
