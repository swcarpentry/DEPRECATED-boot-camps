# Testing - cheat sheet

Thanks to Gordon Webster, the [Digital Biologist](http://www.digitalbiologist.com), for allowing use of his [Python DNA function](http://www.digitalbiologist.com/2011/04/code-tutorial-getting-started-with-python.html).

## Detecting errors

What we know about software development - code reviews work. Fagan (1976) discovered that a rigorous inspection can remove 60-90% of errors before the first test is run. 

What we know about software development - code reviews should be about 60 minutes long. Cohen (2006) discovered that all the value of a code review comes within the first hour, after which reviewers can become exhausted and the issues they find become ever more trivial.

## Runtime tests

[dna.py](python/dna/dna.py)

Dictionary stores molecular weights of 4 standard DNA nucleotides, A, T, C and G

Function takes DNA sequence as input and returns its molecular weight, which is the sum of the weights for each nucelotide in the sequence,

    $ nano dna.py
    weight = calculate_weight('GATGCTGTGGATAA')
    print weight
    print calculate_weight(123)

`TypeError` is an exception, raised, here, by `for...in`.

Runtime tests can make code robust and behave gracefully.

        try:
            for ch in sequence:
                weight += NUCLEOTIDES[ch]
            return weight
        except TypeError:
            print 'The input is not a sequence e.g. a string or list'

Exception is caught by the `except` block.

Pass error to another part of program e.g. UI

    except TypeError:
        raise ValueError('The input is not a sequence e.g. a string or list')

## Correctness tests

Test implementation is correct.

    print "A is ", calculate_weight('A')
    print "G is ", calculate_weight('G')
    print "GA is ", calculate_weight('GA')

Automate checking.

    assert calculate_weight('A') == 131.2
    assert calculate_weight('G') == 329.2
    assert calculate_weight('GA') == 460.4

`assert` checks whether condition is true and, if not, raises an exception.

Define constants in one place only.

    assert calculate_weight('A') == NUCLEOTIDES['A']
    assert calculate_weight('G') == NUCLEOTIDES['G']
    assert calculate_weight('GA') == NUCLEOTIDES['G'] + NUCLEOTIDES['A']

Define test functions for modularity.

    def test_a():
        assert calculate_weight('A') == NUCLEOTIDES['A']
    def test_g():
        assert calculate_weight('G') == NUCLEOTIDES['G']
    def test_ga():
        assert calculate_weight('GA') == NUCLEOTIDES['G'] + NUCLEOTIDES['A']

    test_a()
    test_g()
    test_ga()

Put test functions in separate file for modularity.

    $ nano test_dna.py

Import dictionary and function.

    from dna import calculate_weight
    from dna import NUCLEOTIDES

Run tests.

    $ python test_dna.py

Test function,

* Set up inputs and expected outputs.
* Runs function / component on inputs to get actual outputs.
* Checks actual outputs match expected outputs. 

Verbose, but equivalent, version of `test_a`.

    def test_a():
        expected = NUCLEOTIDES['A']
        actual = calculate_weight('A')                     
        assert expected == actual

## `nose` - a Python test framework

[nose](https://pypi.python.org/pypi/nose/) automatically finds, runs and reports on tests.

[xUnit test framework](http://en.wikipedia.org/wiki/XUnit).

`test_` file and function prefix.

    $ nosetests test_dna.py

`.` denotes successful tests.

    # Remove `test_` function calls.
    $ nosetests test_dna.py

xUnit test report, standard format, convert to HTML, present online.

    $ nosetests --with-xunit test_dna.py
    $ cat nosetests.xml

    $ python
    >>> from nose.tools import *
    >>> expected = 123
    >>> actual = 123
    >>> assert_equal(expected, actual)
    >>> actual = 456
    >>> assert_equal(expected, actual)
    >>> expected = "GATTACCA"
    >>> actual = ["GATC", "GATTACCA"]
    >>> assert_true(expected in actual)
    >>> assert_false(expected in actual)
    >>> assert_true("GTA" in actual, "Expected value was not in the output list")

## Propose some more tests. 

Consider,

* What haven't we tested for so far? 
* Have we covered all the nucleotides? 
* Have we covered all the types of string we can expect?
* In addition to test functions, other types of runtime test could we add to `calculate_weight`?

Examples.

    calculate_weight('T')
    calculate_weight('C')
    calculate_weight('TC')
    calculate_weight(123)

Test for the latter,

    try:
        calculate_weight(123) 
        assert False
    except ValueError:
        assert True

Alternatively,

    from nose.tools import assert_raises

    def test_123():
        assert_raises(ValueError, calculate_weight, 123)

Another run-time test, for `GATCX`,

        ...
    except KeyError:
        raise ValueError('The input is not a sequence of G,T,C,A')

## Testing in practice

Legacy code of 10000s of lines, with many input and output files,

* Run code on set of input files.
* Save output files.
* Refactor code e.g. to optimise it or parallelise it.
* Run code on set of input files.
* Check that outputs match saved outputs. 

EPCC and the Colon Cancer Genetics Group (CCGG) of MRC Human Genetics Unit at the Western General Hospital, Edinburgh - [Oncology](http://www.edikt.org/edikt2/OncologyActivity) project to optimise and parallelise FORTRAN genetics code.

Continuous integration server e.g. [Jenkins](http://jenkins-ci.org/) - detect commit to version control, build, run tests, publish.

[Muon Ion Cooling Experiment](http://www.mice.iit.edu/) (MICE) - Bazaar version control, Python tests, Jenkins, [published online](https://micewww.pp.rl.ac.uk/tab/show/maus).

## When 1 + 1 = 2.0000001

Computers don't do floating point arithmetic too well.

    $ python
    >>> expected = 0
    >>> actual = 0.1 + 0.1 + 0.1 - 0.3
    >>> assert expected == actual
    >>> print actual

Compare to within a threshold, or delta e.g. expected == actual  if expected - actual < 0.0000000000000001.

Thresholds are application-specific. 

Python [decimal](http://docs.python.org/2/library/decimal.html), floating-point arithmetic functions.

    $ python
    >>> from nose.tools import assert_almost_equal
    >>> assert_almost_equal(expected, actual, 0)
    >>> assert_almost_equal(expected, actual, 10)
    >>> assert_almost_equal(expected, actual, 15)
    >>> assert_almost_equal(expected, actual, 16)

`nose.testing` uses absolute tolerance: abs(x, y) <= delta

[Numpy](http://www.numpy.org/)'s `numpy.testing` uses relative tolerance: abs(x, y) <= delta * (max(abs(x), abs(y)). 

`assert_allclose(actual_array, expected_array, relative_tolerance, absolute_tolerance)`

## When should we test?

* Always!
* Early, and not wait till after we've used it to generate data for our important paper, or given it to someone else to use.
* Often, so that we know that any changes we've made to our code, or to things that our code needs (e.g. libraries, configuration files etc.) haven't introduced any bugs.

How much is enough? 

What we know about software development - we can't test everything. It is nearly impossible to test software at the level of 100 percent of its logic paths", fact 32 in R. L. Glass (2002).

No excuse for testing nothing! Learn by experience, like writing a paper.

Review tests, like code, to avoid

* Pass when they should fail, false positives.
* Fail when they should pass, false negatives.
* Don't test anything. 

Example.

    def test_critical_correctness():
        # TODO - will complete this tomorrow!
        pass

## Test-driven development

Complement, cDNA, mapping,

* A => T
* C => G
* T => A
* G => C

DNA strand GTCA, cDNA strand CAGT. 

Antiparallel DNA - calculate inverse of cDNA by reversing it.

[Test driven development](http://www.amazon.com/Test-Driven-Development-By-Example/dp/0321146530) (TDD)

 * Red - write tests based on requirements. They fail as there is no code.
 * Green - write/modify code to get tests to pass.
 * Refactor code - clean it up.

Think about what code should do and test what it should do rather than what it actually does.

YAGNI principle (You Ain't Gonna Need It) avoid developing code for which there is no need.

    $ nano test_dnautils.py
    from dnautils import antiparallel

    $ nosetests test_dnautils.py
    $ nano dnautils.py

    def antiparallel(sequence):
        """
        Calculate the antiparallel of a DNA sequence.
 
        @param sequence: a DNA sequence expressed as an upper-case string.
        @return antiparallel as an upper-case string. 
        """
        pass

    $ nosetests test_dnautils.py

Propose a test and add it.

    $ nosetests test_dnautils.py

Change function to make test pass.

Propose some tests and implement these, but don't update the function.

Discuss!

Update the function.

TDD delivers tests and a function - no risk of us leaving the tests "till later" (never!).

Refactor code with security of tests to flag any changes made introduce a bug.

## Conclusions and further information

Testing

* Saves us time.
* Gives us confidence that our code does what we want and expect it to.
* Promotes trust that our code, and so our research, is correct.

Remember [Geoffrey Chang](http://en.wikipedia.org/wiki/Geoffrey_Chang)

Bruce Eckel, [Thinking in Java, 3rd Edition](http://www.mindview.net/Books/TIJ/),  "If it's not tested, it's broken".

## Find out more...

* [Software Carpentry](http://software-carpentry.org/)'s online [testing](http://software-carpentry.org/4_0/test/index.html) lectures.
* A discussion on [is it worthwhile to write unit tests for scientific research codes?](http://scicomp.stackexchange.com/questions/206/is-it-worthwhile-to-write-unit-tests-for-scientific-research-codes)
* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].
