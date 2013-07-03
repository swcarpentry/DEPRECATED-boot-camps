# Testing - crib sheet

## Detecting errors

What we know about software development - code reviews work. Fagan (1976) discovered that a rigorous inspection can remove 60-90% of errors before the first test is run. M.E., Fagan (1976). [Design and Code inspections to reduce errors in program development](http://www.mfagan.com/pdfs/ibmfagan.pdf). IBM Systems Journal 15 (3): pp. 182-211.

What we know about software development - code reviews should be about 60 minutes long. Cohen (2006) discovered that all the value of a code review comes within the first hour, after which reviewers can become exhausted and the issues they find become ever more trivial. J. Cohen (2006). [Best Kept Secrets of Peer Code Review](http://smartbear.com/SmartBear/media/pdfs/best-kept-secrets-of-peer-code-review.pdf). SmartBear, 2006. ISBN-10: 1599160676. ISBN-13: 978-1599160672.

## Runtime tests

[morse04.py](python/morse/morse04.py)

    $ python morse.py
    encode
    1 + 2 = 3

`KeyError` is an exception.

Traceback shows us Python's exception stack trace.

Runtime tests can make code robust and behave gracefully.

    try:
        print "Encoded is '%s'" % translator.encode(message)
    except KeyError:
        print 'The input should be a string of a-z, A-Z, 0-9 or space'

Exception is caught by the `except` block.

Exception can be converted and passed e.g. if this was deep within a function we would not want to print but to keep UI separate so, can `raise` an exception e.g.

    except KeyError:
        raise ValueError('The input should be a string of a-z, A-Z, 0-9 or space')

## Exercise

Add a runtime test for decoding.

## Correctness tests

Testing manually works but is time-consuming and error prone - might forget to run a test.
Write down set of test steps so won't forget. 
Still time-consuming.

    def test(self):
        print "SOS is ", self.encode('SOS')
        print "...---... is ", self.decode('... --- ...')

Extend UI to invoke:

    while True:

        elif line == "test":
            print "Testing..."
            translator.test()
            break

Automate checking.

    def test(self):
        assert '... --- ...' == self.encode('SOS')
        assert 'sos' == self.decode('... --- ...')
        print "OK"

`assert` checks whether condition is true and, if not, raises an exception.

Put test functions in separate file for modularity.

    $ cp morse.py test_morse.py
    $ nano test_morse.py

    from morse import MorseTranslator

    class TestMorseTranslator:

        def test(self):
            translator = MorseTranslator()
            assert '... --- ...' == translator.encode('SOS')
            assert 'sos' == translator.decode('... --- ...')

    if __name__ == "__main__":    

        test_translator = TestMorseTranslator()
        test_translator.test()
        print "OK"

Remove test code from `MorseTranslator`.

Run tests.

    $ python test_morse.py

Modularise functions.

    def test_encode_sos(self):
        ...
    def test_decode_sos(self):
        ...

    test_translator.test_encode_sos()
    test_translator.test_decode_sos()

Remove duplicated code:

    def __init__(self):
        self.translator = MorseTranslator()

Test function,

* Set up inputs and expected outputs.
* Runs function / component on inputs to get actual outputs.
* Checks actual outputs match expected outputs. 

Verbose, but equivalent, version of `test_encode_sos`.

    def test_encode_sos(self):
        expected = '... --- ...'
        actual = self.translator.encode('SOS')                     
        assert expected == actual

## `nose` - a Python test framework

[nose](https://pypi.python.org/pypi/nose/) automatically finds, runs and reports on tests.

[xUnit test framework](http://en.wikipedia.org/wiki/XUnit).

`test_` file and function prefix.

    $ nosetests test_morse.py

`.` denotes successful tests.

Remove `__main__ code.

    $ nosetests test_morse.py

xUnit test report, standard format, convert to HTML, present online.

    $ nosetests --with-xunit test_dna.py
    $ cat nosetests.xml

## Propose some more tests. 

Consider,

* What haven't we tested for so far? 
* Have we covered all possible strings?
* Have we covered all possible arguments?

Propose examples and add to Etherpad.

    encode('sos')
    encode('')
    decode('')
    encode('1 + 2 = 3')
    decode('...---...')

Implement examples.

Tests for illegal arguments:

    def test_encode_illegal(self):
        try:
            self.translator.encode('1 + 2 = 3')
            assert False
        except KeyError:
            assert True

Alternatively:

    from nose.tools import assert_raises

    def test_encode_illegal(self):
        assert_raises(KeyError, self.translator.encode, '1 + 2 = 3')

Testing components together:

    assert 'sos' == decode(encode('sos'))
    assert '... --- ...' == encode(decode('... --- ...'))

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

## Summary

Testing

* Saves time.
* Gives confidence that code does what we want and expect it to.
* Promotes trust that code, and so research, is correct.

Remember [Geoffrey Chang](http://en.wikipedia.org/wiki/Geoffrey_Chang)

Bruce Eckel, [Thinking in Java, 3rd Edition](http://www.mindview.net/Books/TIJ/),  "If it's not tested, it's broken".

## Links

* [Software Carpentry](http://software-carpentry.org/)'s online [testing](http://software-carpentry.org/4_0/test/index.html) lectures.
* A discussion on [is it worthwhile to write unit tests for scientific research codes?](http://scicomp.stackexchange.com/questions/206/is-it-worthwhile-to-write-unit-tests-for-scientific-research-codes)
* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].
