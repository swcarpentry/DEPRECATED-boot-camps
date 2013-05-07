## Test-driven development

Given a DNA sequence consisting of A, C, T and G, we can create its *complement*, cDNA, which is DNA synthesized from a messenger RNA template by applying a mapping to each nucleotide in turn,

* A => T
* C => G
* T => A
* G => C

For example, given DNA strand GTCA, the cDNA is CAGT. 

We can then create its *antiparallel* by calculating the *inverse* of the sequence, by reversing it. So, the anti-parallel of GTCA is TGAC.

Let's write a function to calculate this antiparallel. Now, before, we had our code and wrote some tests. This time we're going to turn this on it's head and try some test-driven development.

[Test driven development](http://www.amazon.com/Test-Driven-Development-By-Example/dp/0321146530) (TDD), proposed by Kent Beck, is a philosophy and approach where we write code by *writing the tests first*, then write the code to make the tests pass. If a new feature is needed, another test is written and the code is expanded to meet this new use case. This continues until the code does what is needed. This can be summarised as red-green-refactor:

 * Red - write tests based on requirements. They fail as there is no code!
 * Green - write/modify code to get tests to pass.
 * Refactor code - clean it up.

By writing tests first, we're forced to think about what our code should do. In contrast, in writing our code then tests, we risk testing what the code actually *does*, rather than what it *should* do.

TDD operates on the YAGNI principle (You Ain't Gonna Need It) to avoid developing code for which there is no need.

So, back to our example. We'll start by creating a file `test_dnautils.py` and import our function,

    from dnautils import antiparallel

And then run `nosetests`,

    $ nosetests test_dnautils.py

This fails as not only are there no tests, there's no module or function. Let's create a file, `dnautils.py`, and add a function that does nothing,

    def antiparallel(sequence):
        """
        Calculate the antiparallel of a DNA sequence.
 
        @param sequence: a DNA sequence expressed as an upper-case string.
        @return antiparallel as an upper-case string. 
        """
        pass

And let's run the tests to date,

    $ nosetests test_dnautils.py

Zero tests, as we expected! Nnow we need to add some tests. Someone propose a test...

OK we'll add that test...

And let's run the tests to date,

    $ nosetests test_dnautils.py

This fails as our function does nothing. So let's change it to pass...

Now let's add another test...someone give me one...

To get both our tests to pass, we can change our function to be...

Now, add some more tests to `test_dnautils.py` but do not make any more changes to `dnautils.py` or your function, yet.

Let's discuss the tests you've come up with...

Now update `antiparallel` to make your tests pass...

When we're done, not only do we have a working function, we also have a set of tests. There's no risk of us leaving the tests "till later" and then never having time to write them.

We now may want to spend time refactoring our function to clean up our code. We can do this with the security of our tests which allow us to detect if any changes we make introduce a bug.

Previous: [Testing in practice](RealWorld.md) Next: [Conclusions and further information](Conclusion.md)
