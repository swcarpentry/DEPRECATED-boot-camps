## Test Driven Development

Traditionally, we'd write our code, then write the tests. [Test driven development](http://www.amazon.com/Test-Driven-Development-By-Example/dp/0321146530) (TDD), proposed by Kent Beck, is a philosophy that turns this on its head - we write code by *writing the tests first*, then write the code to make the tests pass. If a new feature is needed, another test is written and the code is expanded to meet this new use case. This continues until the code does what is needed. This can be summarised as red-green-refactor:

 * Red - write tests based on requirements. They fail as there is no code!
 * Green - write/modify code to get tests to pass.
 * Refactor code - clean it up.

By writing tests first, we're forced to think about what our code should do. In contrast, in writing our code then tests, we risk testing what the code actually does, rather than what it should do.

TDD operates on the YAGNI principle (You Ain't Gonna Need It) to avoid developing code for which there is no need.

## TDD of a DNA complement function

Given a DNA sequence consisting of A, C, T and G, we can create its complementary DNA, cDNA, by applying a mapping to each nucleotide in turn,

* A => T
* C => G
* T => A
* G => C

For example, given DNA strand GTCA, the cDNA is CAGT. 

So, let's write a `complement` function that creates the cDNA strand, given a DNA strand in a string. We'll use TDD, so to start, let's create a file `test_cdna.py` and add a test,

    from cdna import complement

    def test_complement_a():
        assert_equals complement('A') == 'T'

And let's run the test,

    $ nosetests test_cdna.py

Which fails as we have no function! So, let's create a file `cdna.py`. Our initial function to get the tests to pass could be,

    def complement(sequence):
        return 'T'

This is simplistic, but the test passes. Now let's add another test,

    def test_complement_c():
        assert complement('C') == 'G'

To get both our tests to pass, we can change our function to be,

    def complement(sequence):
        if (sequence == 'A'):
            return 'T'
        else:
            return 'G'

Now, add some more tests. Don't worry about `complement` just now.

Let's discuss the tests you've come up with.

Now update `complement` to make your tests pass. You may want to reuse some of the logic of `calculate_weight`!

When we're done, not only do we have a working function, we also have a set of tests. There's no risk of us leaving the tests "till later" and then never having time to write them.

Previous: [Testing in practice](RealWorld.md) Next: [Conclusions and further information](Conclusion.md)
