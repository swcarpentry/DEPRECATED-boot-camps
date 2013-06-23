# Testing


**Based on materials by Mike Jackson, Katy Huff, Paul Ivanov, Rachel Slaybaugh, and Anthony Scopatz.**

## Some questions for you

Let's open with some questions. How many of you test your code or scripts...

* With a set of input data and check that this gives an expected set of output data?
* With multiple sets of input data?
* With input data you know to be incorrect and check that your code or scripts behave appropriately e.g. giving a warning, or exiting, or anything that doesn't involve producing output data that *seems* OK but is not?
* After every change you've made to fix a bug or optimise your code or to add a new feature?
* Using some form of automation e.g. a set of testing scripts or a test framework?

## What is testing?

Software testing is exercising and confirming expected behaviours and results from code. It allows us to check that,

* Our code behaves as expected and produces valid output data given valid input data.
* Our code does this using *any* set of valid input data.
* Our code fails gracefully if given invalid input data - it does not just crash or behave mysteriously or unpredictably but, for example, exits politely with a message as to why the input data is invalid.
* Our code can handle extreme boundaries of input domains, output ranges, parametric combinations or any other edge cases.
* Our code's existing behaviour is still the same after we've changed it (this is called *regression testing*).

In testing, you'll come across the terms *verification* and *validation*,

* Verification is the process of asking, "Have we built the software correctly?" That is, is the code bug free, precise, accurate, and repeatable?
* Validation is the process of asking, "Have we built the right software?" That is, is the code designed in such a way as to produce the answers we are interested in, data we want, etc.

Testing also gives us the confidence to...

* Add new features.
* Optimise our code.
* Parallelise our code.
* Fix bugs.

...all without introducing bugs. Nothing is worse than fixing a bug only to introduce a new one.

Tests also help us remember what all the parts of our code does. If we are working on a large project over three years and end up with 100s of functions, it may be hard to remember what each function does in detail. If we have a test that checks all of the function's functionality, we can look at the test to remember what it's supposed to do.

## Why isn't testing done?

Do any of these sound familiar?

* "I don't write buggy code", well, we are naturally very protective of our work and may refuse to accept that our code has bugs. Unfortunately, almost all code has bugs.
* "It's too hard", but, if it's hard to write a test for some code, then this is a good sign that the code is not well designed.
* "It's not interesting", sometimes, testing is just viewed as not being interesting which means...
* "It takes too much time and I've research to do"

And, this is a fair point. So, why should we care about testing?

## Why we should do testing?

Testing allows us, and others, to trust our code and trust it enough to answer in the affirmative to at least a few of the following questions:

* Does your code work?
* Always?
* Does it do what we think it does?
* Does it continue to work after changes are made, for example optimisations or bug fixes?
* Does it continue to work after system configurations or libraries are upgraded?
* Does it respond properly for a full range of input parameters?
* Does it handle edge and corner cases?

As a cautionary tale, consider Ariane 5 which used Ariane 4 software. Ariane 5 had new and improved engines which caused the code to produce a buffer overflow...and Ariane 5 blew up! So, some forgotten tests led to millions of pounds down the drain and some very red faces.

Or, consider [Geoffrey Chang](http://en.wikipedia.org/wiki/Geoffrey_Chang) who had to [retract](http://www.sciencemag.org/content/314/5807/1875.2.long) 3 papers from [Science](http://www.sciencemag.org), due to a flipped sign! Or, McKitrick and Michaels' [Climate Research 26(2) 2004](http://www.int-res.com/abstracts/cr/v26/n2/p159-173/) paper, which drew the attention of a blogger Tim Lambert who noted a [problem](http://crookedtimber.org/2004/08/25/mckitrick-mucks-it-up/) which led to their subsequent [erratum](http://www.int-res.com/articles/cr2004/27/c027p265.pdf).

Do this too regularly and people may not trust our research, which could affect our chances for collaborations, publications or funding.

But if this is not compelling, then, if nothing else, writing tests is an investment in time that saves us time in future,

* We can automate the checking of outputs from our software to ensure they're valid.
* We can detect more quickly whether refactoring, optimisation or parallelisation has introduced bugs.
* We can run our tests while doing other, more interesting, things.

## Fixing things before we test...

Before we test our code, it can be very productive to get a colleague to look at it for us...why?

> **What we know about software development - code reviews work** 

> Fagan (1976) discovered that a rigorous inspection can remove 60-90% of errors before the first test is run. 
> M.E., Fagan (1976). [Design and Code inspections to reduce errors in program development](http://www.mfagan.com/pdfs/ibmfagan.pdf). IBM Systems Journal 15 (3): pp. 182-211.

> **What we know about software development - code reviews should be about 60 minutes long** 

> Cohen (2006) discovered that all the value of a code review comes within the first hour, after which reviewers can become exhausted and the issues they find become ever more trivial.
> J. Cohen (2006). [Best Kept Secrets of Peer Code Review](http://smartbear.com/SmartBear/media/pdfs/best-kept-secrets-of-peer-code-review.pdf). SmartBear, 2006. ISBN-10: 1599160676. ISBN-13: 978-1599160672.

## Let's dive in...

* [Let's start writing some tests](Writing.md)
* [Testing in practice](RealWorld.md)
* [Test-driven development](TDD.md)
* [Conclusions and further information](Conclusion.md)

Next: [Let's start writing some tests](Writing.md)
