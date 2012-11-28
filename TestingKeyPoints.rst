
Testing and Automation Key Points
=================================

Mike Jackson, The Software Sustainability Institute.

This work is licensed under the Creative Commons Attribution License. Copyright (c) Software Carpentry and The University of Edinburgh 2010-2012. See http://software-carpentry.org/license.html for more information.

.. Written in reStructuredText, http://docutils.sourceforge.net/rst.html.

Prerequisites
-------------

Python, nose.

Introduction
------------

Use wget or curl e.g.
::

 wget --no-check-certificate https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/test-materials.zip
 curl --insecure https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/test-materials.zip > test-materials.zip

to get test-materials.zip from BitBucket and unzip these.

utilities.py contains functions, test functions and calls to test functions.

Run the tests.
::

 python utilities.py

Modularity is a design principle so separate out functions and test functions. Create test_utilities.py and move test functions and calls to test functions to this file.

Add to test_utilities.py
::

 from utilities import sum_list
 from utilities import calc_mean

Run tests.
::

 python test_utilities.py

Test data is outside functions but validation is within. Put test data into test functions e.g.
::

 numbers = [1, 2, 3, 4]
 def test_calc_mean(numbers):
     assert calc_mean(numbers) == 2.5, "mean of %s is not 2.5" % str(numbers)
 test_calc_mean(numbers)

becomes
::

 def test_calc_mean():
     numbers = [1, 2, 3, 4]
     assert calc_mean(numbers) == 2.5, "mean of %s is not 2.5" % str(numbers)
 test_calc_mean()

Run.
::

 python test_utilities.py

Test what happens if the list is empty.
::

 def test_sum_empty_list():
     numbers = []
     assert calc_mean(numbers) == 0, "sum of %s is not 0" % str(numbers)
 test_sum_empty_list()

Every time we add a test we need to add a call to that test function. Automate!
::

 nosetests test_utilities.py

nosetests looks for functions beginning with "test\_". 

Remove the test function calls.
::

 nosetests test_utilities.py

nosetests reports progress - "."s - time taken and number of tests run.

nosetests reports failures e.g. introduce a bug into calc_mean so it adds 1 to the total.
::

 nosetests test_utilities.py

Correct calc_mean and retest just that function.
::

 nosetests test_utilities.py:test_calc_mean

nosetests can output an "xUnit" test report.
::

 nosetests --with-xunit test_utilities.py
 cat nosetests.xml

xUnit framework by Kent Beck. JUnit, CUnit, fUnit etc. Can present results in different ways.

Automated build-and-test
------------------------

Version control + automated testing => automated build and test.

EPCC oncology project optimized and paralleled medical code. Initial run to get expected results. Create overnight test job to check out code, run code, compare to expected results. Optimize and parallelize in confidence.

VTK test dashboard, built using CDash. http://open.cdash.org/index.php?project=VTK 

Continuous integration tools detect version control commits, check out code, build, run tests, and publish, or run every few minutes and publish.

MICE test dashboard uses Jenkins continuous integration server. Python code and tests, run using nosetests. https://micewww.pp.rl.ac.uk/tab/show/maus. 

Faster you see a failure, faster you can fix it. Public shame is a motivator too!

OGSA-DAI uses Jenkins, Java code and JUnit tests, http://ogsadai-public.epcc.ed.ac.uk:8080/jenkins/

How much testing is enough?
---------------------------

Learn by experience. Analogous to when to finish a proof reading a paper.

If you find bugs when you use your code, you did too little.

Tests, like code, should be reviewed. 

Helps avoid tests that:
 - Pass when they should fail.
 - Fail when they should pass.
 - Don't test anything. For example,

::

 def test_vital_correctness():
     # TODO - will complete this tomorrow!
     # Tomorrow never comes!
     pass

Test driven development
-----------------------

Common to write code then write tests. 

Test-driven development - test first code second. 

Red-green-refactor:
 - Red - write tests based on requirements. They fail as there is no code!
 - Green - write/modify code to get tests to pass.
 - Refactor code - clean it up.

Think about what the code should do, before we write it, not what we know it does.

Conclusion
----------

Cover Testing.ppt.
