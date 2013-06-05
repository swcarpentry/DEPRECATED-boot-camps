# Manual set-up

You should ensure you have the following software and tools available. 

* Web browser
* Bash shell (as accessed from a Linux command-line terminal)
* A text editor e.g. nano, vi or emacs
* [Git](http://git-scm.com/)
* Python 2.6 or 2.7
* Python [pip](https://pypi.python.org/pypi/pip)
* Python [nose](https://nose.readthedocs.org/en/latest/)
* Python [coverage](http://nedbatchelder.com/code/coverage/)
* Python [pytest](http://pytest.org/)
* Python [pytest-cov](https://pypi.python.org/pypi/pytest-cov)
* Python [IPython](https://pypi.python.org/pypi/ipython)

Below are brief instructions for RedHat/Scientific Linux 6 and Ubuntu, and how
to check that the install worked. 

For more detailed installation instructions and those for Windows and Mac OSX, check out the [general boot camp set-up instructions](http://software-carpentry.org/setup/) and the installation instructions specific to each product.

If you have any problems then please e-mail the boot camp mailing list at [soton@lists.software-carpentry.org](mailto:soton@lists.software-carpentry.org).

## To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo yum install nano git python ipython python-nose python-coverage pytest python-pip
    $ sudo pip-python install pytest-cov

## To install under Ubuntu

Ubuntu 10.04 and above already comes with a `bash` shell (via Terminal), `vi`
and `nano` text editors, and Python (version 2.6 if running Ubuntu 10.04,
version 2.7 otherwise). To install the other packages run,

    $ sudo apt-get install git python-nose python-pip ipython
    $ sudo pip install pytest-cov

The second `pip install` command also installs `coverage` and `pytest` as
dependencies.

**Troubleshooting**

If the second command exits with errors talking about invalid versions, it may
be because you already (somehow) have `pytest` and/or `coverage` installed from
packages. In that case, do the command below that matches your Ubuntu version
and then rerun the `pip install` command above:

        [Ubuntu 12.10+]
        $ sudo apt-get remove python-coverage python-pytest
        [Ubuntu 12.04+]
        $ sudo apt-get remove python-coverage python-py
        [Ubuntu 10.04+]
        $ sudo apt-get remove python-coverage python-codespeak-lib

**Note on the Recommended Install Procedure**

Python `pytest` and `coverage` *are* available in Ubuntu packages (rather than
installing them via `pip` as per the instructions here), but there are naming
differences per Ubuntu version and, currently, there is a problem with
dependency versions which means it won't work that way for Ubuntu 10.04. If you
prefer doing it that way, see the per release commands below:

        [Ubuntu 12.10+]
        $ sudo apt-get install git python-nose python-pip ipython python-coverage python-pytest
        $ sudo pip install pytest-cov
        [Ubuntu 12.04+]
        $ sudo apt-get install git python-nose python-pip ipython python-coverage python-py
        $ sudo pip install pytest-cov
        [Ubuntu 10.04+, currently causes errors]
        $ sudo apt-get install git python-nose python-pip ipython python-coverage python-codespeak-lib
        $ sudo pip install pytest-cov

## Checking the Installation

Test that they are installed properly as below:

    $ git --version
    git version 1.7.0.4

    $ nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK

    $ nosetests --with-coverage
    [Lots of snipped output]
    Ran 0 tests in 0.009s

    $ py.test
    ===================== test session starts ===================== 
    platform linux2 -- Python 2.7.3 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    ===================== in 0.00 seconds =====================

    $ py.test --cov .
    ===================== test session starts ======================
    platform linux2 -- Python 2.7.3 -- pytest-2.2.4
    collected 0 items 
    Coverage.py warning: No data was collected.
    ------- coverage: platform linux2, python 2.7.3-final-0 --------
    Name    Stmts   Miss  Cover
    ---------------------------
    =======================  in 0.03 seconds =======================

    $ ipython
    Python 2.7.3 (default, Aug  1 2012, 05:14:39) 
    Type "copyright", "credits" or "license" for more information.
    
    IPython 0.12.1 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object', use 'object??' for extra details.
    
    In [1]: quit()

The version numbers shown will obviously depend on your Linux distribution. On
some machines, commands will take a lot longer than the times shown in the
examples above. `py.test` scans the current directory and all nested ones for
Python tests, so it takes longer for large directory structures and *may* find
tests if you happen to have Python code within your current directory. Run the
`py.test` commands from an empty directory if you want to ensure that the
outputs are roughly as above.
