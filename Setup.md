# Manual set-up

You should ensure you have the following software and tools available. 

* Web browser
* Bash shell
* A text editor e.g. nano, vi or emacs
* Make
* [Git](http://git-scm.com/)
* Python 2.6 or 2.7
* Python [setuptools](https://pypi.python.org/pypi/setuptools)
* Python [pip](https://pypi.python.org/pypi/pip)
* Python [nose](https://nose.readthedocs.org/en/latest/)
* Python [coverage](http://nedbatchelder.com/code/coverage/)
* Python [pytest](http://pytest.org/)
* Python [pytest-cov](https://pypi.python.org/pypi/pytest-cov)

Below there are brief instructions for RedHat/Scientific Linux 6 and Ubuntu. 

For more detailed installation instructions and those for Windows and Mac OSX, check out the [general boot camp set-up instructions](http://software-carpentry.org/setup/) and the installation instructions specific to each product.

## To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    $ yum install nano
    $ yum install git
    $ yum install python
    $ yum install python-setuptools
    $ easy_install nose
    $ nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    $ yum install python-pip
    $ easy_install coverage
    $ pip install pytest
    $ pip-python install pytest-cov
    $ py.test
    ===================== test session starts ======================
    platform linux2 -- Python 2.6.6 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    =======================  in 0.01 seconds =======================

## To install under Ubuntu

Ubuntu 11.04 and above already comes with shell, vi and nano text editors, Python 2.7. To install the other packages run,

    $ sudo su -
    $ apt-get install git
    $ apt-get install python-setuptools
    $ easy_install nose
    $ nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    $ apt-get install python-pip
    $ easy_install coverage
    $ pip install pytest
    $ pip install pytest-cov
    $ py.test
    ===================== test session starts ===================== 
    platform linux2 -- Python 2.7.3 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    ===================== in 0.00 seconds ===================== 

