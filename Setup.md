# Manual set-up

You should ensure you have the following software and tools available. 

* Web browser
* Bash shell - other shells are OK (e.g Cygwin) but please ensure they support the commands `grep`, `find`, `awk`, `cat`, `history`.
* A text editor e.g. nano, vi or emacs
* [Git](http://git-scm.com/)
* Python 2.6 or 2.7
* Python [pip](https://pypi.python.org/pypi/pip)
* Python [nose](https://nose.readthedocs.org/en/latest/)
* Python [coverage](http://nedbatchelder.com/code/coverage/)
* Python [pytest](http://pytest.org/)
* Python [pytest-cov](https://pypi.python.org/pypi/pytest-cov)
* Python [ipython](http://ipython.org)
* Python [numpy](http://numpy.scipy.org/)
* Python [scipy](http://www.scipy.org/)
* Python [matplotlib](http://matplotlib.org/)

Below there are brief instructions for RedHat/Scientific Linux 6 and Ubuntu 12.10. 

For more detailed installation instructions and those for Windows and Mac OSX, check out the [general boot camp set-up instructions](http://software-carpentry.org/setup/) and the installation instructions specific to each product.

If you have any problems then please e-mail the boot camp mailing list at [bath@lists.software-carpentry.org](mailto:bath@lists.software-carpentry.org).

## To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    # yum install nano
    # yum install git
    # yum install python
    # yum install python-nose
    # nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    # yum install python-coverage
    # nosetests --with-coverage
    ...
    # yum install pytest
    # py.test
    ===================== test session starts ======================
    platform linux2 -- Python 2.6.6 -- pytest-2.3.4
    collected 0 items 
    =======================  in 0.01 seconds =======================
    # yum install python-pip
    # pip-python install pytest-cov
    # py.test --cov .
    ===================== test session starts ======================
    platform linux2 -- Python 2.6.6 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    Coverage.py warning: No data was collected.
    ------- coverage: platform linux2, python 2.6.6-final-0 --------
    Name    Stmts   Miss  Cover
    ---------------------------
    =======================  in 0.05 seconds =======================
    $ easy_install ipython
    $ ipython
    Python 2.6.6 (r266:84292, Jun 18 2012, 09:57:52) 
    Type "copyright", "credits" or "license" for more information.

    IPython 0.10 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object'. ?object also works, ?? prints more.

    In [1]: CTRL-D
    $ yum install numpy
    $ python
    >>> import numpy
    >>> print numpy.__version__
    >>> from numpy import *
    >>> from numpy.linalg import *
    >>> a = array([[1.0, 2.0], [3.0, 4.0]])
    >>> eig(a)
    >>> CTRL-D
    $ yum install scipy
    $ python
    >>> import scipy
    >>> print scipy.__version__
    >>> from scipy import *
    >>> a = zeros(1000)
    >>> a[:100]=1
    >>> b=fft(a)
    >>> print b
    >>> CTRL-D
    $ yum install python-matplotlib
    $ python
    >>> import matplotlib
    >>> print matplotlib.__version__
    >>> CTRL-D
    $ ipython -pylab
    In []: plot([1,2,3])
    In []: from scipy import *
    In []: a = zeros(1000)
    In []: a[:100]=1
    In []: b=fft(a)
    In []: plot(abs(b))

## To install under Ubuntu

Ubuntu 12.10 and above already comes with shell, vi and nano text editors, Python 2.7. To install the other packages run,

    $ sudo su -
    # apt-get install nano
    # apt-get install git
    # apt-get install python
    # apt-get install python-nose
    # nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    # apt-get install python-coverage
    # nosetests --with-coverage
    ...

    # apt-get install python-pytest
    # py.test
    ===================== test session starts ===================== 
    platform linux2 -- Python 2.7.3 -- pytest-2.3.4
    collected 0 items 
    ===================== in 0.00 seconds ===================== 
    # apt-get install python-pip
    # pip install pytest-cov
    # py.test --cov .
    ===================== test session starts ======================
    platform linux2 -- Python 2.7.3 -- pytest-2.2.4
    collected 0 items 
    Coverage.py warning: No data was collected.
    ------- coverage: platform linux2, python 2.7.3-final-0 --------
    Name    Stmts   Miss  Cover
    ---------------------------
    =======================  in 0.03 seconds =======================
    $ easy_install ipython
    $ ipython
    Python 2.7.3 (default, Sep 26 2012, 21:53:58) 
    Type "copyright", "credits" or "license" for more information.

    IPython 0.13.1 -- An enhanced Interactive Python.
    ?         -> Introduction and overview of IPython's features.
    %quickref -> Quick reference.
    help      -> Python's own help system.
    object?   -> Details about 'object', use 'object??' for extra details.

    In [1]: CTRL-D
    $ apt-get install python-numpy
    $ python
    >>> import numpy
    >>> print numpy.__version__
    >>> from numpy import *
    >>> from numpy.linalg import *
    >>> a = array([[1.0, 2.0], [3.0, 4.0]])
    >>> eig(a)
    >>> CTRL-D
    $ apt-get install python-scipy
    $ python
    >>> import scipy
    >>> print scipy.__version__
    >>> from scipy import *
    >>> a = zeros(1000)
    >>> a[:100]=1
    >>> b=fft(a)
    >>> print b
    >>> CTRL-D
    $ apt-get install python-matplotlib
    $ python
    >>> import matplotlib
    >>> print matplotlib.__version__
    >>> CTRL-D
    $ ipython -pylab
    In []: plot([1,2,3])
    In []: from scipy import *
    In []: a = zeros(1000)
    In []: a[:100]=1
    In []: b=fft(a)
    In []: plot(abs(b))
