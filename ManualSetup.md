# Manual set-up

If you decide to manually set-up your laptop then you should ensure you have the following software and tools available. 

This page only contains instructions for RedHat/Scientific Linux 6 and Ubuntu. For Windows and Mac OSX, check the installation instructions specific to each product.

## Version control session

You will need,

* Bash shell
* A text editor e.g. nano, vi or emacs
* [Git](http://git-scm.com/)

### To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    # yum install nano
    # yum install git

### To install under Ubuntu

Ubuntu 11.04 and above already comes with shell, vi and nano text editors. To install the other packages run,

    $ sudo su -
    # apt-get install git

## Testing session

You will need,

* Bash shell
* A text editor e.g. nano, vi or emacs
* Python 2.6 or 2.7
* Python [setuptools](https://pypi.python.org/pypi/setuptools)
* Python [pip](https://pypi.python.org/pypi/pip)
* Python [nose](https://nose.readthedocs.org/en/latest/)
* Python [coverage](http://nedbatchelder.com/code/coverage/)
* Python [pytest](http://pytest.org/)
* Python [pytest-cov](https://pypi.python.org/pypi/pytest-cov)

### To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    # yum install nano
    # yum install python
    # yum install python-setuptools
    # easy_install nose
    # nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    # yum install python-pip
    # easy_install coverage
    # pip install pytest
    # pip-python install pytest-cov
    # py.test
    ===================== test session starts ======================
    platform linux2 -- Python 2.6.6 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    =======================  in 0.01 seconds =======================

### To install under Ubuntu

Ubuntu 11.04 and above already comes with shell, vi and nano text editors, Python 2.7. To install the other packages run,

    $ sudo su -
    # apt-get install python-setuptools
    # easy_install nose
    # nosetests
    ------------------------------------------------------------------
    Ran 0 tests in 0.003s
    OK
    # apt-get install python-pip
    # easy_install coverage
    # pip install pytest
    # pip install pytest-cov
    # py.test
    ===================== test session starts ===================== 
    platform linux2 -- Python 2.7.3 -- pytest-2.3.4
    plugins: cov
    collected 0 items 
    ===================== in 0.00 seconds ===================== 

## Data management session

You will need,

* Bash shell
* A text editor e.g. nano, vi or emacs
* Python 2.6 or 2.7
* Python [setuptools](https://pypi.python.org/pypi/setuptools)
* [MongoDB](http://www.mongodb.org/)
* Python [pymongo](http://api.mongodb.org/python/current/)

### To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    # yum install nano
    # yum install python
    # yum install python-setuptools

Edit `/etc/yum.repos.d/10gen.repo` and, for a 32-bit system add the lines,

    [10gen]
    name=10gen Repository
    baseurl=http://downloads-distro.mongodb.org/repo/redhat/os/i686
    gpgcheck=0

For a 64-bit system add the lines,

    [10gen]
    name=10gen Repository
    baseurl=http://downloads-distro.mongodb.org/repo/redhat/os/x86_64
    gpgcheck=0
    enabled=1

    # yum install mongo-10gen
    # yum install mongo-10gen-server
    # /sbin/service mongod start
    Starting mongod: forked process: 4357
                                         [  OK  ]
    all output going to: /var/log/mongo/mongod.log
    # /sbin/service mongod status
    mongod (pid 4357) is running...
    # easy_install pymongo

### To install under Ubuntu

Ubuntu 11.04 and above already comes with shell, vi and nano text editors, Python 2.7. To install the other packages run,

    $ sudo su -
    # apt-get install python-setuptools
    # sudo apt-key adv --keyserver keyserver.ubuntu.com --recv 7F0CEB10

Edit `/etc/apt/sources.list.d/10gen.list` and add the lines,

    deb http://downloads-distro.mongodb.org/repo/ubuntu-upstart dist 10gen

    # apt-get update
    # apt-get install mongodb-10gen
    # service mongodb start
    # easy_install pymongo

### To check MongoDB and pymongo are OK...

Do the following to create a new MongoDB database called "example", collection called "students" and add a document to the collection,

    $ python
    Python 2.7.3 (default, Sep 26 2012, 21:53:58) 
    [GCC 4.7.2] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from pymongo import MongoClient
    >>> connection = MongoClient()
    >>> db = connection["example"]
    >>> print connection.database_names()
    [u'example', u'local']
    >>> collection = db["students"]
    >>> doc={}
    >>> doc["name"]="Fred"
    >>> doc["location"]="Manchester"
    >>> print doc
    {'name': 'Fred', 'location': 'Manchester'}
    >>> doc_id = collection.insert(doc)
    >>> print doc_id
    512e30c41d41c80b6d39dce7
    >>> CTRL-D

Do the following to show that the document is still in the collection in the database, then delete the database,

    $ python
    Python 2.7.3 (default, Sep 26 2012, 21:53:58) 
    [GCC 4.7.2] on linux2
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from pymongo import MongoClient
    >>> connection = MongoClient()
    >>> db = connection["example"]
    >>> collection = db["students"]
    >>> print collection.find_one()
    {u'_id': ObjectId('512e30c41d41c80b6d39dce7'), u'name': u'Fred', u'location': u'Manchester'}
    >>> print connection.database_names()
    [u'example', u'local']
    >>> connection.drop_database("example")
    >>> print connection.database_names()
    [u'local']
    >>> CTRL-D
