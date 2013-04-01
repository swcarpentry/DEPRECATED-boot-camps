# Python in a couple of minutes

Python is interpreted,

    $ python
    >>> print 1+2
    3
    >>> numbers = [1,2,3]
    >>> print numbers
    [1, 2, 3]
    >>> print len(numbers)
    >>> print 1 == 2
    >>> print 1 == 1
    >>> print 'My ticket cost %d pounds!' % 50
    My ticket cost 50 pounds!

Variables are not typed,

    >>> value = 123
    >>> value = 'Hello'
    >>> value = [1, 2, 3, value]
    >>> print value
    [1, 2, 3, 'Hello']

Dictionaries are sets of key-value pairs (equivalent to association lists or hash-tables in other languages),

    >>> data = {}
    >>> data['id'] = 123
    >>> data['name'] = 'Fred'
    >>> data['count'] = 456
    >>> print data 
    {'id': 123, 'name': 'fred', 'count': 456}
    >>> print data['name']
    fred
    >>> print data.keys()
    ['id', 'name', 'count']
    >>> print data.items()
    [('id', 123), ('name', 'fred'), ('count', 456)]

To exit Python, press `CTRL-D` or type,

    exit()

We can write Python scripts,

    $ nano hello.py

    name = 'Fred'
    print 'Hello %s' % name
 
    $ python hello.py
    
Python supports loops and iteration,

    data = [4,5,6]
    for i in data: 
        print i

    data = 'hello'
    for i in data: 
        print i

Indentation is significant - blocks are denoted by indentation (typically 2 or 4 spaces). There are no brackets, braces or other block delimiters, so to end a block, use a blank line.

    $ python hello.py

Functions are declared using the `def` keyword,

    def hello(name):
        print 'Hello %s' % name

    hello('Fred')

    $ python hello.py
