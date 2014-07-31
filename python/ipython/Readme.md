
# iPython Intro
[Up To Schedule](../../README.md) - Back To [Write Code for People](../writing_code_for_people/Readme.md)
* * * * *


**Based on Lecture Materials By: Milad Fatenejad, Katy Huff, Tommy Guy, Joshua 
R. Smith, Will Trimble, and many more**


## iPython
You can run python commands in a handful of ways; you can create executable scripts, you can run the python interpreter, you can run iPython, or you can run iPython notebook.  iPython is an alternative to the built-in Python interpreter with some nice features.  iPython notebook gives you interactive access to the python interpreter from within a browser window, and it allows you to save your commands as a "notebook".

Lets give the built-in interpreter a spin just this once:

```
frodgers@acibootcamp ~ $ python
Enthought Python Distribution -- www.enthought.com
Version: 7.3-2 (64-bit)

Python 2.7.3 |EPD 7.3-2 (64-bit)| (default, Apr 11 2012, 17:52:16) 
[GCC 4.1.2 20080704 (Red Hat 4.1.2-44)] on linux2
Type "credits", "demo" or "enthought" for more information.
>>> print "Hello World"
Hello World
>>> quit()
```

We can also write python commands in a file and execute them from the command line. You will notice that the print command above is located in the file hello.py. Execute the following command at the command line:

```
frodgers@acibootcamp ~ $ python ~/boot-camps/python/ipython/hello.py 
hello world
```

iPython has more useful features for interactive use than the standard python interpreter, so we'll use it from here on out:

```
frodgers@acibootcamp ~ $ ipython
Enthought Python Distribution -- www.enthought.com

Python 2.7.3 |EPD 7.3-2 (64-bit)| (default, Apr 11 2012, 17:52:16) 
Type "copyright", "credits" or "license" for more information.

IPython 0.12.1 -- An enhanced Interactive Python.
?         -> Introduction and overview of IPython's features.
%quickref -> Quick reference.
help      -> Python's own help system.
object?   -> Details about 'object', use 'object??' for extra details.

In [1]: print "Hello World"
Hello World

In [2]: 
```

### Pasting

Unfortunately pasting depends on your operating system and ssh program:

#### Windows

##### Putty
Click with the right mouse button over the window.

##### Bitvise SSH
Click with the right mouse button over the window and then "Paste".

#### Mac OSX
Press ⌘+V.

#### Linux
Click with the right mouse button over the window and then "Paste".

### History

iPython has a history. If you press the up and down keys, you can access the history. Try it now.

### Tab Completion

iPython also has tab completion of previous commands. Try typing "pr" and then hit the tab key. What if you type "pri" followed by tab?

### Getting Help

iPython has some nice help features.
If you wanted to see all the built-in commands available for something, use the dir command.

```
In [2]: dir(__builtin__)
['ArithmeticError',
 'AssertionError',
 'AttributeError',
 'BaseException',
 'BufferError',
 'BytesWarning',
 'DeprecationWarning',
 'EOFError',
 'Ellipsis',
 'EnvironmentError',
 'Exception',
 'False',
 'FloatingPointError',
 'FutureWarning',
 'GeneratorExit',
 'IOError',
 'ImportError',
 'ImportWarning',
 'IndentationError',
 'IndexError',
 'KeyError',
 'KeyboardInterrupt',
 'LookupError',
 'MemoryError',
 'NameError',
 'None',
 'NotImplemented',
 'NotImplementedError',
 'OSError',
 'OverflowError',
 'PendingDeprecationWarning',
 'ReferenceError',
 'RuntimeError',
 'RuntimeWarning',
 'StandardError',
 'StopIteration',
 'SyntaxError',
 'SyntaxWarning',
 'SystemError',
 'SystemExit',
 'TabError',
 'True',
 'TypeError',
 'UnboundLocalError',
 'UnicodeDecodeError',
 'UnicodeEncodeError',
 'UnicodeError',
 'UnicodeTranslateError',
 'UnicodeWarning',
 'UserWarning',
 'ValueError',
 'Warning',
 'WindowsError',
 'ZeroDivisionError',
 '__IPYTHON__',
 '__IPYTHON__active',
 '__debug__',
 '__doc__',
 '__import__',
 '__name__',
 '__package__',
 'abs',
 'all',
 'any',
 'apply',
 'basestring',
 'bin',
 'bool',
 'buffer',
 'bytearray',
 'bytes',
 'callable',
 'chr',
 'classmethod',
 'cmp',
 'coerce',
 'compile',
 'complex',
 'copyright',
 'credits',
 'delattr',
 'dict',
 'dir',
 'divmod',
 'dreload',
 'enumerate',
 'eval',
 'execfile',
 'file',
 'filter',
 'float',
 'format',
 'frozenset',
 'get_ipython',
 'getattr',
 'globals',
 'hasattr',
 'hash',
 'help',
 'hex',
 'id',
 'input',
 'int',
 'intern',
 'isinstance',
 'issubclass',
 'iter',
 'len',
 'license',
 'list',
 'locals',
 'long',
 'map',
 'max',
 'memoryview',
 'min',
 'next',
 'object',
 'oct',
 'open',
 'ord',
 'pow',
 'print',
 'property',
 'range',
 'raw_input',
 'reduce',
 'reload',
 'repr',
 'reversed',
 'round',
 'set',
 'setattr',
 'slice',
 'sorted',
 'staticmethod',
 'str',
 'sum',
 'super',
 'tuple',
 'type',
 'unichr',
 'unicode',
 'vars',
 'xrange',
 'zip']
```

str is the python name for strings.
Let's say we want to know what you can do with strings in python. We can type dir(str)

```
In [6]: dir(str)
['__add__',
 '__class__',
 '__contains__',
 '__delattr__',
 '__doc__',
 '__eq__',
 '__format__',
 '__ge__',
 '__getattribute__',
 '__getitem__',
 '__getnewargs__',
 '__getslice__',
 '__gt__',
 '__hash__',
 '__init__',
 '__le__',
 '__len__',
 '__lt__',
 '__mod__',
 '__mul__',
 '__ne__',
 '__new__',
 '__reduce__',
 '__reduce_ex__',
 '__repr__',
 '__rmod__',
 '__rmul__',
 '__setattr__',
 '__sizeof__',
 '__str__',
 '__subclasshook__',
 '_formatter_field_name_split',
 '_formatter_parser',
 'capitalize',
 'center',
 'count',
 'decode',
 'encode',
 'endswith',
 'expandtabs',
 'find',
 'format',
 'index',
 'isalnum',
 'isalpha',
 'isdigit',
 'islower',
 'isspace',
 'istitle',
 'isupper',
 'join',
 'ljust',
 'lower',
 'lstrip',
 'partition',
 'replace',
 'rfind',
 'rindex',
 'rjust',
 'rpartition',
 'rsplit',
 'rstrip',
 'split',
 'splitlines',
 'startswith',
 'strip',
 'swapcase',
 'title',
 'translate',
 'upper',
 'zfill']
```

Let's look up what some of these functions do.
? displays more information about each datatype/ function. Let's try str.swapcase?.

```
In [7]: str.swapcase?
Type:       method_descriptor
String Form:<method 'swapcase' of 'str' objects>
Namespace:  Python builtin
Docstring:
S.swapcase() -> string
```

Return a copy of the string S with uppercase characters
converted to lowercase and vice versa.

```
In [8]: "Hello world".swapcase()
Out[8]: 'hELLO WORLD'
```

***Excercise***
Can you find with help of the ? which function turns "Hello world" into "HELLO WORLD"?


### Executing code in files

If your code is in a file, you can execute it from the iPython shell with the **%run** command. Execute hello.py like so:

```
In [9] %run ~/boot-camps/python/ipython/hello.py
```

### Clearing iPython

To clear everything from iPython, use the reset command.

```
In [10] reset
Once deleted, variables cannot be recovered. Proceed (y/[n])?
```

Back to [Write Code for People](../writing_code_for_people/Readme.md)

