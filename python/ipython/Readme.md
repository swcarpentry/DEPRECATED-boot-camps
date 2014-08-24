
# iPython Intro
[Up To Schedule](../../README.md) - Return To [Write Code for People](../writing_code_for_people/Readme.md)
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
Press <kbd>⌘</kbd>+<kbd>V</kbd>.

#### Linux
Click with the right mouse button over the window and then "Paste".

### History

iPython has a history. If you press the <kbd>up</kbd> and <kbd>down</kbd> keys, you can access the history. Try it now.

### Tab Completion

iPython also has tab completion of previous commands. Try typing "pr" and then hit the <kbd>tab</kbd> key. What if you type "pri" followed by <kbd>tab</kbd>?

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

## Python basics

### Variables

All programming languages have variables, and python is no different. To create a variable, just name it and set it with the equals sign. One important caveat: variable names can only contain letters, numbers, and the underscore character. Let's set a variable.

```
In [1]: first_name = "Cliff"
In [2]: last_name = "Rodgers"
```

We can glue strings together with the + operator. Computers don't understand context.

```
In [3]: full_name = first_name + last_name
In [4]: print(full_name)
Out[4]: 'CliffRodgers'
```

***Exercise***
Can you add the extra space between my last and first name?

### Conditionals

A conditional (`if` statement) is some statement that in general says :
"When some boolean is true, do the following. Elsewise, do this other
thing."

Many equivalence test statements exist in Python that are similar in
other languages:

```python
i = 1
j = 2
i == j  # i is equal to j : False
i < j   # i is less than j : True
i <= j  # i is less than or equal to j : True
i > j   # i is greater than j : False
i >= j  # i is greater than or equal to j : False
i != j  # i is not equal to j : True
```

However, python has other equivalence test statements that are fairly
unique to python. To check whether an object is contained in a list :

```python 
beatle = "John"
beatles = ["George", "Ringo", "John", "Paul"]
print(beatle in beatles)     # is John one of the beatles? : TRUE
print("Katy" not in beatles) # this is also TRUE. 
```

Conditionals (`if` statements) are also really easy to use in python. Take
a look at the following example:

```python
i = 4
sign = "zero"
if i < 0:
    sign = "negative"

elif i > 0:
    sign = "positive"

else:
    print("Sign must be zero")
    print("Have a nice day")

print(sign)
```

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code.

For example, in MatLab you begin an if statement with the word `if`
and you end it with `end if`. In C/C++ you delimit code blocks with curly
braces.

Python uses *white space* — in this case, *indentation,* to group lines of
code. In this case, there are four spaces at the beginning of the line
following the `if` statement. This is not just to make things look pretty - it
tells Python what the body of the `if`-statement is.

Marking blocks of code is a fundamental part of a language. In C, you'll see
curly braces everywhere. In Python, you'll see indentation everywhere. It's
how you (and the computer) know what your code means.

As an important best practice, other white space in Python (and most other
languages) is only for people. Little things like putting blank lines in
"sane" places, and putting spaces between variables and operators (say, `a +
b` rather than `a+b`) can make your code a lot easier to read.

**Exercise**
Write an if statement that prints whether x is even or odd.

Hint: Try out what the "%" operator. What does 10 % 5 and 10 % 6 return?


- - - - 

Back to [Write Code for People](../writing_code_for_people/Readme.md)

