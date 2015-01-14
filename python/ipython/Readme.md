
# iPython Intro
[Up To Schedule](../../README.md) - Return To [Write Code for People](Readme.md#motivating-example)
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

Unfortunately pasting depends on your operating system and ssh program.  In
general, for multi-line pasting, you should use the `%cpaste` feature of iPython.

#### Windows

##### Git Bash

When in the iPython interpreter, the easiest way to paste is with the right
mouse button over the window, choosing "Paste".

##### Putty
Click with the right mouse button over the window.

#### Mac OSX
Press <kbd>⌘</kbd>+<kbd>V</kbd>.

#### Linux
Click with the right mouse button over the window and then "Paste".

### History

iPython has a history. If you press the <kbd>up</kbd> and <kbd>down</kbd>
keys, you can access the history. Try it now.

### Tab Completion

iPython also has tab completion of previous commands. Try typing "pr" and then
hit the <kbd>tab</kbd> key. What if you type "pri" followed by <kbd>tab</kbd>?

### Getting Help

iPython has some nice help features.

If you wanted to see all the built-in commands available for something, use
the `dir` command.

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

 ....

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

`str` is the python name for strings.
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

 ...
 
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
`?` displays more information about each datatype/ function. Let's try `str.swapcase?`.

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

----

![Exercise](../../common/pics/exercise.jpg) **Exercise**

Can you find with help of the ? which function turns "Hello world" into "HELLO WORLD"?

----

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



- - - - 

Back to [Write Code for People](Readme.md#motivating-example)

