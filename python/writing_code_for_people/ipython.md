
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

----
![Exercise](pics/exercise.jpg) ***Exercise***

Can you add the extra space between my last and first name?

----

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

----
![Exercise](pics/exercise.jpg) **Exercise**
Write an if statement that prints whether x is even or odd.

Hint: Try out what the "%" operator. What does 10 % 5 and 10 % 6 return?

----

### Loops

Most languages offer two main kinds of loops.  Choosing which kind of loop to
use in any given case can also contribute to readability.  In python, these
are the `while` loop and the `for` loop.  The two different kinds of loops can
be forced to behave similarly to each other, but for readability it is best to
use:
* `while` loops when the number of times the loop will execute is not known
* `for` loops when the number of times the loop will execute is known

While Loops
===========

Let's start by looking at `while` loops since they function like while
loops in many other languages. The example below takes a list of
integers and computes the product of each number in the list up to the
-1 element.

A `while` loop will repeat the instructions within itself until the
conditional that defines it is no longer true.

```python
mult = 1
sequence = [1, 5, 7, 9, 3, -1, 5, 3]
while sequence[0] != -1:
    mult = mult * sequence[0]
    sequence.pop(0)

print mult
```

Some new syntax has been introduced in this example.

-   On line 4, we compute the product of the elements just to make this
    more interesting.

-   On line 5, we use list.pop() to remove the first element of
    the list, shifting every element down one.
    Let's verify this with `sequence.pop?`

Notice that the number of times the loop executes depends on the entries of
`sequence` and not on something easily computable such as the length of
`sequence`.

**Watch Out**

Since a `while` loop will continue until its conditional is no longer
true, a **poorly formed** while loop might repeat forever. For example :

```python
i=1
print "Well, there's egg and bacon, egg and spam, egg bacon and"
while i == 1:
  print "spam "
print "or Lobster Thermidor a Crevette with a mornay sauce served in a Provencale manner with shallots..."
```

Since the variable `i` never changes within the while loop, we can expect that
the conditional, `i=1` will remain true forever and the while loop will just
go round and round, as if this restaurant offered nothing but spam. (If you
try this at home, please note that one way to interrupt a non-terminating
process is <kbd>ctrl</kbd>+</kbd>C</kbd> or <kbd>ctrl</kbd>+</kbd>Z</kbd>.)

For Loops
=========

`For` loops in python operate a little differently from other languages.
Let's start with a simple example which prints all of the numbers from 0
to 9:

```python
for i in range(10):
    print i
```

You may be wondering how this works. Start by using `help(range)` to see
what the range function does.

```
range()?
```

    Help on built-in function range in module __builtin__:

    range(...)
        range([start,] stop[, step]) -> list of integers

        Return a list containing an arithmetic progression of integers.
        range(i, j) returns [i, i+1, i+2, ..., j-1]; start (!) defaults to 0.
        When step is given, it specifies the increment (or decrement).
        For example, range(4) returns [0, 1, 2, 3].  The end point is omitted!
        These are exactly the valid indices for a list of 4 elements.

`Range` is a function that returns a list containing a sequence of
integers. So, `range(10)` returns the list [0,1,2,3,4,5,6,7,8,9]. The `for`
loop then simply iterates over that list, setting `i` to each value.

----
![Exercise](pics/exercise.jpg) **Exercise**

Using a loop, calculate the factorial of 6 (the product of all positive integers up to and including 6).

----

break, continue, and else
=========================

A `break` statement cuts off a loop from within an inner loop. It helps
avoid infinite loops by cutting off loops when they're clearly going
nowhere.

```python
reasonable = 5
for n in range(1,10):
    if n == reasonable :
        break
    print n
```

Something you might want to do instead of breaking is to continue to the
next iteration of a loop, giving up on the current one..

```python
reasonable = 5
for n in range(1,10):
    if n == reasonable :
        continue
    print n
```

What is the difference between the output of these two?

We can combine loops and flow control to take actions that are 
more complex, and that depend on the data. First, let us define 
a dictionary with some names and titles, and a list with a 
subset of the names we will treat differently.

```python
favorite = 3
for n in range(1,10):
    if n == favorite :
        print "The number " + str(n) + " is my favorite number"
    else:
        print "The number " + str(n) + " is nothing special"
    print n
```


- - - - 

Back to [Write Code for People](../writing_code_for_people/Readme.md)

