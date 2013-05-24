# Introduction to Python

## Overview
This high-level introduction to the programming language will focus on general purpose development with Python.  There will be a focus on hands-on exercises.

## Topics
 * The Python Interpreter
 * Simple Types (integers, strings, floats)
 * Functions and Indentation
 * Simple Expressions (truthiness, branching, whitespace)
 * Advanced Types (lists, tuples, dictionaries)

# The Python Interpreter

## REPL Functionality
* The Python Interpreter provides a REPL (Read, Evaluate, Print, Loop) interpreter
* Similar to MATLAB or Mathematica
* Good environment for learning, testing, evaluating, and debugging

Here are some simple examples of things you can try with the Python interpreter.  Lines that start with `#` are comments:
<!-- ``` {.python .numberLines} -->

``` {.python}
# Classic Hello World
print "Hello World!"
# Create a string and assign a reference to it
myName = "Aron"
# Concatenate and print strings.
print "Hello from " + myName
```

## Getting more with IPython

IPython, an alternative to the python interpreter developed by Fernando Perez, provides several enhancements over the
default interpreter:

* Tab completion
* `who` and `whos` provide workspace information
* `?` and `??` provide extra information about objects
* Magic functions such as `%run`, `%save`, `%reset`, `%timeit`, and `%paste`
* Improved debugging support
* Excellent history tracking (see `%history`)
* System shell commands with `!`

eWe will be taking advantage of several of these features through the workshop, see the
[iPython documentation](http://ipython.org/documentation.html) for more details.

## Online Python Tutor

 The [Online Python Tutor](http://people.csail.mit.edu/pgbovine/python/), developed by Philip Guo at Stanford
 University, provides line-by-line visualization of memory allocation and variable references on the heap.  This can be
 useful if you are having difficulty understanding how variable and memory references work in Python.

# Simple Types

## The String Container

We have introduced the str (string) object, which contains an immutable (more on this later!) array of characters.  The
string object in Python is very flexible, and supports a number of very useful methods out of the box.

### String Constructors
``` {.python}
# Single-line strings can be delimited with ' ' or " "
i1 = 'Mr. Aron Ahmadia'
i2 = "Mr. Lisandro Dalcin"
# Multi-line strings should be delimited with ''' ''' or """ """
sLong = 
"""Aron's favorite type of string delimiter
extends multiple lines.
"""
# Strings can also be created from numbers
sNum = str(32)
sFloat = str(32.4)
```

### Some Simple Operations

``` {.python}
# length
len(i1)
# forward indexing
i1[0]
# backward indexing 
i1[len(i1)-1]
# backward indexing (simplified)
i1[-1]
```

Note that for direct array access, Python syntax is equivalent to the 0-indexing used in C.

## Slicing 

Slicing is a much more powerful extension of the concept of a index.  Slices are composed of:

* `start` - beginning of the slice, defaults to `0` 
* `stop` - end of the slice, defaults to `len(str)`
* `step` - distance to travel to 'next' element in the slice, defaults to 1

### Visualizing Slicing
One trick to understanding slicing is to think of the indices as vertices, with the data between them as edges.  For
example, `ab` is the set of data between vertices `0` and `2`:  

```{.python}
s = 'abcde'
#  0   1   2   3   4   5
#  +---+---+---+---+---+
#  | a | b | c | d | e |
#  +---+---+---+---+---+
# -5  -4  -3  -2  -1  
```

### Slice Examples

```{.python}
i1 = 'Mr. Aron Ahmadia'
i1[:]
i1[::]
i1[1:]
i1[:-1]
i1[::2]
```

## Negative Slice Steps ( Ugly Truth!)
It is pretty straightforward to understand how slicing works for positive slice steps.  Unfortunately, the default
values of `start` and `stop`, as well as how they are interpreted, changes with negative stepping.  I don't recommend
using negative steps for non-trivial code.

### Negative Steps

* `start` - defaults to `-1` 
* `stop` - defaults to `-len(str) - 1` otherwise

### Visualizing Slicing (Negative steps)

```{.python}
#      0   1   2   3   4
#  +---+---+---+---+---+
#  | a | b | c | d | e |
#  +---+---+---+---+---+
# -6  -5  -4  -3  -2  -1  
```

### Reversing a string

```{.python}
i1[::-1]
```

## Some more String Operations

### Concatenation
```{.python}
theInstructors = i1 + ' and '  + i2
```

### Multiplication
```{.python}
'nah '*8 + 'Batman!'
```

### Substitution
```{.python}
theProfessors = theInstructors.replace('Mr.','Dr.')
```

### Search
```{.python}
theProfessors.find('Lisandro')
```

### Casing
```{.python}
s1 = "the Hitchhiker's guide to the galaxy"
s1.capitalize()
s1.upper()
s1.lower()
s1.title()
```

## Lists 

### Overview

Python gives us a selection of useful container types for manipulating collections of data.  As you may have guessed, a
string object is really a special type of container, one that is only allowed to hold immutable character data
(methods that modify a string object actually creates a new string!).  A very general-purpose container in Python is the
`list`, which contains an ordered, mutable collection of heterogeneous objects.   


### List Constructors
 
Lists can be constructed from other objects using the `list` constructor or directly on the interpeter using the
`[` and `]` bracket operators to demarcate the list, and commas (`,`) to separate each object.

```{.python}
l1 = [1, 2, 3]
l2 = ['a list with just a sentence in it']
l3 = ['a','list','of','strings','and',[],'empty lists',1,2,3,'and numbers']
l4 = []
```

### Length, Indexing and Slicing
```{.python}
len(l1)
l1[0]
l1[0:2]
```

## More List Operations

### Lists are mutable, so we can modify their contents
```{.python}
l1[0] = 'string munches 1'
l1.pop()  # returns the last element in the list and removes it, like a stack
l1.append(4) # places the object as the new last element in the list, like push in a stack 
```

### Concatenation and multiplication work the same way as for strings
```{.python}
l1+l2
l2*2
``` 

## Lists <--> Strings split and join

It is very common to have a string that we would like to `split` into a list of strings, or the reverse, a list of
strings we would like to `join` together into one string.

Python makes this very easy!

In both cases, there is frequently a *separator* that we would like to use either as as a delimiter when we are
separating values, or as a common string between the values in the list.  One thing that trips up newcomers to Python is
that `join` is a member function of string objects which takes a list as an argument, so you need to create the
separator string first, then pass it a list to join. Here are some examples:

```{.python}
l1 = ['This','wants','to','be','a','sentence.']
sep1 = '  '
sep1.join(l1)
' '.join(l1) # this looks funny, but it does the same thing, I promise!
s1 = 'This is a sentence, with a period.'
s1.split() # split with no arguments uses whitespace as a separator
s1.split(',') # but we can really split with anything we want
s1.split('') # except the empty string!
```

### Challenge - Create a list with each element corresponding to a single character from a string.

# Functions, Scope and Indentation

## Scope

Every programming language needs the ability to define and use functions.  Because Python is dynamically typed, you
don't need to declare a function type before using it, though you may need to bring it into the interpreter's *scope*
before you can call the function.  We will talk more about manipulating the interpreter scope later using `import`
statements later.

We **call** functions C-style, by appending a `( )` to the end of the function, with all arguments separated by
commas placed in between.

## Built-in Functions

Python has a number of functions so important that they start off in global scope.  Here are some of the more useful ones:

### type - type inspection of an object
```{.python}
x = 3
s = 'string'
type(x)
type(s)
```

### repr, eval - canonical string representation and scope-constrained expression evaluation
```{.python}
x_repr = repr(x)
s_repr = repr(s)
xx = eval(repr(x))
ss = eval(repr(ss))
```

## More Built-in Functions

### abs, pow, round - simple mathematical operators
```{.python}
abs(-3.2)
pow(9,-0.5)
round(3.2)
```

### bool, int, float, str - simple type constructors 
```{.python}
bool(x)
bool(0.0)
bool(-3.6)
bool(True)
int(3.2)
float(3)
str(int(False))
str(bool(1))
```

### cmp
```{.python}
# comparison between floating point and integral types is legitimate (through coercion)
cmp(2,3)
cmp(3.0,3)
cmp(1/5,0.2) 
# alphabetical comparisons between strings can be made
cmp('alpha','beta')
# 'dictionary' comparisons are legitimate as well
cmp('b','bee')
```

### Challenge - Why does cmp(1/5, 0.2) evaluate to -1?

*Note: cmp has been deprecated in Python3*

## Defining Simple Functions

Defining functions in Python is straightforward, we'll start with an example:

```{.python}
def double(n): return n*2
type(double)
double(3)
double(double(1))
```

Notice that we didn't have to specify the type of `n`, so any object that provides a `*` method
is a valid argument to our function.  For example:

```{.python}
double('Double Me!')
double(double('Quadruple Me!'))
```

Of course, we might want our function to have multiple statements:

```{.python}
def double_twice(n): x = n+n; return x+x
double_twice(1)
```

Note that functions do not have to accept any input arguments or return any parameters:

```{.python}
def do_nothing(): pass;  # pass is a Python keyword that means no operation
do_nothing()
```

## Defining Complex Functions

In many cases it will be necessary for a function to span a multiple line *block* of code. 

In C or C++, the blocks are demarcated by `{` and `}` brackets.  In modern Fortran, the block sentinel names uses an `end` followed by the initial
block statement, e.g. `do ... end do` or `program ... end program`

### Code Blocks in Python (Demarcation by Indentation)

In Python, code blocks are defined by their indentation (using spaces or tabs).  Here is an example:

```{.python}
def factorial(n):
  if n > 1: return n*factorial(n-1)
  else: return 1
  
factorial(1)
factorial(3)
```

Oops, looks like I introduced an `if/else` statement and some recursion there as well!  Note how the indentation of the
function block is 2 spaces, and how the end of the indentation tells the Python interpreter that the function definition
is ended.

We could also have indented the `if` and `else` blocks with a second level of indentation:

```{.python}
def factorial(n):
  if n > 1: 
    return n*factorial(n-1)
  else: 
    return 1
```

The difference in styles is mostly a matter of preference.

# Simple Expressions (truthiness, branching, looping)
 
## Boolean Expressions

Oftentimes in programming we wish to evaluate a statement for truth.   The unsafest (but simplest) was to test an
expression is to directly evaluate on the expression itself.   In Python, this is equivalent to constructing a `bool`
type with the expression as the parameter.   We refer to this method as evaluating the *truthiness* of an expression, as
it occassionally leads to some unexpected consequences.

### Truthiness

The following expressions are truthy:

```{.python}
bool(-0.1) # non-zero number
bool('hi') # non-empty string
bool(True) # boolean True
bool([0]) #non-empty list
```

The following expressions are falsey:

```{.python}
bool(None) # the special None object
bool(0) # zero
bool(False) # boolean False
bool('') # empty string
bool([]) # empty list
```

## is and not
 
It is often safer to explicitly test for equality to Truth, None, or False using the `is`identity test.  .  You can
negate a boolean expression by using the `not` statement, but when it occurs as `is not` the statement should be read as
a single phrase.  You can also use parenthesis to group expressions.  Try to guess the output of each of these before
entering them in the interpreter.    

```{.python}
bool(True is True)
bool(1 is True)
bool(1 is False)
bool(None is (not True))
boolValue = bool(None is not True) 
```

Of course, we don't actually need the `bool` constructor here:

```{.python}
True is True
1 is True
1 is False
None is (not True)
boolValue = None is not True
```

## Combining Boolean Expressions with and, or, any, all

Of course we want to be able to combine pairs and collections of `bool` in more complicated ways if we can:

### Comparison (<,>,==) and Logical Binary Operators (and, or)

```{.python}
a = 5; b = 6; c = 7;
b > a and b < c
b == 6 or b > 10
False or not False
1 <> 2  # don't use this, it has been deprecated
1 != 2 # this form is preferred
```
 
### Collective Binary Operators (any, all)

* `any` returns `True` if **any** object in a collection evaluates "truthy", `False` otherwise
* `all` returns `True` if **all** objects in a collection evaluate "truthy", `False` otherwise

```{.python}
l = ['a list with an empty item',[]]
any(l)
all(l)
```

## Simple Branching with if, elseif, and else

### Ternary Expression Branching
Python allows us to return one of two values based on a Boolean value

```{.python}
test_value = True
response = "It worked!" if test_value else "It failed!"
```

### the assignment isn't necessary if you just want to conditionally call one of two functions

```{.python}
def truth_function(): print "truth_function called!" 
def false_function(): print "false_function called!"
  
truth_function() if test_value else false_function()
```

### now what happens if test_value is False?

```{.python}
test_value = False
truth_function() if test_value else false_function()
```

### Block Expression Branching
You might be more familiar with the following style of branching, which is also legal Python:
 
```{.python}
def my_test(test_value):
  if test_value:
    truth_function()
  elif test_value is None:
    print "test_value is None!"
  else:
    false_function()

my_test(True)
my_test(None)
my_test(False)    
```

## Looping
Sometimes, you just need to write a for loop.  Python makes this easy.  

### For loops

The first thing we need in a for loop is a range to iterate over.  This range can be constructed explicitly using `range`

```{.python}
range(5)
range(5,10)
range(5,10,2)
range(10,5,-1)
# xrange returns an iterator instead of a list, we will discuss this in more detail this afternoon
xrange(5)
```

Here is a simple example of a function that uses a `for` to solve  FizzBuzz problem.  

```{.python}
def FizzBuzz(r):
  """Solves the FizzBuzz problem to r

  (http://www.codinghorror.com/blog/2007/02/fizzbuzz-the-programmers-stairway-to-heaven.html)
  """
  for i in range(1,r+1):
    if i%3 is 0 and i%5 is 0: print "FizzBuzz"
    elif i%3 is 0: print "Fizz"
    elif i%5 is 0: print "Buzz"
    else: print i
  
FizzBuzz(5)
FizzBuzz(20)
```

There is often a more elegant and functional way to express the work performed in a loop (though not necessarily easier to read or
understand!), we will go more into this in the afternoon session.

## Advanced For Loops: Break, and Else

In this example, we iterate directly over the items of a list.  Note how `break` can be used to exit early from a loop.
The `else` clause (notice indentation) is  executed if the `for` loop goes through all iterations (and does not
encounter a`break` statement) 

```{.python}
def drill_oil(test_site):
  """Searches for 'oil' or 'black gold' in the iterable test_site, prints 'struck oil!' 
      and returns immediately if oil found or the current entry otherwise.  If no 
      oil is found, 'drill again!' is also printed.
  """
  for r in test_site:
    if r is 'oil' or r is 'black gold': print 'struck oil!'; break
    else: print "found " + r
  else: print 'drill again!'
  
drill_oil(['shale', 'shale', 'sand'])
drill_oil(['oil', 'oil', 'more oil', 'black gold'])
```

# Unsorted Containers 

## Hashing and Comparisons

Up until now, we have only considered ordered collections.  The order of the collection is preserved across
modifications to the collection, and indices and slices can be taken into it.  Now we consider two important *unordered*
collections, the `set`, a mutable unordered collection of unique objects (no duplicate elements), and the
`dict`, a mutable unordered and collection of unique keys and their associated values.  Note that the keys must be
unique, but the values may be nonunique.

Unordered sets do not directly use the `cmp` function to test for equality, instead they use the `hash` function.

## Sets

Python also possesses the `set` container for holding a mutable, unordered collections with no duplicate elements.
There is no syntax for assembling a set from elements, you must call the constructor explicitly on a list or another
container (simply calling the `set` constructor with no arguments will create an empty set):

```{.python}
s1 = set([1,2,3])
s2 = set([3,4,5])
```

### Basic Manipulation

```{.python}
s2.add(6)
s1.discard(1)
```

### Mathematical Set Operations

```{.python}
s1.union(s2)
s1 - s2
s1.intersection(s2)
s1.difference(s2)
s1.isdisjoint(s2)
```

## Tuples and Dicts

### Tuple Basics

Tuples are immutable, ordered, sequences.  They are very similar to lists, though they cannot be modified.

```{.python}
# tuples are constructed using parenthesis ( ), with commas to separate elements
t1 = (0,'stuff')
# this is also okay, but occasionally dangerous syntax
t2 = 1,'other stuff'
```

Here we are creating tuples of length 2, but we note that tuples may be of arbitrary length.
 The `dict` is one of the most powerful native containers in Python, providing an efficient way to store and index data.
Internally, the `dict` stores each `key` and its associated `value` together in a 2-tuple, and it is often convenient to
express them this way when inserting or extracting data.

### Dict Constructors

Dicts can be directly constructed using angled brackets `{ }`, a colon `:` to separate each key from its associated
value, and a comma `,` to separate the key-value pairs:

```{.python}
d1 = {'roses':'red', 
      'violets':'blue', 
      'the answer to life, the universe, and everything else':42, 
      13:'an unlucky number'}
```

They can also be constructed from a sequence of pairs:

```{.python}
d2 = dict([('roses','red'), ('violets','blue'), ('rhyming','hard')])
```

## More on Dicts

Dicts will accept a wide variety of objects as keys, though the most important restriction is that the object's type
must be `hashable`.  Basic integers, floating point values, strings, and tuples are all hashable, sets and lists are not.

### Testing for Existence of a Given Key
```{.python}
d = dict([('roses','red'), ('violets','blue'), ('rhyming','hard')])
# tests if the string 'roses' is a key in the dictionary
'roses' in d
# tests if the string 'red' is in the dictionary
'red' in d
```

### Accessing, Adding and Removing Data

```{.python}
# you can access any value in the dict by providing the key as an argument using [ ] brackets
d['roses']
# you can also modify the value associated with a given key (or add it if it does not exist) 
d['roses'] = 'green'
d['tulips'] = 'blue'
# use the 'pop' method to remove keys (and their associated value) from the dictionary
# the value associated with the removed key is returned from the function call
d.pop('roses') == 'green'
```

## More on Tuples

### Tuples can be constructed from other containers

```{.python}
t3 = tuple([0,2,3])
t4 = tuple(set(['stuff',0]))
```

### Tuples can be indexed and sliced, but not modified

```{.python}
t5 = t1[0:2]
t2 = t2[0]
```

### You can use tuples to 'pack' and 'unpack' objects

```{.python}
def gimme_five(): return (1,2,3,4,5)

(a,b,c,d,e) = gimme_five()
```

## Summary of Basic Python Types

### Basic Data
* `int` - integer data
* `bool` - boolean True/False values
* `float` - floating-point data

### Containers
* str - ordered, immutable character strings
* list - ordered, mutable sequences of data
* tuple - ordered, immutable sequences of data
* set - unordered, mutable unique data
* dict - unordered, mutable mapping of unique keys to values

## Iterating with in (Afternoon Preview)

One powerful feature of Python is the way the same generic syntax can be used to iterate through containers that are
Iterable.  This inclues all of the collections we have introduced so far.

```{.python}
s = 'string'
l = [0,1,2,3]
t = ('a','b','c','d')
d = dict([('roses','red'), ('violets','blue'), ('rhyming','hard')])

for char in s: print char
for item in l: print item
for item in t: print item
for key in d: print key
```
  
# An Aside on Object Typing

## Static vs. Dynamic Typing

* statically typed - All object types are fixed at compile time.  Examples include C, Fortran, Java, and Go.  
* dynamically typed - The above requirement is relaxed, types are discovered at execution time.  Examples include
  Python, Ruby, Lisp, MATLAB, and Javascript.

Advantages of static typing include compile-time detection of most type errors, at the expense of more verbose code.
Dynamically typed languages are typically more flexible and less verbose, but are susceptible to type errors at run-time.  The
formalism offered by static typing usually lends itself better to business and mission-critical operation, the
flexibility offered by dynamic typing usually lends itself better to prototyping and scientific software development.
Since programs in dynamically typed languages are harder to check statically, programmers usually leverage techniques such as
unit testing and function-level documentation to improve program robustness.

## Strong vs. Weak Typing

* strongly typed - Restrictions exist for mixing objects of different types, and casting (implicit and explicit) are
  disallowed.  Python is strongly typed, you cannot 'cast' data as in C or C++, you must convert the object between
  types.  
* weakly typed - Either implicit casting is supported in the language (Perl, Javascript, Visual Basic), or casting is
  allowed (C, C++).  
  
Weak typing can be very flexible, but extremely dangerous as well, see
[Gary Bernhardt's "Wat" talk](https://www.destroyallsoftware.com/talks/wat) for some good examples of the dangers of
weak typing.  

For a good discussion of other salient details of typing, see cdsmith's ["What to know before debating type systems"](http://blogs.perl.org/users/ovid/2010/08/what-to-know-before-debating-type-systems.html).

