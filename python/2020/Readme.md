# Introductory Python

**Presented By: Jonathan Cooper and Robin Freeman**

**Based on Software Carpentry Lecture Materials and http://scipy-lectures.github.com/intro/language/python_language.html**

What is Python ?
================

Python is an interpreted (byte-compiled) language. Its simple, high
level, human readable interface speeds the programming process in
exchange for some computation time overhead. Python is implemented
mostly in the C programming language, so, as python develops, it is
increasingly possible to do everything in Python that is possible in C.
Python is also free and open source, so if you find a bug or generate a
useful module, the Python Software Foundation will likely be happy to
merge your changes into the language.

During this session you are going to learn about some very basics about
how to execute python code as well as some examples of the built-in
Python data types.

Built-in data types are the basic building blocks of Python programs.
They are really basic things like strings and numbers (either integers,
complex or floating point numbers). There are simple containers like
lists (think of lists as arrays or vectors), tuples and dictionaries.

Useful Python Links
-------------------

### Core Python

 - Main Python Docs
   - [http://docs.python.org](http://docs.python.org)/
 - Global Module Index
   - Built in modules like os, sys, datetime, math, random...
   - [http://docs.python.org/modindex.html](http://docs.python.org/modindex.html)
 - Built-in Functions
   - Built-in, always available functions like open, range, enumerate,
>     zip...
   - [http://docs.python.org/library/functions.html](http://docs.python.org/library/functions.html)
 - String Formatting
   - [http://docs.python.org/library/string.html\#formatstrings](http://docs.python.org/library/string.html#formatstrings)

### Python in Science

 - [NumPy](http://numpy.scipy.org/)
   - Fast arrays, used by almost every scientific Python package
 - [SciPy](http://www.scipy.org/)
   - Minimization, fitting, solvers, statistics, and more
 - [matplotlib](http://matplotlib.sourceforge.net/)
   - 2D and 3D plotting, maps
 - [AstroPy](http://astropy.org) for astronomy
 - [Biopython](http://biopython.org/wiki/Biopython) for bioinformatics
 - [Sage](http://www.sagemath.org/) for mathematic analysis
 - [SymPy](http://sympy.org/en/index.html) for symbolic mathematics
 - [pandas](http://pandas.pydata.org/) data analysis and stats for
    categorical and time series data

Hello World
===========

First, we will use python ''interactively''. This means that we will
type commands directly into iPython. Once we start performing more
complicated tasks we will start writing Python scripts and programs in a
text editor, outside of the interpreter.

To get to the python shell, type **python** into the terminal.

```python
   >>> print "Hello World"
   Hello World
   >>> exit()
```

To get to the interactive python interpreter, a more sophisticated
python shell, type **ipython** into the terminal.

```python
   In [1]: print "Hello World"
   Hello World
   In [2]: exit
```

You can also put the commands in a **.py** file and execute that file in
the terminal by typing **python [filename]**

    $ gedit myfile.py &
    <edit myfile with the hello world program.>
    $ python myfile.py
    Hello World!

Pasting into iPython
====================

**Note:**

To paste text from another application (i.e. the internet browser) into
iPython:

1.  select text from the wiki
2.  copy with **ctrl+c**
3.  in iPython, type **%paste**

The code should paste and execute in iPython.

### History

iPython has a history. If you press the up and down keys, you can access the history.

### Executing code in files

If your code is in a file, you can execute it from the iPython shell with the **%run** command. Execute hello.py like so

```python
In [1] %run hello.py
```

### Clearing iPython

To clear everything from iPython, use the **%reset** command.

```python
In [1] %reset
Once deleted, variables cannot be recovered. Proceed (y/[n])?
```

Variables
=========

Variables are names, while values are the data assigned to those names.

Questions : Variables and Values
--------------------------------

In the code snippet:

```python
    a=2
    b="string"
    c=a
```

 - What is the value of the variable `c`?
 - What is the value of the variable `b` ?
 - What is the name given to the value `2` ?

(The last one is a trick, the value 2 has two names.)

Strings and Numbers
===================

It is really easy to make variables in python. For example, to create a
string, `s`, and print its value, simply type the following into
iPython:

```python
    s = "Hello World"
    print s
```

If you want to see what the type of a variable is, you can use the
built-in python function, `type`. Just enter

```python
   print type(s)
```

into iPython and you should see something like this:

```python
      <type 'str'>
```

This tells us that `s` is of type **str** (i.e. that `s` is a
string). Making numeric variables is equally easy and intuitive. Try
entering the following into IPython. Notice that the \# symbol is used
to start comments so everything after the pound sign is ignored.

```python
   i,r,c = -10, 3.5, 1.0 + 2j  # set i to -10, r to 3.5 and c to 1.0+2j
```

This one line sets the variable `i` to the integer -10 , `r` to the
floating point value 3.5 (a floating point number is just a
real/non-integer number) and `c` to the value 1.0 + 2j (Notice, how
easy and intuitive it is in python to set multiple variables to
something. You'll discover a lot of similar syntax that is designed to
make your life easier). Let's use the built-in type function to determine
the type of each of the three variables we just created:

```python
   print type(i), type(r), type(c)
```

This will give :
```python
    <type 'int'> <type 'float'> <type 'complex'>
```

This tells us that "i" is an integer, "r" is a floating point number,
and "c" is a complex number. As you can see, Python has built-in support
for imaginary numbers!

**Aside: Long integers** Another way python makes our lives easier is by
allowing integers to be arbitrary large. In languages like C/C++ and
FORTRAN integer variables can only store values up to a certain size.
But entering and manipulating the following forty digit number with
iPython is no problem:

```python
   i = 1234567890123456789012345678901234567890
   print i * 6
```

Operations in Python are defined by their type. For instance, look the
difference between these operations:

```python
   In[1]:  1 + 3
     4
   In[2]:  1.0 + 3
     4.0  # This is a float
   In[3]: "Hello " + "world"
     'Hello world'
   In[4]: 1 + "Hello"
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   TypeError: unsupported operand type(s) for +: 'int' and 'str'
```

In the first two cases, addition between numbers meant that 1 was added
to 3 using the standard type rules (float plus int = float). In the
third case, the command was string addition, which concatenates two
strings. The final case broke because an 'int' type can not be added to
a 'str' type. This is because it's unclear how to interpret an int as a
string: should it be the string representation, the ASCII character
code, or something else entirely?

One way to handle this is to explicitly convert the int into a string:

```python
     str(1) + "Hello"
```

Equivalent functions exist for converting to **int**, **float**, and
other types. Let's say you had numerical data in a string.

```python
In [23]: resistanceString = "4.0"

In [24]: resistance = float(resistanceString)

In [25]: resistance
Out[25]: 4.0

In [26]: type(resistance)
Out[26]: <type 'float'>
```

What would happen if you tried to coerce resistanceString to an int? What about coercing resistance to an int? Consider the following:

```python
In [27] resistanceString = "4.0 ohms"
```

Do you think you can coerce that string to a numerical type?


Basic data types in Python have a lot of functionality already built in.
For example, lets say that you are reading names from a file one line at
a time and that sometimes the names have leading and trailing spaces
that we want to strip away. We can just use the `strip` string method
to accomplish this. For example, type the following into iPython:

```python

  In[1]: name = "   Milad    "
  In[2]: print name + "is here"
        Milad     is here
```

Now enter `name.strip()` instead of `name`:

```python
  In[1]: print name.strip() + " is here"
   Milad is here
```

Notice that the extra spaces are gone. We used the `strip()` method,
which removes leading and trailing white space from strings. You can
think of a method as being a function that is attached to a particular
variable. You call methods by typing: `<variable>.<method name>`.

**Aside : Tab Completion**

Maybe you've noticed this already, but check out what happens you begin
typing a variable name (the first two letters of name, for example) and
press tab.

Convenient, right? This is also true of many built in functions.

Dynamic Typing
==============

Importantly, python is a **dynamically typed** language. That is, an
explicit type is not needed when creating a variable. Also, this means
that variables in Python which are initialized to a variable of one type
can be re-assigned to a variable of a different type. Try this:

```python
     sillystring = "What is the airspeed velocity of an unladen swallow?"
     print type(sillystring)
```

You'll see:

```python
     <type 'str'>
```

If you reassign silly string to an integer, what happens? That is, when
you type :

```python
    sillystring = 98    
    print type(sillystring)
```

You should see:

```python
     <type 'int'>
```

This is an interesting feature. Can you think of ways it can be helpful?
Are there ways it might be troublesome?

What is the type of sillystring after this :

```python
    sillystring += 0.1
```

**Aside: In Place Equivalency**

What is the += syntax about? This is an in-place way to write `sillystring =
sillystring + 0.1`. It is common in a number of languages.

Importantly, though we do not explicitly state them, variables always have
exactly one type. The number 98 is an **int**. For the variable holding
this value to be treated as a float, it must be assigned as **98.0**.

Questions : Dynamic Typing
--------------------------

Imagine that I first assign :

```python
    a=2
```

Then, I assign :

```python
    a="Welcome to the ministry of silly walks."
```

What has happened to the memory that was pointing to the number 2??

Getting Help
============

One of the really nice features in Python is that a lot of the help and
documentation is built into the code. Practically, this means that much
of the time you don't have to go digging through some web site to find
help. You can get help in Python using the `help` function. Let's look
at an example - enter

```python
    help(str.strip)
```

into IPython. You should then see documentation for the strip method pop
up. (NOTE: if you don't automatically return to the python interpreter,
just hit "`q`" to exit the help screen). You can also use the question
mark, "`?`", character to display the documentation as well. For
example, enter

```python
    str.strip?
```

into IPython to view the documentation.

Now try entering

```python
    help(str)
```

You should see documentation for the entire string type, including all
of the string methods. This can be useful when you are trying to perform
a specific task, but you don't know the right function to call. For
example, lets say we want to convert the string "cooper" to uppercase,
and we want to know if there is a string method which can do the job for
us. Start by typing "`help(str)`" to pull up the string documentation.
You can scroll through the string methods until you find a method called
"upper" which has documentation that looks like:

    |  upper(...)
    |      S.upper() -> string
    |      |      Return a copy of the string S converted to uppercase.

These lines tell us that the string class has a method called "upper"
which can be used to convert strings to uppercase. Now enter:

```python
    name = "cooper"   print name.upper()
```

At which point, you should see the word "COOPER" printed to the screen.

**Aside: Using Methods Directly on Data**

In the previous example, we first created a string variable, `name`,
assigned it the value "cooper", then used the `upper` string method to
obtain the uppercased version of the string. We didn't have to create a
variable, however. We could simply enter:

```python
    print "cooper".upper()
```

To generate the uppercased version.

As we saw above, the **str** type has a lot of documentation associated
with it, and we had to sift through most of it to find the upper method.
If we had a way to simply print all of the **str** methods, we could
have probably figured out that the `upper` method is what we wanted by
the name and in a lot less time. Luckily, python has a built in
function, "`dir`", for just this situation. The `dir` function takes
a type name and prints all of the methods associated. Try entering
"`print dir(str)`" to see a list of every method and variable
associated with the string class. You can ignore the methods that start
and end with double underscores for now. Try printing the methods
associated with the **int**, and **complex** types.

Finally, there are some really basic functions that are built right into
python that we have been using. For example, we used the "float" function
above to convert a string to a floating point number. You can see a list of
built in functions by entering `dir(__builtins__)`.  If you see something
interesting, such as the `zip` function, you can examine what it does using
help(zip).

Example : Manipulating Basic Data Types
---------------------------------------

Use the basic data types we've learned about along with the `help` and
`dir` functions to figure out how to do the following using either one
function or one method call:

- Take the absolute value of the number -1.4
- Begin with the string "a MaN and His DOG" and create the string "A man
  and his dog"
- Return the position of the character 'e' in the string "my test string"
  (The answer is 4, since `m` is is at position 0 not position 1)


Basic Maths
===========

All the standard operations are available: + - * / ** %

```python
    2 + 5
    2 - 5
    4 * 5
    6 / 3
    5 / 3
    5.0 / 3
    4 ** 2
    4 ** 0.5
    20 % 9
    2 + 3 * 4
    (2 + 3) * 4
```

Compound Data Types
===================

Python would be a fairly useless language if it weren't for the compound data types.
The main two are lists and dictionaries, but I'll mention sets and tuples as well.

Lists
-----

```python
    l = ['zero', 1, 2.0, '3', 4, 5.0, 'six']
    print l[0], l[1] # Indices start at 0 like C, not Matlab
    len(l)
    l[1:4] # Slicing.  Note doesn't include end index
    l[4:]  # Both start and stop are optional.
    l[:4]  # First 4 things
    l[::2] # Can also give a step, which defaults to 1.
    print l[-1], l[-2], l[-2:] # Indexing from the end
    l[1] = l[1] + 99
    l
```

Just like strings have methods, lists do too.

```python
    dir(list)
```

One useful method is append, for adding data to the end of a list.

```python
    l.append(3.0)
    l.append(4.0)
    print l
```

If you have multiple items (i.e. another list!) to add, use extend.

```python
    l.extend([1, 'a', 2.3])
    print l
```

Importantly, lists are **mutable**.  This is particularly important given how Python deals with variable assignments,
When you set a variable equal to another, both variables point to the same thing.
If this thing is mutable, changing the first one ends up changing the second. Be careful about this fact.

```python
In [19]: a = [1,2]

In [20]: b = a

In [21]: a.append(10)

In [22]: b
Out[22]: [1, 2, 10]
```

Tuples
------

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is immutable; once you set the
data in it, the tuple cannot be changed. You define a tuple as follows.

```python
In [1]: tup = ("red", "white", "blue")

In [2]: type(tup)
Out[2]: <type 'tuple'>
```

You can slice and index the tuple exactly like you would a list. Tuples are
used in the inner workings of python, and a tuple can be used as a key in a
dictionary, whereas a list cannot as we will see in a moment.

## Sets

The python set type is similar to the idea of a mathematical set: it is an unordered collection of
unique things. Consider:

```python
In [3]: fruit = set(["apple", "banana", "pear", "banana"]) #You have to use a list to create a set.
```

Since sets contain only unique items, there's only one banana in the set fruit.
You can do things like intersections, unions, etc. on sets just like in maths.

```python
In [1]: firstBowl = set(["apple", "banana", "pear", "peach"])	

In [2]: secondBowl = set(["peach", "watermelon", "orange", "apple"])

In [3]: set.intersection(firstBowl, secondBowl)
Out[3]: set(['apple', 'peach'])

In [4]: firstBowl & secondBowl
Out[4]: set(['apple', 'peach'])

In [5]: firstBowl | secondBowl
Out[5]: set(['apple', 'peach', 'pear', 'watermelon', 'orange', 'banana'])

In [6]: firstBowl ^ secondBowl
Out[6]: set(['orange', 'pear', 'watermelon', 'banana'])
```

Dictionaries
------------

A dictionary is basically an efficient table that maps keys to values. It is an unordered container.


```python
In [7]: dataDict = {"experiment": "current vs. voltage", \
                   "run": 47, \
                   "temperature": 372.756, \
                   "current": [-1.0, -0.5, 0.0, 0.5, 1.0], \
                   "voltage": [-2.0, -1.0, 0.0, 1.0, 2.0]}
```

This model is clearly better than just using lists because you no longer have to remember that
the run number is in the second position of the list, you just refer directly to "run":

```python
In [9]: dataDict["run"]
Out[9]: 47
```

If you wanted the voltage data list:

```python
In [10]: dataDict["voltage"]
Out[10]: [-2.0, -1.0, 0.0, 1.0, 2.0]
```

Or perhaps you wanted the last element of the current data list

```python
In [11]: dataDict["current"][-1]
Out[11]: 1.0
```

Once a dictionary has been created, you can change the values of the data
if you like.

```python
In [12]: dataDict["temperature"] = 3275.39
```

You can also add new keys to the dictionary.

```python
In [13]: dataDict["user"] = "Johann G. von Ulm"
```

Dictionaries, like strings, lists, and all the rest, have built-in methods.
Let's say you wanted all the keys from a particular dictionary.

```python
In [14]: dataDict.keys()
Out[14]: ['run', 'temperature', 'current', 'experiment', 'user', 'voltage']
```

also, values

```python
In [15]: dataDict.values()
Out[15]: 
[47,
 3275.39,
 [-1.0, -0.5, 0.0, 0.5, 1.0],
 'current vs. voltage',
 'Johann G. von Ulm',
 [-2.0, -1.0, 0.0, 1.0, 2.0]]
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

### Exercise

- Make a dictionary and experiment using different types as keys. Can containers be keys?

Since tuples are immutable, they can be used as keys for dictionaries.
Lists are mutable, and therefore cannot.


# Flow Control


Conditionals
============

A conditional (if statement) is some statement that in general says :
"When some boolean is true, do the following. Elsewise, do this other
thing."

Many equivalence test statements exist in Python that are similar in
other languages:

```python
  i=1
  j=2
  i==j # i is equal to j : FALSE
  i<j  # i is less than j
  i<=j # i is less than or equal to j : TRUE
  i>j  # i is greater than j
  i>=j # i is greater than or equal to j : FALSE
  i!=j # i is not equal to j : TRUE
```

However, python has other equivalence test statements that are fairly
unique to python. To check whether an object is contained in a list :

```python
  beatle = "John"
  beatles = ["George", "Ringo", "John", "Paul"]
  print beatle in beatles # is John one of the beatles? : TRUE
  print "Katy" not in beatles # this is also TRUE. 
```

Conditionals (if statements) are also really easy to use in python. Take
a look at the following example:

```python
  i = 4
  sign = "zero"
  if i < 0:
    sign = "negative"
  elif i > 0:
    sign = "positive"
  else:
    print "Sign must be zero"
    print "Have a nice day"
  print sign
```

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code. For example, in MatLab you begin an if
statement with the word "if" and you end it with "end if". In C/C++ you
delimit blocks with curly braces. Python uses ''indentation'' to delimit
code blocks. The **indentation** above is NOT just to make things look
pretty - it tells Python what the body of the if-statement is. This is
true when ever we create any code blocks, such as the bodies of loops,
functions or classes.

**Aside: Compact if-statement:**

Python has an easy to use if-syntax for setting the value of a variable.
Try entering this into IPython:

```python
  i = 5
  sign = "positive" if i > 0 else "negative"
```

While Loops
===========

Let's start by looking at while loops since they function like while
loops in many other languages. The example below takes a list of
integers and computes the product of each number in the list up to the
-1 element.

A while loop will repeat the instructions within itself until the
conditional that defines it is no longer true.

```python
  mult = 1
  sequence = [1, 5, 7, 9, 3, -1, 5, 3]
  while sequence[0] is not -1:
      mult = mult * sequence[0]
      del sequence[0]

  print mult
```

Some new syntax has been introduced in this example.

-   On line 3 We begin the while loop. Notice that instead of using the
    not-equals symbol, !=, we can use "is not". There is a subtle semantic
    difference: 'is' tests object identity, whereas '=' tests value equality.

-   On line 4, we compute the product of the elements just to make this
    more interesting.

-   On line 5, we use the \`del\` keyword to remove the first element of
    the list, shifting every element down one.

**Watch Out**

Since a while loop will continue until its conditional is no longer
true, a **poorly formed** while loop might repeat forever. For example :

```python
  i=1
  print "Well, there's egg and bacon, egg and spam, egg bacon and"
  while i == 1:
    print "spam "
  print "or Lobster Thermidor a Crevette with a mornay sauce served in a Provencale manner with shallots..." 
```

Since the variable **i** never changes within the while loop, we can
expect that the conditional, **i=1** will remain true forever and the
while loop will just go round and round, as if this restaurant offered
nothing but spam. (If you try this at home, please note that one way to
interrupt a non-terminating process is **ctrl+c** or **ctrl+z**.)
An 'infinite' loop (`while True`) is, however, quite useful if you have a 'break'
statement to jump out of it at some point.

To create nested loops, the indentation (preferably four spaces; never tabs) should increase for each looping level.

For Loops
=========

For loops in python operate a little differently from other languages.
Let's start with a simple example which prints all of the numbers from 0
to 9:

```python

  for i in range(10):
      print i
```

You may be wondering how this works. Start by using help(range) to see
what the range function does.

    Help on built-in function range in module __builtin__:

    range(...)
        range([start,] stop[, step]) -> list of integers

        Return a list containing an arithmetic progression of integers.
        range(i, j) returns [i, i+1, i+2, ..., j-1]; start (!) defaults to 0.
        When step is given, it specifies the increment (or decrement).
        For example, range(4) returns [0, 1, 2, 3].  The end point is omitted!
        These are exactly the valid indices for a list of 4 elements.

Range is a function that returns a list containing a sequence of
integers. So, range(10) returns the list [0,1,2,3,4,5,6,7,8,9]. The for
loop then simply iterates over that list, setting i to each value.

### Exercises

1. Using a loop, calculate the factorial of 42 (the product of all integers up to and including 42).

For Loops with Lists and Dictionaries
=====================================

With range, we learned that **for** loops in python are really used to
iterate over sequences of things (they can be used for much more, but
for now this definition will do). Try entering the following to see what
happens:

```python

  for c in ["one", 2, "three", 4, "five"]:
      print c
```

this is equivalent to:

```python

  c = ["one", 2, "three", 4, "five"]
  for i in range(len(c)):
      print c[i]
```

**Warning**: it's generally not safe to modify the sequence you are iterating over!

A common task is to iterate over a sequence while keeping track of the item number.
Python provides the `enumerate` builtin for this:
```python
    words = ('cool', 'powerful', 'readable')
    for index, item in enumerate(words):
        print index, item
```

Use zip to iterate over two lists at once

```python
    fruits = ['apples', 'oranges', 'pears', 'bananas']
    prices = [0.49, 0.99, 1.49, 0.32]
    for fruit, price in zip(fruits, prices):
        print fruit, "cost", price, "each"
```    

With a list, then, it's clear that we can use the **in** keyword to
indicate a list of things. What about nested loops around a list of
lists?

```python
  italy = ["Rome", "Pisa", "Florence", "Venice", "Trieste"]
  argentina = ["Mendoza", "Buenos Aires", "Patagonia"]
  india = ["Ahmedabad","Kolkata", "Chennai", "Jaipur", "Surat"]
  us = ["Chicago", "Austin", "New York", "San Fran"]
  nations = [italy, argentina, india, us]
  nationnames = ["italy","argentina", "india", "us"]
  for nation in nations :
      print nationnames[nations.index(nation)] + ": "
      for city in nation :
          print "  " + city 
```

Of course, this information is better stored in a dictionary, isn't it?
The data makes more sense if the keys were the nation names and the
values were lists of cities. Importantly, python has given us a tool
specifically for dictionary looping.

The syntax for looping through the keys and values of a dictionary is :

    for key, value in dictionary.iteritems():

Importantly, you don't have to use the words key and value. That's just
what will fill those variables. Here, we rewrite the previous loop using
this clever syntax.

```python
  italy = ["Rome", "Pisa", "Florence", "Venice", "Trieste"]
  argentina = ["Mendoza", "Buenos Aires", "Patagonia"]
  india = ["Ahmedabad","Kolkata", "Chennai", "Jaipur", "Surat"]
  us = ["Chicago", "Austin", "New York", "San Fran"]
  nations = {"italy":italy, "argentina":argentina, "india":india, "us":us}
  for nation, cities in nations.iteritems() :
      print nation + " : "
      for city in cities :
          print "  " + city 
```

Note that the order is non-deterministic.  If you care, loop over sorted keys:
```python
    for nation in sorted(nations.keys()):
        print nation
```

break, continue, and else
=========================

A break statement cuts off a loop from within an inner loop. It helps
avoid infinite loops by cutting off loops when they're clearly going
nowhere.

```python
  reasonable = 10
  for n in range(1,2000):
      if n == reasonable :
          break
      print n
```

Something you might want to do instead of breaking is to continue to the
next iteration of a loop, giving up on the current one..

```python
  reasonable = 10
  for n in range(1,2000):
      if n == reasonable :
        continue
      print n
```

What is the difference between the output of these two?

Importantly, Python allows you to use an else statement in a for loop.

That is :

```python
  knights={"Sir Belvedere":"the Wise", "Sir Lancelot":"the Brave", \
          "Sir Galahad":"the Pure", "Sir Robin":"the Brave", "The Black Knight":"John Clease"} 

  favorites=knights.keys()
  favorites.remove("Sir Robin")
  for name, title in knights.iteritems() : 
      string = name + ", "
      for fav in favorites :
          if fav == name :
              string += title
              break
      else:
          string += title + ", but not quite so brave as Sir Lancelot." 
      print string
```

Larger Example
==============

We've seen a lot so far. Let's work through a slightly lengthier example
together. I'll use some of the concepts we already saw and introduce a
few new concepts. To run the example, you'll need to locate a short file
containing phone numbers, called phonenums.txt.
Now we have to move ipython to that directory so it can find the
phonenums.txt file. You navigate within ipython in the same way that you
navigate in the shell, by entering "%cd [path]" .

This example opens a text file containing a list of phone numbers. The
phone numbers are in the format \#\#\#-\#\#\#-\#\#\#\#, one to a line.
The example code loops through each line in the file and counts the
number of times each area code appears. The answer is stored in a
dictionary, where the area code is the key and the number of times it
occurs is the value.

```python

  areacodes = {} # Create an empty dictionary
  f = open("phonenums.txt") # Open the text file
  for line in f: # iterate through the text file, one line at a time (think of the file as a list of lines)
      ac = line.split('-')[0] # Split phone number, first element is the area code
      if not ac in areacodes: # Check if it is already in the dictionary
        areacodes[ac] = 1 # If not, add it to the dictionary
      else:
        areacodes[ac] += 1 # Add one to the dictionary entry
  
  print areacodes # Print the answer
```

Exercise : Iteritems
--------------------

Use the iteritems dictionary method in combination with a for loop to
print the keys/values of the areacodes dictionary one to a line. In
other words, the goal is to write a loop that prints:

    203 4
    800 4
    608 8
    773 3

## List Comprehensions

A very powerful feature is the ability to define lists (and dictionaries and sets in relatively recent versions) by **comprehensions**:

```python
[i**2 for i in range(4)]
```

Functions
=========

```python
def square(x):
    return x * x
print square(2), square(square(2))

def hello(time, name):
    """Print a nice message. Time and name should both be strings.
    Example: hello('morning', 'Software Carpentry')
    """
    print 'Good ' + time + ', ' + name + '!'
hello('afternoon', 'Software Carpentry')
```

**Note:** by default, functions return `None`.

Above we saw mandatory parameters (positional arguments).  
You can also have optional parameters (keyword or named arguments):

```python
def double_it(x=2):
    return x * 2
double_it()
double_it(3)
```

**Warning:**
Default values are evaluated when the function is defined, not when it is called. This can be problematic (or desirable) when using mutable types (e.g. dictionary or list) and modifying them in the function body, since the modifications will be persistent across invocations of the function.

```python
def counter(l=[1]):
    l[0] += 1
    return l[0]
counter()
counter()

my_l = [3]
counter(my_l)
print my_l
```

The latter example illustrates how Python passes arguments by reference to the value. If the value is **immutable**, the function does not modify the caller's variable. If the value is **mutable**, the function may modify the caller's variable in-place.

A more involved example of keyword arguments:
```python
def slicer(seq, start=None, stop=None, step=None):
    """Implement basic python slicing."""
    return seq[start:stop:step]

rhyme = 'one fish, two fish, red fish, blue fish'.split()
['one', 'fish,', 'two', 'fish,', 'red', 'fish,', 'blue', 'fish']
slicer(rhyme)
slicer(rhyme, step=2)
slicer(rhyme, 1, step=2)
slicer(rhyme, step=2, start=1, stop=4)
```

### Exercises

1. Create a function that returns the factorial of a given number


## Variable number of parameters

Special forms of parameters:
-  `*args`: any number of positional arguments packed into a tuple
-  `**kwargs`: any number of keyword arguments packed into a dictionary

```python
def variable_args(*args, **kwargs):
    print 'args is', args
    print 'kwargs is', kwargs

variable_args('one', 'two', x=1, y=2, z=3)
```

## Functions are objects

Functions are first-class objects, which means they can be:
- assigned to a variable
- an item in a list (or any collection)
- passed as an argument to another function (**higher-order** functions).

```python
va = variable_args
va('three', x=1, y=2)
```

# Reading Files

```python
f = open('animals.txt', 'r')
# dir(f)
# help(f)
f.seek(0)
for line in f:
    print line

'2011-06-22 09:23 Moose 2'.split() # numbers come out as strings
print int('2'), float('2'), float('3.14159')

f.seek(0)
date = []
time = []
animal = []
number = []
for line in f:
    d, t, a, n = line.split()
    date.append(d)
    time.append(t)
    animal.append(a)
    number.append(int(n))

date
time
number
```

### Exercise
1. Read the file 'big_animals.txt' and print each line on which more than 10 moose were sighted.
1. Turn the code for #1 into a function and use it on the files 'merida_animals.txt' and 'fergus_animals.txt'.

# Modules

Python has a lot of useful data type and functions built into the language, some of which you have already seen. For a full list, you can type `dir(__builtins__)`. However, there are even more functions stored in modules. An example is the sine function, which is stored in the math module. In order to access mathematical functions, like sin, we need to `import` the math module. Let's take a look at a simple example:

```python
print sin(3) # Error! Python doesn't know what sin is...yet
import math # Import the math module
math.sin(3)
print dir(math) # See a list of everything in the math module
help(math) # Get help information for the math module
```

It is not very difficult to use modules - you just have to know the module name and import it. There are a few variations on the import statement that can be used to make your life easier. Let's take a look at an example:

```python
from math import *  # import everything from math into the global namespace (A BAD IDEA IN GENERAL)
print sin(3)        # notice that we don't need to type math.sin anymore
print tan(3)        # the tangent function was also in math, so we can use that too

%reset # Clear everything from IPython

from math import sin  # Import just sin from the math module. This is a better idea.
print sin(3)          # We can use sin because we just imported it
print tan(3)          # Error: We only imported sin - not tan

%reset                # Clear everything

import math as m      # Same as import math, except we are renaming the module m
print m.sin(3)        # This is really handy if you have module names that are long
```


# Exception handling in Python

It is highly unlikely that you haven't yet raised Exceptions if you have typed all the previous commands of the tutorial. For example, you may have raised an exception if you entered a command with a typo.

Exceptions are raised by different kinds of errors arising when executing Python code. In your own code, you may also catch errors, or define custom error types. You may want to look at the descriptions of the the built-in Exceptions when looking for the right exception type (`dir(__builtins__)`).

## Exceptions

Exceptions are raised by errors in Python:
```python
1/0
1 + 'e'
d = {1:1, 2:2}
d[3]
l = [1, 2, 3]
l[4]
l.foobar
```

As you can see, there are different types of exceptions for different errors.

## Catching exceptions

```python
while True:
    try:
        x = int(raw_input('Please enter a number: '))
        break
    except ValueError:
        print('That was no valid number.  Try again...')
print x
```

```python
try:
    x = int(raw_input('Please enter a number: '))
finally:
    print('Thank you for your input')
```

This version is important for resource management (e.g. closing a file).
(Although see also context managers in Python 2.5+)

## Easier to ask for forgiveness than for permission

```python
def print_sorted(collection):
    try:
        collection.sort()
    except AttributeError:
        pass
    print(collection)

print_sorted([1, 3, 2])
print_sorted(set((1, 3, 2)))
print_sorted('132')
```

## Raising exceptions

Capturing and reraising an exception:

```python
def filter_name(name):
    try:
        name = name.encode('ascii')
    except UnicodeError, e:
        if name == 'Gaël':
            print('OK, Gaël')
        else:
            raise e
    return name

filter_name('Gaël')
filter_name('Stéfan')
```

## Exceptions to pass messages between parts of the code:

```python
def achilles_arrow(x):
    if abs(x - 1) < 1e-3:
        raise StopIteration
    x = 1 - (1-x)/2.
    return x

x = 0
while True:
    try:
        x = achilles_arrow(x)
    except StopIteration:
        break

print x
```

# Object Orientation

```python
class MyClass(object):
    def __init__(self, a):
        self.a = a

    def double(self):
        self.a *= 2

my_obj = MyClass(1)
print type(my_obj)
print my_obj.a
my_obj.a = 2
my_obj.double()
print my_obj.a
```

**Aside: Magic functions**

Methods with leading and trailing double underscores are "magic
functions" in python.

 - Iteration (for x in sequence) uses \_\_next\_\_
 - Slicing ie brackets) (a[1:2]) uses \_\_get\_\_
 - Calling ie parentheses (a(3)) uses \_\_call\_\_
 - Help uses \_\_doc\_\_


# Extra topics

- http://scipy-lectures.github.com/advanced/advanced_python/index.html
    - Iterators, generator expressions and generators
    - Decorators
    - Context managers
- http://scipy-lectures.github.com/advanced/debugging/index.html#using-the-python-debugger
    - Using the Python debugger
- https://github.com/jonc125/boot-camps/tree/2013-05-oxford-dtc/python/debugging#profiling-making-code-fast
    - Profiling
