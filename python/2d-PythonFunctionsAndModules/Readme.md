# Methods and Modules

[Back To Flow
Control](../2c-PythonFlowControl)

* * * * *

**Presented By : Shoaib Sufi**

**Original material by Milad Fatenejad (mostly)**

At this point, we've been over the basics of the python language except for the parts that allow you to reuse code you and others have written. First we will talk about wirting methods (sometimes called functions), then we'll talk about modules. A very important topic that we unfortunately won't have time to cover is classes (aka, object orientation). Dive Into Python has an excellent chapter on this topic that you should read.

## Methods

Lets say you frequently want to calculate a list of values of power given a list of voltages and a list of currents. You could continually drop the for loop we wrote before into your code, but this approach quickly becomes unwieldy and error prone. The solution is to abstract this functionality into a method. Here's how.

```python
def calcpower(voltageList, currentList):
  powerList = [] # Initialize an empty list to be filled later.
  indxList = range(len(currentList))

  for indx in indxList:
    power = voltageList[indx] * currentList[indx]
    powerList.append(power)

  return powerList
```

That's it: just a block containing the code we already wrote plus one extra line. The return statement at the end tells python that when the method is done executing, hand back the value of powerList.

You can call this method at any point in the following way:

```python
calcpower(voltageList, currentList)
```

In fact, we've been using methods almost this entire time:

```python
len()
type()
"string!".upper()
```

and so on.

We just looked at one way to construct a method; there are more. Lets say we wanted to write the ugliest method ever that calculates and returns a value raised to a power. (In fact, this method already exists, so writing a new one is a terrible idea). Here's how it would look:

```python
def topower(base,expt):
  product = 1
  for indx in range(expt):
    product = product * base

  return product
```

Seriously though, don't ever use the above code.

###Default argument value

What if most of the time we were just squaring a number? We'd want a default value of 2 for the exponent.

```python
def topower(base,expt = 2):
  product = 1
  for indx in range(expt):
    product = product * base

  return product
```

In this case, we are setting the default value to 2. We can change that value when we call the command, but we don't always need to set it to 2.

```python
topower(2) # works just fine and returns the square of 2
```

### returning multiple values

Python functions can return multiple values:

```python
def rgb(color_name):

    if color_name is "red":
       return 1,0,0
    if color_name is "green":
       return 1,0,0
    if color_name is "blue":
       return 1,0,0
    if color_name is "yellow":
       return 1,1,0

    print "ERROR: Unknown color"

 r,g,b = rgb("yellow")
 print r,g,b
```

###Exercise: Perimeter ###

Write a function, "\`perimeter\`" that takes one argument - a tuple
containing the lengths of each side of a polygon. Have the function
return the perimeter of the polygon. So, for example to find the area of
a square with side length 3, the function call would be:
perimeter((3,3,3,3)) and the function would return 12. Use the \`help\`
and \`dir\` functions to figure out how to write this function in one
line.

###Variable Number of Arguments ###

In the last example, you wrote a function, called \`perimeter\`, that
you have to call like this:

```python
print perimeter((1,2,3,4,5))
```

You need the extra set of parentheses, because the function takes only
one argument - a single tuple. Wouldn't it be nice if we didn't have to
have the extra parentheses? It turns out that this is easy in python:

```python
def perimeter(*args):
   return sum(args)

print perimeter(1,2,3,4,5)
```

Notice that little star on the first line. That star tells python to
take all of the arguments and bundle them into a single tuple called
args. This feature allows us to create functions that take a variable
number of arguments. The function calls do not require the extra set of
parentheses we saw earlier.

**Note**:  when mixing default value arguments and variable numbers of arguments
in the same function the variable number of arguments need to come last.

```python
def f(a, b=2, *args): # will work
    print a
```

###Exercise: sum floats ###

Write a function that takes a variable number of floating point
arguments and returns their average.

### Keyword Arguments ###

Consider the following function, whose arguments represent the model
year, mileage, and number of accidents that a car has had. The function
attempts to use this information to compute the value of the car.

```python
def carvalue(year, mileage, nacc):
   return 30000 - mileage/10 - nacc*1000 - (2010-year)*1000

print carvalue(2001, 100000, 2)
```

In order to use the \`carvalue\` function, we have to remember that
\`year\` is the first argument, \`mileage\` is the second, and \`nacc\`
is the third. If we accidentally put the wrong argument in the wrong
place then we will compute the wrong answer. Luckily, we can be explicit
when calling functions using the following syntax:

```python
print carvalue(year=2001, mileage = 100000, nacc=2)
print carvalue(mileage= 100000, year = 2001, nacc=2) # Also produces the correct answer!
print carvalue(2001, nacc= 2, mileage = 100000) # Sets the year to 2001, the mileage to 100000, and nacc to 2
print carvalue(year= 2001, mileage = 100000, 2) # ERROR: Keyword arguments must precede non-keyword arguments
```

## Modules

Modules are a way to package code for reuse at a later time. You can buld your own modules, install others' modules, or use built-in modules provided from python.

Lets do the exponents right. It turns out that python ships with a math module containing exactly that functionality. Here's how you bring it in to your environment.

```python
import math
```

Now we can do exponents, say 3^2

```python
math.pow(3,2)
```

How about 2^4?

```python
math.pow(2,4)
```

Writing a method to create this functionality was a terrible idea because the math module already exists. The math module has other methods like sine, cosine, exp, and so on.

There are other ways to pull from modules into your environment. If you only wanted the sine method, you could import it like so

```python
from math import sin
```

and you'd be able to directly use it without a dot.

```python
sin(3)
```

If you want to pull in everything from the math module, but you're too lazy to type out "math" all the time, you can give the module an alias

```python
import math as mth
```

or whatever you want the alias to be.

One final way to import stuff is to use the asterisk notation. This approach is almost always a terrible idea, so you probably shouldn't do it.

```python
from math import *
```

It is a terrible idea because you dump all of the names of data and methods into your local namespace. If you already have data or methods of the same name, you can end up facing some nasty problems as they would have been redefined to those in the math module.

# Using the sys module

The sys module gives you access to the Python interpreter. Some
important objects in this module are: 1) sys.argv is a list of command
line arguments. The first argument is always the name of the file. 2)
sys.path gives a list of all paths on your computer where Python will
look for modules 2) sys.modules gives a dictionary of all currently
loaded modules.

Put the following text in a file output-args.py in your current working
directory.

```python
from sys import argv

for i in xrange(1,len(argv)):
    print argv[i]
```

Then run this from the command line using

    python output-args.py 1 2 3 4

This program should print each input on its own line.

Try changing the first argument of xrange to 0 and re-run to see what happens.

# Fin

So that wraps up the python lessons. Please take a look at the resources mentioned at the top of the first lesson for some really excellent information.

## Exercise

### data.dat to dict ###

Write a program in a file which consists of a function that takes a filename (e.g. data.dat) that will pull the data from the file into a dictionary that looks like the dict we made in lesson 2. This task seems trivial, but its a little harder than it looks. Supply the filename via the command line.

    python data2dict.py data.dat
