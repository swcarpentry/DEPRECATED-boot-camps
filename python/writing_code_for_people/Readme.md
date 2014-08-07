[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)


# Write Code for People: Variables, Data Structures and Conditionals
  

* * * * *


**Based on Lecture Materials By: Milad Fatenejad, Katy Huff, Tommy Guy, Joshua 
R. Smith, Will Trimble, and many more**

This lecture is on basic programming in python. In order to do the examples, we are going to use an environment called iPython.  I expect this lecture to be interactive, so stop me at any point if you have questions. 

This lecture will be structured as follows: I will be teaching the basics of two things: the python programming language (to a greater extent) and the iPython interpreter (to a lesser extent). The iPython interpreter is one of many different ways to implement python code. As far as the python component, I'll shoot for a layered approach: I'll continue building on my previous concepts. It turns out that like any sufficiently complex topic, its not really possible to force the pedagogy into a serial stream. Also, we have a pretty serious time constraint. I'm just going to drop it on you. Because of the brief nature of this tutorial, I've included links to some excellent reference material. Also, if we have time, I'll take questions based on the specific programming needs of this class.

Here is the reference material.

* [Dive into Python](http://www.diveintopython.net/toc/index.html)
* [Software Carpentry's Python Lectures](http://software-carpentry.org/4_0/python/)
* [IPython: A System for Interactive Scientific Computing](http://dx.doi.org/10.1109/MCSE.2007.53)
* [How to Think Like a Computer Scientist](http://www.greenteapress.com/thinkpython/thinkpython.html)

Once we briefly deal with iPython, I'll cover python in the following order:

## What We'll Cover

We'll focus on two overarching concepts that are important to any programming language: 

1) How to write code for people. That is code that is readable and understandable to others in your group and most importantly to your future self 3 months or 3 years down the road.

2) How to not repeat yourself. How to reuse your code with loops and functions. And how to eventually build modules, collections of functions, you and others can use in all your codes. And how to use other people's modules.

### iPython Intro
* iPython

### Lesson 1 (Write Code for People)
* Print statements
* Variables
* Integers
* Floats
* Strings
* Types
* Type Coercion
* Basic Operations: add numbers, concatenate strings, basic data type functionality
* List
* Dictionary 
* Tuple
* Conditional (if) statements

### Lesson 2 (Don't Repeat Yourself)
* For Loops
* While Loops
* Iteration
* Methods
* Reading & Writing Files
* Modules


## iPython

Please follow the link to the [iPython Intro](../ipython/Readme.md).

## Back to Write Code for People

This lesson will introduce the basics of the python programming language while stressing how to make readable code, code for people. As we introduce variables we will discuss how to name them, when to comment and which comments are useful. As we introduce types we will discuss how type choice can be influenced by code readability considerations. The above readability considerations are to a certain degree universal to all programming, regardless of language. Finally, we introduce conditionals, such as `if` statements; we will see what an important role white space plays in python and how it can improve readability.

But first: let's introduce you to python.

## Variables


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

## Naming Your Variables

Which of these variable names are good?

```
v = 4
voltage = 4
TheVoltageInTheCircuitAtPointA = 4 # Camel Case
the_voltage_is = 4 # Pothole
monkey = 4
```
Variable names should be:

* Meaningful (to those who are going to read the code)
* Short enough so you don't misstype them

Variable name choice is important; a well-named variable is self-explanatory without comments and will make your code easier to read as the reader will not have to look up the comments. Remember that context matters a lot: in some cases, you'll want to spell out `voltage` or even `input_voltage`. In other cases, `v` is a shorthand that everyone will understand.

It's also important to choose a naming convention for your whole project, and get the people working on it to agree. Mixing `batteryVoltage` and `capacitor_value` in one place makes code hard to read.

## Writing Comments For People

Above, you may have noticed the `#` character, which denotes a comment in python. Comments should describe meaning but not what the statement is doing.

```
In [15]: voltage = 4 # set the voltage to 4   <- well, duh

In [16]: voltage = 4 # Input to the circuit   <- better; says what this voltage means
```

## Types and Dynamic Typing

Like in most programming languages, things in python are *typed* — the *type* refers to the type of data and what you can do with it. For example, you can do different things with strings and numbers. Numbers can have decimal components or not, and so on. You can see the type of a variable with the `type` command.

```
In [5]: type?
Type:       type
String Form:<type 'type'>
Namespace:  Python builtin
Docstring:
type(object) -> the object's type
type(name, bases, dict) -> a new type

In [6]: type(full_name)
Out[6]: str
In [7]: type(10)
Out[7]: int

```


Python is what is known as a *dynamically typed* language. Dynamic typing means that you don't have to declare the type of a variable when you define it; python just figures it out based on how you are setting the variable. This is in contrast to *statically typed* languages, where you must say up front that a variable is going to be used for strings or numbers or whatever. There are good and bad points to both approaches.

But back to Python. Let's say you set a variable. Sometime later you can just change the type of data assigned to a variable and python is perfectly happy about that. Since it won't be obvious until (possibly much) later why that's important, I'll let you marinate on that idea for a second.

Here's an example of dynamic typing:

```
In [8]: voltage = 2
In [9]: print(voltage)
2

In [10]: type(voltage)
Out[10]: int
```

It's an `int`, which is an integer — a number with no decimal component. Let's assign a value of 2.7 (which has a decimal part) to voltage. What happens to the type?

```
In [11]: voltage = 2.7

In [12]: type(voltage)
Out[12]: float
```

Neat! That's a `float`, which does have a decimal part. You can assign a string to the variable voltage in the same way:

```python
In [13]: voltage = "2.7 volts"

In [14]: type(voltage)
Out[14]: str
```

I'll let you ruminate on the pros and cons of this construction while I change the value of voltage back to an int:

```
In [15]: voltage = 2
```

Choosing an appropriate variable type is not just a practical concern; it can also have an effect on code readability. Is this number used for calculations or only in print statements?


## On Being Precise With floats and ints

Again, the following may seem esoteric and pedantic, but it is very important. So bear with me.

Lets say you had some voltage data that looks like the following

```
0
0.5
1
1.5
2
```

If you just assigned this data individually to a variable, you'd end up with the following types:

```
0   -> int
0.5 -> float
1   -> int
1.5 -> float
2   -> int
```

But what if you wanted all of that data to be floats on its way in? You could assign the variable and then coerce it to type float:

```
In [28]: voltage = float(1)
```

But that's ugly. If you want whats otherwise an integer to be a float, add `.0` to the end:

```python
In [29]: voltage = 1.0

In [30]: type(voltage)
Out[30]: float
```

This point becomes important when we start operating on data in the next section.

## Data Operations

In this section all of the discussion in the previous section becomes important. I don't know if I'd call this stuff fundamental to the language, but it's pretty important and it will zing you if you aren't careful. The takeaway is that you need to be precise with what you are doing. Let's say you want to add some integers.

```
In [31]: a = 1

In [32]: b = 2

In [33]: c = a + b

In [34]: c
Out[34]: 3

In [38]: type(a), type(b), type(c)
Out[38]: (int, int, int)
```

So we got a value of three for the sum, which also happens to be an integer. Any operation between two integers is another integer. Makes sense.

So what about the case where a is an integer and b is a float?

```
In [39]: a = 1

In [40]: b = 2.0

In [41]: c = a + b

In [42]: c
Out[42]: 3.0

In [43]: type(a), type(b), type(c)
Out[43]: (int, float, float)
```

You can do multiplication on numbers as well.

```
In [44]: a = 2

In [45]: b = 3

In [46]: c = a * b

In [47]: c
Out[47]: 6

In [48]: type(a), type(b), type(c)
Out[48]: (int, int, int)
```

Also division.

```
In [49]: a = 1

In [50]: b = 2

In [51]: c = a / b

In [52]: c
Out[52]: 0
```

**ZING!**

Here's why type is important. Dividing two integers returns an integer: this operation calculates the quotient and floors the result to get the answer.

If everything was a float, the division is what you would expect.

```
In [53]: a = 1.0

In [54]: b = 2.0

In [55]: c = a / b

In [56]: c
Out[56]: 0.5

In [57]: type(a), type(b), type(c)
Out[57]: (float, float, float)
```

## Compound Data Types: Lists, Dictionaries, Sets, and Tuples

Python would be a fairly useless language if it weren't for the compound
data types. The main two are lists and dictionaries, but I'll mention sets
and tuples as well. 

## Lists

A list is an ordered, indexable collection of data. Let's say you have
collected some current and voltage data that looks like this:

```
voltage:
-2.0
-1.0
0.0
1.0
2.0

current:
-1.0
-0.5
0.0
0.5
1.0
```

So you could put that data into lists like

```
In [1]: voltage_list = [-2.0, -1.0, 0.0, 1.0, 2.0]

In [2]: current_list = [-1.0, -0.5, 0.0, 0.5, 1.0]
```

obviously voltage_list is of type list:

```
In [3]: type(voltage_list)
Out[3]: list
```

Python lists have the charming (annoying?) feature that they are indexed
from zero. Therefore, to find the value of the first item in voltage_list:

```
In [4]: voltage_list[0]
Out[4]: -2.0
```

And to find the value of the third item

```
In [5]: voltage_list[2]
Out[5]: 0.0
```

Lists can be indexed from the back using a negative index. The last item of
current_list

```
In [6]: current_list[-1]
Out[6]: 1.0
```

and the next-to-last

```
In [7]: current_list[-2]
Out[7]: 0.5
```

You can "slice" items from within a list. Let's say we wanted the second
through fourth items from voltage_list

```
In [8]: voltage_list[1:4]
Out[8]: [-1.0, 0.0, 1.0]
```

Or from the third item to the end

```
In [9]: voltage_list[2:]
Out[9]: [0.0, 1.0, 2.0]
```

and so on.

***Exercise***
What does voltage_list[::2] mean?

### Append and Extend

Just like strings have methods, lists do too.

```
In [10] dir(list)
```

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

```
voltage:
3.0
4.0

current:
1.5
2.0
```

If you want to append items to the end of a list, use the append method.

```
In [11]: voltage_list.append(3.)

In [12]: voltage_list.append(4.)

In [13]: voltage_list
Out[13]: [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
```

You can see how that approach might be tedious in certain cases. If you
want to concatenate a list onto the end of another one, use extend.

```
In [14]: current_list.extend([1.5, 2.0])

In [15]: current_list
Out[15]: [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
```

### Length of Lists

Sometimes you want to know how many items are in a list. Use the len command.

```
In [16]: len(voltage_list)
Out[16]: 7
```

### Heterogeneous Data

Lists can contain hetergeneous data.

```
In [17]: data_list = ["experiment: current vs. voltage", \
                      "run", 47, \
                      "temperature", 372.756, \
                      "current", [-1.0, -0.5, 0.0, 0.5, 1.0], \
                      "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]

```

We've got strings, ints, floats, and even other lists in there. The slashes
are there so we can continue on the next line. They aren't necessary but
they can sometimes make things look better.

## Assigning Variables to Other Variables

Something that might cause you headaches in the future is how python deals
with assignment of one variable to another. When you set a variable equal
to another, both variables point to the same thing. Changing the first one
ends up changing the second. Be careful about this fact.

```
In [19]: a = [1,2]

In [20]: b = a

In [21]: a.append(10)

In [22]: b
Out[22]: [1, 2, 10]
```

There's a ton more to know about lists, but let's press on. Check out Dive
Into Python or the help documentation for more info.

## Tuples

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is immutable; once you set the
data in it, the tuple cannot be changed. You define a tuple as follows.

```
In [1]: tup = ("red", "white", "blue")

In [2]: type(tup)
Out[2]: tuple
```

You can slice and index the tuple exactly like you would a list. 

### Why Use Tuples

Tuples can emphasize intent, although you could make a list and not change it, setting the type to be a tuple makes that intent clear to a reader. 
Tuples are used in the inner workings of python, and a tuple can be used as a key in a
dictionary, whereas a list cannot as we will see in a moment.

***Exercise***
Display the second element of the tuple with two different slices.

## Dictionaries

Recall our variable data_list which contained our current-voltage data
and also some metadata. We were able to store the data as a list, but
clearly the list type is not the optimal choice for a data model. The
dictionary is a much better choice.

A python dictionary is a collection of key, value pairs. The key is a
way to name the data, and the value is the data itself.  Here's a way
to create a dictionary that contains all the data in our data.dat file
in a more sensible way than a list.

```
In [7] data_dict = {"experiment": "current vs. voltage", \
                   "run": 47, \
                   "temperature": 372.756, \
                   "current": [-1.0, -0.5, 0.0, 0.5, 1.0], \
                   "voltage": [-2.0, -1.0, 0.0, 1.0, 2.0]}
```

This model is clearly better because you no longer have to remember that
the run number is in the second position of the list, you just refer
directly to "run":

```
In [9]: data_dict["run"]
Out[9]: 47
```

If you wanted the voltage data list:

```
In [10]: data_dict["voltage"]
Out[10]: [-2.0, -1.0, 0.0, 1.0, 2.0]
```

Or perhaps you wanted the last element of the current data list

```
In [11]: data_dict["current"][-1]
Out[11]: 1.0
```

Once a dictionary has been created, you can change the values of the data
if you like.

```
In [12]: data_dict["temperature"] = 3275.39
```

You can also add new keys to the dictionary.

```
In [13]: data_dict["user"] = "F. C. Rodgers"
```

Dictionaries, like strings, lists, and all the rest, have built-in methods.
Let's say you wanted all the keys from a particular dictionary.

```
In [14]: data_dict.keys()
Out[14]: ['run', 'temperature', 'current', 'experiment', 'user', 'voltage']
```

also, values

```
In [15]: data_dict.values()
Out[15]: 
[47,
 3275.39,
 [-1.0, -0.5, 0.0, 0.5, 1.0],
 'current vs. voltage',
 'F. C. Rodgers',
 [-2.0, -1.0, 0.0, 1.0, 2.0]]
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

Since tuples are immutable, they can be used as keys for dictionaries.
Lists are mutable, and therefore cannot.

When you architect software in python, most data will end up looking either
like a list or a dictionary. These two data types are very important in
python and you'll end up using them all the time.

## Sets

Most introductory python courses do not go over sets this early (or at
all), and in the interest of time we're no different. The python set
type is a useful data type similar to the idea of a mathematical set:
it is an unordered collection of unique things.

Consider the following examples if you're interested in the useful
sorts of things you can do with python sets:

```
In [3] fruit = set(["apple", "banana", "pear", "banana"]) #You have to use a list to create a set.
```

Since sets contain only unique items, there's only one banana in the set
fruit.

You can do things like intersections, unions on sets just like in
math. Here's an example of an intersection of two sets (the common items in
both sets).

```
In [4]: first_bowl = set(["apple", "banana", "pear", "peach"])

In [5]: second_bowl = set(["peach", "watermelon", "orange", "apple"])

In [6]: set.intersection(first_bowl, second_bowl)
Out[6]: set(['apple', 'peach'])
```

You can check out more info using the help docs. The important thing here isn't that sets exist or what exactly they can do, but the main concept: different types of data can do different things easily. Choosing your data types well will both make your code clearer and simpler to write.

## Conditionals

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
and you end it with `end if`. In C/C++ you delimit blocks with curly
braces. Python uses **indentation** to delimit code blocks. The
**indentation** above is NOT just to make things look pretty - it
tells Python what the body of the `if`-statement is. This is true when
ever we create any code blocks, such as the bodies of loops, functions
or classes.

***Exercise***
Write an if statement that prints whether x is even or odd.

Hint: Try out what the "%" operator. What does 10 % 5 and 10 % 6 return?

##Writing Code for People summary

We learned some basics of python and saw that variable type, name, comments and white space affect more than just code functionality, they affect the readability for others and your future self.
Variable names can make a huge difference in code readability and types are important in conveying intent. 

[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)
