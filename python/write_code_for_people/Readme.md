[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)


# Write Code for People: Variables, Data Structures and Conditionals
  

* * * * *


**Based on Lecture Materials By: Nate Vack, Milad Fatenejad, Katy Huff, Tommy Guy, Joshua 
R. Smith, Will Trimble, and many more**

This lecture usess the python scripting language as a mechanism to discuss the best practices in writing software.  Since we can't assume that you all know python, we'll be teaching a little python along the way.  Specifically, we will be using the iPython interpretter, which adds a few bells & whistles to make it easier to get up to speed.  I expect this lecture to be interactive, so stop me at any point if you have questions.

Importantly, this is **not** meant to be a comprehensive lesson in using python or iPython.  Therefore, we've included some reference material where you can get more information.

* [Dive into Python](http://www.diveintopython.net/toc/index.html)
* [Software Carpentry's Python Lectures](http://software-carpentry.org/4_0/python/)
* [IPython: A System for Interactive Scientific Computing](http://dx.doi.org/10.1109/MCSE.2007.53)
* [How to Think Like a Computer Scientist](http://www.greenteapress.com/thinkpython/thinkpython.html)

## What We'll Cover

After a brief introduction to iPython, we'll focus on the following importnat best practices.

1. How to write code for people. That is code that is readable and understandable to others in your group and most importantly to your future self 3 months or 3 years down the road.  This best practice is based largely in making the right choices on a number of topics:

   * Meaningful variable names
   * When to add comments
   * Appropriate data type
   * Effective and efficient compound data types

2. How to not repeat yourself. How to reuse your code with loops and functions. And how to eventually build modules, collections of functions, you and others can use in all your codes. And how to use other people's modules.

   * Using loops to repeat certain tasks locally
   * Conditional statements *PPHW check this for best practice*
   * Writing functions to repeat tasks globally
   * Using other poeple's modules to avoid reinventing the wheel
   * Writing modules to allow others to avoid reinventing the wheel

## iPython

Please follow the link to the [iPython Intro](../ipython/Readme.md).

## Back to Write Code for People

Remember that the purpose of this lesson is to understand how to make readable code, code for people!  It will be easy to think of this as simple an introduction to the basics of the python scripting language, but resist this urge.  Fortunately, you may also learn some python along the way.

## Intro to Variables

All programming languages have variables, and python is no different. To create a variable, just name it and set it with the equals sign. One important caveat: variable names can only contain letters, numbers, and the underscore character. Let's set a variable.

```
In [1]: first_name = "Paul"
In [2]: last_name = "Wilson"
```

We can glue strings together with the + operator. Computers don't understand context.

```
In [3]: full_name = first_name + last_name
In [4]: print(full_name)
Out[4]: 'PaulWilson'
```

***Exercise***
Can you add the extra space between my last and first name?

## Choosing Meaningful Variable Names

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

**Variable name choice is important**; a well-named variable is self-explanatory without comments and will make your code easier to read as the reader will not have to look up the comments. Remember that context matters a lot: in some cases, you'll want to spell out `voltage` or even `input_voltage`. In other cases, `v` is a shorthand that everyone will understand.

It's also important to choose a naming convention for your whole project, and get the people working on it to agree. Mixing `batteryVoltage` (CamelCase) and `capacitor_value` (pothole) in one place makes code hard to read.

## Writing Comments For People

Above, you may have noticed the `#` character, which denotes a comment in python. Comments should describe meaning but not what the statement is doing.

```
In [15]: voltage = 4 # set the voltage to 4   <- well, duh

In [16]: voltage = 4 # Input to the circuit   <- better; says what this voltage means
```

----
![Exercise](../../common/pics/exercise.jpg) **Short Exercise**

Think of a task in your own work where you might use python.  List 4 quantities that would be represented by variables and suggest meaningful variable names for them.

----

## Avoiding Magic Numbers (and Magic Strings)

There are frequently times when a number is needed only once in an equation and it may be convenient to avoid the use of a variable.

```
In [17]: weight = mass * 9.81
```

However, this number may be confusion in the future to you or another reader.  You could add a comment:

```
In [18]: weight = mass * 9.81  # multiply the mass by the acceleration due to gravity
```

This, however, is a nearly trivial comment that makes the code a little messy for readers.  The best choice is to define a variable for this constant:

```
In [19]: gravity = 9.81
In [20]: weight = mass * gravity
```

This option has three important advantages:
1. it is very readable,
2. the `gravity` variable can be reused consistently elsewhere, if necessary, and
3. if you do use it elsewhere it can be changed consistently in only one place.

## Intro to Types and Dynamic Typing

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

It's an `int`, which is an integer — a number with no decimal component. Let's assign a value of 2.7 (which has a decimal part) to `voltage`. What happens to the type?

```
In [11]: voltage = 2.7

In [12]: type(voltage)
Out[12]: float
```

Neat! That's a `float`, which does have a decimal part. You can assign a string to the variable `voltage` in the same way:

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


## Choosing Your Type Explicitly With floats and ints

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

So what about the case where `a` is an integer and `b` is a float?

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

Importantly, this explicit choice of type communicates important information to you and other readers of your code.

----
![Exercise](../../common/pics/exercise.jpg) **Short Exercise**

For the variables you thought of in the last exercise, describe which type they should be.

----


## Introduction to Compound Data Types

Most languages offer compount data types of some kind (arrays, vectors, lists, matrices, maps, sets, etc.) and would not be very useful without them!  Python is no different, and much of its value comes from the ease of using them. The main two are lists and dictionaries, but I'll mention sets and tuples as well. 

## Lists

A `list` is an ordered, indexable collection of data. Let's say you have
collected some current and voltage data that looks like this:

|  voltage | current |
|----------|---------|
|  -2.0    |  -1.0   |
|  -1.0    |  -0.5   |
|   0.0    |   0.0   |
|   1.0    |   0.5   |
|   2.0    |   1.0   |

So you could put that data into lists like

```
In [1]: voltage_list = [-2.0, -1.0, 0.0, 1.0, 2.0]

In [2]: current_list = [-1.0, -0.5, 0.0, 0.5, 1.0]
```

obviously `voltage_list` is of type `list`:

```
In [3]: type(voltage_list)
Out[3]: list
```

Python lists have the charming (annoying?) feature that they are indexed
from zero. Therefore, to find the value of the first item in `voltage_list`:

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
`current_list`

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
through fourth items from `voltage_list`

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

Notice that the second number in a slice should be selected carefully.  You can think of it as the index of the first element that is NOT included in the slice.

----
![Exercise](../../common/pics/exercise.jpg) **Short Exercise**

* What is the difference between the following?
   * `voltage_list[0:0]`
   * `voltage_list[0:1]`
   * `voltage_list[:1]`

* What does `current_list[-3:]` give?

* What does `voltage_list[::2]` mean?

----

### Append and Extend

Just like strings have methods, lists do too.

```
In [10] dir(list)
```

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

|  voltage | current |
|----------|---------|
|   3.0    |   1.5   |
|   4.0    |   2.0   |


If you want to append items to the end of a list, use the `append` method.

```
In [11]: voltage_list.append(3.)

In [12]: voltage_list.append(4.)

In [13]: voltage_list
Out[13]: [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
```

You can see how that approach might be tedious in certain cases. If you
want to concatenate a list onto the end of another one, use `extend`.

```
In [14]: current_list.extend([1.5, 2.0])

In [15]: current_list
Out[15]: [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
```

----
![Exercise](../../common/pics/exercise.jpg) **Short Exercise**

Make a new list:

```
In [16]: exercise_list = [1,1,2,3,5,8]
```

What is the difference between the following?
 * `exercise_list.append([13,21,34])`
 * `exercise_list.extend([13,21,34])`

----


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

If you want to see if that's what's going on in your case, python has a special `is` operation that will test to see if two variables point to the same thing:

```
In [23]: a is b
Out[23]: True
In [24]: c = [1, 2, 10]
In [25]: a is c
Out[25]: False
In [26]: a == c
Out[26]: True
```

There's a ton more to know about lists, but let's press on. [Dive into Python](http://www.diveintopython.net/toc/index.html) or the help documentation for more info.

## Tuples

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is **immutable**; once you set the
data in it, the tuple cannot be changed. You define a tuple as follows.

```
In [1]: tup = ("red", "white", "blue")

In [2]: type(tup)
Out[2]: tuple
```

You can slice and index the tuple exactly like you would a list, but you cannot append or extend it.

----
![Exercise](../../common/pics/exercise.jpg) **Short Exercise**

Display the second element of a tuple with two different slices.

----

### Heterogeneous Data

Lists and Tuples can contain hetergeneous data.

```
In [17]: data_list = ["experiment: current vs. voltage",
                      "run", 47,
                      "temperature", 372.756,
                      "current", [-1.0, -0.5, 0.0, 0.5, 1.0],
                      "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]

```

We've got strings, ints, floats, and even other lists in there. While this is a perfectly valid thing to do in a list, it's often not best practice for a number of reasons:
* the association between the name (e.g. "run", "temperature", etc) and the data (47, 372.756) is implied but not guaranteed,
* lists imply that the order is important which is not the case here, and
* it is common to act on each member of a list with the same operation

Fortunately, python has dictionaries for exactly this purpose!

## Dictionaries

Our variable `data_list` contained our current-voltage data and also some metadata. We were able to store the data as a list, but clearly the list type is not the optimal choice for a data model. The dictionary is a much better choice.

A python dictionary is a collection of key, value pairs. The key is a way to name the data, and the value is the data itself.  Here's a way to create a dictionary that contains all the data in our data.dat file in a more sensible way than a list.

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
In [13]: data_dict["user"] = "P.P.H. Wilson"
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
 'P.P.H. Wilson',
 [-2.0, -1.0, 0.0, 1.0, 2.0]]
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

Since tuples are immutable, they can be used as keys for dictionaries.
Lists are mutable, and therefore cannot.

When you architect software in python, most data will end up looking either like a list or a dictionary. These two data types are very important in python and you'll end up using them all the time.  **Most importantly, choosing the correct data type will make a big difference in the readability of your code!**

## Sets

The python set type is a useful data type similar to the idea of a
mathematical set: it is an unordered collection of unique things.

Consider the following examples if you're interested in the useful
sorts of things you can do with python sets:

```
In [3]: fruit = set(["apple", "banana", "pear", "banana"]) #You have to use a list to create a set.

In [4]: fruit
Out[4]: {'apple', 'banana', 'pear'}
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


### Length of Compound Data Types

Sometimes you want to know how many items are in a compound data type. Use the `len` command.

```
In [16]: len(voltage_list)
Out[16]: 7

In [17]: len(data_dict)
Out[17]: 6

In [18]: len(set.intersectin(first_bowl, second_bowl))
Out[18]: 2
```

## Choosing Effective and Efficient Compound Data Types

Remember that we discussed these compound data types because the choice of data type can make a big difference in the ability for others (and you in 3 months) to understand your code.  You should consider some of the following questions when considering your choice of data type.

1. Does a compound data type make more sense than multiple variables?

2. Is the order of the elements important? (probably need a list)

3. Is the data immutable? (probably want a tuple)

4. Is the data named with keys and values?  (probably want a dictionary)


## Writing Code for People summary

We learned some basics of python and saw that variable type, name, comments and white space affect more than just code functionality, they affect the readability for others and your future self.
Variable names can make a huge difference in code readability and types are important in conveying intent. 

## Preparing for more complicated work

Everything we have done so far relies on single lines of python or short sequences.  As we prepare for longer blocks of python, it may be useful to learn how to [use multiple windows](../dont_repeat_yourself/using_multiple_windows.md) to make life easier.



[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)
