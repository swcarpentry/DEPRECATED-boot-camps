[Up To Schedule](../../README.md) - Return To [Write Code for People](compound_data_types.md)

- - - -

# Compound Data Types: Lists, Dictionaries, Sets, and Tuples

Python would be a fairly useless language if it weren't for the compound data
types. The main two are lists and dictionaries, but I'll mention sets and
tuples as well.

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

Lists have built-in methods, like many other data types.

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

**Try it exercise:**

What happens if you try to add two items to a list using append instead of extend?

```
In [16]: current_list.append([1.5, 2.0])
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
In [17]: data_list = ["experiment: current vs. voltage",
                      "run", 47,
                      "temperature", 372.756,
                      "current", [-1.0, -0.5, 0.0, 0.5, 1.0],
                      "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]

```

We've got strings, ints, floats, and even other lists in there. While this is
a perfectly valid thing to do in a list, it's often not a good
idea. Generally, you want to be able to run the same operation on every
element of a list. If you want to put different kinds of things in a list
structure, you often want to use something called a *tuple,* which we'll talk
about in a minute.

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

If you want to see if that's what's going on in your case, python has a
special `is` operation that will test to see if two variables point to the
same thing:

```
In [23]: a is b
Out[23]: True
In [24]: c = [1, 2, 10]
In [25]: a is c
Out[25]: False
In [26]: a == c
Out[26]: True
```

There's a ton more to know about lists, but let's press
on. [Dive into Python](http://www.diveintopython.net/toc/index.html) or the
help documentation for more info.

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

Tuples can emphasize intent in two ways. First, you can't easily change the
items in a tuple, so if you want to specify that a list shouldn't change, a
tuple is a great way to indicate that. Second, in a tuple, the order of the
items is usually more significant: you might do something like:

```
In [1]: person_data = ('Nate', 35, 'nate@xmail.com')
In [2]: name = person_data[0]
```

... though in practice, you'll more often use a `dictionary` for something
like this.

***Exercise***
Display the second element of a tuple with two different slices.

## Dictionaries

Recall our variable data_list which contained our current-voltage data and
also some metadata. We were able to store the data as a list, but clearly the
list type is not the optimal choice for a data model.  In this example, the
order of the data implied a relationship between some keywords and their data,
e.g. the keyword "run" had value "47".  But since lists are not immutable, it
is easy to accidentally break those connections by deleting one of the list
members.

The dictionary is a much better choice.  It creates a strong binding between
the keyword and its data so that those connections are clear, with intent, and
cannot be accidentally broken.

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

Most introductory python courses do not go over sets this early (or at all),
and in the interest of time we'll just introduce them. The python set type is
a useful data type similar to the idea of a mathematical set: it is an
unordered collection of unique things.

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

You can check out more info using the help docs. The important thing here
isn't that sets exist or what exactly they can do, but the main concept:
different types of data can do different things easily. Choosing your data
types well will both make your code clearer and simpler to write.

## Iterating through compound data types

One of the features of python that makes it very readable is the ease with
which you can iterate over a compound data structure.  All four of these data
structures can use the `for each_item in compound:` pattern:

```
In [7]: for item in voltage_list:
   ...:     print item
In [8]: for item in person_data:
   ...:     print item
In [9]: for item in first_bowl:
   ...:     print item
```

Dictionaries are only a little bit different because they contain both keys and values:

```
In [10]: for item in data_dict.keys():
   ....:     print item + " = " + data_dict[item]
```


- - - -

[Up To Schedule](../../README.md) - Return To [Write Code for People](compound_data_types.md)
