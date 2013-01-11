# Python 2: Lists, Dictionaries, Sets, Tuples

[Back To Python
Variables](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/2a-PythonVariables/)
- [Forward to Python Flow
Control](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/2c-PythonFlowControl/)

* * * * *

**Presented By : Joshua R. Smith**

**Based on lecture materials by Milad Fatenejad**

## Lists

Most languages have some kind of simple syntax for making lists of
things. In python it is extremely easy and intuitive to make a list of
things, for example:

```python
> mylist = [] # Make an empty list
> mylist = [1, 1.0+3j, "aperitivo", True] # Make a list containing four entities
```

Using lists is easy and intuitive. Notice that lists can contain objects
of any data type. Try entering the following lines.

```python
> mylist = [1,2,3,4]
> mylist[2] = 1.0 + 2j # Modify an element
> mylist.append("test") # Add an element to the end of a the list
> print len(mylist) # print the length of mylist (5)

> mylist = [1,9,7,32,2.0]
> del(mylist[2]) # Remove element 2 from the list

> mylist = [1,5,4,2]; mylist.sort() # Sort the list
```

The colon, **:**, provides a syntax for ranges.

```python
> print mylist[1:4] # Prints a list containing elements 1 2 and 3 from mylist 
```

Remember that there is an element 0, so this prints [4, 6, 8]

```python
> print mylist[-2] # Print the second element from the end of the list (8)
```

## Dictionaries

Lists aren't the only compound data type. Another really useful one is a
dictionary (referred to as a map in many other languages). Dictionaries
allow you to set/access elements using a key value relationship. You can
create dictionaries as shown below:

```python
> mydictionary = {} # Make an empty dictionary
> mydictionary = {"one" : 1, "two" : 2, "three" : 3} # Initialize a dictionary 
> with some values

> print type(mydictionary) # Tells you mydictionary is of type "dict"
> print mydictionary["one"] # Prints the number 1
> print mydictionary["two"] # Prints the number 2
> mydictionary["four"] = 4 # Insert an element into the dictionary
> mydictionary["list"] = [1,2,3] # Sets the element "list" to a list containing 
> the numbers 1, 2, and 3
```

### Example: Creating and Sorting Lists

Accomplish the following tasks using Python. Each task should take only
one line. You may need to use the help and dir functions to figure out
parts you don't know:

1.  Create a string and initialize it to "Joshua Katy Milad Anthony"

2. Split the string into a list whose elements are the names Joshua, Katy, Milad, and Anthony. 

3. Sort and print the list 

4. Without deleting anyone, add your own name to the list **so that it comes
first.**

### Example: Manipulating Compound Data

Accomplish the following tasks using Python. Each task should take only
one line. You may need to use the help and dir functions to figure out
parts you don't know:

1.  Create a dictionary containing the key, value pairs:
    -   "Red", 5
    -   "Green", 3
    -   "Purple", 3
    -   "Orange", 1
    -   "Blue", 3
    -   "Teal", 3

2. Extract a list of values from the dictionary (i.e. get a list
containing [3,3,3,3,1,5] from the dictionary, don't make the list on
your own) 

3. Find and use a list method to count the number of times the
value 3 appears (Use the list you produced on step 2, the correct answer
is that the value 3 appears four times)

In a dictionary, the keys must be unique: assigning a second value to a
key overwrites whatever was stored there. What if we want to store a
list of unique items? There are two options using what we know about so
far:

1. Use a list, but every time we add an element, check whether it is
already there. 

2. Use a dictionary to store the object as a key to some dummy value.

## Sets

It turns out there is a third type of container in Python that only
stores unique things: it's called a set.

```python
> s = set()
> s = set([1,1,2,3,4]) # Note that there are 2 1's in the input list.
> print s
set([1, 2, 3, 4])
> 1 in s
True
> 5 in s
False
> s.add(5)
> 5 in s
True
> anotherSet = set([1,2,"hello"])
> s.intersection(anotherSet)
set([1, 2])
```

### Example : Appending/Adding vs Updating

There are two methods to add element(s) to a list: append and update.
Likewise, there are two methods to add element(s) to a set: add and
update. What is the difference?

```python
> myList = [1,2,3]
> myAppendedList = myList.append([4,5])
> myUpdatedList = myList.update([4,5])
```

What is the difference between the appended list and the updated list?
Why did this happen?

1. Try the same thing with the add() and update() functions on a set.
The key is that containers can hold other containers.

## Tuples

There is one other compound data type - the tuple. Think of a tuple as a
list that you can't change. The example below demonstrates how to create
and use tuples:

```python
> mytuple = (1,2,3,4) # Create a four element tuple
> mytuple[2] = 4 # ERROR - tuples can't be modified
> print mytuple[2], len(mytuple)

> myonetuple = ("hello",) # Create a tuple containing only one element (note the trailing comma)
```

You might be asking yourself, why do we need tuples if we have lists?
The answer is that tuples are used internally in Python in a lot of
places. One of the basic differences is that dictionaries cannot use a
list as a key, but they can use a tuple:

```python
> d = {}
> d[(1,2)] = 'numbers'
> d
{(1, 2): 'numbers'}
> d[ [1,2] ] = 'listOnumbers'
 Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
 TypeError: unhashable type: 'list'
```

As you learn more about python you'll see how lists, tuples and
dictionaries are the basic building blocks of the entire language.

## Copy or Reference?

Simple data types like integers and strings behave slightly differently
than more complicated objects. To see one unexpected example, try these
commands:

```python
> list1 = [1, 5, 9, 13]
> list2 = list1
> list2[0] = -1
> print list1, list2
```

What happens? You'll notice that modifying list2 also modifies list1!
This is because line 2 does not copy list1, instead list2 is set to
*reference* the same data as list1. After line 2 is executed, list1 and
list2 refer to the same data. Modifying one list also modifies the
other. This was not the case when we were dealing with simple numbers.
This behavior can be very annoying and can lead to a lot of bugs, so be
careful. We can force python to copy list1 as shown in the example
below:

```python
> list1 = [1, 5, 9, 13]
> list2 = list1[:] # <--- Notice the colon!
> list2[0] = -1
> print list1, list2
```

Why would Python variables be references rather than copied instances?
Let's think about a list file handles or some other complicated data
type. In two lectures, we will talk a lot more about file handles. If
you've ever done file I/O, you know that the file handle has a *state*:
it stores the location of the next readable byte in the file as well as
a lot of other information. When we copy the list, we probably want to
maintain the same set of file handles rather than reinitialize them all.
Python made the decision that the default should be to pass the
reference to the list rather than create a copy of each element.
