# Python 1: Data Types: Lists, Dictionaries, Tuples, and Reading Files

* * * * *

One of the useful features of Python are its compound data types. The main two are lists and dictionaries, but we will mention tuples as well. We will also go over reading text data from files. 

## Lists

A list is an ordered, indexable collection of data. Lets say you have
collected some data that looks like this:

```
Name:
rapA
polB
araD
araA
araB
araC

AT_content:
45.1
44.7
44.7
44.5
41.9
47.4
```

So you could put that data into lists like

```python
names = ['rapA', 'polB', 'araD', 'araA', 'araB', 'araC']
print names
AT_content = [45.1, 44.7, 44.7, 44.5, 41.9, 47.4]
print AT_content
```

We can check the names type (obviously it is of type list):

```python
type(names)
```

Python lists have the charming (annoying?) feature that they are indexed from zero. Therefore, to find the value of the first item in names, we have to check index 0:

```python
names[0]
```

And to find the value of the third item:

```python
names[2]
```

Lists can be indexed from the back using a negative index. The last item of
AT_content

```python
AT_content[-1]
```

and the next-to-last

```python
AT_content[-2]
```

You can "slice" items from within a list. Lets say we wanted the second
through fourth items from names

```python
names[1:4]
```

Or from the third item to the end

```python
names[2:]
```

and so on. Note: if you do not have a number on both sides of the :, the side lacking a number will default to 
the end/beginning of the list


### Append and Extend

Just like strings have methods, lists do too.

```python
dir(list)
```

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

```
extraNames:
yabI
thiQ

extraContents:
46.1
42.9
```

If you want to append items to the end of a list, use the append method.

```python
names.append('yabI')
names.append('thiQ')
print names
```
You can now see these items at the end of the list. Note: append only allows
you to stick one thing onto a list. You can see how that approach might be tedious in certain cases. If you
want to add a list onto the end of another one, we can create a new list and 
stick that list onto our list.

```python
AT_content.extend([46.1, 42.9])
print AT_content
```

**Task:** create a copy of AT\_content (either assign AT\_content to a new 
variable, or just create a new list with the same content). Then, create a new list
with the AT contents of the two new proteins. Try using both append and extend to add the new values.
Can you see what the difference between append and extend is?

### Length of Lists

Sometimes we want to know how many items are in a list. Use the len command.

```python
len(names)
```

### Heterogeneous Data and Lists of Lists

Lists can contain hetergeneous data.

```python
dataList = ["Experiment", "AT content effects", \
          "run", 47, \
          "temperature", 55, \
          "names", ['rapA', 'polB', 'araD', 'araA', 'araB', 'araC', 'yabI', 'thiQ'], \
          "AT_content", [45.1, 44.7, 44.7, 44.5, 41.9, 47.4, 46.1, 42.9]]
```

We've got strings, ints, floats, and even other lists in there. The slashes
are there so we can continue on the next line. They aren't necessary but
they can sometimes make things look better.

Notice that there is actually two lists contained in this list. You can access these using indexing
as before. Let's experiment a bit:

```python
print dataList[7]
print dataList[7][2]
````
You now see you got the name of the third protein in the list (again, python list are zero based).

**Task:** can you figure out what the AT content of this protein is?


## Reading From Files

At this point it is useful to take a detour regarding files. Lets say you
have a file with some data in it:

```
data.dat:

Experiment: AT content effects
run: 47
temperature: 55
names: ['rapA', 'polB', 'araD', 'araA', 'araB', 'araC', 'yabI', 'thiQ']
AT content: [45.1, 44.7, 44.7, 44.5, 41.9, 47.4, 46.1, 42.9]
```

We can read this data into a variable of type list pretty easily.

```python 
fh = open("data.dat")
ptdata = fh.readlines()
fh.close()

print ptdata
```

Your results should look like this:

```
['Experiment: AT content effects\n', 'run: 47\n', 'temperature: 55\n', 
"names: ['rapA', 'polB', 'araD', 'araA', 'araB', 'araC', 'yabI', 'thiQ']\n", 
'AT content: [45.1, 44.7, 44.7, 44.5, 41.9, 47.4, 46.1, 42.9]\n']
```
What happened here, is that each line in your file became one element in the ptdata list. Note, the newlines are also included.

Right now the data in ptdata isn't in a particularly useful format, but you
can imagine that with some additional programming we could straighten it
out. We will eventually do that.

## Tuples

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is immutable; once you put the
data in it, the tuple cannot be changed. We define a tuple as follows.

```python
tup = ("red", "white", "blue")
type(tup)
```

You can slice and index the tuple exactly like you would a list. Tuples are
used in the inner workings of python, and a tuple can be used as a key in a
dictionary (which you will see soon), whereas a list cannot as we will see in a moment.



## Dictionaries

Recall our file data.dat which contained our experiment data and also
some metadata. We were able to import the data as a list, but clearly the
list type is not the optial choice for a data model. The dictionary is a
much better choice. A python dictionary is a collection of key-value
pairs. The key is a way to name the data, and the value is the data itself.
Here's a way to create a dictionary that contains all the data in our
data.dat file in a more sensible way than a list.

```python
dataDict = {'Experiment' : 'AT content effects', 'run': 47,  \
'temperature': 55, \
'names': ['rapA', 'polB', 'araD', 'araA', 'araB', 'araC', 'yabI', 'thiQ'], \
'AT content': [45.1, 44.7, 44.7, 44.5, 41.9, 47.4, 46.1, 42.9]}
```

This model is clearly better because you no longer have to remember that
the run number is in the second position of the list, you just refer
directly to "run":

```python
print dataDict["run"]
```

If you wanted the names data list:

```python
print dataDict["names"]
```

Or perhaps you wanted the last element of the AT content data list

```python
print dataDict["AT content"][-1]
```

Once a dictionary has been created, you can change the values of the data
if you like.

```python
dataDict["temperature"] = 66
print dataDict["temperature"]
```

You can also add new keys to the dictionary.

```python
dataDict["Experimenter"] = "Jane Doe"
print dataDict["Experimenter"]
```

Dictionaries, like strings, lists, and all the rest, have built-in methods.
Lets say you wanted all the keys from a particular dictionary.

```python
dataDict.keys()
```

also, values

```python
dataDict.values()
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

However, the key has to be something that cannot change. This means numbers, strings
and tuples. Lists are mutable, and therefore cannot.

When you architect software in python, most data will end up looking either
like a list or a dictionary. These two data types are very important in
python and you'll end up using them all the time.

Previous: [Variable types](0_Variables_Types.md) Next: [Flow control](2_Flow_Control.md)
