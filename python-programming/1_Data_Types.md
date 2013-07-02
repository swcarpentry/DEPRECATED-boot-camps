# Python 1: Data Types: Lists, Dictionaries, Tuples, and Reading Files

* * * * *

One of the useful features of Python are its compound data types. The main two are lists and dictionaries, but we will mention tuples as well. We will also go over reading text data from files. 

## Lists

A list is an ordered, indexable collection of data. Lets say you have
collected some data that looks like this:

```
Country:
Norway
Sweden
Denmark
Finland
Iceland

lifeExp:
77.67
71.86
70.78
66.55
72.49
```

So you could put that data into lists like this:

```python
countries = ['Norway', 'Sweden', 'Denmark', 'Finland', 'Iceland']
print countries
lifeExp = [77.67, 71.86, 70.78, 66.55, 72.49]
print lifeExp
```

We can check the countries type (obviously it is of type list):

```python
type(countries)
```

Python lists have the charming (annoying?) feature that they are indexed from zero. Therefore, to find the value of the first item in names, we have to check index 0:

```python
countries[0]
```

And to find the value of the third item:

```python
countries[2]
```

Lists can be indexed from the back using a negative index. The last item of
lifeExp

```python
lifeExp[-1]
```

and the next-to-last

```python
lifeExp[-2]
```

You can "slice" items from within a list. Lets say we wanted the second
through fourth items from countries

```python
countries[1:4]
```

Or from the third item to the end

```python
countries[2:]
```

and so on. Note: if you do not have a number on both sides of the :, the side lacking a number will default to 
the end/beginning of the list

### Heterogenous lists

List can contain any kind of data, including other lists. Lets create a new list containing the two lists above and also some other information:

```python
combinedList = ["Information about countries", "1952", countries, lifeExp]
print combinedList
```

How do we access elements of this list? Just like with any other list.

```python
# print the information line:
print combinedList[0]
# print the countries list:
print combinedList[2]
# print out the third country in the countries list
print combinedList[2][2]
```

**Task:** Figure out how to print the year 1952 to the screen.

**Task:** Figure out how to print the life expectancy list on screen.

**Task:** Figure out how to print the life expectancy of Finland to the screen. 


### List methods

Just like strings have methods, lists do too.

```python
dir(list)
```

#### Append and Extend

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

```
extraCountries:
Germany
France

extraLifeExp:
67.5
67.41
```

If you want to append items to the end of a list, use the append method.

```python
countries.append('Germany')
countries.append('France')
print countries
```
You can now see these items at the end of the list. Note: append only allows
you to stick one thing onto a list. You can see how that approach might be tedious in certain cases. If you
want to add a list onto the end of another one, we can create a new list and 
stick that list onto our list. Note, for this use extend instead of append.

```python
extraLifeExp = [67.5, 67.41]
lifeExp.extend(extraLifeExp)
print lifeExp
```

**Task:** Open help for lists - help(list), find extend and append, and figure out the difference between them.


#### Sorting a list

We can also sort lists. 

```python
countries.sort()
print countries
```

Note: we can also reverse the sort by setting reverse as True:

```python
countries.sort(reverse=True)
print countries
```

#### Length of Lists

Sometimes we want to know how many items are in a list. Use the len command.

```python
len(countries)
```


## Assigning Variables to Other Variables

Something that might cause you headaches in the future is how python deals
with assignment of one variable to another. When you set a variable equal
to another, both variables point to the same thing. Changing the first one
ends up changing the second. Be careful about this fact. Note however, this 
only goes for "compound" datatypes, not basic datatypes as numbers and strings. 
Let´s see how this works.

First with basic datatypes:

```python
a = 1
b = a
print a
print b
a = 2
print a
print b
```
Now with a list:

```python
a = [1,2]
b = a
print a
print b
a.append(10)
print a
print b
```


## Reading From Files

At this point it is useful to take a detour regarding files. Lets say you
have a file with some data in it:

```
country	year	pop	continent	lifeExp	gdpPercap
Norway	1952	3327728	Europe	72.67	10095.42172
Sweden	1952	7124673	Europe	71.86	8527.844662
Denmark	1952	4334000	Europe	70.78	9692.385245
Finland	1952	4090500	Europe	66.55	6424.519071
Iceland	1952	147962	Europe	72.49	7267.688428
```

As you can see, some of the data we looked at earlier stems from this file. 

We can read this data into a variable of type list pretty easily.

```python 
fh = open("popdata.txt")
data = fh.readlines()
fh.close()

print data
```

Your results should look like this:

```
['country\tyear\tpop\tcontinent\tlifeExp\tgdpPercap\n', 'Norway\t1952\t3327728\tEurope\t72.67\t10095.42172\n', 'Sweden\t1952\t7124673\tEurope\t71.86\t8527.844662\n', 'Denmark\t1952\t4334000\tEurope\t70.78\t9692.385245\n', 'Finland\t1952\t4090500\tEurope\t66.55\t6424.519071\n', 'Iceland\t1952\t147962\tEurope\t72.49\t7267.688428\n']
```
What happened here, is that each line in your file became one element in the data list. Note, the newlines are also included.

Right now the data isn't in a particularly useful format, but you
can imagine that with some additional programming we could straighten it
out. We will eventually do that.

## Dictionaries

Remember the two lists we worked on earlier? These contained some of the data in the file above. In that case, we worked on each list individually. There was no connection between each country and its life expectancy value. For these kinds of connections a dictionary is a much better choice. 
much better choice. A python dictionary is a collection of key-value
pairs. The key is a way to name the data, and the value is the data itself.
Here's a way to create a dictionary that contains all the data in our
data.dat file in a more sensible way than a list.

```python
dataDict = {"Norway":72.67, "Sweden": 71.86, "Denmark", 70.78, "Finland": 66.55, "Iceland": 72.49}
```

This model is clearly better because you no longer have to remember the position in the list of each country.

```python
print dataDict["Sweden"]
```

Once a dictionary has been created, you can change the values of the data
if you like.

```python
dataDict["Norway"] = 100
print dataDict["Norway"]
```

You can also add new keys to the dictionary.

```python
dataDict["Germany"] = 67.5
print dataDict["Germany"]
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
