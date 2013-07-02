# Python 0: Variable Types (with focus on String and String methods)

* * * * *

There are many different kinds of programming languages. The language we are teaching you today is python, which is an interpreted language. This means that we can write in code and it will be executed right away. 

You can deal with python in two ways, either directly in an interpreter, or you can write bits of code in a file and have that run by the interpreter for you. 

NOTE: a line anywere in python that is preceeded with a # is a comment and is ignored by python!

## The python interpreter

Let's bring up the Python interpreter. This is what it looks like on a mac - linux is pretty similar:
```python
Python 2.7.1 (r271:86832, Jun 16 2011, 16:59:05) 
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2335.15.00)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> 
```
See the three >? You type in your commands after that.

The interpreter can be used as a calculator:

```python
print 2 + 2
print 2*8
```

We can also work with strings in python:

```python
print "Hello world!"
print 'Ada' + ' Lovelace'
```
You can see that these are actually strings and part of the python language in that they have quotation marks around them. Note: to create strings, you can use both single and double quotes.

As you can see, strings can actually be added together.

**Task:** Note that there is a space before Lovelace in the last instruction. Repeat that instruction without this space. Can you see the difference? 

**Task:** type in `print 'Ada', 'Lovelace'`. Can you see what happened here?

## Assigning values to variables

In a programming situation, we usually have data that we want to do something to. We want to be able to work on the data without having to specify it directly, like we did with the numbers and "Hello World" above. Instead of specifying the data directly, we can use variables instead. Variables contains data that we can do things to. Many languages require you to specify or declare what kind the data a variable should contain should be, like strings or numbers. Python does not require this, we can use them directly.

Python creates variables when we assign a value to them:
```python
a = 1
print a
a = 'pigeon'
print a
a = 2 + 7 
print a
```

Assignment is done by having a variable name to the left, then an equal sign, and then the value to the left of the equal sign.

Variables in Python are not returned to you by the interpreter, as you can see. They can be thought of as being just names. However, beware that Python does not assume default variables - you have to assign them before use:
```python
print b
```

We can try to "combine" two Python variables: string and integer:

```python
string = "text"
number = 4
print string * number 
```

Here you see that the word is repeated four times, as many times as number specifies.

However, what happens if we add string and number?

```python
print string + number
```

Python has what we call dynamic and strong typing - dynamic in that you do not have to declare variables, and strong in that it objects if you try to do something to a variable that the data type in that variable does not permit.

## Type conversion

Sometimes a number can be interpreted as a text string when you need it to work as a number. This can be fixed by **converting** it to an integer (a whole number). We can also convert to decimal numbers - floats. These have a decimal point in them.

```python
print int('2') + 3
print '2' + str(3)
print string + str(number)
```

Python supports typical (and popular) arithmetic operations: 
Adding:
```python
print 10 + 3
```

Remainder:
```python
print 10 % 3
```

Integer division:
```python
print 10 / 3
```

Note: python 2x does integer division, i.e. it floors the results. To get a float out, we will need to ensure that at least one of the numbers is a float, by for instance multipylying it with 1.0.
 

## String operations

Python provides a number of very useful methods for string variables.
We can check the string length:

```python
a = 'Oslo'
print len(a)
```

We can also split the string - splitting at space:
```python
a = 'Ada Lovelace'
b = a.split()
print b
```

Or by any other character:
```python
a = 'Ada-Lovelace'
b = a.split('-')
print b
```

We can also count the occurrence of a substring:

```python
a = 'She sells seashells'
ses = a.count("s")
print ses
```

**Task**: did the number of s-es reported correspond to how many you counted?

In order to list all Python built-in methods for string:
```python
dir(str)
```

To get help on methods for a type, do

```python
help(str)
```

**Task:** Can you figure out a python method that will let you remove all whitespace to the right in the string? Hint: search for whitespace. Use this method to remove the trailing whitespace in this string: `"This string has lots of end spaces      "`.


Next: [Data Types](1_Data_Types.md)
