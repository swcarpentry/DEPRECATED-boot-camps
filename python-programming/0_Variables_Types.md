# Python 0: Variable Types (with focus on String and String methods)

* * * * *

**Based on lecture materials by Mike Jackson and Stephen McGough**

**Modified by Karin Lagesen**

There are many different kinds of programming languages. The language we are teaching you today is python, which is an interpreted language. This means that we can write in code and it will be executed right away. 

You can deal with python in two ways, either directly in an interpreter, or you can write bits of code in a file and have that run by the interpreter for you. 

NOTE: a line anywere in python that is preceeded with a # is a comment and is ignored by python!

## The python interpreter

Let's bring up the Python interpreter. This is what it looks like on cygwin - linux is pretty similar:
```python
   Python 2.7.3 (default, Dec 18 2012, 13:50:09)
   [GCC 4.5.3] on cygwin
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
   print 'Ada ' + 'Lovelace'
```
Strings can actually be added together.

Q: what happens if you don't include the space in Thomas?


## Assigning values

In a programming situation, we usually have data that we want to do something to. We want to be able to work on the data without having to specify it directly, like we did with the numbers and "Hello World" above. Instead of specifying the data directly, we can use variables instead. Variables contains data that we can do things to. Many languages require you to specify what kind the data a variable should contain should be, like strings or numbers. Python does not require this, we can use them directly.

Python does not require us to declare variables as they are created when we first use them:
```python
   a = 1
   print a
   a = 'pigeon'
   print a
   a = 2 + 7 
   print a
```

Variables in Python are not typed as you can see. They can be thought of as being just names. However, beware that Python does not assume default variables - you have to assign them before use:
```python
   print b
```

What happens when we try to "combine" two Python variables: string and integer:
```python
   string = "text"
   number = 4
   print string * number 
```
Q: what was the output here?

## Type conversion

Sometimes a number can be interpreted as a text string when you need it to work as a number. This can be fixed by converting it to an integer (a whole number). We can also convert to decimal numbers - floats. These have a decimal point in them.

```python
   print int('2') + 3
   print '2' + str('3')
   print string + str(number)
```

Python supports typical (and popular) arithmetic operations: 
Adding:
```python
   print 10 + 3
```

Integer division:
```python
   print 10 / 3
```

Reminder:
```python
   print 10 % 3
```




## String operations

Python provides a number of very useful methods for string variables.
We can check the string length:

```python
   a = 'Manchester'
   print len(a)
```

We can also split the string - splitting at space:
```python
   a = 'Thomas Mann'
   b = a.split()
   print b
```

Or by any other character:
```python
   a = 'Thomas-Mann'
   b = a.split('-')
   print b
```


In order to list all Python built-in methods for string:
```python
   dir(str)
```



Next: [Data Types](1_Data_Types.md)
