# Python 0: Variable Types (with focus on String and String methods)

* * * * *

**Based on lecture materials by Mike Jackson and Stephen McGough**

Python is an interpreted language. That means that there is no separate compilation step. Let's start and bring up the Python prompt:
```python
   Python 2.7.3 (default, Dec 18 2012, 13:50:09)
   [GCC 4.5.3] on cygwin
   Type "help", "copyright", "credits" or "license" for more information.
   >>>
```

Let's get Python to talk to us:
```python
   print "Hello world!"
   print 2 + 2
   print 'Thomas ' + 'Mann'
```

## Assigning values

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
   print string * number #What is the output?
```

We can convert types as values are typed:
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

Floating point division:
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

In order to list all Python built-in methods for string:
```python
   dir(str)
```

Next: [Data Types](1_Data_Types.md)
