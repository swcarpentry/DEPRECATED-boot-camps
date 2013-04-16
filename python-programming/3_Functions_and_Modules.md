#Python 3 : Functions and Modules
-----------------------

A function is a block of code that performs a specifc task. In this section we
will learn how to write our own functions and modules, and import them. The example used in this section is a simple one. But the main purpose of this part of the tutorial is not only to show you how to write our own functions and modules but `why` do it. The goal is to show you how and why modularisation of the source code is a good programming practice.

Writing our own functions
============



```python
   def get_initials(line):
      names = line.split()
      initial = ''
      i = 0
      while i < len(names):
         initial=initial+(names[i][0])
         i = i+1
      return initial
   
   b = 'Thomas Mann'
   print get_initials(b)
   print get_initials('Anna Karenina')
```


Importing Python modules
============
