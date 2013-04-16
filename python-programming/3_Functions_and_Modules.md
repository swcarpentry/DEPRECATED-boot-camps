#Python 3 : Functions and Modules
-----------------------

A function is a block of code that performs a specifc task. In this section we
will learn how to utilize available Python functions as well as write our own.

# Sections:

* Writing our own functions

```python
   def initials(line):
      names = line.split()
      initial = ''
      i = 0
      while i < len(names):
         initial=initial+(names[i][0])
         i = i+1
      return initial


   b = 'Thomas Mann'

   print initials(b)
```


* Importing Python modules
