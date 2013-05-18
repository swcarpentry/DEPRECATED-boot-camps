
#Python 3 : Functions and Modules
-----------------------

A function is a block of code that performs a specifc task. In this section we
will learn how to write our own functions and modules, and import them. The example used in this section is a simple one. But the main purpose of this part of the tutorial is not only to show you how to write our own functions and modules but also why do it. The goal is to show you how and why modularisation of the source code is a good programming practice.

##Writing our own functions

So far we have used Python prompt and write Python scripts (typing up the source code in a file and then running it from the Python prompt). We will now show how to "package" this source code into functions. This will allow us to reuse the same code easily multiple times (without having to run it each time manually).

Let's say we want to extract people's initials from their full name. The code below should do it.

```python
def get_initials(full_name):
   names = full_name.split() 
   initial = ''
   i = 0
   while i < len(names):
      initial=initial+(names[i][0]) #From each part of the name we get the first letter
      i = i+1
   return initial
   
a = 'Thomas Mann'
b = 'James Bond'
print get_initials(a)
print get_initials(b)
print get_initials('Anna Karenina')
```

The keyword for defining a function is def. After that we have the function name followed by any arguments in paranthesis. After that comes
the code that the function performs. Input to the function is available in the arguments passed to it.

Now instead of rerunning the code for each person, we just had to call our function. Note that our function takes one argument (full_name). 

We can now use our function with the data file [famousauthors.txt](famousauthors.txt).

##Importing Python modules

First, let's save our function in a file named "namehandler.py" (deleting the last 5 lines - beginning from b='Thomas Mann'). Now, we can import our function from our own module:

```python
from namehandler import get_initials

f = open('famousauthors.txt')
for l in f:
   print get_initials(l)
f.close()
```

In some Python code you may notice the following statements `from module import *`. This means that all names defined in the module are imported. However, this is not a good practice as the code becomes more difficult to read. 


**Getting input from the command line**

What if we want to process different files and provide the name of the input file on the command line? We need to pass the file name as one of the Python arguments, and to do that, we have to access a special list that is called sys.argv. This list is available from the sys module, so we have to import that one. The list sys.argv contains everything that is put in on the command line, therefore the first element in this list is actually the name of the script. Anything after the name gets put in position 1, 2 and so forth of this list. Let's try to specify the input to the previous program on the command line:

```python
import sys
from namehandler import get_initials

def main():
  f = open(sys.argv[1])
  for l in f:
     print get_initials(l)
  f.close()

if __name__ == "__main__":
   main()
```

So, the input file is to be found in sys.argv[1]. 

**

Most of the work in this program seems to be done in a function named main. But- how does this one get executed? In any script file python will look for the "if __name__ == "__main__":" line and runs it. If it evaluates to true, i.e. __name__ has the value __main__, it knows that this script is actually being run directly. If we were to use this file as a module somewhere else, this if would not be true. So, in this case the if statement is true, so main() is executed.


##Exercise

You probably noticed that there are two authors whose initials were not correctly extracted. It's Erich von Däniken and Gérard de Villiers. How can we modify our function to address this problem? Hint: check the methods available for string variables ( dir(str) ).


####Solution:

```python
def get_initials(line):
    names = line.split()
    initial = ''
    i = 0
    while i < len(names):
       if not names[i][0].islower():
           initial=initial+(names[i][0])
           i = i+1
       else:
           i = i+1
    return initial
```

Previous: [Flow control](2_Flow_Control.md)
