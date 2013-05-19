
#Python 3 : Functions and Modules
-----------------------

A function is a block of code that performs a specifc task. In this section we
will learn how to write our own functions and modules, and import them. The example 
used in this section is a simple one. But the main purpose of this part of 
the tutorial is not only to show you how to write our own functions and modules 
but also why do it. The goal is to show you how and why modularisation of the source code is a good programming practice.

##Writing our own functions

So far we have used the python interactive Python shell and also written Python scripts (typing up the source code in a file and then running it from the Python prompt). We will now show how to "package" this source code into functions. This will allow us to reuse the same code easily multiple times (without having to run it each time manually).

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

NOTE: any variables defined in the function is invisible outside of it. Any results that is created in it must be returned to the outside
using the return statement.

Now instead of rerunning the code for each person, we just had to call our function. Note that our function takes one argument (full_name). 

##Creating modules

First, let's save our function in a file named "namehandler.py". In the same file, after the function, add the following:

```python

if __name__ == "__main__":
    print get_initials("George Adams")
```

Now run this script on the command line, just like it is. You should get the letters GA printed on your screen. Now you know that your function works.

But, how come this code is run?  What happens is that when you run a script directly
from the command line, a variable that is called __name__ is set to have the value __main__. We then ask python
to test on this variable with an if statement, and if it is true, whatever is inside of it is run. In this case, we test the function. This variable
will not be __name__ if we use this script inside of another script.

We are now going to import this function into a different script:

```python
from namehandler import get_initials

f = open('famousauthors.txt')
for l in f:
   print get_initials(l)
f.close()
```

What happens here is that Python goes into the namehandler file and gets the function that we specified, and uses that inside of this script. We get at the function by using the import statement. 

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


Previous: [Flow control](2_Flow_Control.md) Next: [Final exercise](4_Conflict.md)
