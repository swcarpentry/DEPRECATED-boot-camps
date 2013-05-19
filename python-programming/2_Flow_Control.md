# Python 2 : Flow Control - Loops, Conditionals, etc

**Based on Lecture Materials By: Milad Fatenejad and Katy Huff**

*Modified by Karin Lagesen*

Conditionals
============

A conditional (if statement) is some statement that in general says : 
"When some boolean is true, do the following. Otherwise, do this other
thing."

Many equivalence test statements exist in Python that are similar in
other languages:

```python
i=1
j=2
i==j # i is equal to j : FALSE
i<j  # i is less than j : True
i<=j # i is less than or equal to j : TRUE
i>j  # i is greater than j
i>=j # i is greater than or equal to j : FALSE
i!=j # i is not equal to j : TRUE
```

However, python has other equivalence test statements that are fairly
unique to python. To check whether an object is contained in a list :

```python
beatle="John"
beatles=["George", "Ringo","John", "Paul"]
print beatle in beatles # is John one of the beatles? : TRUE
print "Tom" not in beatles # this is also TRUE. 
```

Conditionals (if statements) are also really easy to use in python. Take
a look at the following example:

```python
i = 4
sign = "zero"
if i < 0:
  sign = "negative"
elif i > 0:
  sign = "positive"
else:
  print "Sign must be zero"
  print "Have a nice day"
```

What do you expect the value of sign to be after it has run?

Now type it in, and run it, and see if the value of sign is what
you expected it to be. Note, you have to press enter twice after 
"Have a nice day", then you can show the value of sign by printing it.

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code. For example, in MatLab you begin an if
statement with the word "if" and you end it with "end if". In C/C++ you
delimit blocks with curly braces. Python uses ''indentation'' to delimit
code blocks. The **indentation** above is NOT just to make things look
pretty - it tells Python what the body of the if-statement is. This is
true whenever we create any code blocks, such as the bodies of loops,
functions or classes. This is also why we had to print enter twice after
the if statement - after the first time, python is offering you to type 
in more things that would be executed in the else, so you press enter to 
tell python that you are done, and that it should evaluate the if.

**Short exercise**

Repeat the if statement above, and test out what the value if i has to
be in order to be able to trigger the else clause.


While Loops
===========

Lets start by looking at while loops since they function like while
loops in many other languages. The example below takes a list of
integers and computes the product of each number in the list up to the
-1 element.

A while loop will repeat the instructions within itself until the
conditional that defines it is no longer true.

```python
counter = 0
sequence = [1, 2, 3, 4, 5, 6, 7, 8]
while counter < len(sequence):
  print sequence[counter] * sequence[counter]
  counter += 1
```

Some new syntax has been introduced in this example.

-   On line 1, we define a counter.
-   On line 3 we begin the while loop. The conditional of the while loop
    is in this case whether the counter is less than the length of the
    sequence we are operating on. Since the last element in the list will
    always be one less than the length of the list, this will become false

-   On line 4, we compute the product of the elements just to make this
    more interesting.

-   On the last line we increment the counter, using the += syntax, this
    means that whatever is after the += sign will be added to the value
    that the variable already has.

**Watch Out**

Since a while loop will continue until its conditional is no longer
true, a **poorly formed** while loop might repeat forever. For example :

```python
i=0
mystring="Well, there's egg and bacon, egg and spam, egg bacon and"
while i < len(mystring):
  fields = mystring.split()
  print fields[i]


print "When will this be printed?!" 
```

Tip: if you run this, you can interrupt it with Ctrl-C!


Since the variable **i** never changes within the while loop, we can
expect that the conditional, **i<1** will remain true forever and the
while loop will just go round and round and continue forever.

**Combining while and if**

To create nested while and if loops, the indentation (preferably two or four
spaces) should increase for each looping level. Let's put the following source code in a file and save it as "script.py".

```python
words =	["first", "second#", "third", "f#ourth", "fifth"]
counter = 0
while counter < len(words):
    if "#" not in words:
       print words[counter]
    counter += 1
```


For Loops
=========

For loops in python operate a little differently from other languages.
Lets start with a simple example which prints all of the numbers from 0
to 9:

```python
for i in range(10):
    print i
```

You may be wondering how this works. Start by using help(range) to see
what the range function does.

    Help on built-in function range in module __builtin__:

    range(...)
        range([start,] stop[, step]) -> list of integers

        Return a list containing an arithmetic progression of integers.
        range(i, j) returns [i, i+1, i+2, ..., j-1]; start (!) defaults to 0.
        When step is given, it specifies the increment (or decrement).
        For example, range(4) returns [0, 1, 2, 3].  The end point is omitted!
        These are exactly the valid indices for a list of 4 elements.

Range is a function that returns a list containing a sequence of
integers. So, range(10) returns the list [0,1,2,3,4,5,6,7,8,9]. The for
loop then simply iterates over that list, setting i to each value.

For Loops with Lists and Dictionaries
=====================================

With range, we learned that **for** loops in python are really used to
iterate over sequences of things (they can be used for much more, but
for now this definition will do). Try entering the following to see what
happens:

```python
for c in ["one", 2, "three", 4, "five"]:
    print c
```

this is equivalent to:

```python
c = ["one", 2, "three", 4, "five"]
for i in range(len(c)):
    print c[i]
```

With a list, then, it's clear that we can use the **in** keyword to
indicate a list of things. What about a nested loops around a list of
lists?

```python
italy = ["Rome", "Pisa", "Florence", "Venice", "Trieste"]
argentina = ["Mendoza", "Buenos Aires", "Patagonia"]
india = ["Ahmedabad","Kolkata", "Chennai", "Jaipur", "Surat"]
us = ["Chicago", "Austin", "New York", "San Fran"]
nations = [italy, argentina, india, us]
nationnames = ["italy","argentina", "india", "us"]
for nation in nations :
    print nationnames[nations.index(nation)] + ": "
    for city in nation :
        print "  " + city 
```

Of course, this information is better stored in a dictionary, isn't it?
The data makes more sense if the keys were the nation names and the
values were lists of cities. Importantly, python has given us a tool
specifically for dictionary looping.

The syntax for looping through the keys and values of a dictionary is :

    for key, value in dictionary.iteritems():

Importantly, you don't have to use the words key and value. That's just
what will fill those variables. Here, we rewrite the previous loop using
this clever syntax.

```python
italy = ["Rome", "Pisa", "Florence", "Venice", "Trieste"]
argentina = ["Mendoza", "Buenos Aires", "Patagonia"]
india = ["Ahmedabad","Kolkata", "Chennai", "Jaipur", "Surat"]
us = ["Chicago", "Austin", "New York", "San Fran"]
nations = {"italy":italy, "argentina":argentina, "india":india, "us":us}
for nation, cities in nations.iteritems() :
    print nation + " : "
    for city in cities :
        print "  " + city 
```

break and continue
=========================

A break statement cuts off a loop from within an inner loop. It helps
avoid infinite loops by cutting off loops when they're clearly going
nowhere.

```python
reasonable = 10
for n in range(1,2000):
    if n == reasonable:
        break
    print n
```

Something you might want to do instead of breaking is to continue to the
next iteration of a loop, giving up on the current one..

```python
reasonable = 10
for n in range(1,2000):
    if n == reasonable:
      continue
    print n
```

What is the difference between the output of these two?

Final Example
=============

We've seen a lot so far. Lets work through a slightly lengthier example
together. I'll use some of the concepts we already saw and introduce a
few new concepts. To run the example, you'll need to locate a short file
containing phone numbers. The file is in your repository is called phonenums.txt.

This example opens a text file containing a list of phone numbers. The
phone numbers are in the format \#\#\#\#-\#\#\#-\#\#\#, one to a line.
The example code loops through each line in the file and counts the
number of times each area code appears. The answer is stored in a
dictionary, where the area code is the key and the number of times it
occurs is the value.

```python
areacodes = {} # Create an empty dictionary
f = open("phonenums.txt") # Open the text file
for line in f: # iterate through the text file, one line at a time (think of the file as a list of lines)
    ac = line.split('-')[0] # Split phone number, first element is the area code
    if not ac in areacodes: # Check if it is already in the dictionary
      areacodes[ac] = 1 # If not, add it to the dictionary
    else:
      areacodes[ac] += 1 # Add one to the dictionary entry

print areacodes # Print the answer
```

**Short exercise**

- Modify the script so that it prints the results nicely, with first the area code, then a 
  tab, and then the number of times it occurs.

- Modify this script so that it completely ignores phone numbers which contain the number 4.



Previous: [Data types](1_Data_Types.md) Next: [Functions and modules](3_Functions_and_Modules.md)
