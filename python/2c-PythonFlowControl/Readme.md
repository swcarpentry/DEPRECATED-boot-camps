# Flow Control

[Back To Python Data
Structures](../2b-PythonDataStructures)
- [Forward To Python Functions and
Modules](../2d-PythonFunctionsAndModules)

* * * * *

**Presented By : Shoaib Sufi**

**Based on Lecture Materials By: Milad Fatenejad (mostly)**

So far we've learned about various data types in python and how to rudimentarily read data out of a file. Now we will learn how to do logical operations with conditional blocks (if statements) and how to iterate using for and while loops. We will also learn to write data to a file.

## For Loops

Iteration is a very useful concept, and it appears all the time in python.

You can use a for loop to print out all the contents of a list using python's compact syntax for iteration.

```python
In [1]: numberList=[1,2,4,5,6]

In [2]: for item in numberList:
   ...:     print item
   ...:     
1
2
4
5
6
```

Using the enumerate() function we can easily print the array index and the value at that index.

```python
In [3]: for item_index, item_value in enumerate(numberList):
   ...:     print item_index, item_value
   ...:     
0 1
1 2
2 4
3 5
4 6
```

Lets now consider our current and voltage data again.


```python
In [1]: voltageList = [-2.0, -1.0, 0.0, 1.0, 2.0]

In [2]: currentList = [-1.0, -0.5, 0.0, 0.5, 1.0]
```

We know that power is the product of current and voltage. We already have the current as a function of voltage, but what if we wanted to generate a list of values of power coresponding to the value of voltage? We use a for loop.

```python
In [3]: powerList = [] # Initialize an empty list to be filled later.
In [4]: indxList = range(len(currentList))

In [5]: for indx in indxList:
   ....:    power = voltageList[indx] * currentList[indx]
   ....:    powerList.append(power)
``` 

There's a lot of stuff that just happened here. Lets break it down.

### Indentation

The indentation is a feature of python that some people hate. Some other programming languages use brackets to denote a command block. Python uses indentation. The amount of indentation doesn't matter, so long as everything in the same block is indented the same amount.

The for loop itself is pretty simple: we take a list (in this case), pull out the presently indexed value, and execute the block below the for command. Once the block has been executed, the for loop increments to the next index and keeps going to the end.

The range command gives us an incremental list of values. We are using it to generate our lists of indices for our for loop iteration.

## Conditional (if) Statements

There are many situations where you want to execute some code based on some condition, say you want a message to print if the value of current is higher than a preset value. In this case, you use a conditional statement. There are many different conditions, here are some of them:

```python
i=1
j=2
i==j # i is equal to j : False
i<j  # i is less than j : True
i<=j # i is less than or equal to j : True
i>j  # i is greater than j : False
i>=j # i is greater than or equal to j : False
i!=j # i is not equal to j : True
```

Consider the following example.

```python
In [11]: for current in currentList:
   ....:     if current > 0.7:
   ....:         print "high"
   ....:     elif current < -0.7:
   ....:         print "low"
   ....:     else:
   ....:         print "ok"
   ....:         
   ....:         
low
ok
ok
ok
high
```

In this example, we are iterating over the list of values of current. We've nested the conditional statement within the for loop statement, and we are looking at each value of current to make a decision. The lowest value causes the code to display "low". The three middle values are "ok" and the last value is "high".


## While Loops

Lets start by looking at while loops since they function like while
loops in many other languages. The example below takes a list of
integers and computes the product of each number in the list up to the
-1 element.

A while loop will repeat the instructions within itself until the
conditional that defines it is no longer true.

```python
mult = 1
sequence = [1, 5, 7, 9, 3, -1, 5, 3]
while sequence[0] != -1:
    mult = mult * sequence[0]
    del sequence[0]

print mult
```

Some new syntax has been introduced in this example.

-   On line 3 We begin the while loop. This while loop will execute
    until sequence[0] == -1 . That is, until it deletes all of the
    entries of the sequence that come before -1.

-   On line 4, we compute the product of the elements just to make this
    more interesting.

-   On line 5, we use the \`del\` keyword to remove the first element of
    the list, shifting every element down one.

**Watch Out**

Since a while loop will continue until its conditional is no longer
true, a **poorly formed** while loop might repeat forever. For example :

```python
i=1
print "Well, there's egg and chips, egg and fries, eggy chips and"
while i == 1:
  print "fries "
print "or Lobster Thermidor a Crevette with a mornay sauce served in a Provencale manner with shallots..." 
```

Since the variable **i** never changes within the while loop, we can
expect that the conditional, **i=1** will remain true forever and the
while loop will just go round and round, as if this restaurant offered
nothing but fries. (If you try this at home, please note that one way to
interrupt a non-terminating process is **ctrl+c** or **ctrl+z**.

## Writing to a File
We could have gone over this topic in the previous lesson, but iteration makes writing data to files a lot less ugly. 

Recall our file data.dat. There is no mention of the name of the user in that file, and I want to insert my own name. I want my name to appear on the second line, immediately below the experiment type. Here's how we do that.

```python
# Create an empty list to contain all of the file data.
dataList = []

# Open the file and iterate over it to read all of the data into the list.
f = open("data.dat")

for line in f:
  dataList.append(line)

f.close()

# Insert my name into the list.
dataList.insert(1, "user: Shoaib Sufi\n")

# Open a new file for writing and write all of the new lines to it.
f = open("newdat.dat", "w")
for line in dataList:
  f.write(line)

f.close()
```

You have to have that newline character or else you're going to get a weird result.

# break, continue, and else

A break statement cuts off a loop from within an inner loop. It helps
avoid infinite loops by cutting off loops when they're clearly going
nowhere.

```python
reasonable = 10
for n in range(1,15):
    if n == reasonable :
        break
    print n
```

Something you might want to do instead of breaking is to continue to the
next iteration of a loop, giving up on the current one..

```python
reasonable = 10
for n in range(1,15):
    if n == reasonable :
      continue
    print n
```

What is the difference between the output of these two?

Importantly, Python allows you to use an else statement in a for loop.

That is :

```python
knights={"Sir Belvedere":"the Wise", "Sir Lancelot":"the Brave", \
        "Sir Galahad":"the Pure", "Sir Robin":"the Brave", "The Black Knight":"John Clease"} 

favorites=knights.keys()
favorites.remove("Sir Robin")
for name, title in knights.iteritems() : 
    string = name + ", "
    for fav in favorites :
        if fav == name :
            string += title
            break
    else:
        string += title + ", but not quite so brave as Sir Lancelot." 
    print string
```

##Exercise##

###Context###
We've seen a lot so far. Lets work through a slightly lengthier example
together. I'll use some of the concepts we already saw and introduce a
few new concepts. To run the example, you'll need to locate a short file
containing phone numbers. The file can be found in your 
repository within the python/2c-PythonFlowControl directory and is called phonenums.txt.
Now we have to move ipython to that directory so it can find the
phonenums.txt file. You navigate within ipython in the same way that you
navigate in the shell, by entering "cd [path]" or you could restart ipython in that directory.

This example opens a text file containing a list of phone numbers. The
phone numbers are in the format \#\#\#-\#\#\#-\#\#\#\#, one to a line.
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

##Iteritems##

Use the iteritems dictionary method in combination with a for loop to
print the keys/values of the areacodes dictionary one to a line. In
other words, the goal is to write a loop that prints:

    203 4
    800 4
    608 8
    773 3

This exercise is a little tricky to figure out, but give it a try.
