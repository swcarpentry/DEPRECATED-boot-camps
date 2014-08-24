[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)


# Write Code for People: Variables, Data Structures and Conditionals

* * * * *

This lecture covers a number of topics under the heading of writing code for people:
* naming variables in meaningful ways
* choosing appropriate data types
* writing effective comments
* don't repeat yourself
* breaking your code/script into small pieces

These concepts are important and apply to any language.  For today's lecture,
we'll use python as the language, but you use a different language, you should
think about how you would use these practices in your language.  In
particular, we will use the [iPython interpreter](../ipython/Readme.md), a
version of python with some added features.

At the same time, this exercise will demonstrate one way to think about
developing and designing code/scripts in a modular way.

This is not designed to teach python, *per se*, but use it as a relatively
accessible language with features that help shine a light on best practices.
For more on python, consider the following reference material.

* [Dive into Python](http://www.diveintopython.net/toc/index.html)
* [Software Carpentry's Python Lectures](http://software-carpentry.org/4_0/python/)
* [IPython: A System for Interactive Scientific Computing](http://dx.doi.org/10.1109/MCSE.2007.53)
* [How to Think Like a Computer Scientist](http://www.greenteapress.com/thinkpython/thinkpython.html)

## Motivating Example

For this lecture, we will return to the simulated cochlear implant data that
we used in the lecture on using the shell.  This data was in many different
files, each distributed over many different directories.  One of the exercises
was to place all of those files in a single directory named `cleandata` and
give them all the same extension.

However, all those data in different files does not make it easy for us use
for additional analysis.  It would be easier if that data was all in a single
file.  We could do that by simply concatenating all the files to make one
single file, but this will make it less obvious which pieces of data belong
together.

Therefore, for this exercise, we will create a single file with all this data,
but with one line per file, with all the data from that file.  The different
pieces of data will appear in the same order in each line, and will be
separated by commas.  This is a so-called *comma-separated values* (or CSV)
file and is a common simple way to represent tables of data.  In fact, this
type of file can be loaded into most spreadsheet software.

We chose this example because it represents a common task: amalgamating data
from a number of files into a single table or database.  By choosing a CSV
file, it may resonate even more strongly since it means you can use such a
file in a spreadsheet.  Imagine the effort of making a spreadsheet from this
set of more than 350 files.  Although spreadsheets are frequently NOT the best
way to process data, they are probably something with which you are familiar.

**Pro-tip:** [Using two windows](using_multiple_windows.md) will make it much
  easier to follow this exercise.

## Step 1: Start at the top

Let's start with the big picture view of what this script needs to do:

* loop through a list of filenames
    * parse each filename into a single data record
    * save that record
* write a header in the correct CSV syntax
* loop through each record
    * write that record in the correct CSV syntax

```python
for filename in filelist:
    data_record = parseFile(filename)
    all_data.append(data_record)

csv_separator = ','

writeCSVHeader(column_labels,csv_separator)

for data_record in all_data:
    writeCSVRow(column_labels,data_record,csv_separator)
```

Let's review the best practices so far:

### Modular development and design

These 8 lines (including 2 blank lines) represent all the main steps that our
script must accomplish.  We have only 6 variables and 4 functions to track in
our heads.

### Writing effective comments

These lines don't really need comments since we have chosen expressive
variable and function names in a short script that performs clear operations.

### Variable and function naming

We have chosen meaningful names for each variable: `filename`, `filelist`,
`data_record`, `all_data`, `column_labels`, `csv_separator`.  It is often
recommended that variables be given names that are nouns, indicating that they
represent things rather than actions.

The length of variable names is also important: too short and they loose some
meaning; too long and they become prone to errors when typing them.

**Short Exercise**

Which of these are bad choices and why?  Which are reasonable alternatives?
* `fn`
* `all_data_from_all_the_files_in_the_list`
* `col_lbls`
* `record`
* `sep`
* `col_labels` 

We have chosen meaningful names for each function: parseFile, writeCVSHeader,
writeCVSRow.  It is often recommended that functions be given names that are
verbs, indiciating that they represent action rather than things.

We have used a consistent style that helps the reader distinguish between
variables, in which words are separated by underscores, and functions, in
which CamelCase is used.

### Choosing appropriate data types

At this point we haven't made very few decisions on data types:

* the `filelist` variable should be some type that allows us to take advantage
  of python's ability to iterate over the members
* the `data_record` variable should be some type that stores all the
  information from a single file
* the `all_data` variable should be some type that allows us to use the
  `append()` method
* the `column_labels` variable should store all the column labels for our CSV
  file

Now is a good time to explore the different
[compound data types in python](compound_data_types.md).  Most languages have
some compound data types and you should learn about them and their features to
make the best choices in your scripts/code in that langauge.

For our problem the following choices are probably wise:

* `filelist` should be a list
  * lists are a common default unless you need other features like key-value
    pairs (dictionaries), immutability (tuples), or uniqueness (sets)
* `data_record` should be a dictionary
  * each record will store the data that is already in key-value pairs in the file
* `all_data` should be a list of those dictionaries
* `column_labels` should be a tuple because we want it to be immutable

### Don't repeat yourself

At this point there are two examples of using this best practice.  In it's
simplest form, we have used loops so that the computer repeats the same task
for many different items. In addition, using the `csv_separator` variable
instead of just typing `','` has two advantages:
1. it provides context for what this means
2. it means we can change it in a single place

### Initialize some of these variables

1. We can initialize the list of files from the command-line arguments.
2. We should also declare the immutable ordered list of column labels.
3. Now that we know `all_data` should be a list, we should initialize it as a
   empty list.

```python
import sys
filelist = sys.argv[1:]
column_labels = ("Subject","Reported","Year/month of birth",
                "Sex","CI type","Volume","Range","Discrimination")
all_data = []
```

## Step 2: Parse a single file

To continue with our script we can now focus on the first funciton,
`parseFile`.  Python functions have the same components of functions in most
languages.

1. a way to declare where the function starts: `def`
2. a function name
3. arguments to the function
4. a way to return some result(s) from a function
5. a way to declare where the function ends: indentation

For our function, we already know that the name is `parseFile` and that it
will take a single argument that is a `filename`.

Just like the last time, we'll consider what the big picture steps are:

* open the file
* loop over all the lines in the file
  * extract the data from that line in a key/value pair
  * add that data to the dictionary
* return the complete dictionary

```python
def parseFile(filename):
    """Read all the lines from a file and return them as a dictionary of key/value pairs"""
    textfile = open(filename, 'r')
    data_record = {}
    for line in textfile:
        # all lines that contain key/value pairs must have at least on colon
        if ':' in line:
            (key,value) = extractData(line)
            data_record[key] = value
    return data_record
```

### Reading from files

At this point we've had to introduce a new concept: reading data from files.
For most purposes, python makes this very easy.  When we `open` a file, we get
an object over which we can iterate, much like a list.  Let's look at the file
`phonenums.txt` using this concept:

**Try it exercise**

```python
In [15]: file_data = open('phonenums.txt')
In [16]: for line in file_data:
   ....:     print line
```

Again we'll review the best practices.
 
### Modular development and design

Once again we have a short method with only 10 lines, including 2 comments,
and only 6 variables.  To keep this method focused, we will rely on another
new function.

### Writing effective comments

In this case there are 2 places where comments are useful.  Every function in
python should have a special comment at the beginning known as a `docstring`.
A docstring has special meaning to tools that, for example, can show help
about a function.

This function also requires a special test to ensure that the line is actually
a data line, using the assumption that all data lines have at least one colon
(:).  This is a useful place for a comment because the conditional clause is
not obvious without it.

### Variable and function naming

Two of the variables used in this function use the same name and same meaning
as the corresponding variables in the main script: `filename` and
`data_record`.  Although this is not necessary, it makes code more readable
when variables in different places that have the same meaning can also have
the same name.  In functions that broader utility, this may be neither
possible nor logical.

The other variables have names that are obvious matches to their purposes:
`textfile`, `line`, `key`, and `value`.  These are all nouns and the only
method that we have declared is a verb: `extractData`.

### Choosing appropriate data types

One of the features of python that we may start to notice is that it relies on
something called [dynamic typing](dynamic_types.md).  Although the type of a
variable can change over its lifetime, too much change can make the script
more difficult to follow.  It is helpful to know what type is implied and
endeavor to use that consistently.

In addition, the type of `textfile` is a built-in object of type `file`,
determined dynamically by the object returned by the call to the `open`
method.

Because `data_record` is a new dictionary each time we call this function, we
need to initialize it as a new dictionary each time.

## Step 3: Extracting the data from each line

The next function that we will add will operate on a small amount of data to
separate a single line of the form:

```
Some keyword phrase: data that is the value
```

into two variables where:

```
key = "Some keyword phrase"
value = "data that is the value"
```

The line of data is a string, so we'll be using some methods of the built-in
string class.  There are many ways to accomplish this task, but we'll rely on
the `split()` method of the string class with the following high level algorithm:

* split the line at the location of the separator into a list of strings
* use the first of this list as the key
* use the remainder of the list as the value

```python
def extractData(line):
    """Assume the first colon (:) separates the key from the value."""
    separator = ':'
    line_data = line.strip().split(separator)
    key = line_data[0]
    # the value may contain a separator so will need to be reassembled from
    # possibly multiple elements
    value = separator.join(line_data[1:])
    return key,value.strip()
```

### Getting help on new functions

One of the features of iPython is it's ability to show you help of a function,
typically in the form of its docstring.  You can learn more about the
`strip()` and `split()` string functions:

```
In [21]: str.strip?
In [22]: str.split?
```

Let's again review the best practices:

### Modular development and design

Not much new to say here: a short method: 9 lines including 3 comments, and
only 5 variables.  At this point, we have used 12 different variables
throughout the script, but no single block of the script uses more than 6 of
them.

### Writing effective comments

In addition to a docstring, we have a comment to remind us that when the
separator happens to appear as part of the value, the `line_data` list will
have more than just 2 elements, and we need to reassemble all but the first
one into the value.

### Variable and function naming

Of the 5 variables used in this function, 3 are based on names defined and
discussed in previous functions.  The `separator` variable is an appropriate
length and expressive in describing its role.  The list `line_data` is a
descriptive name, but other names would be fine.

**Think aloud:** Recommend some alternative variable names for `line_data`.

### Choosing appropriate data types

All the choices in this case were made for us by the requirements of the
string functions we call.


## Step 4: Write the CSV data

**Discussion Exercise** Given the following implementation of the final two
  functions, discuss ways in which they follow the best practices, and ways
  that they could be improved, if any.

```python
def writeCSVHeader(column_labels,csv_separator):
    header = []
    for column in column_labels:
        header.append('"' + column + '"')
    print csv_separator.join(header)

def writeCSVRow(column_labels,data_record,csv_separator):
    row = []
    for column in column_labels:
        row.append('"' + data_record[column] + '"')
    print csv_separator.join(row)
```




Once we briefly deal with iPython, I'll cover python in the following order:

## What We'll Cover

We'll focus on two overarching concepts that are important to any programming language: 

1) How to write code for people. That is code that is readable and
understandable to others in your group and most importantly to your future
self 3 months or 3 years down the road.

2) How to not repeat yourself. How to reuse your code with loops and
functions. And how to eventually build modules, collections of functions, you
and others can use in all your codes. And how to use other people's modules.

### iPython Intro
* iPython

### Lesson 1 (Write Code for People)
* Print statements
* Variables
* Integers
* Floats
* Strings
* Types
* Type Coercion
* Basic Operations: add numbers, concatenate strings, basic data type functionality
* List
* Dictionary 
* Tuple
* Conditional (if) statements

### Lesson 2 (Don't Repeat Yourself)
* For Loops
* While Loops
* Iteration
* Methods
* Reading & Writing Files
* Modules


## iPython

Please follow the link to the [iPython Intro](../ipython/Readme.md).

## Back to Write Code for People

This lesson will introduce the basics of the python programming language while stressing how to make readable code, code for people. As we introduce variables we will discuss how to name them, when to comment and which comments are useful. As we introduce types we will discuss how type choice can be influenced by code readability considerations. The above readability considerations are to a certain degree universal to all programming, regardless of language. Finally, we introduce conditionals, such as `if` statements; we will see what an important role white space plays in python and how it can improve readability.

But first: let's introduce you to python.

## Variables


All programming languages have variables, and python is no different. To create a variable, just name it and set it with the equals sign. One important caveat: variable names can only contain letters, numbers, and the underscore character. Let's set a variable.

```
In [1]: first_name = "Cliff"
In [2]: last_name = "Rodgers"
```

We can glue strings together with the + operator. Computers don't understand context.

```
In [3]: full_name = first_name + last_name
In [4]: print(full_name)
Out[4]: 'CliffRodgers'
```

***Exercise***
Can you add the extra space between my last and first name?

## Naming Your Variables

Which of these variable names are good?

```
v = 4
voltage = 4
TheVoltageInTheCircuitAtPointA = 4 # Camel Case
the_voltage_is = 4 # Pothole
monkey = 4
```
Variable names should be:

* Meaningful (to those who are going to read the code)
* Short enough so you don't misstype them

Variable name choice is important; a well-named variable is self-explanatory without comments and will make your code easier to read as the reader will not have to look up the comments. Remember that context matters a lot: in some cases, you'll want to spell out `voltage` or even `input_voltage`. In other cases, `v` is a shorthand that everyone will understand.

It's also important to choose a naming convention for your whole project, and get the people working on it to agree. Mixing `batteryVoltage` and `capacitor_value` in one place makes code hard to read.

## Writing Comments For People

Above, you may have noticed the `#` character, which denotes a comment in python. Comments should describe meaning but not what the statement is doing.

```
In [15]: voltage = 4 # set the voltage to 4   <- well, duh

In [16]: voltage = 4 # Input to the circuit   <- better; says what this voltage means
```

## Types and Dynamic Typing

Like in most programming languages, things in python are *typed* — the *type* refers to the type of data and what you can do with it. For example, you can do different things with strings and numbers. Numbers can have decimal components or not, and so on. You can see the type of a variable with the `type` command.

```
In [5]: type?
Type:       type
String Form:<type 'type'>
Namespace:  Python builtin
Docstring:
type(object) -> the object's type
type(name, bases, dict) -> a new type

In [6]: type(full_name)
Out[6]: str
In [7]: type(10)
Out[7]: int

```


Python is what is known as a *dynamically typed* language. Dynamic typing means that you don't have to declare the type of a variable when you define it; python just figures it out based on how you are setting the variable. This is in contrast to *statically typed* languages, where you must say up front that a variable is going to be used for strings or numbers or whatever. There are good and bad points to both approaches.

But back to Python. Let's say you set a variable. Sometime later you can just change the type of data assigned to a variable and python is perfectly happy about that. Since it won't be obvious until (possibly much) later why that's important, I'll let you marinate on that idea for a second.

Here's an example of dynamic typing:

```
In [8]: voltage = 2
In [9]: print(voltage)
2

In [10]: type(voltage)
Out[10]: int
```

It's an `int`, which is an integer — a number with no decimal component. Let's assign a value of 2.7 (which has a decimal part) to voltage. What happens to the type?

```
In [11]: voltage = 2.7

In [12]: type(voltage)
Out[12]: float
```

Neat! That's a `float`, which does have a decimal part. You can assign a string to the variable voltage in the same way:

```python
In [13]: voltage = "2.7 volts"

In [14]: type(voltage)
Out[14]: str
```

I'll let you ruminate on the pros and cons of this construction while I change the value of voltage back to an int:

```
In [15]: voltage = 2
```

Choosing an appropriate variable type is not just a practical concern; it can also have an effect on code readability. Is this number used for calculations or only in print statements?


## On Being Precise With floats and ints

Again, the following may seem esoteric and pedantic, but it is very important. So bear with me.

Lets say you had some voltage data that looks like the following

```
0
0.5
1
1.5
2
```

If you just assigned this data individually to a variable, you'd end up with the following types:

```
0   -> int
0.5 -> float
1   -> int
1.5 -> float
2   -> int
```

But what if you wanted all of that data to be floats on its way in? You could assign the variable and then coerce it to type float:

```
In [28]: voltage = float(1)
```

But that's ugly. If you want whats otherwise an integer to be a float, add `.0` to the end:

```python
In [29]: voltage = 1.0

In [30]: type(voltage)
Out[30]: float
```

This point becomes important when we start operating on data in the next section.

## Data Operations

In this section all of the discussion in the previous section becomes
important. I don't know if I'd call this stuff fundamental to the language,
but it's pretty important and it will zing you if you aren't careful. The
takeaway is that you need to be precise with what you are doing. Let's say you
want to add some integers.

```
In [31]: a = 1

In [32]: b = 2

In [33]: c = a + b

In [34]: c
Out[34]: 3

In [38]: type(a), type(b), type(c)
Out[38]: (int, int, int)
```

So we got a value of three for the sum, which also happens to be an
integer. Any operation between two integers is another integer. Makes sense.

So what about the case where a is an integer and b is a float?

```
In [39]: a = 1

In [40]: b = 2.0

In [41]: c = a + b

In [42]: c
Out[42]: 3.0

In [43]: type(a), type(b), type(c)
Out[43]: (int, float, float)
```

You can do multiplication on numbers as well.

```
In [44]: a = 2

In [45]: b = 3

In [46]: c = a * b

In [47]: c
Out[47]: 6

In [48]: type(a), type(b), type(c)
Out[48]: (int, int, int)
```

Also division.

```
In [49]: a = 1

In [50]: b = 2

In [51]: c = a / b

In [52]: c
Out[52]: 0
```

**ZING!**

Here's why type is important. Dividing two integers returns an integer: this operation calculates the quotient and floors the result to get the answer.

If everything was a float, the division is what you would expect.

```
In [53]: a = 1.0

In [54]: b = 2.0

In [55]: c = a / b

In [56]: c
Out[56]: 0.5

In [57]: type(a), type(b), type(c)
Out[57]: (float, float, float)
```

## Compound Data Types: Lists, Dictionaries, Sets, and Tuples

Python would be a fairly useless language if it weren't for the compound
data types. The main two are lists and dictionaries, but I'll mention sets
and tuples as well. 

## Lists

A list is an ordered, indexable collection of data. Let's say you have
collected some current and voltage data that looks like this:

```
voltage:
-2.0
-1.0
0.0
1.0
2.0

current:
-1.0
-0.5
0.0
0.5
1.0
```

So you could put that data into lists like

```
In [1]: voltage_list = [-2.0, -1.0, 0.0, 1.0, 2.0]

In [2]: current_list = [-1.0, -0.5, 0.0, 0.5, 1.0]
```

obviously voltage_list is of type list:

```
In [3]: type(voltage_list)
Out[3]: list
```

Python lists have the charming (annoying?) feature that they are indexed
from zero. Therefore, to find the value of the first item in voltage_list:

```
In [4]: voltage_list[0]
Out[4]: -2.0
```

And to find the value of the third item

```
In [5]: voltage_list[2]
Out[5]: 0.0
```

Lists can be indexed from the back using a negative index. The last item of
current_list

```
In [6]: current_list[-1]
Out[6]: 1.0
```

and the next-to-last

```
In [7]: current_list[-2]
Out[7]: 0.5
```

You can "slice" items from within a list. Let's say we wanted the second
through fourth items from voltage_list

```
In [8]: voltage_list[1:4]
Out[8]: [-1.0, 0.0, 1.0]
```

Or from the third item to the end

```
In [9]: voltage_list[2:]
Out[9]: [0.0, 1.0, 2.0]
```

and so on.

***Exercise***
What does voltage_list[::2] mean?

### Append and Extend

Just like strings have methods, lists do too.

```
In [10] dir(list)
```

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

```
voltage:
3.0
4.0

current:
1.5
2.0
```

If you want to append items to the end of a list, use the append method.

```
In [11]: voltage_list.append(3.)

In [12]: voltage_list.append(4.)

In [13]: voltage_list
Out[13]: [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
```

You can see how that approach might be tedious in certain cases. If you
want to concatenate a list onto the end of another one, use extend.

```
In [14]: current_list.extend([1.5, 2.0])

In [15]: current_list
Out[15]: [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
```

### Length of Lists

Sometimes you want to know how many items are in a list. Use the len command.

```
In [16]: len(voltage_list)
Out[16]: 7
```

### Heterogeneous Data

Lists can contain hetergeneous data.

```
In [17]: data_list = ["experiment: current vs. voltage",
                      "run", 47,
                      "temperature", 372.756,
                      "current", [-1.0, -0.5, 0.0, 0.5, 1.0],
                      "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]

```

We've got strings, ints, floats, and even other lists in there. While this is a perfectly valid thing to do in a list, it's often not a good idea. Generally, you want to be able to run the same operation on every element of a list. If you want to put different kinds of things in a list structure, you often want to use something called a *tuple,* which we'll talk about in a minute.

## Assigning Variables to Other Variables

Something that might cause you headaches in the future is how python deals
with assignment of one variable to another. When you set a variable equal
to another, both variables point to the same thing. Changing the first one
ends up changing the second. Be careful about this fact.

```
In [19]: a = [1,2]

In [20]: b = a

In [21]: a.append(10)

In [22]: b
Out[22]: [1, 2, 10]
```

If you want to see if that's what's going on in your case, python has a special `is` operation that will test to see if two variables point to the same thing:

```
In [23]: a is b
Out[23]: True
In [24]: c = [1, 2, 10]
In [25]: a is c
Out[25]: False
In [26]: a == c
Out[26]: True
```

There's a ton more to know about lists, but let's press on. [Dive into Python](http://www.diveintopython.net/toc/index.html) or the help documentation for more info.

## Tuples

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is immutable; once you set the
data in it, the tuple cannot be changed. You define a tuple as follows.

```
In [1]: tup = ("red", "white", "blue")

In [2]: type(tup)
Out[2]: tuple
```

You can slice and index the tuple exactly like you would a list.

### Why Use Tuples

Tuples can emphasize intent in two ways. First, you can't easily change the items in a tuple, so if you want to specify that a list shouldn't change, a tuple is a great way to indicate that. Second, in a tuple, the order of the items is usually more significant: you might do something like:

```
In [1]: person_data = ('Nate', 35, 'njvack@wisc.edu')
In [2]: name = person_data[0]
```

... though in practice, you'll more often use a `dictionary` for something like this.

***Exercise***
Display the second element of a tuple with two different slices.

## Dictionaries

Recall our variable data_list which contained our current-voltage data
and also some metadata. We were able to store the data as a list, but
clearly the list type is not the optimal choice for a data model. The
dictionary is a much better choice.

A python dictionary is a collection of key, value pairs. The key is a
way to name the data, and the value is the data itself.  Here's a way
to create a dictionary that contains all the data in our data.dat file
in a more sensible way than a list.

```
In [7] data_dict = {"experiment": "current vs. voltage", \
                   "run": 47, \
                   "temperature": 372.756, \
                   "current": [-1.0, -0.5, 0.0, 0.5, 1.0], \
                   "voltage": [-2.0, -1.0, 0.0, 1.0, 2.0]}
```

This model is clearly better because you no longer have to remember that
the run number is in the second position of the list, you just refer
directly to "run":

```
In [9]: data_dict["run"]
Out[9]: 47
```

If you wanted the voltage data list:

```
In [10]: data_dict["voltage"]
Out[10]: [-2.0, -1.0, 0.0, 1.0, 2.0]
```

Or perhaps you wanted the last element of the current data list

```
In [11]: data_dict["current"][-1]
Out[11]: 1.0
```

Once a dictionary has been created, you can change the values of the data
if you like.

```
In [12]: data_dict["temperature"] = 3275.39
```

You can also add new keys to the dictionary.

```
In [13]: data_dict["user"] = "F. C. Rodgers"
```

Dictionaries, like strings, lists, and all the rest, have built-in methods.
Let's say you wanted all the keys from a particular dictionary.

```
In [14]: data_dict.keys()
Out[14]: ['run', 'temperature', 'current', 'experiment', 'user', 'voltage']
```

also, values

```
In [15]: data_dict.values()
Out[15]: 
[47,
 3275.39,
 [-1.0, -0.5, 0.0, 0.5, 1.0],
 'current vs. voltage',
 'F. C. Rodgers',
 [-2.0, -1.0, 0.0, 1.0, 2.0]]
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

Since tuples are immutable, they can be used as keys for dictionaries.
Lists are mutable, and therefore cannot.

When you architect software in python, most data will end up looking either
like a list or a dictionary. These two data types are very important in
python and you'll end up using them all the time.

## Sets

Most introductory python courses do not go over sets this early (or at
all), and in the interest of time we're no different. The python set
type is a useful data type similar to the idea of a mathematical set:
it is an unordered collection of unique things.

Consider the following examples if you're interested in the useful
sorts of things you can do with python sets:

```
In [3] fruit = set(["apple", "banana", "pear", "banana"]) #You have to use a list to create a set.
```

Since sets contain only unique items, there's only one banana in the set
fruit.

You can do things like intersections, unions on sets just like in
math. Here's an example of an intersection of two sets (the common items in
both sets).

```
In [4]: first_bowl = set(["apple", "banana", "pear", "peach"])

In [5]: second_bowl = set(["peach", "watermelon", "orange", "apple"])

In [6]: set.intersection(first_bowl, second_bowl)
Out[6]: set(['apple', 'peach'])
```

You can check out more info using the help docs. The important thing here isn't that sets exist or what exactly they can do, but the main concept: different types of data can do different things easily. Choosing your data types well will both make your code clearer and simpler to write.

## Conditionals

A conditional (`if` statement) is some statement that in general says :
"When some boolean is true, do the following. Elsewise, do this other
thing."

Many equivalence test statements exist in Python that are similar in
other languages:

```python
i = 1
j = 2
i == j  # i is equal to j : False
i < j   # i is less than j : True
i <= j  # i is less than or equal to j : True
i > j   # i is greater than j : False
i >= j  # i is greater than or equal to j : False
i != j  # i is not equal to j : True
```

However, python has other equivalence test statements that are fairly
unique to python. To check whether an object is contained in a list :

```python 
beatle = "John"
beatles = ["George", "Ringo", "John", "Paul"]
print(beatle in beatles)     # is John one of the beatles? : TRUE
print("Katy" not in beatles) # this is also TRUE. 
```

Conditionals (`if` statements) are also really easy to use in python. Take
a look at the following example:

```python
i = 4
sign = "zero"
if i < 0:
    sign = "negative"

elif i > 0:
    sign = "positive"

else:
    print("Sign must be zero")
    print("Have a nice day")

print(sign)
```

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code.

For example, in MatLab you begin an if statement with the word `if`
and you end it with `end if`. In C/C++ you delimit code blocks with curly
braces.

Python uses *white space* — in this case, *indentation,* to group lines of code. In this case, there are four spaces at the beginning of the line following the `if` statement. This is not just to make things look pretty - it
tells Python what the body of the `if`-statement is.

Marking blocks of code is a fundamental part of a language. In C, you'll see curly braces everywhere. In Python, you'll see indentation everywhere. It's how you (and the computer) know what your code means.

Other white space in Python (and most other languages) is only for people. Little things like putting blank lines in "sane" places, and putting spaces between variables and operators (say, `a + b` rather than `a+b`) can make your code a lot easier to read.

**Exercise**
Write an if statement that prints whether x is even or odd.

Hint: Try out what the "%" operator. What does 10 % 5 and 10 % 6 return?

##Writing Code for People summary

We learned some basics of python and saw that variable type, name, comments and white space affect more than just code functionality, they affect the readability for others and your future self.
Variable names can make a huge difference in code readability and types are important in conveying intent. 

[Up To Schedule](../../README.md) - Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) - Forward To [Don't Repeat Yourself](../dont_repeat_yourself/Readme.md)
