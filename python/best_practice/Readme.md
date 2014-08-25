[Up To Schedule](../../README.md) - 
Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) -
Forward To [Make Incremental Changes I](../../version-control/git/local/Readme.md)

- - - -

# Write Code for People: Variables, Data Types, Comments and Modularity

This lecture covers a number of topics under the heading of *"writing code for
people"*:
* naming variables in meaningful ways
* choosing appropriate data types
* writing effective comments
* breaking your code/script into small pieces

In some cases, we'll also look at the best practice of *"don't repeat
yourself"*.

These concepts are important and apply to any language.  For today's lecture,
we'll use python as the language, but you use a different language, you should
think about how you would use these practices in your language.  In
particular, we will use the [iPython interpreter](ipython.md), a
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

Written in python that is:
```python
for filename in filelist:
    data_record = parseFile(filename)
    all_data.append(data_record)

csv_separator = ','

writeCSVHeader(column_labels,csv_separator)

for data_record in all_data:
    writeCSVRow(column_labels,data_record,csv_separator)
```

----
![Exercise](pics/exercise.jpg) **Follow along**

1. Open a new file in your editor called `text2csv.py`
2. Copy & paste the above code into that editor.

---

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

----
![Exercise](pics/exercise.jpg) **Short Exercise**

Which of these are bad choices and why?  Which are reasonable alternatives?
* `fn`
* `all_data_from_all_the_files_in_the_list`
* `col_lbls`
* `record`
* `sep`
* `col_labels` 

----

We have chosen meaningful names for each function: `parseFile`, `writeCSVHeader`,
`writeCSVRow`.  It is often recommended that functions be given names that are
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
for many different items.

In addition, using the `csv_separator` variable instead of just typing `','`
has two advantages:
1. it provides context for what this means
2. it means we can change it in a single place

### Initialize some of these variables

1. We can initialize the list of files from the command-line arguments.
2. We should also declare the immutable ordered list of column labels.
3. Now that we know `all_data` should be a list, we should initialize it as a
   empty list.

----
![Exercise](pics/exercise.jpg) **Follow Along**

Add these lines above the others in your script:
```python
import sys
filelist = sys.argv[1:]
column_labels = ("Subject","Reported","Year/month of birth",
                "Sex","CI type","Volume","Range","Discrimination")
all_data = []
```

----

## Step 2: Parse a single file

To continue with our script we can now focus on the first function,
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

Or written as python:
```python
def parseFile(filename):
    """Read all the lines from a file and return them as a dictionary of key/value pairs"""
    textfile = open(filename, 'r')
    data_record = {}
    for line in textfile:
        # all lines that contain key/value pairs must have at least one colon
        if ':' in line:
            (key,value) = extractData(line)
            data_record[key] = value
    return data_record
```

----
![Exercise](pics/exercise.jpg) **Follow along**

1. Copy & paste the above code into your file at the top.

---

### New concept: Reading from files

At this point we've had to introduce a new concept: reading data from files.
For most purposes, python makes this very easy.  When we `open` a file, we get
an object over which we can iterate, much like a list.  Let's look at the file
`phonenums.txt` using this concept:

----
![Exercise](pics/exercise.jpg) **Try it exercise**

```python
In [15]: file_data = open('phonenums.txt')
In [16]: for line in file_data:
   ....:     print line
```

----

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
    value = separator.join(line_data[1:]).strip()
    return key,value
```

----
![Exercise](pics/exercise.jpg) **Follow along**

1. Copy & paste the above code into your file at the top.

---


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

----
![Exercise](pics/exercise.jpg) **Think aloud:** Recommend some alternative variable names for `line_data`.

----

### Choosing appropriate data types

All the choices in this case were made for us by the requirements of the
string functions we call.


## Step 4: Write the CSV data

![Exercise](pics/exercise.jpg) **Discussion Exercise** Given the following implementation of the final two
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


##Writing Code for People summary

We learned some basics of python and saw that variable type, name, comments
and white space affect more than just code functionality, they affect the
readability for others and your future self.  Variable names can make a huge
difference in code readability and types are important in conveying intent.

* * * * *

[Up To Schedule](../../README.md) -
Back To [Let the Computer Do the Work](../../shell/automation/Readme.md) -
Forward To [Make Incremental Changes I](../../version-control/git/local/Readme.md)
