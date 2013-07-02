# Python 2 : Flow Control - Loops, Conditionals, etc


Conditionals
============

A conditional (if statement) is some statement that in general says: 
"When some boolean is true, do the following. Otherwise, do this other
thing." In a statement we thus need to be able to compare things.

Many equivalence statements exist in Python that are similar to those of
other languages:

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
	print sign
elif i > 0:
	sign = "positive"
	print sign
else:
	print sign 
```

What do you expect the value of sign to be after it has run? 

Please note: the if statement consists of everything from the first if until the end of the else.

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code. For example, in MatLab you begin an if
statement with the word "if" and you end it with "end if". In C/C++ you
delimit blocks with curly braces. Python uses ''indentation'' to delimit
code blocks. The **indentation** above is NOT just to make things look
pretty - it tells Python what the body of the if-statement is. This is
true whenever we create any code blocks, such as the bodies of loops,
functions or classes. This is also why you will have to press enter twice after
the if statement - after the first time, python is offering you to type 
in more things that would be executed in the else, so you press enter to 
tell python that you are done, and that it should evaluate the if.

**Task:** Open an editor (nano for instance), type the script into a text file and save
it as if.py. Run the script like this: 

```
python if.py
```

and see if the value of the sign is what you expected it to be. 

**Task:** Change the value of i so that you trigger all of the parts of the if statement.

For Loops
=========

For loops are something we use a lot in python. Let's try a simple loop:


```python
for i in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
    print i
```

For loops basically iterate over anyting that is iteratable. The things you have seen so far that are iterable are strings, lists and dictionaries. 

##Using for loops to process files

Commonly, when you read in files, it is often done with the readlines method. This method gives us a list of strings. Let´s try to read in a file and print each line:

```python
fh = open("popdata.txt", "r")
lines = fh.readlines()
fh.close()

for line in lines:
	print line
```

Note, we get one extra newline here since each line that we read contains a newline. We can remove the newline by using rstrip():

```python
for line in lines:
	print line.rstrip()
```

Now, let us combine this with some text splitting. This will let us get to each line of the text. The goal in this case is to print out the country and the population: 

```python
for line in lines:
	fields = line.split()
	print fields[0], fields[2]
```

**Task:** Can you figure out how to print the life expectancy instead?

**Task:** Can you figure out how to skip the first line, i.e. leave out the header line?

## Combining for and if

Fairly often, we end up reading in a file and doing something with it. For this we often end up using if statements to figure out what to do with each line. In this case, let´s work with a bigger file, one that contains population data for the year 1952. 

```python
fh = open("1952.txt", "r")
popdat = fh.readlines()
fh.close()
```

So - now we have read in the years and have them as a list of strings in the wholeyear variable. 

Now, let's select only those countries that are in Europe. In this case, that information would be found in field no 3:

```python
for line in popdat:
	fields = line.split()
	if fields[3] == "Europe":
		print line.rstrip()
```

**Task:** Take the code that you have above, and put into a script named europe.py. Run it using python.

**Task:** Figure out how you print out the Country, the Continent and the life expectancy.

## Aggregating data

We will now go a bit further. We can now access each line, process it, and print it out. We will now instead aggregate the life expactancies and calculate the average.

Note: this all takes place in your script file.

```python
fh = open("1952.txt", "r")
popdat = fh.readlines()
fh.close()

#create a list that can keep the life expectancies
aggregate = []
for line in popdat:
	fields = line.split()
	if fields[3] == "Europe":
		lifeexp = float(fields[4])
		aggregate.append(lifeexp)
		
#Now we calculate the average

no_countries = len(aggregate)
total = sum(aggregate)
avg = total/no_countries
print "Average lifespan in Europe:", avg

```

Using dictionaries
----------------------------------
In the previous example we got the ones that corresponded to Europe and put them into a list. We then used this list to calculate the average life expectancy for all of the different continents. We could do that by creating one list for each continent and manually adding each life expectancy to that, but that would be a bit heavy to do by hand. Instead, we can write python code that does this for us. For this, we use dictionaries. We will do this incrementally, so that we can show how this would work.

First, let's get all of the countries and their life expectancies into a dictionary:


```python
fh = open("1952.txt", "r")
popdat = fh.readlines()
fh.close()

country_lifeexp = {}
for line in popdat[1:]:
	fields = line.split()
	country = fields[0]
	lifeexp = float(fields[4])
	country_lifeexp[country] = lifeexp
	
# we can now print it like this:
for key in country_lifeexp:
	print key, country_lifeexp[key] 
```

So, this is how dictionaries work in this case. Now, we change this a bit. We use the continent as the key instead of the country. The value is still the life expectancy, but we now have several values per key, and that always spells using a list. How do we do this? The way to do this is to use this pattern:

```python
#create dictionary
my_dict = {}
for line in popdat[1:]:
    fields = line.split()
    key_field = fields[key_column_number]
    aggregate_value = float(fields[value_column_number])
    # Now, test if you already have it in the dictionary:
    if key_field not in my_dict:
        # have to create a new list to keep the values in
        my_dict[key_field] = []
    # So, now we know we have a key value pair in our dictionary, where
    # the key_value is the key, and we have a list as a value. We can
    # now add to this list
    my_dict[key_field].append(aggregate_value)
	
```

**Task:** use the above pattern to print each continent, plus a list of the life expectancies

**Task:** figure out how to calculate the average for each continent. In this case, you need to have that in the second for loop, where you use the data that you have put into your dictionary.



Previous: [Data types](1_Data_Types.md) Next: [Functions and modules](3_Functions_and_Modules.md)
