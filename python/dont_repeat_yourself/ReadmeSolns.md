#Solutions to exercises in Don't Repeat Yourself

* * * *
##For Loops, Indentation

**Short Exercise**

Using a loop, calculate the factorial of 6 (the product of all positive integers 
up to and including 6).

***Solution***

```python
factorial = 1	#initialize factorial
for i in range(1,7):
    factorial = factorial * i
print factorial
```

* * * *
##Conditionals

**Short Exercise**

Write an if statement that prints whether x is even or odd.

Hint: Try out what the "%" operator. What does 10 % 5 and 10 % 6 return?

***Solution***

```python
x = 5
if x%2 == 0:
	print "Even!"
else:
	print "Odd!"
```

* * * *
##Reading from Files

**Short Exercise**

Use a loop to print the area codes and number of occurences in one line.
Remember how we previously looped through a dictionary using `iteritems`.

Your output should look like this:

    203 4
    800 4
    608 8
    773 3

***Solutions***

```python
for ac, num in areacodes.iteritems():
    print ac,num
```

* * * * 
##Built-in string methods
 
 **Short Exercise**
 
 Because the binding strength of guanine (G) to cytosine (C) is different from
 the binding strength of adenine (A) to thymine (T) (and many other
 differences), it is often useful to know the fraction of a DNA sequence that
 is G's or C's. Go to the [string method
 section](http://docs.python.org/2/library/string.html) of the Python
 documentation and find the string method that will allow you to calculate this
 fraction.

 ```python
 # Calculate the fraction of G's and C's in this DNA sequence
 seq1 = 'ACGTACGTAGCTAGTAGCTACGTAGCTACGTA'
 gc = 
 ```

 Check your work:

 ```python
 round(gc, ndigits = 2) == .47
 ```
 
***Solution***

```python
gc = (seq1.count('G')+seq1.count('C'))/float(len(seq1))
```

* * * *
##Create your own functions

**Short Exercise**

Make a function that calculate the GC content of a given DNA sequence. For the more advanced participants, make your function able to handle sequences of mixed case (see the third test case).

```python
def calculate_gc(x):
    """Calculates the GC content of DNA sequence x.
    x: a string composed only of A's, T's, G's, and C's."""
```

Check your work:

```python
print round(calculate_gc('ATGC'), ndigits = 2) == 0.50
print round(calculate_gc('AGCGTCGTCAGTCGT'), ndigits = 2) == 0.60
print round(calculate_gc('ATaGtTCaAGcTCgATtGaATaGgTAaCt'), ndigits = 2) == 0.34
```

***Solution***

```python
def calculate_gc(x):
    """Calculates the GC content of DNA sequence x.
    x: a string composed only of A's, T's, G's, and C's."""
    x = x.upper()
    gc = (x.count('G')+x.count('C'))/float(len(x))
    return gc
```

* * * *
##Modules

**Short Exercise**

We have written a number of short functions. Collect these in a text file with an extension ".py", for example, "myFunctions.py". Test out the different import methods listed above. You may want to reset the iPython session between imports in the same way as the examples.

Try adding a new function to the module. Note that you need to `reload` the module in python to update it, if the module is already imported. For example:

```python
import myFunctions as myFun
# ... editing myFunctions.py in nano or other text editor...
reload(myFun)
```

***Solution***

Create a file `myFunctions.py` containing the following functions:

```python
def square(x):
    return x * x

def hello(time, name):
    """Print a nice message. Time and name should both be strings.
    
    Example: hello('morning', 'Software Carpentry')
    """
    print 'Good ' + time + ', ' + name + '!'

def calculate_gc(x):
    """Calculates the GC content of DNA sequence x.
    x: a string composed only of A's, T's, G's, and C's."""
    x = x.upper()
    gc = (x.count('G')+x.count('C'))/float(len(x))
    return gc
```

* * * *
##General Problem

**Short Exercise**

One common pattern is to generalize an existing function to work over a wider class of inputs. Try this by generalizing the `calculate_gc` function above to a new function, `calculate_dna_fraction` that computes the fraction for an arbitrary list of DNA bases. Add this to your own module file. Remember to `reload` the module after adding or modifying the python file. (This function will be more complicated than previous functions, so writing it interactively within iPython will not work as well.)

```python
def calculate_dna_fraction(x, bases):
    """Calculate the fraction of DNA sequence x, for a set of input bases.
    x: a string composed only of A's, T's, G's, and C's.
    bases: a string containing the bases of interest (A, T, G, C, or 
       some combination)"""
```

Check your work. Note that since this is a generalization of `calculate_gc`, it should reproduce the same results as that function with the proper input:


```python
test_x = 'AGCGTCGTCAGTCGT'
print calculate_gc(test_x) == calculate_dna_fraction(test_x, 'GC')
print round(calculate_dna_fraction(test_x, 'C'), ndigits = 2) == 0.27
print round(calculate_dna_fraction(test_x, 'TGC'), ndigits = 2) == 0.87
```

Generalization can bring problems, due to "corner cases", and unexpected inputs. You need to keep these in mind while writing the function; this is also where you should think about test cases. For example, what should the results from these calls be?

```python
print calculate_dna_fraction(test_x, 'AA')
print calculate_dna_fraction(test_x, '')
print calculate_dna_fraction(test_x, 2.0)
```

***Solution***

There are many ways to write  a solution to solve this exercise.

```python
def calculate_dna_fraction(x, bases):
    """Calculate the fraction of DNA sequence x, for a set of input bases.
    x: a string composed only of A's, T's, G's, and C's.
    bases: a string containing the bases of interest (A, T, G, C, or 
       some combination)"""

    if type(bases) != str:
        print 'Provided bases are not a string' #Print a helpful message
        return 0.
    x = x.upper()

    if bases.count('A') != 0:
        A = x.count('A')
    else:
        A = 0.
    if bases.count('T') != 0:
        T = x.count('T')
    else:
        T = 0.
    if bases.count('G') != 0:
        G = x.count('G')
    else:
        G = 0.
    if bases.count('C') != 0:
        C = x.count('C')
    else:
        C = 0.
    
    return float(A+T+G+C)/len(x)
```

* * * *
##Reading Cochlear implant data into Python

**Part I**

```python
def view_cochlear(filename):
    """Write your docstring here.
    """
```

Test it out:

```python
view_cochlear('/home/<username>/boot-camps/shell/data/alexander/data_216.DATA')
view_cochlear('/home/<username>/boot-camps/shell/data/Lawrence/Data0525')
```

***Solution***

```python
def view_cochlear(filename):
    """Read in a file and print each line of cochlear implant data.
    filename: a string of the filename for the file to be open and printed
    """
    f = open(filename)
    for line in f:
        print line
	
```

**Part II**

Adapt your function above to exclude the first line using the flow control techniques we learned in the last lesson. The first line is just `#` (but don't forget to remove the `'\n'`).

```python
def view_cochlear(filename):
    """Write your docstring here.
    """
```

Test it out:


```python
view_cochlear('/home/<username>/boot-camps/shell/data/alexander/data_216.DATA')
view_cochlear('/home/<username>/boot-camps/shell/data/Lawrence/Data0525')
```

***Solution***

```python
def view_cochlear(filename):
    """Read in a file and print each line of cochlear implant data.
    filename: a string of the filename for the file to be open and printed
    """
    f = open(filename)
    for line in f:
        line = line.strip()
        if line == '#':
            continue
        print line
```

**Part III**

Adapt your function above to return a dictionary containing the contents of the file. Split each line of the file by a colon followed by a space (': '). The first half of the string should be the key of the dictionary, and the second half should be the value of the dictionary.

```python
def load_cochlear(filename):
    """Write your docstring here.
    """
```

Check your work:

```python
data_216 = load_cochlear("/home/<username>/boot-camps/shell/data/alexander/data_216.DATA")
print data_216["Subject"]

Data0525 = load_cochlear("/home/<username>/boot-camps/shell/data/Lawrence/Data0525")
print Data0525["CI type"]
```

***Solution***

```python
def load_cochlear(filename):
    """Load in a file and return a dictionary containing cochlear implant data.
    filename: a string of the filename for the file to be open and printed
    """
    f = open(filename)
    cochlear_data = {}
    for line in f:
        line = line.strip()
        if line == '#':
            continue
        line = line.split(':')
        cochlear_data[line[0]] = line[1]
    return cochlear_data
```

* * * *
##Bonus Exercise: Transcribe DNA to RNA

During transcription, an enzyme called RNA Polymerase reads the DNA sequence and creates a complementary RNA sequence. Furthermore, RNA has the nucleotide uracil (U) instead of thymine (T). 

**Exercise**

Write a function that mimics transcription. The input argument is a string that contains the letters A, T, G, and C. Create a new string following these rules: 

* Convert A to U
* Convert T to A
* Convert G to C
* Convert C to G
 
Hint: You can iterate through a string using a `for` loop similarly to how you loop through a list.

```python
def transcribe(seq):
    """Write your docstring here.
    """
```

Check your work:

```python
print transcribe('ATGC') == 'UACG'
print transcribe('ATGCAGTCAGTGCAGTCAGT') == 'UACGUCAGUCACGUCAGUCA'
```

***Solution***
```python
def transcribe(seq):
    """Transcribe DRA sequence to RNA. Returns a string containing RNA sequence.
    seq: a string composed only of A's, T's, G's, and C's
    """
    seq = seq.upper()
    rna = ''
    for s in seq:
        if s == 'A':
            rna += 'U'
        elif s == 'T':
            rna += 'A'
        elif s == 'G':
            rna += 'C'
        else:
            rna += 'G'
    return rna
```