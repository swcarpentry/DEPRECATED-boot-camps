# Solutions to exercises in Introduction to Shell

* * * *
##Arguments

**Short Exercises**

*1*. Use the manual page for `ls` to guess what you would expect from
using the arguments `-l`, `-t`, and `-r` at the same time.

***Solution***

This will produce a long listing `-l`, sorted by the time that the
file last changed `-t`, and in the reverse order `-r`.

* * * *
##Full vs. Relative Paths

**Short Exercise**

Now, list the contents of the /bin directory. Do you see anything
familiar in there?

***Solution***
```
ls /bin
```

In this listing, you might recognize the different colors indicating
different file types, especially executables.  You might specifically
recognize that we have programs like 'ls' and 'cp' in this directory.

* * * *
##Wild Cards

**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that contain the letter `a`

    `ls /bin/*a*`

2.  List all of the files in `/bin` that contain the letter `a` or the letter `b`

    `ls /bin/*a* /bin/*b*`

* * * *
## Command history

**Short Exercise**

1. Find the line number in your history for the last exercise (listing
directories in `~/boot-camps/shell/data`) and reissue that command.


* * * *
## Working with Files

**Short Exercises**

*1*.  Print out the contents of the `~/boot-camps/shell/dictionary.txt`
    file. What does this file contain?

***Solution***

```
cat ~/boot-camps/shell/dictionary.txt
```

This file contains a list of english language words in alphabetical order.

* * * *
*2*.  Without changing directories, (you should still be in `shell`),
    use one short command to print the contents of all of the files in
    the `/home/<username>/boot-camps/shell/data/THOMAS` directory.

***Solution***

```
cat data/THOMAS/*
```

* * * *
##Redirection

**Short Exercise**

Use `>>`, to append the contents of all of the files whose names
contain the number 4 in the directory `/home/<username>/boot-camps/shell/data/gerdal` 
to the existing `all_data` file. Thus, when you are done, `all_data`
should contain all of the experiment data from Bert and any
experimental data file from gerdal with filenames that contain the
number 4.

***Solution***

This exercise starts after a learners have already created a file
named `all_data` as follows:

```
cat Bert/* > all_data
```

To now add the data from the files with a 4 in the name in the
`gerdal` directory:

```
cat gerdal/*4* >> all_data
```

* * * * 
## Creating, moving, copying, and removing

**Short Exercise**

Do the following:

1.  Rename the `all_data_IMPORTANT` file to `all_data`.
2.  Create a directory in the `data` directory called `foo`
3.  Then, copy the `all_data` file into `foo`

***Solution***
This assumes that a file named `all_data_IMPORTANT` has already been created.
Then the following set of commands should work for 1-3 above:

```
mv all_data_IMPORTANT all_data
mkdir foo
cp all_data foo/.
```

* * * * 
## Finding files

**Exercises**

Navigate to the `data` directory. Use one `find` command to perform each
of the operations listed below (except number 2, which does not
require a `find` command):

*1*.  Find any file whose name is "NOTES" within `data` and delete it 

***Solution***
One solution is to use find with its -exec argument:

```
find . -name "*NOTES*" -exec rm {} \;
```

This will find all files in the hierarchy beginning with the current
directory (.), with a name containing the word NOTES, and execute the
remove (rm) command on each {}.

Another solution is to use find and pipe the results to xargs:

```
find . -name "*NOTES*" | xargs rm
```

This is similar but instead of using find to run the rm command, it
passes all of the results of find as if they were command line
arguments for rm.

* * * *
*2*.  Create a new directory called `cleaneddata`

***Solution***

```
mkdir cleaneddata
```

* * * *
*3*.  Move all of the files within `data` to the `cleaneddata` directory

***Solution***

```
mv */* cleaneddata
```

Will move all the files (*) in each directory (*/) into the director cleaneddata.


* * * *

