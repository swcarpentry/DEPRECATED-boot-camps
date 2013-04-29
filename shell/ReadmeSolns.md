# Solutions to execises in Introduction to Shell

* * * *
**Short Exercise**

1. Use the manual page for `ls` to guess what you would expect from
using the arguments `-l`, '-t', '-r' at the same time.

``` 
This will produce a long listing `-l`, sorted by the time that the
file last changed `-t`, and in the reverse order `-r`.
```

2. Try the following and see if you can figure out what they do, either by examining the results or consulting the manual page.
   * `ls -lS` (equivalent to `ls -l -S`)
   * `ls -lt` (equivalent to `ls -l -t`)
   * `ls -1`  (that's the number one, not a letter 'ell')

```
In this order, these will produce;
* `ls -lS` lists files in a long listing sorted by their size
* `ls -lt` lists files in a long listing sorted by its time
* `ls -1`  lists files in a brief listing, but 1 per line
```

* * * *
**Short Exercise**

Now, list the contents of the /bin directory. Do you see anything
familiar in there?

```
ls /bin

In this listing, you might recognize the different colors indicating
different file types, especially executables.  You might specifically
recognize that we have programs like 'ls' and 'cp' in this directory.
```

* * * *
**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that contain the letter `a`

    `ls /bin/*a*`

2.  List all of the files in `/bin` that contain the letter `a` or the letter `b`

    `ls /bin/*a* /bin/*b*`

3.  List all of the files in `/bin` that contain the letter `a` AND the letter `b`

    `ls /bin/*a*b* /bin/*b*a*`

* * * *

**Short Exercises**

1.  Print out the contents of the `~/boot-camps/shell/dictionary.txt`
    file. What does this file contain?

```
cat ~/boot-camps/shell/dictionary.txt

This file contains a list of english language words in alphabetical order.
```

2.  Without changing directories, (you should still be in `shell`),
    use one short command to print the contents of all of the files in
    the `/home/<username>/boot-camps/shell/data/THOMAS` directory.

```
cat data/THOMAS/*
```

* * * *
**Short Exercise**

Use the commands we've learned so far to figure out how to search
in reverse while using `less`.

```
Once you have started a search with `/` you can continue it forward with `n` and in reverse with `N`.

You can start a search in the reverse direction with `?`.
```

* * * * 

**Short Exercise**

Use `>>`, to append the contents of all of the files whose names
contain the number 4 in the directory:

    /home/<username>/boot-camps/shell/data/gerdal

to the existing `all_data` file. Thus, when you are done `all_data`
should contain all of the experiment data from Bert and any
experimental data file from gerdal with filenames that contain the
number 4.

``` 
This exercise starts after a learners have already created a file
named `all_data` as follows:

cat Bert/* > all_data

To now add the data from the files with a 4 in the name in the
`gerdal` directory:

cat gerdal/*4* >> all_data
```

* * * * 
**Short Exercise**

Do the following:

1.  Rename the `all_data_IMPORTANT` file to `all_data`.
2.  Create a directory in the `data` directory called `foo`
3.  Then, copy the `all_data` file into `foo`

```

This assmes that a file named `all_data_IMPORTANT` has already been created.
Then the following commands:

    mv all_data_IMPORTANT all_data
    mkdir foo
    cp all_data foo/.
```

**Short Exercise**

Navigate to the `data` directory. Use one `find` command to perform each
of the operations listed below (except number 2, which does not
require a `find` command):

1.  Find any file whose name is "NOTES" within `data` and delete it 

```
One solution is to use find with its -exec argument:

find . -name "*NOTES*" -exec rm {} \;

This will find all files in the hierarchy beginning with the current
directory (.), with a name containing the word NOTES, and execute the
remove (rm) command on each {}.

Another solution is to use find and pipe the results to xargs:

find . -name "*NOTES*" | xargs rm

This is similar but instead of using find to run the rm command, it
passes all of the results of find as if they were command line
arguments for rm.
```

2.  Create a new directory called `cleaneddata`

```
mkdir cleaneddata
```

3.  Move all of the files within `data` to the `cleaneddata` directory

```
mv */* cleaneddata

Will move all the files (*) in each directory (*/) into the director cleaneddata.
```

4.  Rename all of the files to ensure that they end in `.txt` (note:
    it is ok for the file name to end in `.txt.txt`
```
find cleaneddata -type f -exec mv {} {}.txt \;

This will find only files and not directories (-type f) in the
hierarchy starting with cleaneddata and for each one {} will move it
to the same name with .txt added.
```

**BONUS**

Redo exercise 4, except rename only the files which do not already end
in `.txt`. You will have to use the `man` command to figure out how to
search for files which do not match a certain name. 

```
find cleaneddata -type f -not -name "*.txt" -exec mv {} {}.txt \;

This adds the condition that the name does not end in ".txt".
```

* * * *

