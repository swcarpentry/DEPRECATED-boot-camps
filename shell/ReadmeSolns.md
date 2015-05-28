# Solutions to exercises in Introduction to Shell

* * * *
##Arguments

**Short Exercises**

Use the manual page for `ls` to guess what you would expect from
using the arguments `-l`, `-t`, and `-r` at the same time.

***Solution***

This will produce a long listing (`-l`), sorted by the time that the
file last changed (`-t`), and in the reverse order (`-r`).

* * * *
##Full vs. Relative Paths

**Short Exercise**

Now, list the contents of a directory of your own (one you know the 
location of by using the full path (and without using `cd`), so that 
you are running the command from your current location).

***Solution***

Each student will have a different answer, but you'll use `ls` to view 
the contents of a directory location that begins with `/`.

* * * *
##Wild Cards

**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that contain the letter `a` _or_ the letter `b` 
(including files that may contain both).

    `ls /bin/*a* /bin/*e*`

2.  List all of the files in `/bin` that contain the letter `a` _and_ the letter `b`.

    `ls /bin/*a*b* /bin/*b*a*`

* * * *
## Command history

**Short Exercise**

Find the line number in your history for the last exercise (listing
directories in `/bin`) and reissue that command.

**Solution**

You'll use your own history command number, and re-issue that command with `!` 
followed by the command number (as indicated in the example just before this exercise). 
So, if the command you want is number 181, you'd reissue it with:
```
!181
```


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
    the `~/boot-camps/shell/data/THOMAS` directory.

***Solution***

```
cat data/THOMAS/*
```

* * * *
##Redirection

**Short Exercise**

Use `>>`, to append the contents of all of the files whose names
contain the number 4 in the directory `/home/<username>/boot-camps/shell/data/h_jackson` 
to the existing `all_data` file. Thus, when you are done, `all_data`
should contain all of the experiment data from `Bert` *AND* any
experimental data file from `h_jackson` with filenames that contain the
number 4.

***Solution***

This exercise starts after learners have already created a file
named `all_data` as follows:

```
cat Bert/* > all_data
```

To now add the data from the files with a 4 in the name in the
`h_jackson` directory:

```
cat h_jackson/*4* >> all_data
```

* * * * 
## Creating, moving, copying, and removing

**Short Exercise**

Do the following:

1.  Rename the `all_data_IMPORTANT` file to `all_data`.
2.  Create a directory in the `data` directory called `foo`.
3.  Then, *copy* the `all_data` file into `foo`.

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

Navigate to the `data` directory. Use a single 
`find` command to perform each of the below exercises (aside from exercise
#1, which does not require `find`):

*1*. Create a new directory called `cleaneddata`.
*2*. Copy all of the files within the subdirectories of `data` into `cleaneddata`.
(Hint: Remember the wildcard. If you mess up, you can just delete the 
contents of cleaneddata, and try again.)
*3*. Find any files in `cleaneddata` containing "NOTES" in the name, and then delete them.
*4*. Rename all of the files to ensure that they end in `.txt`. 
(Note: it is okay for certain files to end in `.txt.txt`, as some 
already end with `.txt`.)


*1*.  Create a new directory called `cleaneddata`.

***Solution***

```
mkdir cleaneddata
```

* * * *
*2*.  Copy all of the files (only) within the subdirectories of `data` into `AllClean` directory.

***Solution***

```
find . -type f -exec cp {} cleaneddata \;
```

You can first just run the `find . -type f` portion, which will let you know you're selecting the 
correct files. Note that piping to `xargs` in this case woudl not work. Can you think of why you 
get output like the below:
```
cp: cleaneddata/data_103.DATA and ./cleaneddata/data_103.DATA are identical (not copied)
```
Hint: It has to do with the order by which find processes output through the commands following `-exec`.

* * * *
*3*. Find any files in `cleaneddata` containing "NOTES" in the name, and then delete them.

***Solution***

There are two options, here, You can either use `-exec`:
```
find cleaneddata -name *NOTES* -exec rm {} \;
```
or pipe to `xargs`:
```
find cleaneddata -name *NOTES* | xargs rm
```

* * * *
*4*. Rename all of the files to ensure that they end in `.txt`. 
(Note: it is okay for certain files to end in `.txt.txt`, as some 
already end with `.txt`.)

***Solution***
```
find cleaneddata -type f -exec mv {} {}.txt \;
```

* * * *
