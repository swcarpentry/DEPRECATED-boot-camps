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
