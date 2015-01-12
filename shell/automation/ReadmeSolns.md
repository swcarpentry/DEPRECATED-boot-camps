# Solutions to execises in Automating Workflows

* * * *
**Short Exercise**

From the `data` directory, use a single command to view the `Subject` 
line of every file in `THOMAS` and `jamesm`.

**Solution**
```
grep Subject THOMAS/* jamesm/*
```

* * * *

**Short Exercise**

Figure out how to get `wc` to print the length of the longest line in
`all_data`.

**Solution**
```
wc -L all_data
```

* * * *

**Short Exercise**

Use the `echo` command and the append operator, `>>`, to append your
name to the file and then sort the output into a new file called `Sorted`.

**Solution**

This assumes that there is already a file named `toBeSorted` that we
made with an editor (nano) that contains:

```
Bob
Alice
Diane
Charles
```

To add your name to the file:

```
echo Lauren >> toBeSorted
```

The file can be sorted into a new file with:

```
sort toBeSorted > Sorted
```

* * * *
**Short Exercise**

1. Use the `man` command to find out how to sort the output from `wc` in
reverse order.

**Solution**
To reverse the order of the `sort` command you add the `-r` argument.

2. Combine the `wc`, `sort`, `head` and `tail` commands so that only the
`wc` information for the largest file is listed

Hint: To print the smallest file, use:

    wc Bert/* | sort -k 3 -n | head -n 1

**Solution**

It is important to realize that when we are using `wc` on a set of files
`Bert/*` it will always conclude with a line that is the sum.  No
matter how we sort this, the largest single file is always either the
2nd line or the 2nd last line.

The above line prints the info for the smallest file because it is the
first line.  If we wanted the last line:

    wc Bert/* | sort -k 3 -n | tail -n 1

but, again, this will be the sum of all the files.

We can get the last 2 lines with:

    wc Bert/* | sort -k 3 -n | tail -n 2

and then we can take the first of those lines:

    wc Bert/* | sort -k 3 -n | tail -n 2 | head -n 1

* * * * 
**Short Exercise**

Create an executable script called `smallestrange` in the `data`
directory, that is similar to the `smallest` script, but prints the
file containing the file with the smallest Range. Use the commands
`grep`, `sort`, and `tail` to do this.

**Solution**

We can experiment on the command line to find the right commands.
Starting with the command to just get all the "Range" data:

    grep Range *

will search for the word "Range" in all the files (*).

Next we can sort this numerically based on the numerical value of the
range, found in the 2nd column:

    grep Range * | sort -n -k 2

The smallest range will be the first item in that list:

    grep Range * | sort -n -k 2 | head -n 1

Now that we've found this line, we can put it in a script to help us
reuse this with less typing and be less error prone.

We create a new file by opening it with nano

    nano smallestrange

and putting in the following lines:

    #!/bin/bash
    grep Range * | sort -n -k 2 | head -n 1

The first line tells the shell what program to use to run this script.

After exiting the editor (Control-X) we can make it executable:

    chmod a+x smallestrange

* * * * 
