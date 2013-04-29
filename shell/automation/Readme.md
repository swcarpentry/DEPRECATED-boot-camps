[Up To Schedule](../../README.md) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People I](../../python/variables_and_types/Readme.md)

# Automating Workflows

**Material by Paul Wilson, Milad Fatenejad, Sasha Wood, and Radhika Khetani**

## Count the words

The `wc` program (word count) counts the number of lines, words, and
characters in one or more files. Make sure you are in the `data`
directory, then enter the following command:

    wc Bert/* gerdal/*4*

For each of the files indicated, `wc` has printed a line with three
numbers. The first is the number of lines in that file. The second is
the number of words. Finally, the total number of characters is
indicated. The final line contains this information summed over all of
the files. Thus, there were 10445 characters in total. 

Remember that the `Bert/*` and `gerdal/*4*` files were merged
into the `all_data` file. So, we should see that `all_data` contains
the same number of characters:

    wc all_data

Every character in the file takes up one byte of disk space. Thus, the
size of the file in bytes should also be 10445. Let's confirm this:

    ls -l all_data

Remember that `ls -l` prints out detailed information about a file and
that the fifth column is the size of the file in bytes.

* * * *
**Short Exercise**

Figure out how to get `wc` to print the length of the longest line in
`all_data`.

* * * *

## The awesome power of the Pipe

Suppose I wanted to only see the total number of character, words, and
lines across the files `Bert/*` and `gerdal/*4*`. I don't want to
see the individual counts, just the total. Of course, I could just do:

    wc all_data

Since this file is a concatenation of the smaller files. Sure, this
works, but I had to create the `all_data` file to do this. Thus, I
have wasted a precious 10445 bytes of hard disk space. We can do this
*without* creating a temporary file, but first I have to show you two
more commands: `head` and `tail`. These commands print the first few,
or last few, lines of a file, respectively. Try them out on
`all_data`:

    head all_data
    tail all_data

The `-n` option to either of these commands can be used to print the
first or last `n` lines of a file. To print the first/last line of the
file use:

    head -n 1 all_data
    tail -n 1 all_data

Let's turn back to the problem of printing only the total number of
lines in a set of files without creating any temporary files. To do
this, we want to tell the shell to take the output of the `wc Bert/*
gerdal/*4*` and send it into the `tail -n 1` command. The `|`
character (called pipe) is used for this purpose. Enter the following
command:

    wc Bert/* gerdal/Data0559 | tail -n 1

This will print only the total number of lines, characters, and words
across all of these files. What is happening here? Well, `tail`, like
many command line programs will read from the *standard input* when it
is not given any files to operate on. In this case, it will just sit
there waiting for input. That input can come from the user's keyboard
*or from another program*. Try this:

    tail -n 2

Notice that your cursor just sits there blinking. Tail is waiting for
data to come in. Now type:

    French
    fries
    are
    good

then CONTROL+d. You should is the lines:

    are
    good

printed back at you. The CONTROL+d keyboard shortcut inserts an
*end-of-file* character. It is sort of the standard way of telling the
program "I'm done entering data". The `|` character is replaces the
data from the keyboard with data from another command. You can string
all sorts of commands together using the pipe. 

The philosophy behind these command line programs is that none of them
really do anything all that impressive. BUT when you start chaining
them together, you can do some really powerful things really
efficiently. If you want to be proficient at using the shell, you must
learn to become proficient with the pipe and redirection operators:
`|`, `>`, `>>`.


### A sorting example

Let's create a file with some words to sort for the next example. We
want to create a file which contains the following names:

    Bob
    Alice
    Diane
    Charles

To do this, we need a program which allows us to create text
files. There are many such programs, the easiest one which is
installed on almost all systems is called `nano`. Navigate to
`/tmp/<username>` and enter the following command:

    nano toBeSorted

Now enter the four names as shown above. When you are done, press
CONTROL+O to write out the file. Press enter to use the file name
`toBeSorted`. Then press CONTROL+x to exit `nano`.

When you are back to the command line, enter the command:

    sort toBeSorted

Notice that the names are now printed in alphabetical order.

* * * *
**Short Exercise**

Use the `echo` command and the append operator, `>>`, to append your
name to the file, then sort it and make a new file called Sorted.

* * * *

Let's navigate back to `~/boot-camps/shell/data`. Enter the following command:

    wc Bert/* | sort -k 3 -n

We are already familiar with what the first of these two commands
does: it creates a list containing the number of characters, words,
and lines in each file in the `Bert` directory. This list is then
piped into the `sort` command, so that it can be sorted. Notice there
are two options given to sort:

1.  `-k 3`: Sort based on the third column
2.  `-n`: Sort in numerical order as opposed to alphabetical order

Notice that the files are sorted by the number of characters.

* * * *
**Short Exercise**

1. Use the `man` command to find out how to sort the output from `wc` in
reverse order.

2. Combine the `wc`, `sort`, `head` and `tail` commands so that only the
`wc` information for the largest file is listed

Hint: To print the smallest file, use:

    wc Bert/* | sort -k 3 -n | head -n 1

* * * * 

Printing the smallest file seems pretty useful. We don't want to type
out that long command often. Let's create a simple script, a simple
program, to run this command. The program will look at all of the
files in the current directory and print the information about the
smallest one. Let's call the script `smallest`. We'll use `nano` to
create this file. Navigate to the `data` directory, then:

    nano smallest

Then enter the following text:

    #!/bin/bash
    wc * | sort -k 3 -n | head -n 1

Now, `cd` into the `Bert` directory and enter the command
`../smallest`. Notice that it says permission denied. This happens
because we haven't told the shell that this is an executable
file. If you do `ls -l ../smallest`, it will show you the permissions on 
the left of the listing.

Enter the following commands:

    chmod a+x ../smallest
    ../smallest

The `chmod` command is used to modify the permissions of a file. This
particular command modifies the file `../smallest` by giving all users
(notice the `a`) permission to execute (notice the `x`) the file. If
you enter:

    ls -l ../smallest

You will see that the file name is green and the permissions have changed. 
Congratulations, you just created your first shell script!

# Searching files

You can search the contents of a file using the command `grep`. The
`grep` program is very powerful and useful especially when combined
with other commands by using the pipe. Navigate to the `Bert`
directory. Every data file in this directory has a line which says
"Range". The range represents the smallest frequency range that can be
discriminated. Lets list all of the ranges from the tests that Bert
conducted:

    grep Range *

* * * * 
**Short Exercise**

Create an executable script called `smallestrange` in the `data`
directory, that is similar to the `smallest` script, but prints the
file containing the file with the smallest Range. Use the commands
`grep`, `sort`, and `tail` to do this.

* * * * 

# Variables

We were earlier introduced to the `$PATH` variable.  This is a
variable that the shell expects to be able to function properly.
There are other variables defined and used by default.  To see a list
of variables just enter:

    set

or better yet, pipe that through `less` like:

   set | less

Now you will see a long list of variables that are already set.  Some
of these are built-in to the bash shell, others are set by the system
administrator for all users.



# Aliases

Another way to avoid having to retype commands is to use an **alias**.
Most shells allow you to define an alias after which it is available
to use in that shell.  Let's define an alias to perform the same
funtion as the `smallest` script that we made earlier.  If you enter:

    alias my_smallest='wc * | sort -k 3 -n | head -n 1'

you will have created a new command `my_smallest` that will be
available in any directory.

Some people also like to guard against accidentally deleting a file
and will create an alias for the `rm` command:

    alias rm='rm -i'

that double checks their intent every time they delete a file.  After
creating the above alias, try:

    touch testfile
    rm testfile

Some system administrators will make this a default, and some users
find it annoying.  You can find a list of aliases by just entering:

    alias

If you would like to remove this constant questioning, you can just
enter:

    unalias rm


# Startup Scripts

One of the most common scripts is one that is invoked every time you
login.  When using the bash shell it is called the `.bashrc` file.


## Bonus:

**backtick, xargs**: Example find all files with certain text

**alias** -> rm -i

**variables** -> use a path example

**.bashrc**

**du**

**ln**

**ssh and scp**

**Regular Expressions**

**Permissions**

**Chaining commands together**
