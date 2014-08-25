[Up To Schedule](../../README.md) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People](../../python/best_practice/Readme.md)

# Let the Computer Do the Work: Automating Workflows

**Material by Paul Wilson, Lauren Michael, Milad Fatenejad, Sasha Wood, and Radhika Khetani**

#Working with Files (cont...)

## Searching Within Files

We've searched for *files* with `find`, but you can also search for 
specific *contents* in one or more files using the command `grep` 
(and without laboriously searching through each file with something like `less`). 
The `grep` program is very powerful and useful, especially when combined
with other commands by using the pipe. Navigate to the `Bert`
directory. Every data file in this directory has a line which says
"Range". The range represents the smallest frequency range that can be
discriminated. Lets list all of the ranges from the tests that Bert
conducted:

    grep Range *

Notice that this is a much shorter command than the `find` command we 
used to view the "Volume" line of each file.

* * * *
**Short Exercise**

From the `data` directory, use a single command to view the `Subject` 
line of every file in `THOMAS` and `jamesm`.

* * * *

## Count the Words

What if you want to understand the total size of files you have in 
particular locations? The `wc` program (word count) counts the number 
of lines, words, and characters in one or more files. Make sure you 
are in the `data`
directory, then enter the following command:

    wc Bert/* gerdal/*4*

For each of the files indicated, `wc` has printed a line with three
numbers. The first is the number of lines in that file. The second is
the number of words. Finally, the total number of characters is
indicated. The final line contains this information summed over all of
the files. Thus, there were 10449 characters in total. 

Remember that the `Bert/*` and `gerdal/*4*` files were merged
into the `all_data` file. So, we should see that `all_data` contains
the same number of characters:

    wc all_data

Every character in the file takes up one byte of disk space. Thus, the
size of the file in bytes should also be 10449. Let's confirm this:

    ls -l all_data

Remember that `ls -l` prints out detailed information about a file and
that the fifth column is the size of the file in bytes.

* * * *
**Short Exercise**

Figure out how to get `wc` to print the length of the longest line in
`all_data`. (Hint: Where will the list of options for `wc` be?)

* * * *

# The awesome power of the Pipe

Suppose I wanted to only see the total number of lines, words, and 
characters across the files `Bert/*` and `gerdal/*4*`. I don't want to
see the individual counts, just the total. Because the `all_data` file 
is just a concatenation of this information, I could just do 
the following

    wc all_data

Sure, this works, but I had to create the `all_data` file to do this. Thus, I
have wasted a precious 10449 bytes of hard disk space that is already used 
by the number of individual files. We can do the same count
*without* creating a temporary file, but first I have to show you two
more commands: `head` and `tail`. These commands print the first few,
or last few, lines of data from a file or set of text, respectively. Try them out on
`dictionary.txt`:

    head ../dictionary.txt
    tail ../dictionary.txt

The `-n` option to either of these commands can be used to print the
first or last `n` number of lines. To print the first/last line of the
file use:

    head -n 1 all_data
    tail -n 1 all_data

Let's turn back to the problem of printing only the total number of
lines in a set of files without creating any temporary files. To do
this, we want to tell the shell to take the output of  
`wc Bert/* gerdal/*4*` and send it into the `tail -n 1` command. The `|`
character (called "pipe") is used for this purpose. Enter the following
command:

    wc Bert/* gerdal/*4* | tail -n 1

This will print only the total number of lines, words, and characters 
across all of these files. What is happening here? Well, `tail`, like
many command line programs will read from the *standard input* when it
is not given any files to operate on. In this case, it will just sit
there waiting for input. That input can come from the user's keyboard, 
from a file, *or from another program* via the pipe. Try this:

    tail -n 2

Notice that your cursor just sits there blinking. Tail is waiting for
data to come in. Now type:

    French
    fries
    are
    good

then <kbd>CONTROL</kbd>+<kbd>d</kbd>. You should see the lines:

    are
    good

printed back at you. The <kbd>CONTROL</kbd>+<kbd>d</kbd> keyboard shortcut inserts an
*end-of-file* character. It is sort of the standard way of telling the
program "I'm done entering data". 

The `|` character feeds output from the first program (to the left of 
the `|`) as input to the second program on the right. Therefor, you can 
string all sorts of commands together using the pipe. 

The philosophy behind these command line programs is that none of them
really do anything all that impressive. BUT when you start *chaining them together*, 
you can do some really powerful things really
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
files. There are many such programs, the easiest one of which is
installed on almost all systems and is called `nano`. Navigate to
the `shell` directory and enter the following command to create and 
add text to a new file called `toBeSorted`:

    nano toBeSorted

Now enter the four names as shown above. When you are done, press
<kbd>CONTROL</kbd>+<kbd>O</kbd> to write out these changes to the file. 
Press <kbd>ENTER</kbd> to use the file name
`toBeSorted`. Then press <kbd>CONTROL</kbd>+<kbd>X</kbd> to exit `nano`.

When you are back to the command line, enter the command:

    sort toBeSorted

Notice that the names are now printed in alphabetical order.

* * * *
**Short Exercise**

Use the `echo` command and the append operator, `>>`, to append your
name to the file, then sort it and make a new file called `Sorted`.

* * * *

Let's navigate back to `~/boot-camps/shell/data`. Enter the following command:

    wc Bert/* | sort -k 3 -n

We are already familiar with what the first of these two commands
does: it creates a list containing the number of lines, words, and 
characters in each file in the `Bert` directory. This list is then
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
`wc` information for the largest file is listed. Remember that the last 
and largest line of the `wc` output will always be the sum (which is not 
what we want).

Hint: To print the smallest file, use:

    wc Bert/* | sort -k 3 -n | head -n 1

* * * * 

#Shell Scripts

Printing the smallest file seems pretty useful. We don't want to type
out that long command often. Let's create a simple script, a simple
program, to run this command. The program will look at all of the
files in the current directory and print the information about the
smallest one. Let's call the script `smallest`. We'll use `nano` to
create this file. Navigate to the `data` directory, then type:

    nano smallest

Enter the following text and save your changes before exiting:

    #!/bin/bash
    wc * | sort -k 3 -n | head -n 1

Now, `cd` into the `Bert` directory and enter the command
`../smallest`. Mac and Linux users will likely get a "permission denied" 
error. (Windows users can likely execute without problem.). This happens
on unix-based operating systems because we haven't told the shell that this 
file can be *executed*. 

*Executables* are just files with executable permissions. 
If you do `ls -l ../smallest`, it will show you the permissions on 
the left of the listing.

To add executable permissions on Mac or Linux, enter the following commands:

    chmod a+x ../smallest
    ../smallest

The `chmod` command is used to modify the permissions of a file. This
particular command modifies the file `../smallest` by giving all users
(notice the `a`) permission to execute the file (notice the `x`). If
you enter the following:

    ls -l ../smallest

you will see that the permissions have changed. 

Congratulations, you just created your first shell script! You can now execute it
by entering the path location of `smallest` (absolute or relative) from within any 
directory you'd like to analyze.

* * * * 
**Short Exercise**

Create an executable script called `smallestrange` in the `data`
directory, that is similar to the `smallest` script, but prints the
line with the smallest "Range" value. Use the commands
`grep`, `sort`, and `head` to do this.

* * * * 

## Variables

We were earlier introduced to the `$PATH` variable.  This is a
variable that the shell expects to be able to function properly.
There are other variables defined and used by default.  To see a list
of variables just enter:

    set

You may recognize that the output includes some scripted portions 
(as your admins have created 
somethings specifically useful to your group), but the default "set" 
listing on Linux servers is a _set_ of variables defined for 
your shell session. 

To work with the information, pipe it through 
`less` like this:

    set | less

Now you will see a long list of variables that are already set. Some
of these are built-in to the bash shell, others are set by the system
administrator for all users.

Some important variables you can expect to see:

| Variable | What it stores |
| -------- | -------------- |
| HOME     | the full path to your home directory |
| HOSTNAME | the name of this computer |
| SHELL    | the full path to your current shell command |
| USER     | your user name |

You can use these variables in your commands.  For example, this is
yet another way to go to your home directory:

    cd $HOME

You can also define your own variables.  If you are always needing to
go to the same directory for your work, you could store it's path as a
variable.  Change to your home directory and try:

    DATADIR="~/boot-camps/shell/data"
    cd $DATADIR

You can also use one variable as part of the definition to another:

    DATADIR="$HOME/boot-camps/shell/data"
    cd $DATADIR

Let's make a directory for your data in the `$DATADIR`:

    mkdir $DATADIR/$USER

and then add another variable:

    MYDATADIR="$DATADIR/$USER"

##Example: The bash prompt
In `bash`, the prompt that you see at the beginning of each line is
just another variable, `$PS1`.  Examine the prompt:

    echo $PS1

and then update it to show a longer path relative to your home:

    PS1="[\u@\h \w]\$ "

When you define a variable, it is only available in this current shell
process.  If you spawn a new shell from this one, the variables will be
back to the default set.  After changing your prompt, as above, start
a new shell:

    bash

and notice that the prompt is now the same as the original one.  Let's
exit that shell and see how to keep those changes.

    exit

You can make variables available to processes that are spawned from
this shell with the `export` command:

    export PS1="[\u@\h \w]\$ "

Now launch a new shell and confirm that the prompt has changed, and
exit that shell again.

# Aliases

Another way to avoid having to retype commands is to use an **alias**.
Most shells allow you to define an alias so that it is available
to use in that shell.  Let's define an alias to perform the same
function as the `smallest` script that we made earlier.  If you enter:

    alias my_smallest='wc * | sort -k 3 -n | head -n 1'

you will have created a new command `my_smallest` that will be
available in any directory.


# Startup Scripts

Alas, none of those changes to variables or aliases will be there when
you logout and come back.  Thankfully, we can turn to a script to make
that happen.  The most common script is one that is invoked every time
you login.  When using the bash shell it is called the `.bashrc` file.

Go to your home directory and open the `.bashrc` file:

    nano .bashrc

At the end of this file, add the new variables and prompts that we
defined above:

    export DATADIR="$HOME/boot-camps/shell/data"
    export MYDATADIR="$DATADIR/$USER"
    export PS1="[\u@\h \w]\$ "
    alias my_smallest='wc * | sort -k 3 -n | head -n 1'

These variables should then be available to you in the future. On a server, 
these variables would be available even after you log out and log back in. Many 
experienced users will gradually build up a long list of
important variables and aliases for tasks they do frequently.


* * * * 
Jump to look at the [solutions to all the Automation exercises.](ReadmeSolns.md)


## Bonus topics:
You may wish to look into the below topics/tools to become even more advanced in using the shell:

[**ssh and scp**](http://software-carpentry.org/v4/shell/ssh.html)

[**permissions**](http://software-carpentry.org/v4/shell/perm.html)

**backticks and variables**

**loops**

**du**

**ln**

**regular expressions**

----

[Up To Schedule](../../README.md) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People](../../python/best_practice/Readme.md)
