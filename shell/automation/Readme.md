[Up To Schedule](../../README.md#schedule) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People](../../python/best_practice/Readme.md)

# Let the Computer Do the Work: Automating Workflows

**Material by Paul Wilson, Lauren Michael, Milad Fatenejad, Sasha Wood, and Radhika Khetani**

#5. Working with Files (cont...)

## Searching Within Files

We've searched for *files* with `find`, but you can also search for 
specific *contents* in one or more files using the command `grep` 
(and without laboriously searching through each file with something like `less`). 
The `grep` program is very powerful and useful, especially when combined
with other commands by using the pipe. Navigate to the `Bert`
directory. Every data file in this directory has a line which says
"Age". Lets list all of the ages from the survey data that Bert collected:

    grep Age *

Notice that this is a much shorter command than the `find` command we 
used to view the "Income" line of each file.

* * * *
**Short Exercise**

From the `data` directory, use a single command to view the `Income` 
line of every file in `THOMAS` and `jamesm`.

* * * *

## Count the Words

What if you want to understand the total size of files you have in 
particular locations? The `wc` program (word count) counts the number 
of lines, words, and characters in one or more files. Make sure you 
are in the `data`
directory, then enter the following command:

    wc Bert/* THOMAS/*4*

For each of the files indicated, `wc` has printed a line with three
numbers. The first is the number of lines in that file. The second is
the number of words. Finally, the total number of characters is
indicated. The final line contains this information summed over all of
the files. Thus, there were 4754 characters in total (~5060 on Windows,
which counts another character for each line ending). 

Remember that the `Bert/*` and `THOMAS/*4*` files were merged
into the `all_data` file. So, we should see that `all_data` contains
the same number of characters:

    wc all_data

Every character in the file takes up one byte of disk space. Thus, the
size of the file in bytes should also be 4754 (Windows: 5060). Let's confirm this:

    ls -l all_data

Remember that `ls -l` prints out detailed information about a file and
that the fifth column is the size of the file in bytes.

## The awesome power of the Pipe

Suppose I wanted to only see the total number of lines, words, and 
characters across the files `Bert/*` and `THOMAS/*4*`. I don't want to
see the individual counts, just the total. Because the `all_data` file 
is just a concatenation of this information, I could just do 
the following

    wc all_data

Sure, this works, but I had to create the `all_data` file to do this. Thus, I
have wasted a precious 4754 bytes (5060 on Windows) of hard disk space that is already used 
in the individual files. We can do the same count
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
`wc Bert/* THOMAS/*4*` and send it into the `tail -n 1` command. The `|`
character (called "pipe") is used for this purpose. Enter the following
command:

    wc Bert/* THOMAS/*4* | tail -n 1

This will print only the total number of lines, words, and characters 
across all of these files. What is happening here?

Well, the `|` character feeds output from the first program (to the left of 
the `|`) as input to the second program on the right. Therefor, you can 
string all sorts of commands together using the pipe. 

The philosophy behind the command line programs we're learning in the shell
is that none of them really do anything all that impressive. BUT when you 
start *chaining them together*, you can do some really powerful things really
efficiently. If you want to be proficient at using the shell, you must
learn to become proficient with the pipe and redirection operators:
`|`, `>`, `>>`.


## A sorting example

There is another useful program for working with files called `sort`. 
Let's create a file with some words to *sort* for the next example. We
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
name to the file and then sort the contents into a new file called `Sorted`.

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

# 6. Shell Scripts

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

Create an executable script called `smallestage` in the `data`
directory, that is similar to the `smallest` script, but prints the
line with the smallest "Age" value. Use the commands
`grep`, `sort`, and `head` to do this.

* * * * 

# 7. Automation via Shell Configuration

## Variables

We were earlier introduced to the `$PATH` variable.  This is a
variable that the shell expects to be able to function properly.
There are other variables defined and used by default.  To see a list
of variables just enter:

    set

To work with the information, pipe it through 
`less` like this:

    set | less

Now you will see a long list of variables that are already set. Some
of these are built-in to the bash shell, others are set for all users.

Some important variables you can expect to see:

| Variable | What it stores |
| -------- | -------------- |
| HOME     | the full path to your home directory |
| HOSTNAME | the name of this computer |
| PATH     | where programs are searched for |
| SHELL    | the full path to your current shell command |
| USER     | your user name (Linux, Mac) |
| USERNAME | your user name (Windows) |

You can use these variables in your commands.  For example, this is
yet another way to go to your home directory:

    cd $HOME

You can also define your own variables.  If you are always needing to
go to the same directory for your work, you could store it's path as a
variable so that you don't have to type out (and remember!) the whole 
path every time. Note that you can use one variable as part of the
definition of another. For example, `cd` to your home directory and try:

    DATADIR="$HOME/boot-camps/shell/data"
    cd $DATADIR

Let's make a directory for your data in the `$DATADIR`.

For Mac/Linux:

    mkdir $DATADIR/$USER

For Windows:

    mkdir $DATADIR/$USERNAME

##Example: The bash prompt
In `bash`, the prompt that you see at the beginning of each line is
just another variable, `$PS1`.  Examine the prompt:

    echo $PS1

and then update it to show a longer path relative to your home:

    PS1='[\u@\h \w]\$ '

Or, you can switch back, if you want, by using the output from the
`echo $PS1` command, making sure to leave a space after `$`.

An important note: When you define a variable on the command line, as 
we did just above, it is only available in this current shell session. 
If you later spawn a new shell from this one, or close this shell 
and open a new one later, you won't have the changes you made in this 
shell process. (Feel free to test this yourself by opening a second 
shell process!)

Later, we'll cover one way to have the default variable values set to 
something else (like our new value for PS1, above), whenever you open
a new shell process on *this* computer.

## Aliases

Another way to avoid having to retype paths and commands is to use an 
**alias**. Most shells allow you to define an alias so that it is available
to use in the current shell process.  Let's define an alias to perform the same
function as the `smallest` script that we made earlier.  If you enter:

    alias my_smallest='wc * | sort -k 3 -n | head -n 1'

you will have created a new command `my_smallest` that will be
available in any directory.

Want to know how to *color* your `ls` output? (The below works for 
Windows Git Bash, Mac, and *most* 'flavors' of Linux, but you can also 
search online for other color options pertaining to your computer's 
operating system and shell type.)

    alias ls='ls -G'

## Startup Scripts

Alas, none of the changes to variables or aliases will be there *when
you logout and come back*, as we described above.  Thankfully, we can 
turn to a script to make that happen.  The most common script is one 
that is already invoked every time you login (if it already exists). 
When using the bash shell it is called the `.bashrc` file. 

(On a Mac, it's called `.bash_profile`. Actually, the truth is a bit
more complicated; see [this discussion](http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html).)

Go to your home directory and open the `.bashrc` file (on a Mac, use
`.bash_profile`, which you'll have to create if it doesn't already exist):

    nano .bashrc

At the end of this file, add any of the new variables and prompts that we
defined above, if you WANT them on your computer.

    export DATADIR="$HOME/boot-camps/shell/data"
    export PS1="[\u@\h \w]\$ "
    alias my_smallest='wc * | sort -k 3 -n | head -n 1'
    alias ls='ls -G'

After exiting the shell and then opening a new session, these variables
should be available to you in the future. Many experienced users will 
gradually build up a long list of important variables and aliases for 
tasks they do frequently.

You can *always* regain the defaults by just removing these lines 
from the same startup script (.bashrc; or .bash_profile for Mac), which 
is always in your home directory, and then opening a new shell session.


* * * * 
Jump to look at the [solutions to all the Automation exercises.](ReadmeSolns.md)


## Bonus topics:
You may wish to look into the below topics/tools to become even more advanced in using the shell.
There are numerous learning resources online for each topic. We've linked to Software Carpentry's online
tutorials for some of them.:

[**ssh and scp**](http://software-carpentry.org/v4/shell/ssh.html) for connecting to remote servers

[**permissions**](http://software-carpentry.org/v4/shell/perm.html)

**backticks and variables**

**loops**

**du**

**regular expressions**

----

[Up To Schedule](../../README.md#schedule) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People](../../python/best_practice/Readme.md)
