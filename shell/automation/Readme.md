[Up To Schedule](../../README.md) -
Back To [Introduction to the Shell](../Readme.md) - Forward To [Write Code for People](../../python/writing_code_for_people/Readme.md)

# Let the Computer Do the Work: Automating Workflows

**Material by Paul Wilson, Lauren Michael, Milad Fatenejad, Sasha Wood, and Radhika Khetani**

## Count the Words

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
`all_data`. (Hint: Where will the list of options for `wc` be?)

* * * *

## The awesome power of the Pipe

Suppose I wanted to only see the total number of lines, words, and 
characters across the files `Bert/*` and `gerdal/*4*`. I don't want to
see the individual counts, just the total. Because the `all_data` file 
is just a concatenation of this information, I could just do 
the following

    wc all_data

Sure, this works, but I had to create the `all_data` file to do this. Thus, I
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
this, we want to tell the shell to take the output of  
`wc Bert/* gerdal/*4*` and send it into the `tail -n 1` command. The `|`
character (called pipe) is used for this purpose. Enter the following
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
`/tmp/<username>` and enter the following command to create and 
add text to a new file called `toBeSorted`:

    nano toBeSorted

Now enter the four names as shown above. When you are done, press
<CONTROL>+O to write out these changes to the file. Press ENTER to use the file name
`toBeSorted`. Then press <CONTROL>+X to exit `nano`.

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
`../smallest`. Notice that it says "permission denied". This happens
because we haven't told the shell that this is an executable
file. If you do `ls -l ../smallest`, it will show you the permissions on 
the left of the listing.

Enter the following commands:

    chmod a+x ../smallest
    ../smallest

The `chmod` command is used to modify the permissions of a file. This
particular command modifies the file `../smallest` by giving all users
(notice the `a`) permission to execute the file (notice the `x`). If
you enter the following:

    ls -l ../smallest

you will see that the file name is now green and the permissions have changed. 
Congratulations, you just created your first shell script! You can now execute it
by entering the path location of `smallest` (absolute or relative) from within the 
directory you'd like to analyze.

## Searching Files

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
file containing the line with the smallest "Range" value. Use the commands
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

    DATADIR="/home/<username>/boot-camps/shell/data"
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

## Backticks to Capture Output

While redirection is useful to capture output in a file, you may want
to sometimes capture output in a variable.  The bash shell lets you do
that using so-called *backticks*.  These are the backwards apostrophes
that appear on a US English keyboard at the upperr left, under the `~`
(tilde).

A command inside a pair of these ticks is substituted with the
results of that command.  Since `pwd` told us the current directory,
we can capture that in a variable:

    THIS_DIR=`pwd`
    echo $THIS_DIR
    cd
    echo $THIS_DIR

# Aliases

Another way to avoid having to retype commands is to use an **alias**.
Most shells allow you to define an alias so that it is available
to use in that shell.  Let's define an alias to perform the same
funtion as the `smallest` script that we made earlier.  If you enter:

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

Now logout and login again, and confirm that those variables and
aliases are set with the `set` command and `alias` command,
respectively.

Many experienced users will gradually build up a long list of
important variables and aliases.

# Workflow Scripts

Perhaps we'd like to rearrange and combine the data we have in all of our data files so
that we can more easily import it into a spreadsheet.  One simple
format is known as a *comma separated values* or "CSV" format, which can be 
opened in Excel and other spreadsheet applications.  It takes all the
data from one record an puts it on a single line, separated by commas.

Let's start a new script in our data directory callaed `data2csv`:

    nano data2csv

To make this a shell script, the first line will be

    #!/bin/bash

Since we'll be able to capture any output using redirection, we'll
just write everything we want out to the standard output.

Our data is fixed so we'll write a standard header:

    echo '"Reported","Subject","Year/month of birth","Sex","CI type","Volume","Range","Discrimination"'

If this script is intended to process a single file, we'll need a way
to indicate which one.  When you run a shell script, you can provide arguments 
to the shell script that may be used as variables.  The first argument is used 
by the shell script as `$1`,
the second is `$2`, and so on.  Our script will be told which file to use 
from the first argument, and we'll save this filename in a variable to
remind us what it's for.

    datafile=$1

Now let's capture the data in other variables. 

    reported=`grep Reported $datafile | sed -e 's/Reported: \(.*\)$/\1/'`
    subject=`grep Subject $datafile | sed -e 's/Subject: \(.*\)$/\1/'`
    year_month=`grep Year $datafile | sed -e 's/Year\/month of birth: \(.*\)$/\1/'`
    sex=`grep Sex $datafile | sed -e 's/Sex: \(.*\)$/\1/'`
    type=`grep type $datafile | sed -e 's/CI type: \(.*\)$/\1/'`
    volume=`grep Volume $datafile | sed -e 's/Volume: \(.*\)$/\1/'`
    range=`grep Range $datafile | sed -e 's/Range: \(.*\)$/\1/'`
    discrimination=`grep Discrimination $datafile | sed -e 's/Discrimination: \(.*\)$/\1/'`

The `sed` program is one of many powerful tools to process text on the
command line.  We won't be going into much detail about sed today, but it 
is used for this script to combine the data from our data files into appropriate 
fields within a csv file (that can later be treated as columns in Excel).  
The later Python sessions will introduce other
tools for this kind of text processing.

We now need to write out a CSV line for this file:

    echo \"$reported\", \"$subject\", \"$year_month\", $sex, $type, $volume, $range, $discrimination

The escape character `\` before the quotation marks ensures that the quotation 
marks will appear in the file.

Now save the file and exit nano before making the file executable:

    chmod u+x data2csv

Now let's try it for the `0213` file in `THOMAS`:

    cat THOMAS/0213
    ./data2csv THOMAS/0213

## Loops in Shell Scripts

Let's add a loop so that we can make a single CSV file from many data
files.  The general form of a loop is:

   for var in list
   do
      act on $var
   done

In our case, the `list` will be a list of filenames that will come from
the arguments to the script, known in bash as `$@`.  The action for
each file will be the parts necessary to generate a CSV line.

Open your script again:

    nano data2csv

Replace the line that set the `$datafile` variable with

    for datafile in $@
    do

And after the last `echo` command add:

    done

After we save and exit, we can now run that script on the entire set of 
files in `THOMAS`:

    ./data2csv THOMAS/*


* * * * 
Jump to look at the [solutions to all the Automation exercises.](ReadmeSolns.md)


## Bonus topics:
You may wish to look into the below topics/tools to become even more advanced in using the shell:


[**ssh and scp**](http://software-carpentry.org/v4/shell/ssh.html)

[**Permissions**](http://software-carpentry.org/v4/shell/perm.html)

**du**

**ln**

**Regular Expressions**
