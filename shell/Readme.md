# The Shell

**Material by Milad Fatenejad, Sasha Wood, and Radhika Khetani**
* modified by Tracy Teal *

# The Shell

# What is the shell and how do I access the shell?

The *shell* is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

Use a browser to open the tutorial on github, located at:
    https://github.com/swcarpentry/boot-camps/blob/2013-05-umass/

Click on the directory named `shell`.

A *terminal* is a program you run that gives you access to the
shell. There are many different terminal programs that vary across
operating systems.

To get started, you will need the data and notes for the bootcamp.  Just 
enter the command:

      git clone -b 2013-05-ucdavis https://github.com/swcarpentry/boot-camps.git

More materials from SWC on the shell
http://software-carpentry.org/4_0/shell/index.html

Shell cheat sheet
https://github.com/swcarpentry/boot-camps/blob/2013-05-umass/shell/shell_cheatsheet.md



Topics in the shell:

1.  The hierarchical file structure
    'Getting around at the command line'

2.  Redirection
    'Input and output'


# The heirarchical file structure

Unix has a hierarchical file structure.

## Knowing where you are

Show where you are in the file system.  
print working directory

   pwd

If you just started up the terminal you're in the 'home' directory.
Type 'pwd' and see what your home directory is.

## Knowing what is in your directory

To see what is in your current directory, you list the files

   ls

This will show all the files and directories where you are.

Type 'ls' and you should see a bunch of things, including the 'boot-camps' folder.

However, when we look at all these things, we don't know if they're files or
folders or programs.  We can give arguments to 'ls' to tell it to give us more 
information.  

   ls -F

modifies the results so that a slash is placed at the
end of directories. If the file is *executable* meaning that it can be
run like a program, then a star will be placed at the end of of the
file name.

   ls -l 

gives a lot more
information too, such as the size of the file and information about
the owner. If the entry is a directory, then the first letter will be
a "d". The fifth column shows you the size of the entries in
bytes.

You can combine them, to see everything at once.

   ls -lF

If you ever wonder what types of options there are for a command you can
use Unix's built-in help and type:

   man ls

This will open the manual page for `ls`. Use the space key to go
forward and b to go backwards. When you are done reading, just hit `q`
to exit.  This will work for any command.


Exercise:  What does ls -alF give you?  Are there files you didn't see
with ls -lF?


** Getting around the file system  **

** Changing directories **

Go in to the boot camps directory by changing directories

   cd boot-camps

See what's in this directory using 'ls'.  Now go in to the 'shell' directory 
and see what's in there.


you can also
give `ls` the names of other directories to view. Go to your home directory
by typing 'cd'.  Then enter the
command:

    ls boot-camps

This will list the contents of the `boot-camps` directory without
you having to navigate there. Now enter:

    ls boot-camps/shell

This prints the contents of `shell`. The `cd` command works in a
similar way. Try entering:

    cd boot-camps/shell

and you will jump directly to `shell` without having to go through
the intermediate directory.


## Full vs. Relative Paths

So far we've just been using relative paths.  This means that we've been using the information about where we are to get somewhere else.  You can also 
use absolute paths.

Go back to your home directory by typing 'cd'

Type

  cd boot-camps/shell

Then type

  cd /HOMEDIRECTORY/boot-camps/shell

Are you in the same place?  Do you see why?

**Shortcuts**

There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
`~`, is a shortcut for your home directory. Navigate to the `shell`
directory, then enter the command:

    ls ~

This prints the contents of your home directory, without you having to
type the full path. The shortcut `..` always refers to the directory
above your current directory. Thus:

    ls ..

prints the contents of the `/home/swc/boot-camps`. You can chain
these together, so:

    ls ../../

prints the contents of `/home/swc` which is your home
directory. Finally, the special directory `.` always refers to your
current directory, so:

    ls .

print what is in your current directory.


# The Example: Manipulating Experimental Data Files


**Our data set: Cochlear Implants**

A cochlear implant is a small electronic device that is surgically
implanted in the inner ear to give deaf people a sense of
hearing. More than a quarter of a million people have them, but there
is still no widely-accepted benchmark to measure their effectiveness.
In order to establish a baseline for such a benchmark, our supervisor
got teenagers with CIs to listen to audio files on their computer and
report:

1.  the quietest sound they could hear
2.  the lowest and highest tones they could hear
3.  the narrowest range of frequencies they could discriminate

To participate, subjects attended our laboratory and one of our lab
techs played an audio sample, and recorded their data - when they
first heard the sound, or first heard a difference in the sound.  Each
set of test results were written out to a text file, one set per file.
Each participant has a unique subject ID, and a made-up subject name.
Each experiment has a unique experiment ID. The experiment has
collected 351 files so far.

The data is a bit of a mess! There are inconsistent file names, there
are extraneous "NOTES" files that we'd like to get rid of, and the
data is spread across many directories. We are going to use shell
commands to get this data into shape. By the end we would like to:

1.  Put all of the data into one directory called "alldata"

2.  Have all of the data files in there, and ensure that every file
    has a ".txt" extension

3.  Get rid of the extraneous "NOTES" files

If we can get through this example in the available time, we will move
onto more advanced shell topics...

**Wild cards**

Navigate to the `~/boot-camps/shell/data/THOMAS` directory. This
directory contains our hearing test data for THOMAS. If we type `ls`,
we will see that there are a bunch of files which are just four digit
numbers. By default, `ls` lists all of the files in a given
directory. The `*` character is a shortcut for "everything". Thus, if
you enter `ls *`, you will see all of the contents of a given
directory. Now try this command:

    ls *1

This lists every file that ends with a `1`. This command:

    ls /usr/bin/*.sh

Lists every file in `/usr/bin` that ends in the characters `.sh`. And
this command:

    ls *4*1

lists every file in the current directory which contains the number
`4`, and ends with the number `1`. There are four such files: `0241`,
`0341`, `0431`, and `0481`. 

So how does this actually work? Well...when the shell (bash) sees a
word that contains the `*` character, it automatically looks for files
that match the given pattern. In this case, it identified four such
files. Then, it replaced the `*4*1` with the list of files, separated
by spaces. In other the two commands:

    ls *4*1
    ls 0241 0341 0431 0481

are exactly identical. The `ls` command cannot tell the difference
between these two things.

* * * *
**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that start with the letter 'c'
2.  List all of the files in '/bin' that contain the letter 'a'
3.  List all of the files in '/bin' that end in 'o'

BONUS:
List all of the files in `/bin` that contain the letter `a` or the letter `t`


* * * *

**Tab Completion**

Navigate to the home directory. Typing out directory names can waste a
lot of time. When you start typing out the name of a directory, then
hit the tab key, the shell will try to fill in the rest of the
directory name. For example, enter:

    cd b<tab>

The shell will fill in the rest of the directory name for
`boot-camps`. Now enter:

    ls s<tab><tab>

When you hit the first tab, nothing happens. The reason is that there
are multiple directories in the home directory which start with
s. Thus, the shell does not know which one to fill in. When you hit
tab again, the shell will list the possible choices. 

Tab completion can also fill in the names of programs. For example,
enter `e<tab><tab>`. You will see the name of every program that
starts with an `e`. One of those is `echo`. If you enter `ech<tab>` you
will see that tab completion works.

**Command History**

You can easily access previous commands.  Hit the up arrow.  
Hit it again.  You can step backwards through your command history. 
The down arrow takes your forwards in the command history.  

^-C will cancel the command you are writing, and give you a fresh prompt.

^-R will do a reverse-search through your command history.  This
is very useful.

## Which program? ##

Commands like `ls`, `rm`, `echo`, and `cd` are just ordinary programs
on the computer. A program is just a file that you can *execute*. The
program `which` tells you the location of a particular program. For
example:

    which ls

Will return "/bin/ls". Thus, we can see that `ls` is a program that
sits inside of the `/bin` directory. Now enter:

    which find

You will see that `find` is a program that sits inside of the
`/usr/bin` directory.

So ... when we enter a program name, like `ls`, and hit enter, how
does the shell know where to look for that program? How does it know
to run `/bin/ls` when we enter `ls`. The answer is that when we enter
a program name and hit enter, there are a few standard places that the
shell automatically looks. If it can't find the program in any of
those places, it will print an error saying "command not found". Enter
the command:

    echo $PATH

This will print out the value of the `PATH` environment variable. More
on environment variables later. Notice that a list of directories,
separated by colon characters, is listed. These are the places the
shell looks for programs to run. If your program is not in this list,
then an error is printed. The shell ONLY checks in the places listed
in the `PATH` environment variable. 

Navigate to the `shell` directory and list the contents. You will
notice that there is a program (executable file) called `hello` in
this directory. Now, try to run the program by entering:

    hello

You should get an error saying that hello cannot be found. That is
because the directory `/home/swc/boot-camps/shell` is not in the
`PATH`. You can run the `hello` program by entering:

    ./hello

Remember that `.` is a shortcut for the current working
directory. This tells the shell to run the `hello` program which is
located right here. So, you can run any program by entering the path
to that program. You can run `hello` equally well by specifying:

    /home/swc/boot-camps/shell/hello

Or by entering:

    ../shell/hello

When there are no `/` characters, the shell assumes you want to look
in one of the default places for the program.


## Examining Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Enter the following command:

    cat ex_data.txt

This prints out the contents of the `ex_data.txt` file. If you enter:

    cat ex_data.txt ex_data.txt

It will print out the contents of `ex_data.txt` twice. `cat` just
takes a list of file names and writes them out one after another (this
is where the name comes from, `cat` is short for concatenate). 

* * * *
**Short Exercises**

1.  Print out the contents of the `~/boot-camps/shell/dictionary.txt`
    file. What does this file contain?

2.  Without changing directories, (you should still be in `shell`),
    use one short command to print the contents of all of the files in
    the `/home/swc/boot-camps/shell/data/THOMAS` directory.

* * * *

`cat` is a terrific program, but when the file is really big, it can
be annoying to use. The program, `less`, is useful for this
case. Enter the following command:

    less ~/boot-camps/shell/dictionary.txt

`less` opens the file, and lets you navigate through it. The commands
are identical to the `man` program. Use "space" to go forward and hit
the "b" key to go backwards. The "g" key goes to the beginning of the
file and "G" goes to the end. Finally, hit "q" to quit.

`less` also gives you a way of searching through files. Just hit the
"/" key to begin a search. Enter the name of the word you would like
to search for and hit enter. It will jump to the next location where
that word is found. Try searching the `dictionary.txt` file for the
word "cat". If you hit "/" then "enter", `less` will just repeat
the previous search. `less` searches from the current location and
works its way forward. If you are at the end of the file and search
for the word "cat", `less` will not find it. You need to go to the
beginning of the file and search.

Remember, the `man` program uses the same commands, so you can search
documentation using "/" as well!


* * * * 


## Redirection

Let's turn to the experimental data from the hearing tests that we
began with. This data is located in the `~/boot-camps/shell/data`
directory. Each subdirectory corresponds to a particular participant
in the study. Navigate to the `Bert` subdirectory in `data`.  There
are a bunch of text files which contain experimental data
results. Lets print them all:

    cat au*

Now enter the following command:

    cat au* > ../all_data

This tells the shell to take the output from the `cat au*` command and
dump it into a new file called `../all_data`. To verify that this
worked, examine the `all_data` file. If `all_data` had already
existed, we would overwritten it. So the `>` character tells the shell
to take the output from what ever is on the left and dump it into the
file on the right. The `>>` characters do almost the same thing,
except that they will append the output to the file if it already
exists.

* * * *
**Short Exercise**

Use `>>`, to append the contents of all of the files which contain the
number 4 in the directory:

    /home/swc/boot-camps/shell/data/gerdal

to the existing `all_data` file. Thus, when you are done `all_data`
should contain all of the experiment data from Bert and any
experimental data file from gerdal that contains the number 4.

* * * *


## Creating, moving, copying, and removing

We've created a file called `all_data` using the redirection operator
`>`. This file is critical - it's our analysis results - so we want to
make copies so that the data is backed up.
Lets copy the file using the `cp` command. The `cp`
command backs up the file. Navigate to the `data` directory and enter:

    cp all_data all_data_backup

Now `all_data_backup` has been created as a copy of `all_data`. We can
move files around using the command `mv`. Enter this command:

    mv all_data_backup /tmp/

This moves `all_data_backup` into the directory `/tmp`. The directory
`/tmp` is a special directory that all users can write to. It is a
temporary place for storing files. Data stored in `/tmp` is
automatically deleted when the computer shuts down.

The `mv` command is also how you rename files. Since this file is so
important, let's rename it:

    mv all_data all_data_IMPORTANT

Now the file name has been changed to all_data_IMPORTANT. Let's delete
the backup file now:

    rm /tmp/all_data_backup

The `mkdir` command is used to create a directory. Just enter `mkdir`
followed by a space, then the directory name. 

* * * *
**Short Exercise**

Do the following:

1.  Rename the `all_data_IMPORTANT` file to `all_data`.
2.  Create a directory in the `data` directory called `foo`
3.  Then, copy the `all_data` file into `foo`

* * * *

By default, `rm`, will NOT delete directories. You can tell `rm` to
delete a directory using the `-r` option. Enter the following command:

    rm -r foo


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

## The awesome power of the Pipe

Suppose I wanted to only see the total number of character, words, and
lines across the files `Bert/*` and `gerdal/*4*`. I don't want to
see the individual counts, just the total. Of course, I could just do:

    wc all_data

Since this file is a concatenation of the smaller files. Sure, this
works, but I had to create the `all_data` file to do this. Thus, I
have wasted a precious 7062 bytes of hard disk space. We can do this
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


**A sorting example**

Let's create a file with some words to sort for the next example. We
want to create a file which contains the following names:

    Bob
    Alice
    Diane
    Charles

To do this, we need a program which allows us to create text
files. There are many such programs, the easiest one which is
installed on almost all systems is called `nano`. Navigate to `/tmp`
and enter the following command:

    vi toBeSorted

Vi is a very useful text editor to know for the shell.  

    vi        opens the editor
    i         allows you to edit in the file
    Esc :wq   quits and saves the file
    Esc :q!   quits without saving

Now enter the four names as shown above. When you are done, 
write out the file with wq. 


When you are back to the command line, enter the command:

    sort toBeSorted

Notice that the names are now printed in alphabetical order.

* * * *
**Short Exercise**

Use the `echo` command and the append operator, `>>`, to append your
name to the file, then sort it and make a new file called Sorted.

* * * *


* * * * 
**Short Exercise**

Combine the `wc`, `sort`, `head` and `tail` commands so that only the
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


# For Future Reference

# Finding files

The `find` program can be used to find files based on arbitrary
criteria. Navigate to the `data` directory and enter the following
command:

    find . -print

This prints the name of every file or directory, recursively, starting
from the current directory. Let's exclude all of the directories:

    find . -type f -print

This tells `find` to locate only files. Now try these commands:

    find . -type f -name "*1*"
    find . -type f -name "*1*" -or -name "*2*" -print
    find . -type f -name "*1*" -and -name "*2*" -print

The `find` command can acquire a list of files and perform some
operation on each file. Try this command out:

    find . -type f -exec grep Volume {} \;

This command finds every file starting from `.`. Then it searches each
file for a line which contains the word "Volume". The `{}` refers to
the name of each file. The trailing `\;` is used to terminate the
command.  This command is slow, because it is calling a new instance
of `grep` for each item the `find` returns.

A faster way to do this is to use the `xargs` command:

    find . -type f -print | xargs grep Volume

`find` generates a list of all the files we are interested in, 
then we pipe them to `xargs`.  `xargs` takes the items given to it 
and passes them as arguments to `grep`.  `xargs` generally only creates
a single instance of `grep` (or whatever program it is running).

* * * * 
**Short Exercise**

Navigate to the `data` directory. Use one `find` command to perform each
of the operations listed below (except number 2, which does not
require a `find` command):

1.  Find any file whose name is "NOTES" within `data` and delete it 

2.  Create a new directory called `cleaneddata`

3.  Move all of the files within `data` to the `cleaneddata` directory

4.  Rename all of the files to ensure that they end in `.txt` (note:
    it is ok for the file name to end in `.txt.txt`

Hint: If you make a mistake and need to start over just do the
following:

1.  Navigate to the `shell` directory

2.  Delete the `data` directory

3.  Enter the command: `git checkout -- data` You should see that the
    data directory has reappeared in its original state

**BONUS**

Redo exercise 4, except rename only the files which do not already end
in `.txt`. You will have to use the `man` command to figure out how to
search for files which do not match a certain name. 

* * * * 



















