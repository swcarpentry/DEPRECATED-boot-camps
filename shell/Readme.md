[Up To Schedule](../README.md) -
Forward To [Let the Computer Do the Work](automation/Readme.md)

# Introduction to The Shell

**Material by Paul Wilson, Lauren Michael, Milad Fatenejad, Sasha Wood, and Radhika Khetani**

# What is the shell how do I access the shell?

The *shell* is a program that presents a command line interface
which allows you to control your computer. Tasks are accomplished 
by entering commands with a keyboard instead of controlling graphical user interfaces
(GUIs or application "windows") that use a mouse/keyboard combination.

A *terminal* is a program you run that gives you access to the
shell. There are many different terminal programs that vary across
operating systems.
	 
There are many reasons to learn about the shell. In our opinion, the
most important reasons are that: 

1.  It is very common to encounter the shell and
    command-line-interfaces in scientific computing, so you will
    probably have to learn it eventually 

2.  The shell is a really powerful way of interacting with your
    computer. GUIs (clickable window interfaces) and the shell are 
    complementary - by knowing both
    you will greatly expand the range of tasks you can accomplish with
    your computer. You will also be able to perform many tasks more
    efficiently.

The shell is just a program and there are many different shell
programs that have been developed. The most common shell (and the one
we will use) is called the Bourne-Again SHell (bash). Even if bash is
not the default shell, it is usually installed on most systems unix-based 
systems (Mac or Linux operating systems), but Windows users can use something 
like Git Bash. Most commands, especially a
lot of the basic ones, work across the various shells but some things may 
be different. We recommend sticking with bash and learning it well.
([Here is a link for more information](http://en.wikipedia.org/wiki/Bash_(Unix_shell)))

To open a terminal on a Mac or Linux computer, just single click on the 
"Terminal" icon on the Desktop (or in Applications). For Windows, open 
the Git Bash program you installed for the bootcamp.

## The Example: Manipulating Experimental Data Files

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data from a hearing test. To get
the data for this test, you will need internet access. Just enter the
command:

    git clone -b 2014-08-25 https://github.com/UW-Madison-ACI/boot-camps.git

This command will grab all of the data needed for this workshop from our central GitHub repository. 
(We will talk about the `git` command later in the workshop.)

# Interacting with the Shell

One very basic command is `echo`. This command just prints text to
the terminal. Try the command:

    echo Hello, World

Then press enter. You should see the text "Hello, World" printed back
to you as "standard output". The echo command is useful for printing from 
a shell script, for displaying variables, and for generating known values 
to pass to other programs.

# 1. Navigating the file system

Let's learn how to move around your file system using command line
programs. This is really easy to do using a GUI (you just click on
things). Once you learn the basic commands, you'll see that it is
really easy to do in the shell too. 

First we have to know where we are. The program `pwd` (print working
directory) tells you where you are sitting in the directory tree. The
command `ls` will list the files and other directories in the current
directory. Directories are often called "folders" because of how they
are represented in GUIs. Directories are just organizational groupings 
of files. They can contain files or other directories.

Whenever you start up a terminal, you will start in a special
directory called the *home* directory. You can always
find out your user name by entering the command `whoami`.

## File Types

When you enter the `ls` command, you'll get back a list of contents of 
the current directory. There are several items in the home directory, 
notice that they are all colored blue (but not colored in Windows). These directories and files 
should be generally familiar to you.

Lets create an empty file using the `touch` command. Enter the
command:

    touch testfile

Then list the contents of the directory again. You should see that a
new entry, called `testfile`, exists.

Some terminals will not color the directory entries (though yours might, 
already). In those terminals, use `ls -F` instead of `ls`. The
`-F` argument modifies the results so that a slash is placed at the
end of directories. If the file is *executable* meaning that it can be
run like a program, then a star will be placed at the end of of the
file name.

You can also use the command `ls -l` to see whether items in a
directory are files or directories. `ls -l` gives a lot of other
information too, such as the size of the file and information about
the owner. If the entry is a directory, then the first letter will be
a "d". The fifth column shows you the size of the entries in
bytes. Notice that `testfile` has a size of zero.

Now, let's get rid of `testfile`. To remove a file, just enter the
command:

    rm testfile

The `rm` command can be used to remove files. If you enter `ls` again,
you will see that `testfile` is gone.


## Changing Directories

Now, let's move to a different directory. The command `cd` (change
directory) is used to move around. Let's move into the `boot-camps`
directory. Enter the following command:

    cd boot-camps

Use the `ls` command to see what is inside this directory.  This
directory contains all of the material for this boot camp. Now move to
the directory containing the data for the shell tutorial:

    cd shell

Now use the `ls` command to see what is inside this directory. Do you 
see files of different colors?

If you enter the `cd` command by itself, you will return to the home
directory. Try this, and then navigate back to the `shell`
directory.

## Arguments

Most programs take additional arguments that control their exact
behavior. For example, `-F` and `-l` are arguments to `ls`.  The `ls`
program, like many programs, can take a variety of arguments. But how do we
know what the options are to particular commands?

Most commonly-used shell programs have a manual. Unix-based (Mac and Linux) 
users can access the manual using the `man` program. Windows users may not 
have the `man` program, but can use Linux man pages online: http://www.linuxmanpages.com/ 
On Mac or Linux, try entering:

    man ls

This will open the manual page for `ls`. Use the <kbd>space</kbd> key to go
forward and <kbd>b</kbd> to go backwards. When you are done reading, just hit <kbd>q</kbd>
to quit.

Programs that are run from the shell can get extremely complicated. To
see an example, open up the manual page for the `find` program,
which we will use later this session. No one can possibly learn all of
these arguments, of course. So you will probably find yourself
referring back to the manual frequently.

* * * *
**Short Exercise**

Use the manual for `ls` to guess what you would expect from
using the arguments `-l`, `-t`, `-r` at the same time.

* * * *


## Examining the contents of other directories

By default, the `ls` commands lists the contents of the working
directory (i.e. the directory you are in). You can always find the
directory you are in using the `pwd` command. However, you can also
give `ls` the names of other directories to view. Navigate to the
home directory if you are not already there. Then enter the
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

The `cd` command takes an argument which is the directory
name. Directories can be specified using either a *relative* path or 
a *full* path. The directories on the computer are arranged into a
hierarchy. The full path tells you where a directory is in that
hierarchy. Navigate to your "home directory" with

    cd ~

Now, enter the `pwd` command and you should see:

the full name of *your* home directory. This tells you that you
are in a directory called `<username>`, and indicates the *full path* 
for that directory, starting with the top of the directory structure, 
which is indicated by `/`.

For example, you can enter a command to `cd` into `boot-camps/shell`, but 
using the *full* path based upon your output from `pwd`. How is the full 
path different for Windows versus unix-based computers?

Now go back to the home directory. We saw earlier that the
command:

    cd boot-camps/shell

had the same effect - it took us to the `shell` directory. But,
instead of specifying the full path, we specified a *relative path*. In
other words, we specified the path relative to our current
directory. A full path always starts with a `/`. A relative path does
not. You can usually use either a full path or a relative path
depending on what is most convenient for you. If we are in the home directory,
it is more convenient to just enter the relative path since it
involves less typing.

Over time, it will become easier for you to keep a mental note of the
structure of the directories that you are using and how to quickly
navigate amongst them.

* * * *
**Short Exercise**

Now, list the contents of a directory of your own, by using the full
path (and without using `cd`, so that you are running the command 
from your current location).

* * * * 

# Saving time with shortcuts, wild cards, and tab completion

## Shortcuts

There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
`~`, is a shortcut for your home directory, as we used it with `cd` before. 
Navigate to the `shell` directory, then enter the command:

    ls ~

This prints the contents of your home directory, without you having to
type the full path. The shortcut `..` always refers to the directory
above your current directory. Thus: 

    ls ..

prints the contents of the `boot-camps` folder when you are in `shell`. 
You can chain these together, so:

    ls ../../

prints the contents of your home directory from `shell`. 

Finally, the special directory `.` always refers to your
current directory. So, `ls`, `ls .`, and `ls ././././.` all do the
same thing, they print the contents of the current directory. This may
seem like a useless shortcut right now, but we'll see when it is
needed in a little while.

To summarize, while you are in the `shell` directory, the commands
`ls ~`, `ls ~/.`, and `ls ../../` all do exactly the
same thing. *Or* you could use the full path, too.

## Introducing our data set: Cochlear Implants

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
set of test results were written out to a text file, with one set per file.
Each participant has a unique subject ID, and a made-up subject name.
Each experiment has a unique experiment ID. The experiment has
collected 351 files so far.

The data is a bit of a mess! There are inconsistent file names, there
are extraneous "NOTES" files that we'd like to get rid of, and the
data is spread across many directories. We are going to use shell
commands to get this data into shape. By the end we would like to:

1.  Put all of the data into one directory called "cleaneddata"

2.  Have all of the data files in there, and ensure that every file
    has a ".txt" extension

3.  Get rid of the extraneous "NOTES" files

If we can get through this example in the available time, we will move
on to more advanced shell topics...

## Wild cards

Navigate to the `~/boot-camps/shell/data/THOMAS` directory. This
directory contains our hearing test data for THOMAS. If we type `ls`,
we will see that there are a bunch of files which are just four digit
numbers. By default, `ls` lists all of the files in a given
directory. The `*` character is a shortcut for "everything". Thus, if
you enter `ls *`, you will see all of the contents of a given
directory. Now try this command:

    ls *1

This lists every file that ends with a `1`. This command:

    ls *4*1

lists every file in the current directory whose name contains the
number `4`, and ends with the number `1` (so it has also 
established that `4` must be before `1` in the filename). There are four such files:
`0241`, `0341`, `0431`, and `0481`.

So how does this actually work? Well...when the shell (bash) sees a
word that contains the `*` character, it automatically looks for filenames
that match the given pattern. Using the wildcard does *NOT* 
exclude cases where there are no characters between `4` and `1`, as the 
wildcard includes cases where there are *any* number of characters
in place of the `*` character, even if there are zero characters. 

In this case, it identified four such
files. Then, it replaced the `*4*1` with the list of files, separated
by spaces. In other words, the two commands:

    ls *4*1
    ls 0241 0341 0431 0481

are exactly identical. The `ls` command cannot tell the difference
between these two things.

* * * *
**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the directories in `~/boot-camps/shell/data` that contain the letter `a` or the letter `e` (including files that may contain both). Hint: `ls -d` will list directories and not files.
2.  List all of the directories in `~/boot-camps/shell/data` that contain the letter `a` *AND* the letter `e`.

* * * *

## Tab Completion

Navigate to the home directory. Typing out directory names can waste a
lot of time. Instead, you can just start typing out the name of a directory, 
and then hit the tab key; the shell will try to fill in the rest of the
directory name you! For example, enter:

    cd b<tab>

The shell will fill in the rest of the directory name for
`boot-camps`. Now enter:

    ls s<tab><tab>

When you hit the first tab, nothing happens. The reason is that there
are *multiple* directories in the home directory which start with
`s`. Thus, the shell does not know which one to fill in. When you hit
tab again, the shell will list the possible choices. 

# Command History

You can easily access previous commands.  Hit the up arrow.  
Hit it again.  You can step backwards through your command history. 
The down arrow takes you forward in the command history.  

<kbd>CONTROL</kbd>+<kbd>c</kbd> will cancel the command you are writing, and give you a fresh prompt.

<kbd>CONTROL</kbd>+<kbd>r</kbd> will do a reverse-search through your command history. You can also 
use the up and down arrows to navigate through your command history, which can 
be usefule for easily repeating a recent command. Useful, right?

You can also review all of your recent commands at once with the `history` command.  Just enter:

    history

to see a numbered list of recent commands, including this just-issued
`history` command.

If your history looked like this:

    259  ls ../../
    260  cd ~/boot-camps/shell/data/THOMAS
    261  ls *1

then you could repeat command #260 by simply entering:

    !260

(that's an exclamation mark).

* * * * 
**Short Exercise**

1. Find the line number in your history for the last exercise (listing
directories in `~/boot-camps/shell/data`) and reissue that command.

* * * * 

# Which program?

Commands like `ls`, `rm`, `echo`, and `cd` are just ordinary programs
that exist on the computer (or that come along with Git Bash, for Windows 
users). A program is just a file that you can *execute*, otherwise 
known as an *executable*. If you want to find the location of a program, you 
can use `which`. For example:

    which rm

will return `/bin/rm` on unix-based computers. Thus, we can see that `rm` is a program that sits inside of the `/bin` directory. Now enter:

    which find

You will see that `find` is a program that sits inside of the
`/usr/bin` directory, if you're on a unix-based computer.

So when we enter a program name, like `ls`, and hit enter, how does
the shell know where to look for that program? How does it know to run
`/bin/ls` when we enter `ls`. The answer is that when we enter a
program name and hit enter, there are a few standard places that the
shell automatically looks. If it can't find the program in any of
those places, it will print an error saying "command not found". 

These standard places are stored in your *path*. To see your path, 
enter the following command:

    echo $PATH

This will print out the value of the `PATH` environment variable (more
on environment variables later...). Notice that a list of directories 
is displayed, separated by colon characters. These are the places the
shell looks for programs to run. 

If your program is not in this list, then an error is printed. For 
example, Navigate to the `shell` directory and list the contents. You will
notice that there is a program (executable file) called `hello` in
this directory. Now, try to run the program by entering:

    hello

You should get an error saying that hello cannot be found. That is
because the directory `~/boot-camps/shell` is not in the
`PATH`. Instead, you can run the `hello` program by entering:

    ./hello

Remember that `.` is a shortcut for the current working
directory. This tells the shell to run the `hello` program that is
located in your current location. So, you can run any program by entering the path
to that program. You can run `hello` equally well by specifying:

    ~/boot-camps/shell/hello

# Working with Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Enter the following command:

    cat ex_data.txt

This prints out the contents of the `ex_data.txt` file to the terminal 
as standard output. If you enter:

    cat ex_data.txt ex_data.txt

It will print out the contents of `ex_data.txt` twice. `cat` just
takes a list of file names and writes out their contents one after another (this
is where the name comes from, `cat` is short for "concatenate"). 

* * * *
**Short Exercises**

1.  Print out the contents of the `~/boot-camps/shell/dictionary.txt`
    file. What does this file contain?

2.  Without changing directories (you should still be in `shell`),
    use one short command to print the contents of all of the files in
    the `~/boot-camps/shell/data/THOMAS` directory.

* * * *

##Viewing file contents with `less`

`cat` is a terrific program, but when the file is really big, it can
be annoying to use. The program, `less`, is useful for this
case. Enter the following command:

    less ~/boot-camps/shell/dictionary.txt

`less` opens the file, and lets you navigate through it. The commands
are identical to the `man` program. 

**Some commands in `less`**

| key     | action |
| ------- | ---------- | 
| "space" | to go forward |
|  "b"    | to go backward |
|  "g"    | to go to the beginning |
|  "G"    | to go to the end |
|  "q"    | to quit |

`less` also gives you a way of searching through files. Just hit the
"/" key to begin a search, type the name of the word you would like
to search for, and hit `<enter>`. It will jump to the next location where
that word is found. Try searching the `dictionary.txt` file for the
word "cat". If you hit "/" then "enter", `less` will just repeat
the previous search. `less` only searches forward from the current 
location. Note: If you are at the end of the file and search
for the word "cat", `less` will not find it. You need to go to the
beginning of the file (by typing `g`) and then search. Now quit the less 
program.

Pro-tip: The `man` program actually uses `less` to show you the contents 
of the manual file for each program, so you can search program manuals 
using "/" and use all other less commands as well!

* * * *


## Redirection

Let's turn to the experimental data from the hearing tests that we
began with. This data is located in the `~/boot-camps/shell/data`
directory. Each subdirectory corresponds to a particular participant
in the study. Navigate to the `Bert` subdirectory in `data`. There
are a bunch of text files which contain experimental data
results. Lets print them all:

    cat au*

Now enter the following command:

    cat au* > ../all_data

This tells the shell to take the output from the `cat au*` command and
redirect it into a new file called `../all_data`. To verify that this
worked, examine the `all_data` file. If `all_data` had already
existed, we would have overwritten it because the `>` character tells the shell
to take the output from the command on the left and turn that into the contents
of the file on the right. The `>>` characters do almost the same thing,
except that they will *append* the output to the file if it already
exists.

* * * *
**Short Exercise**

Use `>>`, to append the contents of all of the files whose names
contain the number 4 in the directory `~/boot-camps/shell/data/gerdal` 
to the existing `all_data` file. Thus, when you are done, `all_data`
should contain all of the experiment data from Bert AND any
experimental data file from gerdal with filenames that contain the
number 4.

* * * *


## Creating, moving, copying, and removing

We first created a file called `all_data` using the redirection operator
`>`. This file is critical - it's our analysis results - so we want to
make copies so that the data is backed up.
Let's copy the file using the `cp` command. The `cp`
command backs up the file. Navigate to the `data` directory and enter:

    cp all_data all_data_backup

Let's make a temporary directory to store that file.  The `mkdir`
command is used to make a directory. Just enter `mkdir` followed by a
space, then the directory name:

    mkdir backup

This makes a directory with your username in the directory `backup`.

Now `all_data_backup` has been created as a copy of `all_data`. We can
move files and directories around using the command `mv`. Enter this command:

    mv all_data_backup backup/

This moves `all_data_backup` into your directory within `backup`. 

The `mv` command is also how you can rename files and directories. Since this file is so
important, let's rename it:

    mv all_data all_data_IMPORTANT

Now the file name has been changed to `all_data_IMPORTANT`. Let's delete
the backup file now:

    rm backup/all_data_backup

* * * *
**Short Exercise**

Do the following:

1.  Rename the `all_data_IMPORTANT` file back to `all_data`.
2.  Create a directory in the `data` directory called `foo`
3.  Then, *copy* the `all_data` file into `foo`

* * * *

By default, `rm`, will NOT delete directories. You can tell `rm` to
delete a directory and all of its contents using the `-r` option. 
Enter the following command:

    rm -r foo


## Finding files

The `find` program can be used to find files based on arbitrary
criteria. Navigate to the `shell` directory and enter the following
command:

    find . -print

This prints the name of every file or directory, recursively, starting
from the current directory. Let's exclude all of the directories:

    find . -type f -print

This tells `find` to locate only files. Now try these commands:

    find . -type f -name "*1*"
    find . -maxdepth 1 -type f -print
    find . -mindepth 2 -type f -print
    find . -type f -name "*1*" -print -o -name "*2*" -print
    find . -type f -name "*1*" -and -name "*2*" -print

The `find` command can acquire a list of files and perform some
operation on each file. Try this command out:

    find . -type f -exec grep Volume {} \;

This command finds every file within and below the current directory, 
and then searches each file for a line which contains the word "Volume". 
How does it do this? The `-exec` argument allow `find` to run to the 
program `grep`, multiple times, such that each file name is inserted 
whenever the `{}` occurs (as an argument to `grep`. The trailing `\;` is used to terminate the
command, in order to end the task run by `-exec`. 

We'll talk a bit more 
about grep after the break, but you'll use the `{}` trick in the 
exercises below.

## BONUS Topic: Using xargs to pass information to another program

The above command is slow, because it is calling a new instance
of `grep` for each item the `find` returns. A faster way to repeat 
the same task is to use the `xargs` command:

    find . -type f -print | xargs grep Volume

`find` generates a list of all the files we are interested in, 
then we *pipe* (`|`) them to `xargs`. `xargs` takes the items given to it 
and passes them as arguments to `grep`. `xargs` generally only creates
a single instance of `grep` (or whatever program it is running). We'll 
also use the *pipe* (`|`) more after the break.

* * * * 
**Exercises**

Let's clean up this data! Navigate to the `data` directory. Use a single 
`find` command to perform the below exercises. Number 2 does not require `find`:

1.  Find any file whose name is "NOTES" within `data` and its subdirectories, and delete it 

2.  Create a new directory called `cleaneddata`

3.  Copy all of the files (only) within the subdirectories of `data` into `cleaneddata`. (Hint: remember the wildcard. If you mess up, you can just delete the contents of cleaneddata, and try again.)

4.  Rename all of the files to ensure that they end in `.txt` (note:
    it is okay for certain files to end in `.txt.txt`, as some already end with `.txt`.)

**BONUS**

Redo exercise 4, except rename only the files which do not already end
in `.txt`. You will have to use the `man` command to figure out how to
search for files which do not match a certain name.

* * * * 

Jump to look at the [solutions to all the exercises.](ReadmeSolns.md)
