# The Shell

# What is the shell? How do I access the shell?

A *shell* is a program that allows a user to interact with a computer by
running programs, sending input to those programs (via the keyboard, mouse,
etc.) and display output from those programs. This typically features programs
for browsing the filesystem and viewing files and directory contents (programs
themselves are files with the special feature that they contain code understood
by the operating system as executable instructions).

There are many different shell programs.  On Windows and OSX the GUI-based
desktop/file browser inteface typically presented users *is* a kind of
shell--it allows you to double-click on icons representing programs in order to
execute those programs, and view the files on your computer.

More traditionally-speaking, however, a shell is a text-based program
implementing what you've probably heard called a "command line".  The "command
line" really just refers to a single line on the screen (typically at the
bottom), in which the user types the name of a program they want to run,
usually with some options (such as the name of a file they want to open with
that program).  The program may or may not output some text (or open a window
if it's a GUI program), then present the user with a new command line.  On most
modern UNIX shells the command line is represented by a dollar sign:

    $

In most tutorials or help on shells, this dollar sign typically signifies that
one should type everything *after* the dollar sign into the command line.

A shell is just a program among many programs on a computer.  In fact, on older
machines the first program that runs when connecting a terminal is *not* the
shell, but rather the login prompt--the screen that asks you for your username
and password is just another program too.  Upon successful login the login
prompt then runs that user's preferred shell program.

Bash (which stands for the Bourne Again SHell) is one of many text-based shell
programs, and is the most popular in use today.  It is the default shell on Mac
OSX and on most modern Linux systems.  It is also the shell bundled with Git on
Windows, since it is more effective for working with programs like Git that
were designed for UNIX-like systems than Windows' own text shell which is based
on DOS.  Since Bash is just another program, even if it is not your default
shell you can always run it by entering

    $ bash

at your command line.  This will replace your previous shell with a bash shell.
If you type 'exit' that will quit bash and return you to your previous shell.
You can even run a bash shell inside a bash shell (this can be useful for many
reasons).


# Shell vs. Terminal

Sometimes you'll hear the words "shell" and "terminal" used interchangeably.
Indeed, on Mac OSX and most Linux systems the program you run to get a command
line shell is called "Terminal".  But in the olden days a terminal was
literally a machine you sat at:

<img alt="A man sitting at a DECwriter teletype terminal" src="https://raw.github.com/swcarpentry/boot-camps/2013-07-cmu/shell/images/decwriter-la36.jpg" />

The earliest computers were interacted with through primitive devices
such as switches and lights and oscilloscope.  But very soon it became common
to interact with mainframe with teletype machines, at which one typed commands
on a typewriter that converted the keypresses to signals to the computer.  The
output from the computer then controlled the typewriter to display.  By the
late 70s these started to be replaced with actual CRT screens, but those too
were still just dumb terminals that sent input to, and displayed output from
much bigger multi-user computers over a serial port.

Today use of multi-user computers is less common, and most computers meant for
humans to do their work on boot into GUI interfaces.  (This is opposed to
servers and supercomputers that spend most of their time working without human
interaction and may still just have a text-based shell.)  Since the text-based
shells we still use today were developed for users of terminals they still
carry some of that history with them.  So to get to a shell from our GUI
interfaces we run a Terminal application that *emulates* the old text-based
terminals of yore and runs the shell program inside them.


# Why use a shell?

So why are we teaching you a computing technology that was fashionable when
bell bottoms were still in?

There are many reasons to learn about the shell. In my opinion, the most
important reasons are that:

1.  It is very common to encounter the shell and command-line-interfaces in
    scientific computing.  This is because writing GUI applications is hard,
    and takes a lot of effort that usually is not germane to producing
    scientific results.  Much of the software written by and for scientists
    is text-only and run at a command line.  So it's good to get comfortable
    with it.

2.  The shell is a really powerful way of interacting with your computer. GUIs
    and the shell are complementary - by knowing both you will greatly expand
    the range of tasks you can accomplish with your computer. You will also be
    able to perform many tasks more efficiently.  The shell is especially good
    for performing repetitive tasks that are hard to program a GUI to do.

3.  Because of the first two reasons most software written by and for
    developers (especially in the open source community) are command
    line-based.  Some of this is driven by the preference of programmers, who
    are already efficient users of the shell.  And part of it is out of
    practicality--not much is gained by writing a complicated GUI wrapper to
    your C compiler itself (though there *are* advantages to GUI-based code
    editors with built-in interfaces to the compiler).

# Let's get started

One very basic command is `echo`. This command just prints text to
the terminal. Try the command:

    echo Hello, World

Then press enter. You should see the text "Hello, World" printed back
to you. The echo command is useful for printing from a shell script,
for displaying variables, and for generating known values to pass
to other programs.

## Moving around the file system

Let's learn how to move around the file system using command line
programs. This is really easy to do using a GUI (just click on
things). Once you learn the basic commands, you'll see that it is
really easy to do in the shell too.

First we have to know where we are. The program `pwd` (print working
directory) tells you where you are sitting in the directory tree. The
command `ls` will list the files in files in the current
directory. Directories are often called "folders" because of how they
are represented in GUIs. Directories are just listings of files. They
can contain other files or directories.

Whenever you start up a terminal, it will start your shell in a special
directory called the *home* directory. Every user has their own home
directory where they have full access to do whatever they want. This
is more relevant on multi-user machines used for research ("back in
the day" all UNIX-like computers were multi-users).  Users can normally
only read and write files in their home directory--this gives them
privacy from other users.

* On Linux your home directory is usually `/home/<your_username>`
* On Mac OSX the home directory is usually `/Users/<your_username>`
* On Windows the situation is a little more complicated: The closest
  thing Windows has natively to a home directory is
  `/C/Documents and Settings/<your_username>`, but if you run a bash
  shell with Git Bash it also sets up a UNIX-like home directory for
  you. The exact location of this varies, but it is typically
  something like `/C/Program Files/Git/home/<your_username>`

You can always find out your user name by entering the
command `whoami`.

## File Types

When you enter the `ls` command lists the contents of the current
directory. There are several items in the home directory, notice that
they are all colored blue. This tells us that all of these items are
directories as opposed to files.

Lets create an empty file using the `touch` command. Enter the
command:

    touch testfile

Then list the contents of the directory again. You should see that a
new entry, called `testfile`, exists. It is colored white meaning that
it is a file, as opposed to a directory. The `touch` command just
creates an empty file.

Some terminals will not color the directory entries in this very
convenient way. In those terminals, use `ls -F` instead of `ls`. The
`-F` argument modifies the results so that a slash is placed at the
end of directories. If the file is *executable* meaning that it can be
run like a program, then a star will be placed at the end of of the
file name.

You can also use the command `ls -l` to see whether items in a
directory are files or directories. `ls -l` gives a lot more
information too, such as the size of the file and information about
the owner. If the entry is a directory, then the first letter will be
a "d". The fifth column shows you the size of the entries in
bytes. Notice that `testfile` has a size of zero.

Now, let's get rid of `testfile`. To remove a file, just enter the
command:

    rm testfile

The `rm` command can be used to remove files. If you enter `ls` again,
you will see that `testfile` is gone.


## Experimental Data Files

Let's download a small set of sample data files to use in some of our
excersizes.  Are sample files are text descriptions of molecules in
the Protein Data Bank (PDB) format--it's a simple format and one we
can read easily-enough. But it's also an example of the kind of thing
scientists have to deal with regularly--specially made-up data formats
specialized to one domain.  We need to read files like these all the
time, but when making up your own formats please try not to do this--
if possible find a way to model your data with a standard data format
like XML or JSON.

Enough preaching: To get the sample files we can try using a command
line tool like `curl` or `wget` to download files from the web to your
current working directory.

First run `cd` to get back to your home directory.

* On Mac OSX run:
      
    curl -O https://github.com/swcarpentry/boot-camps/raw/2013-07-cmu/shell/shell-data.zip

* On Linux and on Windows (in Git Bash) run:

    wget https://github.com/swcarpentry/boot-camps/raw/2013-07-cmu/shell/shell-data.zip
      
OSX systems should have `curl`.  Linux will usually have `wget` but might also have `curl`.
If neither of those work, try downloading the file by entering the URL in your web browser.
But for doing the exercises make sure to save the file in your home directory!

Now unzip the shell-data.zip file by running:

    unzip shell-data.zip
    
(Note: Most versions of the `unzip` command are smart enough that if you have a file called
`shell-data.zip` you can just write `unzip shell-data` and it will fill in the `.zip` extension
automatically.)


## Changing Directories

Unzipping the `shell-data.zip` file should have created a new directory under
our home directory called `shell-data/`.  Use the `ls -F` command from earlier
to check this.

Now, let's move to a different directory. The command `cd` (change
directory) is used to move around. Let's move into the `shell-data`
directory. Enter the following command:

    cd shell-data

Use the `ls` command to see what is inside this directory. Before we learn
version control you will want to get comfortable using the `cd` and `ls`
commands almost without thinking to move around and explore the file system.


## Arguments

Most programs take additional arguments that control their exact
behavior. For example, `-F` and `-l` are arguments to `ls`.  The `ls`
program, like many programs, take a lot of arguments. But how do we
know what the options are to particular commands?

The first way to get a quick overview is to use the `--help` argument
to a command. Not all commands support this, but it's standard
practice for well-written command-line programs (something to keep in
mind if/when you write your own command-line programs).  Usually
`--help` will give a brief overview of what the commonly used arguments
do.  But to get a more detailed explanation of how to use a command
one typically finds more in the "manpage" (short for "manual page").

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

    man ls

This will open the manual page for `ls`. Use the space key to go
forward and b to go backwards. When you are done reading, just hit `q`
to quit.

Note to Window users: Git Bash does *not* come with man pages by default,
so the `man` command won't work. But you can always google the man 
pages for common commands (for example, just google "man ls").

Programs that are run from the shell can get extremely complicated. To
see an example, open up the manual page for the `find` program,
which we will (time-permitting) use later this session. No one can
possibly learn all of these arguments, of course. So you will 
probably find yourself referring back to the manual page frequently.

* * * *
**Short Exercise**

1. Use the manual page for `ls` to guess what you would expect from
using the arguments `-l`, `-t`, `-r` at the same time (remember: Windows
users just google it or look on your neighbor's screen).

2. Try the following and see if you can figure out what they do, either by examining the results or consulting the manual page.
   * `ls -lS` (equivalent to `ls -l -S`)
   * `ls -lt` (equivalent to `ls -l -t`)
   * `ls -1`  (that's the number one, not a letter 'ell')

* * * *


## Examining the contents of other directories

By default, the `ls` commands lists the contents of the working
directory (i.e. the directory you are in). You can always find the
directory you are in using the `pwd` command. However, you can also
give `ls` the names of other directories to view. Navigate to the
home directory if you are not already there. Then enter the
command:

    ls shell-data

This will list the contents of the `molecules` directory without
you having to navigate there. Now enter:

    ls shell-data/molecules

This prints the contents of `molecules/`. The `cd` command works in a
similar way. Try entering:

    cd shell-data/molecules

and you will jump directly to `molecules/` without having to go through
the intermediate directory.

## Full vs. Relative Paths

The `cd` command takes an argument which is the directory
name. Directories can be specified using either a *relative* path a
full *path*. The directories on the computer are arranged into a
hierarchy. The full path tells you where a directory is in that
hierarchy. Navigate to the home directory. Now, enter the `pwd`
command and you should see something like

    /home/<your_username>

which is the full name of your home directory, as we saw before.
This tells you that you are in a directory called `<your_username>`,
which sits
inside a directory called `home` (or `Users` on OSX, etc.) which sits
inside the very top directory in the hierarchy. The
very top of the hierarchy is a directory called `/` which is usually
referred to as the *root directory*. So, to summarize: `<your_username>`
is a directory in `home` which is a directory in `/`.

Now enter a command like:

    cd /home/<your_username>/shell-data/molecules
    
The `/home/<your_username>/` part will differ for your computer--replace
it with whatever was returned by the `pwd` command.

This jumps to `molecules`. Now go back to the home directory. We saw
earlier that the command:

    cd shell-data/molecules

had the same effect - it took us to the `shell` directory. But,
instead of specifying the full path
(`/home/<your_username>/boot-camps/shell`), we specified a *relative path*. In
other words, we specified the path relative to our current
directory. A full path always starts with a `/`. A relative path does
not. You can usually use either a full path or a relative path
depending on what is most convenient. If we are in the home directory,
it is more convenient to just enter the relative path since it
involves less typing.

Over time, it will become easier for you to keep a mental note of the
structure of the directories that you are using hand how to quickly
navigate amongst them.

* * * *
**Short Exercise**

Now, list the contents of the /bin directory. Do you see anything
familiar in there?

* * * *

## Saving time with shortcuts, wild cards, and tab completion

### Shortcuts

There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
`~`, is a shortcut for your home directory. Navigate to the
`shell-data/molecules` directory, then enter the command:

    ls ~

This prints the contents of your home directory, without you having to
type the full path. The shortcut `..` always refers to the directory
above your current directory. Thus:

    ls ..

prints the contents of the `/home/<your_username>/shell-data`. You can chain
these together, so:

    ls ../../

prints the contents of your home
directory. Finally, the special directory `.` always refers to your
current directory. So, `ls`, `ls .`, and `ls ././././.` all do the
same thing, they print the contents of the current directory. This may
seem like a useless shortcut right now, but we'll see when it is
needed in a little while.

To summarize, while you are in the `molecules/` directory, the commands
`ls ~`, `ls ~/.`, `ls ../../`, and `ls /home/<your_username>` all do exactly the
same thing. These shortcuts are not necessary, they are provided for
your convenience.


### Wild cards

Navigate to the `~/shell-data/molecules/` directory. If we type `ls`,
we will see that there are a bunch of files with `.pdb` extensions.
By default, `ls` lists all of the files in a given
directory. The `*` character is a shortcut for "everything". Thus, if
you enter `ls *`, you will see all of the contents of a given
directory. Now try this command:

    ls *.pdb

This lists every file that ends with the string `.pdb`. This command:

    ls /usr/bin/*.sh

Lists every file in `/usr/bin` that ends in the characters `.sh` 
(depending on your system there might not be any). And
this command:

    ls meth*

lists every file in the current directory whose name begins with `meth`.
The star can go anywhere.  For example `ls p*.pdb` would list all
`.pdb` files whose name begins with "p" (as it turns out we only have
one, but there could be many).

So how does this actually work? Well...when the shell (bash) sees a
word that contains the `*` character, it automatically looks for filenames
that match the given pattern. In this case, it identified two such
files. Then, it replaced the `p*.pdb` with the list of files, separated
by spaces. In other words, the two commands:

    ls p*.pdb
    ls pentane.pdb propane.pdb

are exactly identical. The `ls` command cannot tell the difference
between these two things.

* * * *
**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that contain the letter `a`
2.  List all of the files in `/bin` that contain the letter `a` or the letter `b`
3.  List all of the files in `/bin` that contain the letter `a` AND the letter `b`

* * * *

### Tab Completion

Navigate to the home directory. Typing out directory names can waste a
lot of time. When you start typing out the name of a directory, then
hit the tab key, the shell will try to fill in the rest of the
directory name. For example, enter:

    cd s<tab>

The shell will fill in the rest of the directory name for
`shell-data`. Now enter:

    ls shell-data/molecules/p<tab><tab>

When you hit the first tab, nothing happens. The reason is that there
are multiple molecules which start with
`p`. Thus, the shell does not know which one to fill in. When you hit
tab again, the shell will list the possible choices.

Tab completion can also fill in the names of programs. For example,
enter `e<tab><tab>`. You will see the name of every program that
starts with an `e`. One of those is `echo`. If you enter `ec<tab>` you
will see that tab completion works.

## Command History

You can easily access previous commands.  Hit the up arrow.
Hit it again.  You can step backwards through your command history.
The down arrow takes your forwards in the command history.

^-C will cancel the command you are writing, and give you a fresh prompt.

^-R will do a reverse-search through your command history.  This
is very useful.

You can also review your recent commands with the `history` command.  Just enter:

    history

to see a numbered list of recent commands, including this just issues
`history` command.  You can reuse one of these commands directly by
referring to the number of that command.

If your history looked like this:

    259  ls *!
    260  ls /usr/bin/*.sh
    261  ls p*.pdb

then you could repeat command #260 by simply entering:

    !260

(that's an exclamation mark).

* * * *
**Short Exercise**

1. Find the line number in your history for the last exercise (listing
files in /bin) and reissue that command.

* * * *

## Which program?

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

Navigate to the `shell-data` directory and list the contents. You will
notice that there is a program (executable file) called `hello` in
this directory. Now, try to run the program by entering:

    hello

You should get an error saying that hello cannot be found. That is
because the directory `~/shell-data` is not in the
`PATH`. You can run the `hello` program by entering:

    ./hello

Remember that `.` is a shortcut for the current working
directory. This tells the shell to run the `hello` program which is
located right here. So, you can run any program by entering the path
to that program. You can run `hello` equally well by specifying:

    ~/shell-data/hello

Or by entering:

    ../shell-data/hello

When there are no `/` characters, the shell assumes you want to look
in one of the default places for the program.


## Examining Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Change directories to `~/shell-data/molecules`
and run:

    cat propane.pdb

This prints out the contents of the `propane.pdb` file. If you enter:

    cat propane.pdb ethane.pdb

It will print out the contents of `propane.pdb` followed by the
contents of `ethane.pdb` (you can check this by running just
`cat ethane.pdb`). `cat` just takes a list of file names and writes them
to the screen out one after another (this is where the name comes from,
`cat` is short for concatenate).

* * * *
**Short Exercises**

Use the `cat` command with a wildcard to print the contents of all the
files in `~/shell-data/modules` with a single command.

* * * *

`cat` is a terrific program, but when the file is really big, it can
be annoying to use. Try `cat ~/shell-data/dictionary.txt`--you'll just
see the end of the file.  Most of it will go shooting off screen.  You
can probably scroll your terminal window up, but only so far for a very
large file.  The program, `less`, is useful for this
case. Enter the following command:

    less ~/shell-data/dictionary.txt

`less` opens the file, and lets you navigate through it. The commands
are identical to the `man` program (in fact `man` is using `less` as
its "pager" on most systems).

**Some commands in `less`**

| key     | action |
| ------- | ---------- |
| "space" | to go forward |
|  "b"    | to go backwarsd |
|  "g"    | to go to the beginning |
|  "G"    | to go to the end |
|  "q"    | to quit |

`less` also gives you a way of searching through files. Just hit the
"/" key to begin a search. Enter the name of the word you would like
to search for and hit enter. It will jump to the next location where
that word is found. Try searching the `dictionary.txt` file for the
word "cat". If you hit "/" then "enter", `less` will just repeat
the previous search. `less` searches from the current location and
works its way forward. If you are at the end of the file and search
for the word "cat", `less` will not find it. You need to go to the
beginning of the file and search.

Remember, the `man` program actually uses `less` internally and
therefore uses the same commands, so you can search documentation
using "/" as well!

* * * *
**Short Exercise**

Use the commands we've learned so far to figure out how to search
in reverse while using `less`.

* * * *


## Redirection

Let's say we want to dump all molecules into a single file. We can do
this with a single command (from within the `molecules/` directory):

    cat *.pdb > ../all_data

This tells the shell to take the output from the `cat *.pdb` command and
dump it into a new file called `../all_data`. To verify that this
worked, examine the `all_data` file. If `all_data` had already
existed, we would overwritten it. So the `>` character tells the shell
to take the output from what ever is on the left and dump it into the
file on the right. The `>>` characters do almost the same thing,
except that they will append the output to the end of the file if
it already exists.

* * * *
**Exercise**

The `wc` command returns the number of lines, words, and characters
in a text file (it stands for "word count"). For example:

    $ wc ~/shell-data/molecules/propane.pdb
    15 111 825 shell-data/molecules/propane.pdb
    
Say you just want the number of lines in the file. The `-l` option
(that's the lower-case letter "l") does this:

    $ wc -l ~/shell-data/molecules/propane.pdb
    15 shell-data/molecules/propane.pdb

Use a wildcard with the `wc -l` command to see the lengths of all
the `.pdb` files in the `molecules/` directory.

Now run the command again and use the `>` redirection to write
the output to a file called `lengths`

* * * *


## Pipes

Pipes are one of the most important concepts in the shell. It's really
what brings command-line programs together.  As a quick motivating example,
try concatenating all the molecule files to the screen like we did before:

    cat ~/shell-data/molecules/*.pdb
    
That was a fair amount of output and probably scrolled off screen. But we
have a good program for quickly viewing longer files: `less`.  In a previous
exercise we already wrote all the molecule files to a single file.  In case
you didn't get it, here it is again:

    cat ~/shell-data/molecules/*.pdb > ~/shell-data/all_data
    
Now we can view the output with `less`:

    less ~/shell-data/all_data
    
But what if we just wanted to view the output of a command in `less` without
making all these intermediate files that we may not care to keep (it is
duplicating data after all)?  That's where power of the pipe comes in.  We
can "pipe" the output of `cat` *directly* to `less` without making an
intermediate file:

    cat ~/shell-data/molecules/*.pdb | less
    
Note the use of the "pipe" symbol `|`.  It's similar to the `>` we saw
before, but rather than outputting to file it passes the data directly
to the `less` command, which is designed to be able to read data from
the pipe.

### A bit of UNIX philosophy

We've already seen a handful of simple command-line programs (and we'll
see a few more time-permitting).  Most of these traditional UNIX commands
are designed to do one thing and one thing only. Many of them don't even
produce output (like `cp`, which we'll see later) if everything worked as
expected.  Rather than a few very complicated programs that do everything,
UNIX-like programs are *designed* to be chained together with pipes into a
sort of pipeline to accomplish more complex functionality.  To make an
analogy to mathematics, think of this sort of as composing functions, but
here we're composing programs.

<img alt="A wizard performing tricks with the shell" src="https://raw.github.com/swcarpentry/boot-camps/2013-07-cmu/shell/images/unix-magic-overacre-poster.jpg" />

* * * *

**Exercise**

It was easy to see with just a few molecule files which is the shortest
and which is the longest. But what if we have 6000 of these files and we
want a quick look?  The `sort` command can do this for us.  Sort the output
of the `wc -l *.pdb` command we saved in the `lengths` file:

    sort lengths
    
Now saw we wanted the 


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

The `mkdir` command is used to make a directory. Just enter `mkdir`
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




## Writing a Shell Script

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
