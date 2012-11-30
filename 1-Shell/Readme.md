# The Shell

* * * * *

**Updated and presented by : Erik Bray**

**This presentation originally developed by: Milad Fatenejad**

# What is the shell how do I access the shell?

The *shell* is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

Use the GUI to open the tutorial on github.  Single click on the "Firefox
Web Browser".  Type in the URL:
    github.com/swcarpentry/2012-11-uh

A *terminal* is a program you run that gives you access to the
shell. There are many different terminal programs that vary across
operating systems.

There are many reasons to learn about the shell. In my opinion, the
most important reasons are that:

1.  It is very common to encounter the shell and
    command-line-interfaces in scientific computing, so you will
    probably have to learn it eventually

2.  The shell is a really powerful way of interacting with your
    computer. GUIs and the shell are complementary - by knowing both
    you will greatly expand the range of tasks you can accomplish with
    your computer. You will also be able to perform many tasks more
    efficiently.

The shell is just a program and there are many different shell
programs that have been developed. The most common shell (and the one
we will use) is called the Bourne-Again SHell (bash). Even if bash is
not the default shell, it usually installed on most systems and can be
started by typing `bash` in the terminal. Many commands, especially a
lot of the basic ones, work across the various shells but many things
are different. I recommend sticking with bash and learning it well.

To open a terminal, just single click on the "Terminal" icon on the
Desktop.

# The Example: Manipulating Experimental Data Files

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data from a hearing tests. To get
the data for this test, you will need internet access. Just enter the
command:

    git clone -b shell-data git://github.com/swcarpentry/2012-11-uh.git

This will grab all of the data needed for the shell exercises from the
internet (specifically from our version control repository). The details of
this will be explained tomorrow.  For now just think of it as a way to quickly
download an entire directory of files.

If you're having trouble with git or github you can download the data as a
zip file and extract it. You can enter the following URL into a web browser:

    https://github.com/swcarpentry/2012-11-uh/archive/shell-data.zip

Then just unzip it using your favorite archive tool.  Just make sure to unzip
to `/Users/(Your Name)/` on OSX, or to your home directory on Linux.


If you want to get your feet wet you can do also this on the command line using either
the `wget` or `curl` command. `curl` is usually installed on OSX, where as
Linux users are more likely to have `wget`.  First change to your home
directory:

    cd ~

Then use `wget` *or* `curl`:

    wget https://github.com/swcarpentry/2012-11-uh/archive/shell-data.zip

or

    curl -OL https://github.com/swcarpentry/2012-11-uh/archive/shell-data.zip

Then unzip, and rename the extracted directory to something a little shorter:

    unzip shell-data.zip
    mv 2012-11-uh-shell-data 2012-11-uh

# Let's get started

One very basic command is `echo`. This command is just prints text to
the terminal. Try entering the command:

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

Whenever you start up a terminal, you will start in a special
directory called the *home* directory. Every user has their own home
directory where they have full access to do whatever they want. In
this case, the `pwd` command tells us that we are in the `/home/swc`
directory. This is the home directory for the `swc` user. That is our
user name. You can always find out your user name by entering the
command `whoami`.

**File Types**

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
run like a program, then a star fill be placed of the file name.

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


**Changing Directories**

Now, let's move to a different directory. The command `cd` (change
directory) is used to move around. Let's move into the `2012-11-uh`
directory. Enter the following command:

    cd 2012-11-uh

Now use the `ls` command to see what is inside this directory. You
will see that there is an entry which is green. This means that this
is an executable. If you use `ls -F` you will see that this file ends
with a star.

If you enter the `cd` command by itself, you will return to the home
directory. Try this, and then navigate back to the `2012-11-uh`
directory.

## Arguments

Most programs take additional arguments that control their exact
behavior. For example, `-F` and `-l` are arguments to `ls`.  The `ls`
program, like many programs, take a lot of arguments. But how do we
know what the options are to particular commands?

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

    man ls

This will open the manual page for `ls`. Use the space key to go
forward and b to go backwards. When you are done reading, just hit `q`
to exit.

Programs that are run from the shell can get extremely complicated. To
see an example, open up the manual page for the `grep` program,
which we may use later this session. No one can possibly learn all of
these arguments, of course. So you will probably find yourself
referring back to the manual page frequently.

**Examining the contents of other directories**

By default, the `ls` commands lists the contents of the working
directory (i.e. the directory you are in). You can always find the
directory you are in using the `pwd` command. However, you can also
give `ls` the names of other directories to view. Navigate to the
home directory if you are not already there. Then enter the
command:

    ls 2012-11-uh

This will list the contents of the `2012-11-uh` directory without
you having to navigate there.

    ls 2012-11-uh/fits_data

This prints the contents of `fits_data/`. The `cd` command works in a
similar way. Try entering:

    cd 2012-11-uh/fits_data

and you will jump directly to `fits_data` without having to go through
the intermediate directory.

## Full vs. Relative Paths

The `cd` command takes an argument which is the directory
name. Directories can be specified using either a *relative* path a
full *path*. The directories on the computer are arranged into a
hierarchy. The full path tells you where a directory is in that
hierarchy. Navigate to the home directory. Now, enter the `pwd`
command and you should see something like

    /home/<your_username>

    or (if on OSX)

    /Users/Your Name
    
    on Windows this might be something like:
    
    C:\Documents and Settings\Your Name\
    
    but that depends on what Windows version you are using, and whether or not
    you are using a UNIX-like environment for Windows such as Cygwin or MinGW32

which is the full name of your home directory. This tells you that you are in a
directory with the same name as your login name (typically), which sits inside
a directory called `home` which sits inside the very top directory in the
hierarchy. The very top of the hierarchy is a directory called `/` which is
usually referred to as the *root directory*. So, to summarize: your home
directory is a directory in `home` or `Users` which is a directory in `/`.

Now enter the following command:

    cd /home/your_username/2012-11-uh/

or (for OSX):

    cd /Users/Your Name/2012-11-uh/

This jumps to `2012-11-uh`. Now go back to the home directory. We saw
earlier that the command:

    cd 2012-11-uh/

had the same effect - it took us to the `2012-11-uh` directory. But,
instead of specifying the full path (`/home/swc/2012-11-uh/`), we
specified a *relative path*. In other words, we specified the path relative to
our current directory. A full path always starts with a `/`. A relative path
does not. You can usually use either a full path or a relative path depending
on what is most convenient. If we are in the home directory, it is more
convenient to just enter the relative path since it involves less typing.

Now, list the contents of the /bin directory. Do you see anything
familiar in there?


## Saving time with shortcuts, wild cards, and tab completion

**Shortcuts**

There are some shortcuts which you should know about. Dealing with the home
directory is very common. So, in the shell the tilde character, `~`, is a
shortcut for your home directory. Navigate to the `2012-11-uh/fits_data`
directory, then enter the command:

    ls ~

This prints the contents of your home directory, without you having to
type the full path. The shortcut `..` always refers to the directory
above your current directory. Thus:

    ls ..

prints the contents of the `/home/swc/2012-11-uh`. You can chain
these together, so:

    ls ../../

prints the contents of `/home/swc` which is your home directory. Finally, the
special directory `.` always refers to your current directory. So, `ls`, `ls
.`, and `ls ././././.` all do the same thing, they print the contents of the
current directory. This may seem like a useless shortcut right now, but we'll
see when it is needed in a little while.

To summarize, the commands `ls ~`, `ls ~/.`, `ls ../../`, and `ls
/home/swc` all do exactly the same thing. These shortcuts are not
necessary, they are provided for your convenience.

**An example set: Erik's test FITS files**

This is a real world example from straight off Erik's filesystem.
FITS stands for Flexible Image Transport System.  It's a file format
that was developed in the 80s primarily for Astronomy data, especially
images, though it can store non-image data too, such as spectra, data
cubes, and tables.  It has found use in other scientific fields too.
By 1980s standards it was a very Flexible system, but by today it has
grown a number of warts, and requires complex specialized software for
reading and writing them.

Erik works on one such piece of software written primarily in Python
called PyFITS.  One of the goals of PyFITS is to make it possible to at
least be able to *read* every single FITS file that was ever written to
within some approximation of the actual FITS standard (there are many
FITS files out there written by software that takes the meaning of the
word "standard" very loosely).

When Erik receives a bug report on PyFITS he asks users to send a test
file that reproduces the problem.  Eventually this bug report may become
a regression test (which we'll learn about on day 2) so it's important to
keep the sample data, and it's important to keep it organized.  All of
Erik's collected sample data is in the `fits_data` directory.  But this is
a mess!  There are a mixture of FITS files and non-FITS files.  There are
some FITS files that don't use the standard `.fits` filename extension,
and there's no way of knowing what issues each file is associated with.

If we have time to get through everything this morning we'll be able to 
get this data organized. And then perhaps we will move on to some more
advanced shell topics.


**Wild cards**

Navigate to the `~/2012-11-uh/fits_data` directory. If we type `ls`,
we will see that there are a bunch of files with `.fits` extensions.
By default, `ls` lists all of the files in a given directory. The 
`*` character is a shortcut for "everything". Thus, if you enter
`ls *`, you will see all of the contents of a given directory.
Now try this command:

    ls *.fits

This lists every file that ends with a `.fits`. This command:

    ls /usr/bin/*.sh

Lists every file in `/usr/bin` that ends in the characters `.sh`. And
this command:

    ls *many*.fits

lists every file in the current directory which contains the word
`many`, and ends with `.fits`.

So how does this actually work? Well...when the shell (bash) sees a
word that contains the `*` character, it automatically looks for files
that match the given pattern. In this case, it identified four such
files. Then, it replaced the `*many*.fits` with the list of files, separated
by spaces. In other the two commands:

    ls *many*.fits
    ls manyhdus2.fits manyhdus.fits

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

**Tab Completion**

Navigate to the home directory. Typing out directory names can waste a
lot of time. When you start typing out the name of a directory, then
hit the tab key, the shell will try to fill in the rest of the
directory name. For example, enter:

    cd 2012<tab>

The shell will fill in the rest of the directory name for
`2012-11-uh`. Now from within the `2012-11-uh/fits_data` directory
enter:

    ls j<tab><tab>

When you hit the first tab, nothing happens. The reason is that there
are multiple directories in the home directory which start with
j. Thus, the shell does not know which one to fill in. When you hit
tab again, the shell will list the possible choices.

Tab completion can also fill in the names of programs. For example,
enter `e<tab><tab>`. You will see the name of every program that
starts with an `e`. One of those is `echo`. If you enter `ec<tab>` you
will see that tab completion works.

** Command History**

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

Navigate to the `2012-11-uh` directory and list the contents. You will
notice that there is a program (executable file) called `hello` in this
directory. Now, try to run the program by entering:

    hello

You should get an error saying that hello cannot be found. That is
because the directory `/home/swc/2012-11-uh/` is not in the
`PATH`. You can run the `hello` program by entering:

    ./hello

Remember that `.` is a shortcut for the current working directory. This tells
the shell to run the `hello` program which is located right here. So, you can
run any program by entering the path to that program. You can run `hello`
equally well by specifying:

    /home/swc/2012-11-uh/hello

Or by entering:

    ./hello

When there are no `/` characters, the shell assumes you want to look
in one of the default places for the program.


## The EDITOR environment variable

There are many plain text console-based text editors available on UNIX-like
systems.  Many programs (notably git, which we will use tomorrow for version
control) will sometimes drop the user into a text editor if they require some
written input from the user.

On many systems the default editor is a user-unfriendly program called vi, or
its more "modern" cousin vim.  Though vi/vim have many power users it is not
an easy to learn program for beginners.  For that reason we recommend setting
your EDITOR environment variable to something more user-friendly like `nano` by
typing:

    EDITOR=nano
    
at the command line, followed by:

    export EDITOR
    
Unforunately this change is not permanent--when you open a new terminal the default
is restored.  However there are several ways to make it permanent, which unforunately
may depend a lot on what OS you are using.  If you would like to make a permanent change
just ask for help later after class or during a break.

Finally, there is no requirement that the `$EDITOR` environment variable point to a
console-based editor.  If you have a GUI editor of choice, such as TextMate, it's possible
to use that too so long as the correct name of the command is used.


## Examining Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Enter the following command:

    cat fits_data/filenames.txt

This prints out the contents of the `filenames.txt` file. If you enter:

    cat filenames.txt filenames.txt

It will print out the contents of `filenames.txt` twice. `cat` just
takes a list of file names and writes them out one after another (this
is where the name comes from, `cat` is short for concatenate).

`cat` is a terrific program, but when the file is really big, it can
be annoying to use. The program, `less`, is useful for this
case. Enter the following command:

    less ~/2012-11-uh/dictionary.txt

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
**Short Exercise**

Use the commands we've learned so far to figure out how to search
in reverse while using `less`.

* * * *


## Redirection

Let's return to the FITS sample files. This data is located in the
`~/2012-11-uh/fits_data` directory.  Recall the `filenames.txt` file
we were looking at earlier.  Now we'll learn how we can make files like
that by saving the output of one command to a file.  Recall the `ls`
command.  Normally it outputs multiple columns, but upon checking the
man page for `ls` we see that there's a `-1` option to output one column.
Let's output all the `.fits` files in a single column:

    ls -1 *.fits

Now let's run the command again, but instead save the output to a file like
so:

    ls -1 *.fits > filenames.txt

This tells the shell to take the output from the `ls -1 *.fits` command and
dump it into a new file called `filenames.txt`. To verify that this
worked, examine the `filenames.txt` file. If `filenames.txt` file had already
existed (which it did), we would overwritten it. So the `>` character tells
the shell to take the output from what ever is on the left and dump it into the
file on the right. The `>>` characters do almost the same thing,
except that they will append the output to the file if it already
exists.


## Moving, copying, and removing

We've created a file called `filenames.txt` using the redirection operator
`>`. Let's copy the file using the `cp` command. The `cp`
command backs up the file. From within the `fits_data` directory enter:

    cp filenames.txt filenames.txt.bak

Now `filenames.txt.bak` has been created as a copy of `filenames.txt`. We can
move files around using the command `mv`. Enter this command:

    mv filenames.txt.bak /tmp/

This moves `filenames.txt.bak` into the directory `/tmp`. The directory
`/tmp` is a special directory that all users can write to. It is a
temporary place for storing files. Data stored in `/tmp` is
automatically deleted when the computer shuts down.

The `mv` command is also how you rename files. Since this file is so
important, let's rename it:

    mv filenames.txt fits_files.txt

Now the file name has been changed to fits_files.txt. Let's delete
the backup file now:

    rm /tmp/filenames.txt.bak

The `mkdir` command is used to create a directory. Just enter `mkdir`
followed by a space, then the directory name.

* * * *
**Short Exercise**

Do the following:

1.  Rename the `fits_files.txt` file back to `filenames.txt`.
2.  Create a directory in the `data` directory called `issue-001`
3.  Then, copy the `image1.fits` file into `issue-001`

* * * *

By default, `rm`, will NOT delete directories. You can tell `rm` to
delete a directory using the `-r` option. Enter the following command:

    rm -r issue-001

* * * *

## The awesome power of the Pipe

![UNIX magic](http://swcarpentry.github.com/2012-11-uh/images/learning_unix_is_not_magic.jpg)

Let's say we've already gone through and added issue numbers to each
file that they're assocated with in filenames.txt.  This has already
been done and the result is in `fits_data/issues_to_file.txt`.

What if we wanted to list just the issue numbers listed in this file
and not the associated filenames.  A program called `cut` can do that
for us.  `cut` has many options that you can read about in its man
page.  We won't go much detail about it, but in short if you use it like

    cut -c 1-5 issues_to_file.txt

it will return just the first 5 characters of each line in that file.
If we look at `issues_to_file.txt` we can see that it's already
formatted so that each issue actually takes up 9 chacters--6 for `issue-`
and 3 for the number (with zero padding).  So now let's try

    cut -c 1-9 issues_to_file.txt

Great.  This is still a long output (it was a long-ish file to begin with)
so what if we want to view its output with the `less` command?  One
way we could do that is to write this to a temporary file and then open
that with `less`:

    cut -c 1-9 issues_to_file.txt > /tmp/issues.txt
    less /tmp/issues.txt
    rm /tmp/issues.txt

But that seems like almost more hassle than its worth given that in most
modern terminal emulators we can scroll the output anyways.  But this is
where the "magic" of pipes comes in.  Give this a try:

    cut -c 1-9 issues_to_file.txt | less

As you can see, the function of the `|` command is to take the output of
one program and "pipe" it to another program.  The program to the right
of the pipe takes the output of the program on the left of the pipe as its
input.  Let's see another example.  Say we wanted to go through this list
of issues in order by issue number.  There's a program called `sort` that
will take lines of input and return those lines sorted by some criterion
(which can be controlled by command-line arguments).  Try this:

    cut -c 1-9 issues_to_file.txt | sort

This looks like it worked, but it still gave us long output like before.
Never fear: As you would hope pipes can be *chained* in a single command-
line operation:

    cut -c 1-9 issues_to_file.txt | sort | less

This can be read simply left to right: Cut out the first 9 columns of
the file, sort them, and then read the result in `less`.

One more pipe example: You may have noticed that several issue numbers are
repeated. This is because some issues required looking at multiple test files.
It might be easier if we want to look at all the issues if we could cut out
the duplicates.  This is exactly what the `uniq` command is for:

    cut -c 1-9 issues_to_file.txt | uniq | sort

A much shorter list than before.

The philosophy behind these command line programs is that none of them
really do anything all that impressive. BUT when you start chaining
them together, you can do some really powerful things really
efficiently. If you want to be proficient at using the shell, you must
learn to become proficient with the pipe and redirection operators:
`|`, `>`, `>>`.  We haven't even scratched the surface of how powerful
these tools really can be.


* * * *
**Exercise**

Now that we have pipes in our toolkit we're getting close to doing what we
originally set out to do: Organize these files sensibly.  We already have
issues_to_file.txt telling us which files belong to which issue.

Since each issue may have one or more files associated with it, and because
we may want to keep the original filenames, let's create *directories* for
each issue in which to place the associated files.

We already learned how to get a unique list of issue numbers out of the
`issues_to_file.txt` file.  Now using pipes (`|`), `cut`, `uniq`, `xargs`
and `mkdir` create a single directory under `2012-11-uh/fits_data` for each
issue.  The meaning and usage of the `xargs` command will be explained
in class (or you can read the `man` page).

* * * *

# SSH

SSH is a protocol for securely sending data across an insecure network. An
important feature of ssh is the ability to perform public key cryptography.
We are about to take a little digression, but it is important because ssh is
built in to git at a very low level.

Public key cryptography consists of a public key, a private key, and an
algorithm called a cypher. Information can be combined with the public key
using hte cypher in such a way that it appears like nonsense to anyone not
holding the private key. The only way to recover the originial information is
to use the cypher and the private key; using the cypher and the public key only
generates more nonsense-looking data. In this way, a person can use a public
key, encrypt data using the cypher + public key, and send the encrypted data
over the network without fear that someone will intercept the information.

In addition to secure communication, public key cryptography gives
authentication: you can know that a message sent by someone hasn't been
altered. A person uses their private key and the cypher on some data to create
a hash, then sends the data and the resulting hash to someone holding the
public key. The person on the other end can take the public key + hash and
verify the data wasn't changed in transit.

## SSHing to another machine

SSH stands for Secure Shell, because it was originally designed in part as a
protocol for encrypted access to a remote machine.  This was meant as a
replacement for earlier programs such as rlogin that performed the same
function.  The idea is that when you "ssh to" a remote machine, it allows you
to remotely log into that machine and start a shell session, which is then
piped back to your local shell over an encrypted network connection.  You are,
in effect, communicating with a remote shell process wrapped inside your local
shell.

To ssh to a remote machine, simply enter:

    $ ssh example.com

This will then prompt you for your login password on the remote machine.

Note that because you have to log in to the remote machine to start a shell
session, you must have an account (and generally a home directory) on that
machine.  By default, ssh will assume that your username on the remote 
machine is the same as your local username.  If that is not the case you may
specify a different username using the following syntax:

    $ ssh username@example.com

If you are doing lots of work on remote machines, it becomes inconvenient to
type in your password every time.  However, if you have a private SSH key, and
the remote machine has a copy of your public key, we can use public key
cryptography to securely (and effortlessly on our part) identify ourselves to
the remote machine.

## Generating an SSH key pair

    $ cd ~/.ssh

It will likely say "no such file or directory."

    $ ssh-keygen -t rsa -C "your_email@youremail.com"
    Generating public/private rsa key pair.
    Enter file in which to save the key (/home/swc/.ssh/id_rsa):  <press enter>

The path that it provides will be to this home directory. This is okay. **Press
enter.** You may enter a passphrase. You'll see something like this :

    Created directory '/home/swc/.ssh'.
    Enter passphrase (empty for no passphrase):
    Enter same passphrase again:
    Your identification has been saved in /home/swc/.ssh/id_rsa.
    Your public key has been saved in /home/swc/.ssh/id_rsa.pub.
    The key fingerprint is:
    09:06:c6:0f:24:b7:84:ef:22:74:de:95:f0:99:64:5d your_email@youremail.com
    The key's randomart image is:
    +--[ RSA 2048]----+
    |  .+*   . .E     |
    |  .=o+ o .       |
    |   ..oB +        |
    | . ....B .       |
    |. o.. . S        |
    |. ....           |
    | . .             |
    |                 |
    |                 |
    +-----------------+

Now if we copy the content of our public key (id_rsa.pub) into a special file
called .ssh/authorized_keys in our home directory on the remote machine, we can
log in remotely without entering a password, so long as we have access to the
matching private key locally.

Having this set up will save us considerable hassle when we start wanting to
push to remote version control repositories later.



* * * *

## Bonus:

**Permissions**

**grep**

**find**

**shell scripts**

**backtick**

**alias** -> rm -i

**variables** -> use a path example

**.bashrc**

**du**

**scp**

**Regular Expressions**
