[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself/Readme.md) - Forward To [Plan for Mistakes](../../../python/testing/Readme.md)

# Make Incremental Changes: Use Version Control

**Based on materials by Katy Huff, Anthony Scopatz, Joshua R. Smith, Sri
Hari Krishna Narayanan, and Matthew Gidden**

## Motivation

See [motivation.pdf](motivation.pdf).

## Introduction

We will cover [git](http://git-scm.com) (and
[GitHub](http://github.com)) in several sessions.  There exist
numerous graphical user interfaces for git (including
[GitHub for Mac](https://mac.github.com/) and
[GitHub for Windows](https://windows.github.com/)), but we will focus
on using git from the command-line in the bash shell.

## Initial setup: `git config`

We first need to setup git with our user name and email address, and
tell it what editor to use.

     $ git config --global user.name "YOUR NAME"
     $ git config --global user.email "YOUR EMAIL"
     $ git config --global core.editor nano

These commands modify the file `~/.gitconfig`. Type
`less ~/.gitconfig` to see what it did.

## git init : Creating a Local Repository

To keep track of numerous versions of your work without saving numerous
copies, you can make a local repository for it on your computer. What git
does is to save the first version, then for each subsequent version it
saves only the changes. That is, git only records the difference between
the new version and the one before it. With this compact information,
git is able to recreate any version on demand by adding the changes to
the original in order up to the version of interest.

To create your own local (on your own machine) repository, you must
initialize the repository with the infrastructure git needs in order to
keep a record of things within the repository that you're concerned
about. The command to do this is **git init** .

### Exercise : Create a Local Repository

Step 1 : Initialize your repository.

    $ cd
    $ mkdir simplestats
    $ cd simplestats
    $ git init
    Initialized empty Git repository in /Users/swc/simplestats/.git/

Step 2 : Browse the directory's hidden files to see what happened here.
Open directories, browse file contents. Learn what you can in a minute.

    $ ls -A
    .git
    $ cd .git
    $ ls -A
    HEAD        config      description hooks       info        objects     refs

Step 3 : Use what you've learned. You may have noticed the file called
description. You can describe your repository by opening the description
file and replacing the text with a name for the repository.  We will be
creating a module with some simple statistical methods, so mine will be
called "Some simple methods for statistical analysis". You may call yours anything you like.

    $ nano description

Step 4 : Applications sometimes create files that are not needed. For
example, emacs creates a temporary file called 'filename~' when you edit
the file 'filename'. You can ask git to ignore such files by editing
the file '.git/info/exclude'. Edit the file to ignore files that end with '~'.

     git ls-files --others --exclude-from=.git/info/exclude
    # Lines that start with '#' are comments.
    # For a project mostly in C, the following would be a good set of
    # exclude patterns (uncomment them if you want to use them):
    # *.[oa]
    # *~

## git add : Adding a File To Version Control

For the git repository to know which files within this directory you
would like to keep track of, you must add them. First, you'll need to
create one, then we'll learn the **git add** command.

### Exercise : Add a File to Your Local Repository

Step 1 : Create a file to add to your repository.

    $ cd ~/simplestats/
    $ touch README.md

Step 2 : Add some text to the readme.

    $ nano README.md

Step 3 : Inform git that you would like to keep track of future changes
in this file.

    $ git add README.md

## git status : Checking the Status of Your Local Copy

The files you've created on your machine are your local "working" copy.
The changes your make in this local copy aren't stored in the repository
automatically. Until you commit them, the changes you make are local
changes. When you change anything, your set of files becomes different
from the files in the most recent official repository copy, known
as the repository HEAD. To find out what's different about them in the
terminal, try:

    $ git status
    # On branch master
    #
    # Initial commit
    #
    # Changes to be committed:

    #   (use "git rm --cached <file>..." to unstage)
    #
    #       new file:   README.md
    #

This result indicates that the current difference
between the repository HEAD (which, so far, is empty) and your
`simplestats` directory is this new README.md file.

## git commit : Saving a Snapshot

In order to save a snapshot of the current state (revision) of the
repository, we use the commit command. This command is always associated
with a message describing the changes since the last commit and
indicating their purpose. Informative commit messages will serve you
well someday, so make a habit of never committing changes without at
least a full sentence description.

**ADVICE: Commit often**

In the same way that it is wise to often save a document that you are
working on, so too is it wise to save numerous revisions of your code.
More frequent commits increase the granularity of your **undo** button.

**ADVICE: Good commit messages**

There are no hard and fast rules, but good commits are atomic: they are the smallest change that remain meaningful. A good commit message usually contains a one-line description followed by a longer explanation if necessary.

[Our repo](https://github.com/UW-Madison-ACI/boot-camps/commits/2013-08-uwmadison) has some good commit messages.

### Exercise : Commit Your Changes

Step 1 : Commit the file you've added to your repository.

    $ git commit -m "This is the first commit. It adds a readme file."
    [master (root-commit) 1863aef] This is the first commit. It adds a readme file.
     1 files changed, 1 insertions(+), 0 deletions(-)
     create mode 100644 README.md

Step 2 : Admire your work.

    $ git status
    # On branch master
    nothing to commit (working directory clean)

## git diff : Viewing the Differences

There are many diff tools.

If you have a favorite you can set your default git diff tool to execute
that one. Git, however, comes with its own diff system.

Let's recall the behavior of the diff command on the command line.
Choosing two files that are similar, the command :

    $ diff file1 file2

will output the lines that differ between the two files. This
information can be saved as what's known as a patch, but we won't go
deeply into that just now.

The only difference between the command line diff tool and git's diff
tool is that the git tool is aware of all of the revisions in your
repository, allowing each revision of each file to be treated as a full
file.

Thus, git diff will output the changes in your working directory that
are not yet staged for a commit. To see how this works, make a change in
your README.md file, but don't yet commit it.

    $ git diff

A summarized version of this output can be output with the --stat flag :

    $ git diff --stat

To see only the differences in a certain path, try:

    $ git diff HEAD -- [path]

To see what IS staged for commit (that is, what will be committed if you
type git commit without the -a flag), you can try :

    $ git diff --cached

## git log : Viewing the History

A log of the commit messages is kept by the repository and can be
reviewed with the log command.

    $ git log
    commit 1863aefd7db752f58226264e5f4282bda641ddb3
    Author: Joshua Smith <joshua.r.smith@gmail.com>
    Date:   Wed Feb 8 16:08:08 2012 -0600

        This is the first commit. It adds a readme file.

There are some useful flags for this command, such as

    -p
    -3
    --stat
    --oneline
    --graph
    --pretty=short/full/fuller/oneline
    --since=X.minutes/hours/days/weeks/months/years or YY-MM-DD-HH:MM
    --until=X.minutes/hours/days/weeks/months/years or YY-MM-DD-HH:MM
    --author=<pattern>

## git reset : Unstaging a staged file

There are a number of ways that you may accidentally stage a file that
you don't want to commit.  Create a file called `temp_notes` that
describes what you had for breakfast, and then add that file to your
repo.  Check with `status` to see that it is added but not committed.

You can now unstage that file with:

    git reset temp_notes

Check with `status`.

## git checkout : Discarding unstaged modifications (git checkout has other purposes)

Perhaps you have made a number of changes that you realize are not
going anywhere.  Add a line to `README.md` that describes your dinner
last night.  Check with `status` to see that the file is changed and
ready to be added.

You can now return to previous checked in version with:

    git checkout -- README.md

Check with `status` and take a look at the file.

## git rm : Removing files

There are a variety of reasons you way want to remove a file from the
repository after it has been committed.  Create a file called
`READYOU.md` with the first names of all your immediate family
members, and add/commit it to the repository.

You can now remove the file from the repository with:

    git rm READYOU.md

List the directory to see that you have no file named `READYOU.md`.
Use `status` to determine if you need any additional steps.

What if you delete a file in the shell without `git rm`? Try deleting
`README.md`

     rm README.md

What does `git status` say?  Oops! How can you recover this important
file?

     git checkout -- README.md

## Resources

* [git book](http://git-scm.com/book)

[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself) - Forward To [Plan for Mistakes](../../../python/testing)
