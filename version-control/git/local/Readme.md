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
tell it what editor we use. The last command tells git to provided output in color.

    $ git config --global user.name "YOUR NAME"
    $ git config --global user.email "YOUR EMAIL"
    $ git config --global core.editor nano
    $ git config --global color.ui true

These commands modify the file `~/.gitconfig`. Type
`less ~/.gitconfig` to see what it did.

## Create a Local Repository: `git init`

A git repository is a directory on your computer. You will use git to
track the files and subdirectories in that directory. To begin to use
git to track a directory, you must first initialize the repository
with `git init`.

Starting from scratch, you first create a directory for your
repository, change to that directory, and type `git init`.

### <a class=exercise>Exercise : Create a Local Repository</a>

Initialize your repository: change to your home directory, create a
new directory, change to that directory, and type `git init`.

    $ cd
    $ mkdir simplestats
    $ cd simplestats
    $ git init
    Initialized empty Git repository in /Users/swc/simplestats/.git/

A `.git` subdirectory is created; git will store all of its material (configuration
information and the entire history) here.

    $ ls -A
    .git
    $ ls -A .git
    HEAD        config      description hooks       info        objects     refs



## Routine use of git

You use `git init` just once for a project, to initialize the git
repository.

Day-to-day, the basic use of git is the following:

* Change some files

* See what you've changed

  ```
  git status  
  git diff  
  git log  
  ```

* Indicate what changes to save

  ```
  git add
  ```

* Commit to those changes

  ```
  git commit
  ```

## git add : Adding a File To Version Control

For the git repository to know which files within this directory you
would like to keep track of, you must add them. First, you'll need to
create one, then we'll learn the **git add** command.

### <a class=exercise>Exercise : Add a File to Your Local Repository</a>

Step 1 : Create a file to add to your repository.

    $ cd ~/simplestats/
    $ touch README.md

Step 2 : Add some text to the Readme file.

    $ nano README.md

Step 3 : Inform git that you would like to keep track of future changes
in this file.

    $ git add README.md



## git status : Checking the Status of Your Local repository

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

[The repository for this course](https://github.com/UW-Madison-ACI/boot-camps/commits/2014-08-25) has some good commit messages.

### <a class=exercise>Exercise : Commit Your Changes</a>

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


<style>
a.exercise {
  background: #ffcfff;
}
</style>
