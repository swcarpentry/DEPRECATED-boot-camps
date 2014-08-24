[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself/Readme.md) - Forward To [Don't Repeat Others](./broken)

# Make Incremental Changes I: Use Version Control

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

## Create a local repository: `git init`

A git repository is a directory on your computer. You will use git to
track the files and subdirectories in that directory. To begin to use
git to track a directory, you must first initialize the repository
with `git init`.

Starting from scratch, you first create a directory for your
repository, change to that directory, and type `git init`.

### ![Exercise](pics/exercise.jpg) Exercise: Create a Local Repository

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

- Change some files
- See what you've changed, with `git status`, `git diff`, and `git log`.
- Indicate what changes to save, with `git add`.
- Commit to those changes, with `git commit`.

## Check the status of your repository: `git status` 

Use `git status` to check the current status of things in your repository

```
$ git status
On branch master

Initial commit

nothing to commit (create/copy files and use "git add" to track)
```

This says that nothing much is going on.
We will talk about _branches_ tomorrow.


## Adding a File To Version Control: `git add`

You need to _specify_ which files within the directory you would like to
keep track of. To do so, you must add them to the repository with `git
add`. But first you need to create a file.

Typically the first file to add would be a ReadMe file describing the
project. This can be a plain text file or a [Markdown](http://daringfireball.net/projects/markdown/)
file (with the `.md` extension). Markdown is simple system for adding
some light mark-up, like bold and italics, to indicate sections,
and hyperlinks.

### ![Exercise](pics/exercise.jpg) Exercise: Add a File to Your Local Repository

**Step 1**: Create a file to add to your repository.

```
$ cd ~/simplestats/
$ touch README.md
```

**Step 2**: Add some text to the file.

```
$ nano README.md
```

**Step 3**: Use `git status` to check the status of the repository.

```
$ git status
```

**Step 4**: Tell git that you want to keep track of this file.

```
$ git add README.md
```

**Step 5**: Check the status again.

```
$ git status
```

## Commit your changes: `git commit`

Committing changes to your repository involves two steps: indicating
the changes to be committed with `git add` (known as "staging the changes",
which we have just done), and then actually _committing_ those
changes with `git commit`.

If you type just `git commit`, an editor will open for you to add a
comment describing the changes. Alternatively, use can use the `-m`
flag followed by the comment in quotes.

### ![Exercise](pics/exercise.jpg) Exercise: Commit Your Changes

**Step 1**: Commit the file you just added to your repository.

```
$ git commit -m "This is the first commit. It adds a readme file."
```

**Step 2**: Admire your work.

```
$ git status
```

**ADVICE: Commit often**

In the same way that it is wise to often save a document that you are
working on, so too is it wise to save numerous revisions of your code.
More frequent commits increase the granularity of your **undo** button.

**ADVICE: Good commit messages**

There are no hard and fast rules, but good commits are atomic: they
are the smallest change that remain meaningful. A good commit message
usually contains a one-line overview followed by a longer
explanation if necessary.

[The repository for this course](https://github.com/UW-Madison-ACI/boot-camps/commits/2014-08-25) has some good commit messages.

## Viewing the differences: `git diff`

`git diff` is similiar to the shell command `diff`, but rather than
comparing two files, it is used to show the historical changes in a 
repository.

If you type `git diff` alone, it will show all changes in your working
directory that have not yet been staged for a commit.

### ![Exercise](pics/exercise.jpg) Exercise: Try `git diff`

**Step 1**: Try out `git diff` without having made any changes

```
$ git diff
```

**Step 2**: Make a change to the `README.md` file.

```
$ nano README.md
```

**Step 3**: Use git diff to view the changes.

```
$ git diff
```

A summarized version of this output can be output with the `--stat` flag:

    $ git diff --stat

To see only the differences in a certain file or subdirectory, use
`git diff [path]`. For example:

    $ git diff README.md

To see the changes that **are** staged for commit, use

    $ git diff --cached

### ![Exercise](pics/exercise.jpg) Exercise: Use `git diff` with staged changes

**Step 1**: Stage the change you made to `README.md`.

```
$ git add README.md
```

**Step 2**: Try `git diff` on its own.

```
$ git diff
```

**Step 3**: Use `git diff` to see the staged changes.

```
$ git diff --cached
```

**Step 4**: Commit your change

```
$ git commit -m "Small change to README.md"
```

## Viewing the history: `git log`

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

## Unstaging a staged file: `git reset`

There are a number of ways that you may accidentally stage a file that
you don't want to commit.  Use `git reset` to unstage a file.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git reset`

**Step 1**: Make a change to the `README.md` file and stage the change.

```
$ nano README.md
$ git add README.md
```

**Step 2**: Check the status of the repository, to see that the file
  has been staged.

```
$ git status
```

**Step 3**: Unstage the file with `git reset`. `HEAD` refers to the
  most commit to the repository.

```
$ git reset HEAD README.md
```

**Step 4**: Check the status again.

```
$ git status
```

## Discarding unstaged modifications: `git checkout`

If you've made changes to a file and want to just scrap those changes
and go back to the last committed version of the file, use `git
checkout`.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git checkout`

**Step 1**: Check the status of the repository, and look at your
  unstaged changes.
  
```
$ git status
$ git diff
```

**Step 2**: Discard the changes.

```
$ git checkout README.md
```

**Step 3**: Look at the status of things again.

```
$ git status
$ git diff
```

## Removing files: `git rm`

If you want to remove a file from your repository, use `git rm`.
Note that past versions of the file will remain in the repository history.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git rm`

**Step 1**: Create a file, and add and commit it to the repository.

```
$ nano READYOU.md
$ git add READYOU.md
$ git commit -m "Add READYOU.md file"
```

**Step 2**: Remove the file from the repository.

```
$ git rm READYOU.md
```

**Step 3**: Check the status.

```
$ git status
```

**Step 4**: Commit the change (removing the file).

```
$ git commit -m "Remove READYOU.md"
```

**Step 5**: Check the status again.

```
$ git status
```

What happens if you delete a file in the shell without `git rm`? Try deleting
`README.md` with `rm` rather than `git rm`.

```
$ rm README.md
```

What does `git status` say?  Oops! How can you recover this important
file?

```
$ git checkout README.md
```

Note that, just as you should use `git rm` rather than `rm` for
removing files, you should use `git mv` rather than `mv` for moving or
renaming files.

## Resources

* [git book](http://git-scm.com/book)

[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself/Readme.md) - Forward To [Don't Repeat Others](./broken)
