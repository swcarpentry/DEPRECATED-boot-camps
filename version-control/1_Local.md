## Tracking your changes with a local repository

Version control is centred round the notion of a *repository* which holds your directories and files. We'll start by looking at a local repository. The local repository is set up in a directory in your local filesystem (local machine).

### Create a new repository with Git

First, create a directory:

    $ mkdir papers
    $ cd papers

Now, we need to set up this directory up to be a Git repository (or "initiate the repository"):

    $ git init
    Initialized empty Git repository in /home/user/papers/.git/

The directory "papers" is now our *working directory*. 

If we look in this directory, we'll find a `.git` directory:

    $ ls -a
    .git:
    branches  config  description  HEAD  hooks  info  objects  refs

`.git` directory contains Git's configuration files. Be careful not to accidentally delete this directory! 

### Tell Git who we are

As part of the information about changes made to files Git records who made those changes. In teamwork this information is often crucial (do you want to know who rewrote your 'Conclusions' section?). So, we need to tell Git about who we are,

    $ git config --global user.name "Your Name"
    $ git config --global user.email "yourname@yourplace.org"

### Set a default editor

When working with Git we will often need to provide some short but useful information. In order to enter this information we need an editor. We'll now tell Git which editor we want to be the default one (i.e. Git will always bring it up whenever it wants us to provide some information).

You can choose any available in you system editor. For the purpose of this session we'll use *nano*:

    $ git config --global core.editor nano

To set up vi as the default editor:

    $ git config --global core.editor vi

To set up xemacs as the default editor:

    $ git config --global core.editor xemacs
    
### Git's global configuration

We can now preview (and edit, if necessary) Git's global configuration (such as our name and the default editor which we just set up). If we look in our home directory, we'll see a `.gitconfig` file,

    $ cat ~/.gitconfig
    [user]
  name = Your Name
	email = yourname@yourplace.org
    [core]
	editor = nano

This file holds global configuration that is applied to any Git repository in your file system.

### Add a file to the repository

Now, we'll create a file. Let's say we're going to write a journal paper:

    $ nano journal.txt

and add headings for Title, Author, Introduction, Conclusion and References, and save the file.

`git status` allows us to find out about the current status of files in the repository. So, we can run,

    $ git status journal.txt

Information about what Git knows about the file is displayed. For now, the important bit of information is that our file is listed as `Untracked` which means it's in our working directory but Git is not tracking it - that is, any changes made to  this file will not be recorded by Git. To tell Git about the file, we first need to add it to Git's *staging area*, which is also known as the *index* or the *cache*.

    $ git add journal.txt
    $ git status journal.txt

Now, our file is now listed as one of some `Changes to be committed`. The *staging area* can be viewed as a "loading dock", a place to hold files  we've added, or changed, until we're ready to tell Git to record those  changes in the repository. 

In order to tell Git to record our change, our new file, into the repository, we need to  *commit* it:

    $ git commit

Our default editor will now pop up. Why? Well, Git can automatically figure out that directories and files are committed, and who by (thanks to the information we provided before) and even, what changes were made, but it cannot figure out why. So we need to provide this in a *commit message*. So let's type in a message.

> **Initial structure and headings for the journal paper.**

> Ideally, commit messages should have meaning to others who may read them - or you 6 months from now. Messages like "made a change" or "added changes" or "commit 5" aren't that helpful (in fact, they're redundant!).  A good commit message usually contains a one-line description followed by a longer explanation, if necessary.

If we save our commit message, Git will now commit our file.

     [master (root-commit) c381e68] This is my journal paper.
     1 file changed, 9 insertions(+)
     create mode 100644 journal.txt

This output shows the number of files changed and the number of lines inserted or deleted across all those files. Here, we've changed (by adding) 1 file and inserted 8 lines.

Now, if we look at its status,

    $ git status journal.txt
    # On branch master
    nothing to commit (working directory clean)

`nothing to commit` means that our file is now in the repository, our working directory is up-to-date and we have no uncommitted changes in our staging area.

To see the history of changes that we made to our repository (the most recent changes will be displayed at the top):

    $ git log

The output shows: the commit identifier (also called revision number) which uniquely identifies the changes made in this commit, author, date, and your comment. *git log* command has many options to print information in various ways, for example:

    $ git log --relative-date

Now let's make some more changes to our journal.txt file . If we now run,

    $ git status journal.txt

we see a `Changes not staged for commit` section and our file is marked as `modified`. This means that a file Git knows about has been modified by us but has not yet been committed. So we can add it to the staging area and then commit the changes:

    $ git add journal.txt
    $ git commit

It can sometimes be quicker to provide our commit messages at the command-line by doing:

    $ git add journal.txt
    $ git commit -m "Added subsection headings." 

> **Top tip: When to commit changes**

> There are no hard and fast rules, but good commits are atomic - they are the smallest change that remain meaningful. 
> For code, it's useful to commit changes that can be reviewed by someone else in under an hour. 

> **What we know about software development - code reviews work** 

> Fagan (1976) discovered that a rigorous inspection can remove 60-90% of errors before the first test is run. 
> M.E., Fagan (1976). [Design and Code inspections to reduce errors in program development](http://www.mfagan.com/pdfs/ibmfagan.pdf). IBM Systems Journal 15 (3): pp. 182-211.

> **What we know about software development - code reviews should be about 60 minutes long** 

> Cohen (2006) discovered that all the value of a code review comes within the first hour, after which reviewers can become exhausted and the issues they find become ever more trivial.
> J. Cohen (2006). [Best Kept Secrets of Peer Code Review](http://smartbear.com/SmartBear/media/pdfs/best-kept-secrets-of-peer-code-review.pdf). SmartBear, 2006. ISBN-10: 1599160676. ISBN-13: 978-1599160672.

Let's add a directory, `common` and a file `references.txt` for references we may want to reuse:

    $ mkdir common
    $ nano common/references.txt

And add Fagan and Cohen to this:

    $ git add common
    $ git commit -m "Added common directory and references file with Cohen and Fagan"

and Git will add, then commit, both the directory and the file.

> **Top tip: Commit anything that cannot be automatically recreated**

> Typically we use version control to save anything that we create manually e.g. source code, scripts, notes, plain-text documents, LaTeX documents. Anything that we create using a compiler or a tool e.g. object files (`.o`, `.a`, `.class`, `.pdf`, `.dvi` etc), binaries (`exe` files), libraries (`dll` or `jar` files) we don't save as we can recreate it from the source. Adopting this approach also means there's no risk of the auto-generated files becoming out of synch with the manual ones.

In order to add all tracked files to the staging area (which may be very useful when you edited, let's say, 10 files and now you want to commit all of them):

    $ git commit -a


### Discarding changes

Let us suppose we've made a change to our file and not yet committed it. We can see the changes that we've made:

    $ git diff journal.txt

This shows the difference between the latest copy in the repository and the changes we've made. 

* `-` means a line was deleted. 
* `+` means a line was added. 
* Note that, a line that has been edited is shown as a removal of the old line and an addition of the updated line.

Maybe we made our change just to see how something looks, or, for code, to quickly try something out. But we may be unhappy with our changes. If so, we can just throw them away and return our file to the most recent version we committed to the repository by using:

    $ git checkout -- journal.txt

and we can see that our file has *reverted* to being the most up-to-date one in the repository:

    $ git status journal.txt

### Looking at our history

As we saw above, we can see the history of changes to our repository:

    $ git log

Or to a specific file in our repository:

    $ git log journal.txt

Git automatically assigns an identifier (*COMMITID*) to each commit made to the repository. In order to see the changes made between any earlier commit and our current version, we can use  `git diff`  providing the commit identifier of the earlier commit:

    $ git diff COMMITID

And, to see changes between two commits:

    $ git diff OLDER_COMMITID NEWER_COMMITID

Using our commit identifiers we can set our working directory to contain the state of the repository as it was at any commit. So, let's go back to the very first commit we made,

    $ git log
    $ git checkout COMMITID

If we look at `journal.txt` we'll see it's our very first version. And if we look at our directory,

    $ ls
    journal.txt

then we see that our `common` directory is gone. But, rest easy, while it's gone from our working directory, it's still in our repository. We can jump back to the latest commit by doing:

    4 git checkout master

And `common` will be there once more,

    $ ls
    common journal.txt

So we can get any version of our files from any point in time. In other words, we can set up our working directory back to any stage it was when we made a commit.

> **Top tip: Commit often**

> In the same way that it is wise to frequently save a document that you are working on, so too is it wise to save numerous revisions of your files. More frequent commits increase the granularity of your "undo" button.

While DropBox and GoogleDrive also preserve every version, they delete old versions after 30 days, or, for GoogleDrive, 100 revisions. DropBox allows for old versions to be stored for longer but you have to pay for this. Using revision control the only bound is how much space you have!

### Using tags as nicknames for commit identifiers

Commit identifiers are long and cryptic. Git allows us to create tags, which act as easy-to-remember nicknames for commit identifiers.

For example,

    $ git tag VER_REVIEWED_BY_JOHN

We can list tags by doing:

    $ git tag

Now if we change our file,

    $ git add journal.txt
    $ git commit -m "..." journal.txt

We can checkout our previous version using our tag instead of a commit identifier.

    $ git checkout VER_REVIEWED_BY_JOHN

And return to the latest checkout,

    $ git checkout master

> **Top tip: tag significant "events"**

> When do you tag? Well, whenever you might want to get back to the exact version you've been working on. For a paper this, might be a version that has been submitted to an internal review, or has been submitted to a conference. For code this might be when it's been submitted to review, or has been released.

### What is a branch?

You might have noticed the term `branch` in status messages,

    $ git status journal.txt
    # On branch master
    nothing to commit (working directory clean)

and when we wanted to get back to our most recent version of the repository, we used,

    $ git checkout master

Not only can our repository store the changes made to files and directories, it can store multiple sets of these, which we can use and edit and update in parallel. Each of these sets, or parallel instances, is termed a *branch* and `master` is Git's default branch. 

A new branch can be created from any commit. Branches can also be *merged* together. 

Why is this useful? 
Suppose we've developed some software and released this to our users.
We then rewrite some functionality and add some more.
Suppose a user finds a bug and we want to fix it and get this user a bug-fixed version.
Our rewrites may not be complete and we have the choice of either releasing a half-complete, rewrites-in-progress version, or telling our user to wait until the next release comes out - which may cause them to lose interest if the bug is preventing them from using our software.

Branches are a solution to this. When we release our software we can create a new branch:

    -o---o---o                               master
              \
               o                             release1.0

We can then continue developing our software in our default, or master, branch,

    -o---o---o---o---o---o                   master
              \
               o                             release

And, if a user reports a bug, we can fix it in the release branch, and immediately release the bug-fixed version,

    -o---o---o---o---o                       master
              \
               o---o---o                     release

We can then merge the changes we've made to the release back into our master branch:

    -o---o---o---o---o---M                   master
              \         /
               o---o---o                     release

And continue developing,

    -o---o---o---o---o---M---o---o           master
              \         /
               o---o---o                     release

One popular model is to have,

* A release branch, representing a released version of the code.
* A master branch, representing the most up-to-date stable version of the code.
* Various feature and/or developer-specific branches representing work-in-progress, new features etc.

For example,

             0.1      0.2        0.3
              o---------o------o------    release
             /         /      /
     o---o---o--o--o--o--o---o---o---    master
     |            /     /
     o---o---o---o     /                 fred
     |                /
     o---o---o---o---o                   kate

If a bug is found by a user, a bug fix can be applied to the release branch, and then merged with the master branch. When a feature or developer-specific branch, is stable and has been reviewed and tested it can be merged with the master branch. When the master branch has been reviewed and tested and is ready for release, a new release branch can be created from it.

If time permits, we'll look at branches but for now we'll be working solely with the `master` branch.

## The story so far...

So far, we have seen how to use Git to,

* Keep track of changes like a lab notebook for code and documents.
* Roll back changes to any point in the history of changes to our files - "undo" and "redo" for files.

But, we still have some problems. What might these be?

* If we delete our repository not only have we lost our files we've lost all our changes!
* Suppose we're away from our usual computer, for example we've taken our laptop to a conference and are far from our workstation, how do we get access to our repository then?

The answers to these questions can be found in the next section.

Previous: [Version control and Git](README.md) Next: [Working from multiple locations with a remote repository](2_Remote.md)
