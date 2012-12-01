Version Control Collaboratively

----

**Presented By Erik Bray**
**Based on material by Katy Huff**

## What is Version Control ?

Very briefly, version control is a way to keep a backup of changing
files, to store a history of those changes, and most importantly to
allow many people in a collaboration to make changes to the same files
concurrently. There are a lot of verson control systems. Wikipedia
provides both a nice vocabulary list and a fairly complete table of some
popular version control systems and their equivalent commands.

Today, we'll be using git. Git is an example of a distributed version
control system, distinct from centralized verson control systems. I'll
make the distinction clear later, but for now, the table below will
suffice.

Version Control System Tool Options

- **Distributed**
  - Decentralized CVS (dcvs)
  - mercurial (hg)
  - git (git)
  - bazaar (bzr)
- **Centralized**
  - concurrent versions system (cvs)
  - subversion (svn)

## git --help : Getting Help

The first thing I like to know about any tool is how to get help. From
the command line type

    $ man git

The manual entry for the git version control system will appear before
you. You may scroll through it using arrows, or you can search for
keywords by typing **/** followed by the search term. I'm interested in
help, so I type **/help** and then hit enter. It looks like the syntax
for getting help with git is **git --help**.

To exit the manual page, type q.

Let's see what happens when we type :

    $ git --help

Excellent, it gives a list of commands it is able to help with, as well
as their descriptions.

    $ git --help
    usage: git [--version] [--exec-path[=<path>]] [--html-path]
               [-p|--paginate|--no-pager] [--no-replace-objects]
               [--bare] [--git-dir=<path>] [--work-tree=<path>]
               [-c name=value] [--help]
               <command> [<args>]

    The most commonly used git commands are:
       add        Add file contents to the index
       bisect     Find by binary search the change that introduced a bug
       branch     List, create, or delete branches
       checkout   Checkout a branch or paths to the working tree
       clone      Clone a repository into a new directory
       commit     Record changes to the repository
       diff       Show changes between commits, commit and working tree, etc
       fetch      Download objects and refs from another repository
       grep       Print lines matching a pattern
       init       Create an empty git repository or reinitialize an existing one
       log        Show commit logs
       merge      Join two or more development histories together
       mv         Move or rename a file, a directory, or a symlink
       pull       Fetch from and merge with another repository or a local branch
       push       Update remote refs along with associated objects
       rebase     Forward-port local commits to the updated upstream head
       reset      Reset current HEAD to the specified state
       rm         Remove files from the working tree and from the index
       show       Show various types of objects
       status     Show the working tree status
       tag        Create, list, delete or verify a tag object signed with GPG

    See 'git help <command>' for more information on a specific command.


## github.com?

GitHub is a site where many people store their open (and closed) source
code repositories. It provides tools for browsing, collaborating on and
documenting code. Your home institution may have a repository hosting system
of it's own. To find out, ask your system administrator.
GitHub, much like other forge hosting services
( [launchpad](https://launchpad.net),
[bitbucket](https://bitbucket.org), [googlecode](http://code.google.com),
[sourceforge](http://sourceforge.net) etc.)
provides :

-   landing page support
-   wiki support
-   network graphs and time histories of commits
-   code browser with syntax highlighting
-   issue (ticket) tracking
-   user downloads
-   varying permissions for various groups of users
-   commit triggered mailing lists
-   other service hooks (twitter, etc.)

**NOTE** Public repos have public licences **by default**. If you don't want to
share (in the most liberal sense) your stuff with the world, pay github money
for private repos, or host your own.

## github password

Setting up github at first requires a github user name and password.
Please take a moment to [create a free one](https://github.com/signup/free)
(if you want to start paying, you can add that to your account some
other day).

## github ssh keys

It will help you to set up automatic authorization, so that github can handshake
with your computer (in this case, your virtual machine).
There are [some setup instructions](http://help.github.com/set-up-git-redirect)
on the website, but I'll do this along with you at the front of the room as
well.

In this way people can be certain who made writes to a git repository.
Furthermore, a git server can be certain that the person who is pushing changes
to a particular repository actually has write access to that repository.

Remember, your public/private SSH keys from yesterday?  We can use those to
authenticate to github so that we don't have to enter our password every time
we push to or fetch from the remote repository.

## git config : Configuring your git environment

Once you've set up your rsa keys, you need to tell github who you are.
Crack open a terminal.

    $ git config --global user.name "Firstname Lastname"
    $ git config --global user.email "your_email@youremail.com"

Unless your name is Firstname Lastname, please don't copy the above
lines verbatim. Make the appropriate substitutions.

If you did this properly, you'll have a file in your home **(\~)**
directory that's called **.gitconfig** . It's contents should look like
:

    [user]
          name = Joshua Smith
          email = joshua.r.smith@gmail.com

This configuration step allows github to properly credit the authorship
of changes you make in your repository. For projects with numerous
authors, this is essential.

Another configuration step for some will be to set their favorite text
editor as git's text editor of choice. This is optional, since vi is
usually the default, but can be done with the following command (if you
like **nano** for example):

    $ git config --global core.editor nano


## git clone : Copying a Repository

Today, we'll check out a git type repository at
https://github.com/swcarpentry/uh-sandbox

When you clone the Original repository, the one that is created on your
local machine is a copy, and will behave as a fully fledged local
repository locally. However, with the right configuration, it will be
able to pull changes from collaborators to your local machine and push
your changes to the Original repository. We'll get to that soon, but for
now, let's **clone** the repository from GitHub.

### Exercise : Cloning a Repository from GitHub

Step 1 : Point your web browser to https://github.com/swcarpentry/uh-sandbox

Step 2 : Look for the row of buttons like "ZIP", "SSH", "Git Read-Only".
Click on "SSH" then copy the URL in text box to the right of it.

Step 3 : In your shell, enter

    $ git clone https://github.com/swcarpentry/uh-sandbox.git

This uses the same URL you just copied, so you can type "git clone " and
then paste.  You'll see something like:

    $ git clone https://github.com/swcarpentry/uh-sandbox.git
    Cloning into uh-sandbox...
    remote: Counting objects: 24, done.
    remote: Compressing objects: 100% (21/21), done.
    remote: Total 24 (delta 7), reused 17 (delta 1)
    Receiving objects: 100% (24/24), 74.36 KiB, done.
    Resolving deltas: 100% (7/7), done.

Step 4 : Let's make sure it worked. Change directories to the source
code and list the contents.  You shouldn't see anything.  There is
actually one hidden file called .gitignore but don't worry about that
right now.  It's just a stub.

    $ cd uh-sandbox

## git add : Adding a File To Version Control

For the git repository to know which files within this directory you
would like to keep track of, you must add them. First, you'll need to
create one, then we'll learn the **git add** command.

### Exercise : Add a File to Your Local Repository

Step 1 : Create a file to add to your repository. Give it your first
name and put a little one or two line bio for yourself in it.

    $ nano erik.txt

Step 2 : Inform git that you would like to keep track of future changes
in this file.

    $ git add erik.txt

## Checking The Status of Your Local Copy

The files you've created on your machine are your local "working" copy.
The changes your make in this local copy aren't backed up online
automatically. Until you commit them, the changes you make are local
changes. When you change anything, your set of files becomes different
from the files in the official repository copy. To find out what's
different about them in the terminal, try:

    $ git status
    # On branch master
    #
    # Initial commit
    #
    # Changes to be committed:

    #   (use "git rm --cached <file>..." to unstage)
    #
    #       new file:   readme.rst
    #

The null result means that you're up to date with the current version of
the repository online. This result indicates that the current difference
between the repository HEAD (which, so far, is empty).

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

There are no hard and fast rules, but good commits are atomic: they are the
smallest change that remain meaningful. A good commit message usually
contains a one-line description followed by a longer explanation if
necessary.

[This guy's repo](https://github.com/iguananaut/PyFITS/commits/master) has some good commit messages.

### Exercise : Commit Your Changes

Step 1 : Commit the file you've added to your repository.

    $ git commit -am "This is the first commit. It adds a readme file."
    [master (root-commit) 1863aef] This is the first commit. It adds a readme file.
     1 files changed, 2 insertions(+), 0 deletions(-)
     create mode 100644 readme.rst

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
your readme.rst file, but don't yet commit it.

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
    --short
    --full
    --pretty
    --since
    --until
    --author

## git pull : Pull = Fetch + Merge

The command **git pull** is the same as executing **git fetch** followed
by **git merge**. Though it is not recommened for cases in which there
are many branches to consider, the pull command is shorter and simpler
than fetching and merging as it automates the branch matching.
Specificially, to perform the same task as we did in the previous
exercise, the pull command would be :

    $ git pull upstream
    Already up-to-date.

When there have been remote changes, the pull will apply those changes
to your local branch, unless there are conflicts with your local
changes.


## git push : Sending Your Commits to Remote Repositories

The **git push** command pushes commits in a local working copy to a
remote repository. The syntax is git push [remote] [local branch].
Before pushing, a developer should always pull (or fetch + merge), so
that there is an opportunity to resolve conflicts before pushing to the
remote.

We'll talk about conflicts later, but first, since we have no conflicts
and are up to date, we can make a minor change and send our changes to
your fork, the "origin."

    $ git push origin master

If you have permission to push to the upstream repository, sending
commits to that remote is exactly analagous.

    $ git push upstream master


## git branch : Listing, Creating, and Deleting Branches

Branches are parallel instances of a repository that can be edited and
version controlled in parallel. They are useful for pursuing various
implementations experimentally or maintaining a stable core while
developing separate sections of a code base.

Without an argument, the **branch** command lists the branches that
exist in your repository.

    $ git branch
    * master

The master branch is created when the repository is initialized. With an
argument, the **branch** command creates a new branch with the given
name.

    $ git branch experimental
    $ git branch
    * master
      experimental

To delete a branch, use the **-d** flag.

    $ git branch -d experimental
    $ git branch
    * master


## git checkout : Switching Between Branches, Abandoning Local Changes

The **git checkout** command allows context switching between branches
as well as abandoning local changes.

To switch between branches, try

    $ git branch newbranch
    $ git checkout newbranch
    $ git branch

How can you tell we've switched between branches? When we used the
branch command before there was an asterisk next to the master branch.
That's because the asterisk indicates which branch you're currently in.


## git merge : Merging Branches

At some point, the experimental branch may be ready to become part of
the core or two testing branches may be ready to be combined for further
integration testing. The method for combining the changes in two
parallel branches is the **merge** command.

### Exercise : Create and Merge Branches

Step 1 : Create two new branches and list them

    $ git branch first
    $ git branch second

Step 2 : Make changes in each new branch and commit them.

    $ git checkout first
    Switched to branch 'first'
    $ touch firstnewfile
    $ git add firstnewfile
    $ git commit -am "Added firstnewfile to the first branch."
    [first 68eba44] Added firstnewfile to first branch.
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 firstnewfile
    $ git checkout second
    Switched to branch 'second'
    $ touch secondnewfile
    $ git add secondnewfile
    $ git commit -am "Added secondnewfile to the second branch."
    [second 45dd34c] Added secondnewfile to the second branch.
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 secondnewfile

Step 3 : Merge the two branches into the core

    $ git checkout first
    Switched to branch 'first'
    $ git merge second
    Merge made by recursive.
     0 files changed, 0 insertions(+), 0 deletions(-)
      create mode 100644 secondnewfile
    $ git checkout master
    Switched to branch 'master'
    $ git merge first
    Updating 1863aef..ce7e4b5
    Fast-forward
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 firstnewfile
     create mode 100644 secondnewfile


## git merge : Conflicts

This is the trickiest part of version control, so let's take it very
carefully.

In the uh-sandbox code, you'll find a file called Readme.md. This is a
standard documentation file that appears rendered on the landing page
for the repository in github. To see the rendered version, visit your
fork on github, (https://github.com/username/uh-sandbox/).

For illustration, let's imagine that, suddenly, each of the developers
on the 2012-11-uh code would like to welcome visitors in a language other
than English. Since we're all from so many different places and speak
so many languages, there will certainly be disagreements about what to
say instead of "Welcome."

I, for example, am from North Carolina, so I'll push (to the upstream repository) my own version of Welcome on line 2 of Readme.md.

You may speak another language, however, and may want to replace the
english word Welcome with an equivalent word that you prefer (willkommen,
bienvenido, benvenuti, etc.).

You'll want to start a new branch for development. It's a good
convention to think of your master branch as the "production branch,"
typically by keeping that branch clean of your local edits until they
are ready for release. Developers typically use the master branch of
their local fork to track other developers changes in the remote
repository until their own local development branch changes are ready
for production.


### Exercise : Experience a Conflict

Step 1 : Make a new branch, edit the readme file in that branch, and
commit your changes.

    $ git branch development
    $ git checkout development
    Switched to branch 'development'
    $ nano Readme.md
    <edit the readme file and exit nano>
    $ git commit -am "Changed the welcome message to ... "

Step 2 : Mirror the remote upstream repository in your master branch by
pulling down my changes

    $ git checkout master
    Switched to branch 'master'
    $ git fetch upstream
    $ git merge upstream/master
    Updating 43844ea..3b36a87
    Fast-forward
     README.rst |   2 +-
     1 files changed, 1 insertions(+), 1 deletions(-)

Step 3 : You want to push it to the internet eventually, so you pull
updates from the upstream repository, but will experience a conflict.

    $ git merge development
    Auto-merging Readme.md
    CONFLICT (content): Merge conflict in Readme.md
    Automatic merge failed; fix conflicts and then commit the result.


## git resolve : Resolving Conflicts

Now what?

Git has paused the merge. You can see this with the **git status**
command.

    # On branch master
    # Unmerged paths:
    #   (use "git add/rm <file>..." as appropriate to mark resolution)
    #
    #       unmerged:      Readme.md
    #
    no changes added to commit (use "git add" and/or "git commit -a")

The only thing that has changed is the Readme.md file. Opening it,
you'll see something like this at the beginning of the file.

    =====================
    <<<<<<< HEAD
    Howdy
    =======
    Willkommen
    >>>>>>> development
    =====================

The intent is for you to edit the file, knowing now that I wanted the
Welcome to say Howdy. If you want it to say Willkommen, you should
delete the other lines. However, if you want to be inclusive, you may
want to change it to read Howdy and Willkommen. Decisions such as this
one must be made by a human, and why conflict resolution is not handled
more automatically by the version control system.

    =====================
    Howdy and Willkommen
    =====================

This results in a status To alert git that you have made appropriate
alterations,

    $ git add Readme.md
    $ git commit
    Merge branch 'development'

    Conflicts:
      Readme.md
    #
    # It looks like you may be committing a MERGE.
    # If this is not correct, please remove the file
    # .git/MERGE_HEAD
    # and try again.
    #
    $ git push origin master
    Counting objects: 10, done.
    Delta compression using up to 2 threads.
    Compressing objects: 100% (6/6), done.
    Writing objects: 100% (6/6), 762 bytes, done.
    Total 6 (delta 2), reused 0 (delta 0)
    To git@github.com:swcarpentry/uh-sandbox.git

## git remote : Steps for Forking a Repository

A key step to interacting with an online repository that you have forked
is adding the original as a remote repository. By adding the remote
repository, you inform git of a new option for fetching updates and
pushing commits.

The **git remote** command allows you to add, name, rename, list, and
delete repositories such as the original one **upstream** from your
fork, others that may be **parallel** to your fork, and so on.

### Exercise : Fork Our GitHub Repository

While you probably already have a copy of the uh-sandbox repository, GitHub
doesn't know about it yet. You'll need to tell github you want to have an
official fork of this repository.

Step 1 : Go to our
[repository](https://github.com/swcarpentry/uh-sandbox) from your
browser, and click on the Fork button. Choose to fork it to your username
rather than any organizations.

Step 2 : Clone it. From your terminal :

    $ git clone git@github.com:username/uh-sandbox.git
    $ cd uh-sandbox

Step 3 :

    $ git remote add upstream git://github.com:swcarpentry/uh-sandbox.git
    $ git remote -v
    origin  git@github.com:username/uh-sandbox.git (fetch)
    origin  git@github.com:username/uh-sandbox.git (push)
    upstream        git://github.com/swcarpentry/uh-sandbox.git (fetch)
    upstream        git://github.com/swcarpentry/uh-sandbox.git (push)

All repositories that are clones begin with a remote called origin.

## gitolite

[Gitolite](https://github.com/sitaramc/gitolite) is a way for you to host your own multi-user git repositories. I'm not going to go into details here, but all you need is a machine with some drive space and network access. You can install [minimal ubuntu](https://help.ubuntu.com/community/Installation/MinimalCD), then sudo apt-get install gitolite will pull in everything you need. At that point, your collaborators will only need to send you their public ssh keys for you to configure pull and push access to the repos.
