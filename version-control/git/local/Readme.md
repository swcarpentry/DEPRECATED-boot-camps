[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself/Readme.md) - Forward To [Plan for Mistakes](../../../python/testing/Readme.md)

# Make Incremental Changes: Use Version Control

**Based on materials by Katy Huff, Anthony Scopatz, Joshua R. Smith, Sri 
Hari Krishna Narayanan, and Matthew Gidden**

# Motivation

From a recent [tweet](https://twitter.com/kcranstn/statuses/370914072511791104)

```
@mtholder motivating git: You mostly collaborate with yourself, and
me-from-two-months-ago never responds to email. @swcarpentry
```

## git : What is Version Control ?

Very briefly, version control is a way to keep a backup of changing
files, to store a history of those changes, and most importantly to
allow many people in a collaboration to make changes to the same files
concurrently. There are a lot of version control systems. Wikipedia
provides both a nice vocabulary list and a fairly complete table of some
popular version control systems and their equivalent commands.

Today, we'll be using git. Git is an example of a distributed version
control system, distinct from centralized versing control systems. I'll
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

## git config : Controls the behavior of git

     $ git config --global user.name "YOUR NAME"
     $ git config --global user.email "YOUR EMAIL"
     $ git config --global core.editor nano
     
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

    $ touch README.md

Step 2 : Inform git that you would like to keep track of future changes
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

    $ git commit -am "This is the first commit. It adds a readme file."
    [master (root-commit) 1863aef] This is the first commit. It adds a readme file.
     1 files changed, 2 insertions(+), 0 deletions(-)
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

## git revert : the promised "undo" button

It is possible that after many commits, you decide that you really
want to "rollback" a set of commits and start over.  It is easy to
revert your code to a previous version.

You can use `git log` and `git diff` to explore your history and
determine which version you are interested in.  Choose a version and
note the *hash* for that version. (Let's assume `abc456`)

     git revert abc456

**Importantly,** this will not erase the intervening commits.  This
will create a new commit that is changed from the previous commit by a
change that will recreate the desired version.  This retains a
complete provenance of your software, and be compared to the
prohibition in removing pages from a lab notebook.

### Exercise :

1. Create 5 files in your directory with one line of content in each
   file.
2. Commit the files to the repository.
3. Change 2 of the 5 files and commit them.
4. Undo the changes in step 3.
5. Print out the last entry in the log.
    

# Branching in Version Control
    
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

    $ git branch add_stats
    $ git checkout add_stats
    $ git branch

How can you tell we've switched between branches? When we used the
branch command before there was an asterisk next to the master branch.
That's because the asterisk indicates which branch you're currently in.

### Exercise : Copy files into your repo

Let's make sure we have a good copy of `stats.py` and `test_stats.py`.

```
$ cd ~/simplestats
$ cp ~/boot-camps/python/testing/stats.py .
$ cp ~/boot-camps/python/testing/test_stats.py .
```

Now let's add them to our repo, but in the current branch.

```
$ git add *stats.py
$ git commit -m "Adding a first version of the files for mean."
```

### Exercise : Add an additional test for std() and commit the changes.

1. Write an additional test for std().  *(Ask us for a tip if necessary)*
2. Improve std() to pass this test.
3. Commit the changed files to your repo.

## git merge : Merging Branches

At some point, the `add_stats` branch may be ready to become part of
the `master` branch.  In real life, we might do a lot more testing and
development.  For now, let's assume that our mean function is ready
and merge this back to the master.  One method for combining the
changes in two parallel branches is the **merge** command.

```
$ git checkout master
$ git merge add_stats
```

## Aside: Make your Prompt Pretty

In the next section, we'll get into the gritty details of remotes and branches
as we head toward web-based storage of your repositories. It turns out that some
folks have created a way to make this kind of navigation more convenient,
showing you what branch you're on using your bash prompt. Some super nice
properties also include color-coding when you've got changed files or when your
branch is fresh.

### Exercise : Update your prompt

Step 1 : Copy the following lines into your ~/.bashrc file (taken from a
combination of [two](http://stackoverflow.com/a/6086978)
[sources](https://gist.github.com/woods/31967)).

```
function color_my_prompt {
    local __user_and_host="\[\033[01;32m\]\u@\h"
    local __cur_location="\[\033[01;34m\]\w"
    local __git_branch='`git branch 2> /dev/null | grep -e ^* | sed -E  s/^\\\\\*\ \(.+\)$/\(\\\\\1\)\ /`'
    local __prompt_tail="\[\033[35m\]$"
    local __last_color="\[\033[00m\]"

    RED="\[\033[0;31m\]"
    YELLOW="\[\033[0;33m\]"
    GREEN="\[\033[0;32m\]"

    # Capture the output of the "git status" command.                                                                                               
    git_status="$(git status 2> /dev/null)"

    # Set color based on clean/staged/dirty.                                                                                                           
    if [[ ${git_status} =~ "working directory clean" ]]; then
        state="${GREEN}"
    elif [[ ${git_status} =~ "Changes to be committed" ]]; then
        state="${YELLOW}"
    else
        state="${RED}"
    fi

    export PS1="$__user_and_host $__cur_location ${state}$__git_branch$__prompt_tail$__last_color "
}

# Tell bash to execute this function just before displaying its prompt.                                                                              
PROMPT_COMMAND=color_my_prompt
```

Step 2 : Source your bashrc (it'll change immediately)

    $ source ~/.bashrc

Step 3 : Play around with it.

## Resources

* [git book](http://git-scm.com/book)
* [git game](http://pcottle.github.io/learnGitBranching/index.html)

[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself](../../../python/dont_repeat_yourself) - Forward To [Plan for Mistakes](../../../python/testing)
