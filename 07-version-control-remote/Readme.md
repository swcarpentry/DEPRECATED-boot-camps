Version Control

----

**Based on material by Katy Huff, Anthony Scopatz, and Sri Hari Krishna
Narayanan**

## github.com?

GitHub is a site where many people store their open (and closed) source
code repositories. It provides tools for browsing, collaborating on and
documenting code. Your home institution may have a repository hosting
system of it's own. To find out, ask your system administrator.  GitHub,
much like other forge hosting services (
[launchpad](https://launchpad.net), [bitbucket](https://bitbucket.org),
[googlecode](http://code.google.com), [sourceforge](http://sourceforge.net)
etc.) provides :

-   landing page support 
-   wiki support
-   network graphs and time histories of commits
-   code browser with syntax highlighting
-   issue (ticket) tracking
-   user downloads
-   varying permissions for various groups of users
-   commit triggered mailing lists
-   other service hooks (twitter, etc.)

## github login 

Setting up github requires a github user name and password.  Please take a
moment to [create a free github account](https://github.com/signup/free) (if you
want to start paying, you can add that to your account some other day).

## git remote : Steps for Forking a Repository

A key step to interacting with an online repository that you have forked
is adding the original as a remote repository. By adding the remote
repository, you inform git of a new option for fetching updates and
pushing commits.

The **git remote** command allows you to add, name, rename, list, and
delete repositories such as the original one **upstream** from your
fork, others that may be **parallel** to your fork, and so on.

### Exercise : Fork Our GitHub Repository

In step 1, you will make a copy "fork" of our test repository testrepo-2013-06-chicago on github.

In step 2, you will make a copy of **your** fork of the repository on your hard drive.

In step 3, you will let git know that in addition to your local copy and your fork on github, there is another github repostitory (called "upstream") that you might want to get updates from.

Step 1 : Go to our
[repository](https://github.com/wltrimbl/testrepo-2013-06-chicago)
from your browser, and click on the Fork button. Choose to fork it to your
username rather than any organizations.

Step 2 : Clone it. From your terminal :

    $ git clone https://github.com/YOU/testrepo-2013-06-chicago.git 
    $ cd testrepo-2013-06-chicago 
Note: YOU is a placeholder for YOUR GITHUB USERNAME.  If git asks you for a password, it really means there isn't any repository at the url that you provided.

Step 3 : 

    $ git remote add upstream https://github.com/USERNAME/testrepo-2013-06-chicago.git
    $ git remote -v
    origin  https://github.com/YOU/testrepo-2013-06-chicago.git (fetch)
    origin  https://github.com/YOU/testrepo-2013-06-chicago.git (push)
    upstream        https://github.com/wltrimbl/testrepo-2013-06-chicago.git (fetch)
    upstream        https://github.com/wltrimbl/testrepo-2013-06-chicago.git (push)

All repositories that are clones begin with a remote called origin.
### What's going on here?
The **git remote add** merely defines a nickname and a location that git will be able to communicate with for making copies of your repository.  "origin" and "upstream" are nicknames for your fork of testrepo and the "original" testrepo, respectively.

## git fetch : Fetching the contents of a remote

Now that you have alerted your repository to the presence of others, it
is able to pull in updates from those repositories. In this case, if you
want your master branch to track updates in the original wltrimbl/testrepo-2013-06-chicago 
repository, you simply **git fetch** that repository into the master
branch of your current repository.

The fetch command alone merely pulls down information recent changes
from the original master (upstream) repository. By itself, the fetch
command does not change your local working copy. To update your local
working copy to include recent changes in the original (upstream)
repository, it is necessary to also merge.

## git merge : Merging the contents of a remote

To incorporate upstream changes from the original master repository (in
this case USERNAME/testrepo-2013-06-chicago) into your local working copy, you
must both fetch and merge. The process of merging may result in
conflicts, so pay attention. This is where version control is both at
its most powerful and its most complicated.

### Exercise : Fetch and Merge the Contents of Our GitHub Repository

Step 1 : Fetch the recent remote repository history

    $ git fetch upstream

Step 2 : Make certain you are in the YYYY-MM-PLACE branch and merge the
upstream YYYY-MM-PLACE branch into your YYYY-MM-PLACE branch

    $ git checkout YYYY-MM-PLACE
    $ git merge upstream\YYYY-MM-PLACE

Step 3 : Check out what happened by browsing the directory.

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

### Exercise : Push a change to github
We'll talk about conflicts later, but first, since we have no conflicts
and are up to date, we can make a minor change and send our changes to
your fork, the "origin."

    $ git push origin master

This will update your github fork with any changes you've committed.

If you have permission to push to the upstream repository, sending
commits to that remote is exactly analagous.

    $ git push upstream master 

In the case of the upstream push, new developer accounts will not allow
this push to succeed. You're welcome to try it though.

## git merge : Conflicts

This is the trickiest part of version control, so let's take it very
carefully.

Note, your fork contains a file called Readme.md.  This is a
standard documentation file that appears rendered on the landing page
for the repository in github. To see the rendered version, visit your
fork on github, (https://github.com/YOU/testrepo-2013-06-chicago/tree/master/README.md).

In the testrepo-2013-06-chicago, you'll find a directory called githubids.  
We'd like each of you to create a file whose filename is your github id and 
whose content is a message you might want to send to the bootcamp instructors.
Use ```git add``` to add your file to the index, and ```git commit``` to 
commit it.

You'll want to start a new branch for development, make your changes there,
and then merge these changes into your main branch. It's a good convention
to think of your master branch (in this case your YYYY-MM-PLACE branch) as
the "production branch," typically by keeping that branch clean of your
local edits until they are ready for release. Developers typically use the
master branch of their local fork to track other developers changes in the
remote repository until their own local development branch changes are
ready for production.

### Exercise : Experience a Conflict

Step 1 : Make a new branch, edit the readme file in that branch, and
commit your changes.  Let's change the greeting "Vanakkam" to some other
greeting.

    $ git branch development
    $ git checkout development
    Switched to branch 'development'
    $ kate Readme.md &
    <edit the readme file and exit kate>
    $ git commit -am "Changed the welcome message to ... "

Step 2 : Mirror the remote upstream repository in your master branch (in
this case your YYYY-MM-PLACE branch) by pulling down my changes

    $ git checkout YYYY-MM-PLACE
    Switched to branch 'YYYY-MM-PLACE'
    $ git fetch upstream
    $ git merge upstream/YYYY-MM-PLACE
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

    # On branch YYYY-MM-PLACE
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
    Vanakkam
    =======
    Willkommen
    >>>>>>> development
    =====================

The intent is for you to edit the file, knowing now that I wanted the
Welcome to say Vanakkam. If you want it to say Willkommen, you should
delete the other lines. However, if you want to be inclusive, you may
want to change it to read Vanakkam and Willkommen. Decisions such as this
one must be made by a human, and why conflict resolution is not handled
more automatically by the version control system.

    Vanakkam and Willkommen

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
    To git@github.com:username/testrepo-2013-06-chicago.git

## gitolite

[Gitolite](https://github.com/sitaramc/gitolite) is a way for you to host
your own multi-user git repositories. I'm not going to go into details
here, but all you need is a machine with some drive space and network
access. You can install [minimal
ubuntu](https://help.ubuntu.com/community/Installation/MinimalCD), then
sudo apt-get install gitolite will pull in everything you need. At that
point, your collaborators will only need to send you their public ssh keys
for you to configure pull and push access to the repos.
