[Up To Schedule](../../../README.md) - Back To
[Plan for Mistakes](../../../python/testing/Readme.md) - Forward To
[Mobility: Using Version Control at Work and Home](../mobility/Readme.md)

# Github and Remote Version Control
----

**Based on material by Katy Huff, Anthony Scopatz, Sri Hari Krishna
Narayanan, and Matt Gidden**

## GitHub.com?

GitHub is a site where many people store their open (and closed) source
code repositories. It provides tools for browsing, collaborating on and
documenting code. Your home institution may have a repository hosting
system of it's own. To find out, ask your system administrator.  GitHub,
much like other forge hosting services ([launchpad](https://launchpad.net), 
[bitbucket](https://bitbucket.org),
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

**NOTE** Public repos have public licenses **by default**. If you don't want
to share (in the most liberal sense) your stuff with the world, you can make a
private repo.  While that often costs money on Github, they now
have [education discounts](https://github.com/blog/1775-github-goes-to-school).

## GitHub password 

Setting up GitHub requires a GitHub user name and password.  Please take a
moment to [create a free GitHub account](https://github.com/signup/free) (if
you want an
[education discount](https://education.github.com/discount_requests/new) or to
[start paying](https://github.com/pricing), you can add that to your account some other day).

## git remote : Steps for Forking a Repository


<!--  Define forking if this is the first mention of it.
For example from https://help.github.com/articles/fork-a-repo/

A fork is a copy of a repository. Forking a repository allows you to freely experiment with changes without affecting the original project.

or 

From https://github.com/octocat/Spoon-Knife

Creating a fork is producing a personal copy of someone else's project. A fork acts as a bridge between the original repository and your personal copy. You can submit Pull Requests to help make other people's projects better by offering your changes to the original project. Forking is at the core of social coding at GitHub.

-->


A key step to interacting with an online repository that you have forked
is adding the original as a remote repository. 

<!--

This isn't clear.  The reader does not yet know what "adding the original as a remote repository" means. 

I think this needs to be explicit, like: 
1. Forking creates a copy of the repo on the servers.
2. Cloning makes a local copy of that repo on your computer.

That local copy is connected to your fork.  Now you need to connect that fork 
to the original repository.  Here are the steps...

-->

By adding the remote
repository, you inform git of a new option for fetching updates and
pushing commits.

The **git remote** command allows you to add, name, rename, list, and
delete repositories such as the original one **upstream** from your
fork, others that may be **parallel** to your fork, and so on.

We'll be continuing our testing exercises using GitHub as the online repository,
so you'll need to start off by getting a copy of that repository to work on!

### Exercise : Fork Our GitHub Repository

Step 0 : Clean Up Your Local Space

We'll be interacting with remote repositories now, so let's clean up the
simplestats folder on your machine.

    $ cd
    $ rm -r simplestats

Or if you'd like to keep it around

    $ cd
    $ mv simplestats old-simplestats

Step 1 : Go to our
[repository](https://github.com/UW-Madison-ACI/simplestats)
from your browser, and click on the Fork button. Choose to fork it to your
user name rather than any organizations.

Step 2 : Clone it. From your terminal :

    $ git clone https://github.com/YOU/simplestats
    $ cd simplestats
    $ git remote -v
    origin  https://github.com/YOU/simplestats (fetch)
    origin  https://github.com/YOU/simplestats (push)

Your local repository is now connected to the remote repository using the
alias `origin`.

Step 3 : Add a connection to the common repository :

    $ git remote add upstream https://github.com/UW-Madison-ACI/simplestats
    $ git remote -v
    origin  https://github.com/YOU/simplestats (fetch)
    origin  https://github.com/YOU/simplestats (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats (push)

<!-- This is an aside.  Maybe it should be labeled "Naming convention" and in 
smaller font.
-->
All repositories that are clones begin with a remote called `origin`.  The most
common convention is clone from your own fork (`origin`) and add a remote to the
common repository as `upstream`.

## git fetch : Fetching the contents of a remote

Now that you have connected your repository to the "upstream" original, it
is able to pull in updates from that repository. In this case, if you
want your master branch to track updates in the original simplestats
repository, you simply **git fetch** that repository into the master
branch of your current repository.

    $ git fetch upstream

The fetch command alone merely pulls down information recent changes
from the original master (`upstream`) repository. By itself, the fetch
command does not change your local working copy. To update your local
working copy to include recent changes in the original (`upstream`)
repository, it is necessary to also merge.

## git diff : Examine the differences

Now that you have fetched the `upstream` repo, you can look at the differences
between that and your local copy.  To see a summary of which files have
changed and by how much:

    $ git diff --stat upstream/master

To explore the actual changes:

    $ git diff upstream/master


## git rebase vs git merge: Insert changes that have happened on the remote

To incorporate upstream changes from the original master repository (in this
case UW-Madison-ACI/simplestats) into your local working copy, you must do
more than simply `fetch` the changes.  After fetching the changes, your local
repo know about the upstream changes, but hasn't combined them with any local
changes you may have already made.  There are two mechanisms for doing this,
with slightly different behavior.

The role of git is to keep track of little bundles of change (each commit).
In theory, it doesn't matter in what order these changes are applied, it
should end up with the same version of the files.  In practice, however, you
may want to take some control of this.  In particular, when you are combining
upstream changes into a branch where you are making local changes, it is
almost always better to **insert** all the upstream changes before your local
changes, using `rebase`.  This takes each of your commits, since the point at
which the two branches began to differ, and replays them at the end of the
upstream branch.  If there are conflicts, you will be notified and asked to
review them manually.

By contract, `merge` takes each commit from the upstream branch, since the
point at which the two branches began to differ, and replays them at the end
of your branch.  Again, if there are conflicts, you will be notified and asked
to review them manually.

There are lots of details to consider when choosing between `rebase` and
`merge`, but the simplest guidelines are:

* **rebase** when incorporating changes from an authoritative upstream
  repository
* **merge** when incorporating changes from a feature branch or collaborator

The process of rebasing/merging may result in conflicts, so pay
attention. This is where version control is both at its most powerful and its
most complicated.

### Exercise : Fetch and Rebase the Contents of Our GitHub Repository

This exercise is meant to represent the general work flow you should use to
update your fork. Let's say that you come in and sit down in the morning, you've
gotten your coffee (or tea) and you're ready to get started. However, someone
from your research group has added something to the project you're working on,
and you need to add it into your work to keep up to date. I'll add a comment,
then let's get started!

Step 1 : Fetch the recent remote repository history

    $ git fetch upstream

Step 2 : Merge the master branch

    $ git checkout master
    $ git rebase upstream/master

Step 3 : Check out what happened by browsing the directory.

## git pull : Pull = Fetch + Merge

The command **git pull** is the same as executing **git fetch** followed
by **git merge**. Though it is not recommend for cases in which there
are many branches to consider, the pull command is shorter and simpler
than fetching and merging as it automates the branch matching.
Specifically, to perform the same task as we did in the previous
exercise, the pull command would be :

    $ git pull origin
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

## Exercise: Update your remote to an upstream change

Assume that your lab group collectively works on a project (like `simplestats`),
and someone has updated the `master` branch (we can simulate that by a helper
doing an update -- helpers?).

It is now your job to: 

* get the upstream changes
* apply them to your local repository
* apply them to your fork

----

[Up To Schedule](../../../README.md) - Back To
[Plan for Mistakes](../../../python/testing/Readme.md) - Forward To
[Mobility: Using Version Control at Work and Home](../mobility/Readme.md)
