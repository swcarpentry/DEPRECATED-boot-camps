[Up To Schedule](../../../README.md) - Back To [Plan for Mistakes](../../../python/testing/Readme.md)

# Collaborate: Remote Version Control
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

A key step to interacting with an online repository that you have forked
is adding the original as a remote repository. By adding the remote
repository, you inform git of a new option for fetching updates and
pushing commits.

The **git remote** command allows you to add, name, rename, list, and
delete repositories such as the original one **upstream** from your
fork, others that may be **parallel** to your fork, and so on.

We'll be continuing our testing exercises using GitHub as the online repository,
so you'll need to start off by getting a copy of that repository to work on!

### Exercise : Fork Our GitHub Repository

Step 1 : Go to our
[repository](https://github.com/UW-Madison-ACI/simplestats)
from your browser, and click on the Fork button. Choose to fork it to your
user name rather than any organizations.

Step 2 : Clone it. From your terminal :

    $ git clone https://github.com/YOU/simplestats.git
    $ cd simplestats
    $ git remote -v
    origin  https://github.com/YOU/simplestats.git (fetch)
    origin  https://github.com/YOU/simplestats.git (push)

Your local repository is now connected to the remote repository using the
alias `origin`.

Step 3 : Add a connection to the common repository :

    $ git remote add upstream https://github.com/UW-Madison-ACI/simplestats.git
    $ git remote -v
    origin  https://github.com/YOU/simplestats.git (fetch)
    origin  https://github.com/YOU/simplestats.git (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (push)

All repositories that are clones begin with a remote called `origin`.  The most
common convention is clone from your own fork (`origin`) and add a remote to the
common repository as `upstream`.

## git fetch : Fetching the contents of a remote

Now that you have alerted your repository to the presence of others, it
is able to pull in updates from those repositories. In this case, if you
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
* **merge** when incorporating changes from a feature branch or collabortor

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

# Collaboration : An exercise in GitHub and Testing

The remainder of this section will outline an exercise to get your feet wet in
using some of GitHub's features. We'll be continuing our work on testing as an
example.

For the rest of this section, I'll assume that there are two collaborators,
Alpha and Beta. I'll assume that they have super-easy GitHub names, and that
their repositories are at github.com/alpha and github.com/beta.

### Exercise : Get set up

Step 1 : Group up in pairs

Step 2 : Add your collaborator as a remote and check to make sure you're
connected, e.g., Beta would type the following

    $ git remote add alpha https://github.com/alpha/simplestats.git
    $ git remote -v
    origin  https://github.com/YOU/simplestats.git (fetch)
    origin  https://github.com/YOU/simplestats.git (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (push)
    alpha           https://github.com/alpha/simplestats.git (fetch)
    alpha           https://github.com/alpha/simplestats.git (push)
    $ git fetch alpha

and Alpha would type

    $ git remote add beta https://github.com/beta/simplestats.git
    $ git remote -v
    origin  https://github.com/YOU/simplestats.git (fetch)
    origin  https://github.com/YOU/simplestats.git (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats.git (push)
    beta            https://github.com/beta/simplestats.git (fetch)
    beta            https://github.com/beta/simplestats.git (push)
    $ git fetch beta

Let's say that Beta is interested in adding a feature to the code that Beta and
Alpha are working on. Previously, we worked on a mean function, so let's add a
median function.

```python
def median(numlist):
    numlist.sort()
    length = len(numlist)
    index = length/2
    if length % 2 == 0:
       return mean([numlist[index], numlist[index - 1]])
    else:
       return numlist[index]
```

## Pull Requests : Sending Your Collaborators an Update 

From GitHub's [website](https://help.github.com/articles/using-pull-requests), a
pull request

> lets you tell others about changes you've pushed to a GitHub repository. Once
a pull request is sent, interested parties can review the set of changes,
discuss potential modifications, and even push follow-up commits if necessary.

### Exercise : Issue a Pull Request and Review it

For Beta:

Step 1 : Start a new feature branch, named median (you could do this in single
```git checkout -b median``` command)

    $ git branch median
    $ git checkout median

Step 2 : Modify the stats.py module to add the median function (and maybe a test
if you're feeling up to it!)

Step 3 : Commit your changes

    $ git add stats.py
    $ git commit -m "I added a median function!"

Step 4 : Update your remote

    $ git push origin median

Step 5 : Issue a Pull Request

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as **alpha/simplestats**, the base as **master**, the 
    head fork as **beta/simplestats**, and the compare as **median**
  - write a descriptive message and send it off!

For Alpha:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments. Beta can always 
    update their version of commits during a pull request!

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository.  At this point, all the changes exist
**only** on the remote repository.

    $ git checkout master 
    $ git fetch origin
    $ git rebase origin/master

For Beta:

Step 5 : Update your local repository

    $ git checkout master 
    $ git fetch alpha
    $ git rebase alpha/master

### Exercise : Swap Roles

Ok, so we've successfully issued a pull request and merged the updated code
base. Let's swap the roles of pull requester and reviewer. This time, Alpha will
add some tests to the median function.

For Alpha:

Step 1 : Start a new feature branch, named median-tests (you could do this in
single ```git checkout -b median-tests``` command)

    $ git branch median-tests
    $ git checkout median-tests

Step 2 : Modify the test_stats.py module to add tests for the median
function.

Now continue the exercise as was done previously with roles swapped.

Step 3 : Commit your changes

    $ git add test_stats.py
    $ git commit -m "I added tests to the median function!"

Step 4 : Update your remote

    $ git push origin median-tests

Step 5 : Issue a Pull Request

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as **beta/simplestats**, the base as **master**, the 
    head fork as **alpha/simplestats**, and the compare as **median-tests**
  - write a descriptive message and send it off!

For Beta:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments. Alpha can always 
    update their version of commits during a pull request!

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository

    $ git checkout master 
    $ git fetch origin
    $ git rebase origin/master

For Alpha:

Step 5 : Update your local repository

    $ git checkout master 
    $ git fetch beta
    $ git rebase beta/master

## git rebase/merge : Conflicts

This is the trickiest part of version control, so let's take it very carefully.

Remember that there are actually three remotes that have a relationship in this
example: upstream, alpha, and beta. To put this in more realistic terms, imagine
that the upstream branch is managed by your PI or another manager and the alpha
and beta branches are students working on a project. All of you have a copy of
stats.py, but Alpha and Beta have made changes to that file in sync with each
other. What happens if the PI (upstream) also makes changes on the same lines? A
dreaded conflict...

Now, I'll assume the roll of PI. Let's say that I know there's a series of
functions we want to add to our simplestats module. Instead of waiting around
for my grad students to finish their work, I've chosen to add some basic
function signatures, e.g., 

```python
def median(numlist):
    pass
```

I'll add this to stats.py and push it to the upstream repository. Sadly,
this addition overlaps with your recent median addition. It is standard in using
version control for the person or group who is working on the *feature* to
remain up-to-date with the upstream branch. With git, this is easy to do (and is
one of its strengths vs. centralized version control systems like SVN).

### Exercise : Experience a Conflict

Step 1 : Experience the Conflict

    $ git fetch upstream
    $ git rebase upstream/master
    First, rewinding head to replay your work on top of it...
    Applying: added a conflicting change
    Using index info to reconstruct a base tree...
    M	stats.py
    Falling back to patching base and 3-way merge...
    Auto-merging stats.py
    CONFLICT (content): Merge conflict in stats.py
    Failed to merge in the changes.
    Patch failed at 0001 added a conflicting change
    The copy of the patch that failed is found in:
       /home/YOU/simplestats/.git/rebase-apply/patch

    When you have resolved this problem, run "git rebase --continue".
    If you prefer to skip this patch, run "git rebase --skip" instead.
    To check out the original branch and stop rebasing, run "git rebase --abort".

## Resolving Conflicts

Now what?

Git has paused the rebase. You can see this with the ```git status`` command.

    # HEAD detached at c23f1e4
    # You are currently rebasing branch 'test_change2' on 'c23f1e4'.
    #   (fix conflicts and then run "git rebase --continue")
    #   (use "git rebase --skip" to skip this patch)
    #   (use "git rebase --abort" to check out the original branch)
    #
    # Unmerged paths:
    #   (use "git reset HEAD <file>..." to unstage)
    #   (use "git add <file>..." to mark resolution)
    #
    #	both modified:      stats.py
    no changes added to commit (use "git add" and/or "git commit -a")

If you open your stats.py file, you'll notice that git has added some strange
characters to it. Specifically, you'll see something like:

    <<<<<<< HEAD:stats.py
    ** your version of the code **
    =======
    ** upstream's version of the code **
    >>>>>>> upstream:stats.py

Now, your job is to determine how the code *should* look. For this example, that
means you should replace the PI's ```median``` function with yours, and keep the
PI's ```median``` placeholder below it.

### Exercise : Resolve a Conflict

Step 1 : Resolve the conflict by editing your stats.py file. It should
run as expected and should look exactly like your version, but with the
PI's changes included.

Step 2 : Add the updated version and commit

    $ git add stats.py
    $ git rebase --continue

## A GitHub Tour

Let's take a look at Issues and Milestones, both of which are great project
planning tools.

## Extra Credit

Repeat the median function exercise with a mode function. You might find the
[defaultdict](http://docs.python.org/2/library/collections.html#collections.defaultdict)
container useful -- it provides default values for key-value pairs! Here's an
example of its use.

```
In [1]: from collections import defaultdict
In [2]: number_frequencies = defaultdict(int)
In [3]: number_found = 42
In [4]: number_frequencies[number_found] += 1
In [5]: number_frequencies[number_found]
Out[5]: 1
```

You might also ask how to get the maximum value in a python dictionary. Here's
one way.

```
In [6]: max_counts = max(number_frequencies, key = number_frequencies.get)
In [7]: max_counts
Out[7]: 42
```

It works great, right? Maybe we should add a test for bimodal distributions...

# Extra Information

## Moving from Home to Work

[Here](../mobility/Readme.md)'s a tutorial about how to get yourself set up for
keeping your work at home and the office in sync.

## gitolite

[Gitolite](https://github.com/sitaramc/gitolite) is a way for you to host
your own multi-user git repositories. I'm not going to go into details
here, but all you need is a machine with some drive space and network
access. You can install [minimal
ubuntu](https://help.ubuntu.com/community/Installation/MinimalCD), then
sudo apt-get install gitolite will pull in everything you need. At that
point, your collaborators will only need to send you their public ssh keys
for you to configure pull and push access to the repos.

## Test your version control skills!

Feel up to testing all of your skills? Check out
[this](http://pcottle.github.com/learnGitBranching/) excellent website. We
haven't taught you all the things you'll need to progress through the entire
exercise, but feel free to take a look and try it out!

## A little nudge

Feel free to read
[this](http://blogs.biomedcentral.com/bmcblog/2013/02/28/version-control-for-scientific-research/)
blog post, which talks about the usefulness of version control in scientific
work. Furthermore, how many people use Google Drive to get work done, especially
in a collaborative mode? It turns out that the Drive app [includes version
control](http://support.google.com/drive/bin/answer.py?hl=en&answer=190843) as
one of its features, albeit in a limited mode (you can't really control how
often it commits, nor can you give messages to your commits). Here's an
[example](https://docs.google.com/document/d/115sweom_AnEdTvDjipsSK6Ml1HMnIKuQwz2ys52A2Bo/edit?usp=sharing).

## Why should you share your code with the world?

Take the time to do a little [background
reading](http://www.siam.org/news/news.php?id=2064&goback=.gde_112393_member_232769759).
There are arguments for and against sharing your code. Why would you,
personally, choose to do so or not?
