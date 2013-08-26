[Up To Schedule](../../../README.md) - Back To [Use Version Control](../local/Readme.md)

# Collaborate: Remote Version Control
----

**Based on material by Katy Huff, Anthony Scopatz, Sri Hari Krishna
Narayanan, and Matt Gidden**

## GitHub.com?

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

**NOTE** Public repos have public licenses **by default**. If you don't
want to share (in the most liberal sense) your stuff with the world, pay
girths money for private repos, or host your own.

## GitHub password 

Setting up GitHub requires a GitHub user name and password.  Please take a
moment to [create a free GitHub account](https://github.com/signup/free) (if you
want to start paying, you can add that to your account some other day).

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

Step 3 : 

    $ git remote add upstream https://github.com/UW-Madison-ACI/boot-camps.git
    $ git remote -v
    origin  https://github.com/YOU/boot-camps.git (fetch)
    origin  https://github.com/YOU/boot-camps.git (push)
    upstream        https://github.com/UW-Madison-ACI/boot-camps.git (fetch)
    upstream        https://github.com/UW-Madison-ACI/boot-camps.git (push)

All repositories that are clones begin with a remote called origin.

## git fetch : Fetching the contents of a remote

Now that you have alerted your repository to the presence of others, it
is able to pull in updates from those repositories. In this case, if you
want your master branch to track updates in the original simplestats
repository, you simply **git fetch** that repository into the master
branch of your current repository.

The fetch command alone merely pulls down information recent changes
from the original master (upstream) repository. By itself, the fetch
command does not change your local working copy. To update your local
working copy to include recent changes in the original (upstream)
repository, it is necessary to also merge.

## git merge : Merging the contents of a remote

To incorporate upstream changes from the original master repository (in
this case UW-Madison-ACI/simplestats) into your local working copy, you
must both fetch and merge. The process of merging may result in
conflicts, so pay attention. This is where version control is both at
its most powerful and its most complicated.

### Exercise : Fetch and Merge the Contents of Our GitHub Repository

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
    $ git merge upstream/master

Step 3 : Check out what happened by browsing the directory.

## git pull : Pull = Fetch + Merge

The command **git pull** is the same as executing **git fetch** followed
by **git merge**. Though it is not recommend for cases in which there
are many branches to consider, the pull command is shorter and simpler
than fetching and merging as it automates the branch matching.
Specifically, to perform the same task as we did in the previous
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
    $ git fetch alpha

and alpha would type

    $ git remote add beta https://github.com/beta/simplestats.git
    $ git remote -v
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

Step 2 : Modify the mean.py module to add the median function (and maybe a test
if you're feeling up to it!)

Step 3 : Update your remote

    $ git add mean.py
    $ git commit -m "I added a median function!"
    $ git push origin master

Step 4 : Issue a Pull Request

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as alpha/simplestats, the base as master, the head fork 
    as beta/simplestats, and the compare as median
  - write a descriptive message and send it off!

For Alpha:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments. Beta can always 
    update their version of commits during a pull request!

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository

    $ git checkout master 
    $ git fetch origin
    $ git merge origin/master

For Beta:

Step 5 : Update your local repository

    $ git checkout master 
    $ git fetch alpha
    $ git merge alpha/master

### Exercise : Swap Roles

Ok, so we've successfully issued a pull request and merged the updated code
base. Let's swap the roles of pull requester and reviewer. This time, Alpha will
add some tests to the median function.

For Alpha:

Step 1 : Start a new feature branch, named median-tests (you could do this in
single ```git checkout -b median-tests``` command)

    $ git branch median-tests
    $ git checkout median-tests

Step 2 : Modify the mean.py module to add tests for the median function

Now continue the exercise as was done previously with roles swapped.

Step 3 : Update your remote

    $ git add mean.py
    $ git commit -m "I added tests to the median function!"
    $ git push origin master

Step 4 : Issue a Pull Request

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as beta/simplestats, the base as master, the head fork 
    as alpha/simplestats, and the compare as median-tests
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
    $ git merge origin/master

For Alpha:

Step 5 : Update your local repository

    $ git checkout master 
    $ git fetch beta
    $ git merge beta/master

## git merge : Conflicts

This is the trickiest part of version control, so let's take it very carefully.

Remember that there are actually three remotes that have a relationship in this
example: upstream, alpha, and beta. To put this in more realistic terms, imagine
that the upstream branch is managed by your PI or another manager and the alpha
and beta branches are students working on a project. All of you have a copy of
mean.py, but Alpha and Beta have made changes to that file in sync with each
other. What happens if the PI (upstream) also makes changes on the same lines? A
dreaded conflict...

Now, I'll assume the roll of PI, and I want to add an additional test to the
mean function. I'll add the test and push it to the upstream repository. Sadly,
this test addition overlaps with your recent median addition. It is standard in
using version control for the person or group who is working on the *feature* to
remain up-to-date with the upstream branch. With git, this is easy to do (and is
one of its strengths vs. centralized version control systems like SVN). 

### Exercise : Experience a Conflict

Step 1 : Experience the Conflict

    $ git fetch upstream
    $ git merge upstream/master
    Auto-merging mean.py
    CONFLICT (content): Merge conflict in mean.py
    Automatic merge failed; fix conflicts and then commit the result.

## Resolving Conflicts

Now what?

Git has paused the merge. You can see this with the ```git status`` command.

    # On branch master
    # Unmerged paths:
    #   (use "git add/rm <file>..." as appropriate to mark resolution)
    #
    #       unmerged:      mean.py
    #
    no changes added to commit (use "git add" and/or "git commit -a")

If you open your mean.py file, you'll notice that git has added some strange
characters to it. Specifically, you'll see something like:

    <<<<<<< HEAD:mean.py
    ** your version of the code **
    =======
    ** upstream's version of the code **
    >>>>>>> upstream:mean.py

Now, your job is to determine how the code *should* look. For this example, that
means there should be the changes the PI added and your median function should
fit underneath it. 

### Exercise : Resolve a Conflict

Step 1 : Resolve the conflict by editing your mean.py file. It should look run
as expected and should look exactly like your version, but with the PI's changes
included.

Step 2 : Add the updated version and commit

    $ git add mean.py
    $ git commit -m "merged from upstream"
    $ git push origin master

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
In [5]: print number_frequencies[number_found]
Out[5]: 1
```

You might also ask how to get the maximum value in a python dictionary. Here's
one way.

```
In [6]: max_counts = max(number_frequencies, key = number_frequencies.get)
In [5]: print max_counts
Out[5]: 42
```

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
