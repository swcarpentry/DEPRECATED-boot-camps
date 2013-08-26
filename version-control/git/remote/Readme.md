[Up To Schedule](../../../README.md) - Back To [Use Version Control](../local/Readme.md)

# Collaborate: Remote Version Control
----

**Based on material by Katy Huff, Anthony Scopatz, Sri Hari Krishna
Narayanan, and Matt Gidden**

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

**NOTE** Public repos have public licenses **by default**. If you don't
want to share (in the most liberal sense) your stuff with the world, pay
girths money for private repos, or host your own.

## github password 

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

We'll be continuing our testing exercises using GitHub as the online repository,
so you'll need to start off by getting a copy of that repository to work on!

### Exercise : Fork Our GitHub Repository

Step 1 : Go to our
[repository](https://github.com/UW-Madison-ACI/REPO_NAME)
from your browser, and click on the Fork button. Choose to fork it to your
username rather than any organizations.

Step 2 : Clone it. From your terminal :

    $ git clone https://github.com/YOU/REPO_NAME.git
    $ cd REPO_NAME

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
want your master branch to track updates in the original REPO_NAME
repository, you simply **git fetch** that repository into the master
branch of your current repository.

The fetch command alone merely pulls down information recent changes
from the original master (upstream) repository. By itself, the fetch
command does not change your local working copy. To update your local
working copy to include recent changes in the original (upstream)
repository, it is necessary to also merge.

## git merge : Merging the contents of a remote

To incorporate upstream changes from the original master repository (in
this case UW-Madison-ACI/REPO_NAME) into your local working copy, you
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
Alpha and Beta. I'll assume that they have super-easy github names, and that
their repositories are at github.com/alpha and github.com/beta.

### Exercise : Get set up

Step 1 : Group up in pairs

Step 2 : Add your collaborator as a remote and check to make sure you're
connected, e.g., Beta would type the following

    $ git remote add alpha https://github.com/alpha/REPO_NAME.git
    $ git remote -v
    $ git fetch alpha

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

  - Go to your remote's page (github.com/beta/REPO_NAME)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as alpha/REPO_NAME, the base as master, the head fork 
    as beta/REPO_NAME, and the compare as median
  - write a descriptive message and send it off!

For Alpha:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments.

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository

    $ git checkout master 
    $ git fetch origin
    $ git git merge origin/master

For Beta:

Step 5 : Update your local repository

    $ git checkout master 
    $ git fetch alpha
    $ git git merge alpha/master


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
