[Up To Schedule](../../../README.md) - Back To [Mobility: Using Version Control at Work and Home](../mobility/Readme.md)

# Collaboration : An exercise in GitHub and Testing
----

**Based on material by Katy Huff, Anthony Scopatz, Sri Hari Krishna
Narayanan, and Matt Gidden**

This section will outline an exercise to get your feet wet in using some of
GitHub's features. We'll be continuing our work on testing as an example.

For the rest of this section, I'll assume that there are two collaborators,
Alpha and Beta. I'll assume that they have super-easy GitHub names, and that
their repositories are at github.com/alpha and github.com/beta.

Let's start off by relocating back to the original simplestats repository.

    $ cd ~/simplestats

To put this in more realistic terms, imagine that the upstream repository
(UW-Madison-ACI) is managed by your PI and the alpha and beta forks are students
working on a project, tasked with implementing some stats functions. Like good
SWC followers, we'll be working in a branch, called `median`, which I will now
create. Once I have, update your local copies and remotes:

    $ git fetch upstream
    $ git checkout median
    $ git push origin median

### Exercise : Get set up

Step 1 : Group up in pairs

Step 2 : Add your collaborator as a remote and check to make sure you're
connected, e.g., Beta would type the following

    $ git remote add alpha https://github.com/alpha/simplestats
    $ git remote -v
    origin  https://github.com/YOU/simplestats (fetch)
    origin  https://github.com/YOU/simplestats (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats (push)
    alpha           https://github.com/alpha/simplestats (fetch)
    alpha           https://github.com/alpha/simplestats (push)
    $ git fetch alpha

and Alpha would type

    $ git remote add beta https://github.com/beta/simplestats
    $ git remote -v
    origin  https://github.com/YOU/simplestats (fetch)
    origin  https://github.com/YOU/simplestats (push)
    upstream        https://github.com/UW-Madison-ACI/simplestats (fetch)
    upstream        https://github.com/UW-Madison-ACI/simplestats (push)
    beta            https://github.com/beta/simplestats (fetch)
    beta            https://github.com/beta/simplestats (push)
    $ git fetch beta

## Pull Requests : Sending Your Collaborators an Update 

From GitHub's [website](https://help.github.com/articles/using-pull-requests), a
pull request

> lets you tell others about changes you've pushed to a GitHub repository. Once
a pull request is sent, interested parties can review the set of changes,
discuss potential modifications, and even push follow-up commits if necessary.

### Exercise : Issue a Pull Request and Review it

For Beta:

Step 1 : Modify the stats.py module to add the median function (shown below).

```python
def median(vals):
    vals.sort()
    z = len(vals)
    index = z / 2
    if z % 2 == 0:
       return mean([vals[index], vals[index - 1]])
    else:
       return vals[index]
```

Step 2 : Commit your changes

    $ git add stats.py
    $ git commit -m "I added a median function."

Step 3 : Update your remote

    $ git push origin median

Step 4 : Issue a Pull Request to Alpha's `median` branch

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as **alpha/simplestats**, the base branch as **median**, the 
    head fork as **beta/simplestats**, and the compare branch as **median**
  - write a descriptive message and send it off.

For Alpha:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments. Beta can always 
    update their version of commits during a pull request.

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository.  At this point, all the changes exist
**only** on the remote repository.

    $ git checkout median 
    $ git fetch origin
    $ git merge origin/median

### Exercise : Swap Roles

Ok, so we've successfully issued a pull request and merged the updated code
base. Let's swap the roles of pull requester and reviewer. This time, Alpha will
add some tests to the median function.

For Alpha:

Step 1 : Modify the test_stats.py module to add tests for the median
function.

Now continue the exercise as was done previously with roles swapped.

Step 2 : Commit your changes

    $ git add test_stats.py
    $ git commit -m "I added tests to the median function."

Step 3 : Update your remote

    $ git push origin median

Step 4 : Issue a Pull Request

  - Go to your remote's page (github.com/beta/simplestats)
  - Click Pull Requests (on the right menu) -> New Pull Request -> Edit
  - choose the base fork as **beta/simplestats**, the base as **median**, the 
    head fork as **alpha/simplestats**, and the compare as **median**
  - write a descriptive message and send it off.

For Beta:

Step 1 : Review the pull request

  - Is the code clear? Does it need comments? Is it correct? Does something 
    need clarifying? Feel free to provide in-line comments. Alpha can always 
    update their version of commits during a pull request.

Step 2 : Merge the pull request using the merge button

Step 3 : Update your local repository

    $ git checkout median
    $ git fetch origin
    $ git rebase origin/median

## git rebase/merge : Conflicts

This is the trickiest part of version control, so let's take it very carefully.

Alpha and Beta have made changes to that file in sync with each other. What
happens if the PI (upstream) also makes changes on the same lines? A dreaded
conflict... Now, I'll assume the roll of PI. Let's say that I know there's a
series of functions we want to add to our simplestats module. Instead of waiting
around for my grad students to finish their work, I've chosen to add some basic
function signatures, e.g.,

```python
def median(vals):
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

Git has paused the rebase. You can see this with the ``git status`` command.

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
In [4]: number_frequencies[number_found]
Out[4]: 0
In [5]: number_frequencies[number_found] += 1
In [6]: number_frequencies[number_found]
Out[6]: 1
```

You might also ask how to get the maximum value in a python dictionary. Here's
one way.

```
In [6]: max_counts = max(number_frequencies, key=number_frequencies.get)
In [7]: max_counts
Out[7]: 42
```

It works great, right? Maybe we should add a test for bimodal distributions...

# Extra Information

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

----

[Up To Schedule](../../../README.md) - Back To [Mobility: Using Version Control at Work and Home](../mobility/Readme.md)
