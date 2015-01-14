[Up To Schedule](../../../README.md) - Back To [Github and Remote Version Control](../github/Readme.md) - Forward To [Collaborate](../collaborate/Readme.md)

# Mobility: Using Version Control at Work and Home

**Based on material by Matt Gidden**

## Overview

One of the powerful ways to use version control is to maintain your workflow
between the office, home, and wherever else you find yourself being
productive. This type of workflow can be used extensively with both research
work (i.e., coding project collaboration) and thesis writing (i.e., "personal"
activities).

The workflow in this section describes three repository locations - a server, a
work computer, and a home (laptop) computer. The server will host the "base"
repository, and the server can live anywhere you have a connection to. For
example, you could use your GitHub repository (described in the [Github and
Remote Version Control](../git-and-github/Readme.md) section) as the server. If
you have access to a server on campus (e.g., server.uni.edu), you can host your
repository there (and it's private!).

For the purposes of this exercise, all of the repositories will be represented
by different folders in order to provide you with the "flavor" of how such a
workflow would work. For example, we'll use ~/work to represent your work
station. You should assume that ~/work is effectively your work station's home
directory.

## Exercise: Adding a Report to Simplestats

You've been working on your stats module for a while, and you'd like to add a
report to it. You've heard about [Latex](http://www.latex-project.org/) and that
it works really well with version control systems, so you wanted to try it
out. You're spending your time writing, so you'd like to be able to move between
work, some coffee shops, and home.

Let's start by making a laptop and work directory.

    $ cd
    $ mkdir work
    $ mkdir laptop

### Setting Up the "Work" Repository


Ok, so now we have an (empty) bare repository. 

<!-- the git init --bare command is in the aside below.  So we don't actually have a bare repo.  It looks like  "git clone" command superceded an earlier 
version that created a bare repo. -->

Let's clone it on our "work
computer"! You'll find the repository is named properly.

    $ cd ~/work
    $ git clone https://github.com/YOU/simplestats.git
    $ ls
    simplestats

Go ahead and ```cd``` into simplestats and look around at the files. Also, take a
gander at that remote.

    $ cd simplestats
    $ ls -a
    .  ..  .git  README.md  stats.py  test_stats.py
    $ git remote -v
    origin  https://github.com/YOU/simplestats.git (fetch)
    origin  https://github.com/YOU/simplestats.git (push)

Let's do some work in a branch.

    $ git branch report
    $ git checkout report

Go ahead an add a file and commit it.

    $ touch report.tex
    $ git add report.tex
    $ git commit -m "added the base tex file for my report"
    $ git push origin report

Note that if you're working on your own repository's master branch, that last
command would look like ```git push origin master```.

### Setting Up the "Laptop" Repository

Ok, you've spent a long day at work. Maybe you still have a little more to do,
but you'd really rather go home and cook dinner first. Let's set up that "laptop
repository" so you can pick up exactly where you left off.

    $ cd ~/laptop
    $ git clone https://github.com/YOU/simplestats.git
    $ ls
    simplestats

Let's investigate what's inside.

    $ cd simplestats
    $ git checkout report
    $ ls -a
    .  ..  .git  README.md  report.tex stats.py  test_stats.py

Ok, awesome! We were able to checkout the updated version of the repository.

Let's try making one set of changes. We'll add some content to the report

    $ echo "this is one fancy report" >> report.tex
    $ git add report.tex
    $ git commit -m "added some content to the report"
    $ git push origin report

So now the home machine is synced with the repository that's on the
server. Any time you're doing work, as long as you **commit and push** the work
you're doing, it will be available to you anywhere you have access to the
internet. In fact, you don't even **need** access to the internet. Once your
repository is up to date, you can do all your editing and committing without
being online. Once you have a connection again, you can push your changes.

And now let's update our work machine to also be synced against our github
repository.

    $ cd ~/work/simplestats
    $ git pull origin report
    $ tail report.tex
    this is one fancy report

Perfect! Work and laptop are synced again!

At this point, you're fully set up to work in a best-practice, version-control
work flow. Experience shows that it's best to work in branches

<!-- what branch?  master or create a new branch?  -->
 
to make sure the
master branch stays up-to-date with your server's (origin's) master
branch.

<!-- This is confusing.  Is this suggesting that you keep the master branches synced?   How does this advice relate to the previous clause in the sentence about working in branches?  -->


	
 This stuff may not be intuitive when you're first starting out, though,
so just play around and get used to the general work flow for now. You'll get
better at it over time.

#### Aside: Latex and the Limits of the Version Control Workflow

Have you ever struggled with formatting Word's equations, chapters,
bibliography, etc.? [Latex](http://www.latex-project.org/) works wonders with
that. Here's a great graph taken from John Cook's
[website](http://www.johndcook.com/blog/2008/04/03/microsoft-word-and-latex/)
that explains the difference.

<!-- Actually, John Cook copied the graphic from Marko Pinteric 
	http://www.pinteric.com/miktex.html
-->

![wordvlatex](https://raw.github.com/gidden/boot-camps/mobility/version-control/git/mobility/wordvslatex.gif "Word vs. Latex")

With the advent of Google Drive, it's often as easy to use that tool if a
document is simple enough, i.e., on the left side of the curve (where Word is
easier than Latex). Note that Google Docs is version controlled as well!

Furthermore, simply **imagine** having to write something as complicated as a
prelim or thesis using Word. You'd spend as much time formatting the thing as
you do actually writing the content. In other words, it's worth the (smallish)
headache of getting used to Latex in order to use it for bigger
documents. There's even a [Wisconsin Thesis
Template](https://github.com/willb/wi-thesis-template)! That's right, you'd have
to do **0** work to correctly format your thesis. How ridiculous is that!

Finally, and this is pure aesthetics, Latex looks **good**. Have you ever read a
paper and thought "wow, those equations look great"? It's likely written in
Latex. Plus, once you write your first paper, you have all the infrastructure to
write the next one. You can literally copy the files into a different directory
and rewrite content. Super simple!

Latex works great with the workflow described here because it's text-based. You
are literally altering text files, so there's **nothing else** going on behind
the scenes. Word files, etc., have lot's going on under the hood, and so are
poor candidates for version control. 

#### Aside: Setting Up a "Base" Repository

If you want to use GitHub as your repository host, you can safely skip this. If
you want to use a university server as your repository host, you'll have to go
through these steps.

We'll start off in the home directory and create a bare repository in a new
directory.

    $ mkdir server
    $ cd ~/server
    $ git init --bare myrepo.git

You'll see a new directory in ~/server named myrepo.git. If you ```cd``` into
myrepo.git and give an ```ls``` command, you'll see the contents of the .git
directory you saw earlier in the [Use Version Control](../local/Readme.md)
section.

    $ cd myrepo.git
    $ ls
    branches  config  description  HEAD  hooks  info  objects  refs

This is git's way of storing your repository's information. You shouldn't touch
this, and you can safely ignore it.

#### Aside: Bare Repositories

A bare repository is meant to simply **store** your files. It actually stores
the contents of the .git directory that you see in all normal repositories. It's
generally not meant to be touched by a human's hands, and is designed to
communicate through git with other non-bare repositories. 

In fact, when you initialize a new repository on GitHub, GitHub's version is a
bare repository! 

Why use a bare repository? The answer is that non-bare repositories don't always
play nice together, and it turns out it helps to have a single, base repository
that's "always right". You can get a more detailed answer
[here](http://gitolite.com/concepts/bare.html).
<!-- broken link -->

#### More than One Way to Clone

If the repository is served on an external (e.g., university)
server, you'll likely have to use git's ssh cloning
protocol. [Here](http://git-scm.com/book/en/Git-on-the-Server-The-Protocols#The-SSH-Protocol)'s
a great, short explanation of how to do that. GitHub's cloning protocol is
pretty simple, and described on the [Fork help
page](https://help.github.com/articles/fork-a-repo#step-2-clone-your-fork).

----

[Up To Schedule](../../../README.md) - Back To [Github and Remote Version Control](../github/Readme.md) - Forward To [Collaborate](../collaborate/Readme.md)
