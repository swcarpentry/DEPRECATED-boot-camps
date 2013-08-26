[Up To Schedule](../../../README.md) - Back To [Use Version Control](../local/Readme.md) - Forward To [Collaborate](../remote/Readme.md)

# Mobility: Using Version Control at Work and Home

**Based on material by Matt Gidden**

## Overview

One of the powerful ways to use version control is to maintain your workflow
between the office, home, and wherever else you find yourself being
productive. I use this type of workflow extensively with both my research work
(i.e., coding project collaboration) and thesis writing (i.e., "personal"
activities).

The workflow in this section describes three repository locations - a server, a
work computer, and a home (laptop) computer. The server will host the "base"
repository, and the server can live anywhere you have a connection to. For
example, you could use your GitHub repository (described in the
[Collaborate](../remote/Readme.md) section as the server. If you have access to
a server on campus (e.g., I have user space on CAE's server
@best-tux.cae.wisc.edu), you can host your repository there.

For the purposes of this exercise, we will have all of the repositories
represented by different folders on your ACI user space in order to provide you
with the "flavor" of how such a workflow would work. For example, we'll use
~/work to represent your work station. You should assume that ~/work is
effectively your work station's home directory.

Let's start by making each "home directory".

    $ cd
    $ mkdir server
    $ mkdir work
    $ mkdir home

## Setting Up the "Base" Repository

If you want to use GitHub as your repository host, you can safely skip this
step. If you want to use a university server as your repository host, you'll
have to go through these steps.

We'll start off in the home directory and create a bare repository in a new
directory.

    $ cd server
    $ git init --bare myrepo.git

You'll see a new directory in ~/server named myrepo. If you ```cd``` into myrepo
and give an ```ls``` command, you'll see a bunch of weird stuff.

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

## Setting Up the "Work" Repository

Ok, so now we have an (empty) bare repository. Let's clone it on our work
computer! You'll find the repository is named properly.

    $ cd ~/work
    $ git clone ~/server/myrepo.git
    $ ls
    myrepo

Go ahead an cd into myrepo and look around at the files, you'll find that it's
empty (except for the .git folder)! Also, take a gander at that remote.

    $ cd myrepo
    $ ls -a
    . .. .git
    $ git remote -v
    origin	 ~/../server/myrepo.git (fetch)
    origin	 ~/../server/myrepo.git (push)

Sweet, it automatically added our bare repository as the origin remote and we're
set up to get some work done! Let's go ahead an add a file and commit it.

    $ touch report.tex
    $ git add report.tex
    $ git commit -m "added the base tex file for my report"
    $ git push origin master

## Setting Up the "Home" Repository

Ok, you've spent a long day at work. Maybe you still have a little more to do,
but you'd really rather go home and cook dinner first. Let's set up that home
repository so you can pick up exactly where you left off.

    $ cd ~/home
    $ git clone ~/home/myrepo.git
    $ ls
    myrepo

Let's investigate what's inside.

    $ cd myrepo
    $ ls -a 
    . .. .git report.tex

Ok, awesome! So now our work and home machines are synced against the repository
that's on the server. Any time you're doing work, as long as you **commit and
push** the work you're doing, it will be available to you anywhere you have
access to the internet. In fact, you don't even **need** access to the
internet. Once your repository is up to date, you can do all your editing and
committing without being online. Once you have a connection again, you can push
your changes.

At this point, you're fully set up to work in a best-practice, version-control
work flow. I always work in branches to make sure I don't mess anything up on my
master branch, which is another best practice. This stuff isn't intuitive when
you're first starting out, though, so just play around and get used to the
general work flow for now. You'll get better at it over time.

#### Aside: Latex and the Limits of the Version Control Workflow

Have you ever struggled with formatting Word's equations, chapters,
bibliography, etc.? [Latex](http://www.latex-project.org/) works wonders with
that. Here's a great graph taken from John Cook's
[website](http://www.johndcook.com/blog/2008/04/03/microsoft-word-and-latex/)
that explains the difference.

![wordvlatex](https://raw.github.com/gidden/boot-camps/mobility/version-control/git/mobility/wordvslatex.gif "Word vs. Latex")

I'll add that, from my experience, the types of documents needed for the left
side of the curve (where Word is easier than Latex) are simple enough to just
use Google Docs (which is version controlled as well!). 

Furthermore, I simply **can't imagine** having to write something as complicated
as a prelim or thesis using Word. You'd spend as much time formatting the thing
as you do actually writing the content. In other words, it's worth the
(smallish) headache of getting used to Latex in order to use it for bigger
documents. There's even a [Wisconsin Thesis
Template](https://github.com/willb/wi-thesis-template)! That's right, you'd have
to do **0** work to correctly format your thesis. How ridiculous is that!

Finally, and this is pure aesthetics, Latex looks **good**. I mean, **real
good**. Have you ever read a paper and thought "wow, those equations look
great"? I'd be willing to give 10-to-1 odds it was written in Latex. Plus, once
you write your first paper, you have all the infrastructure to write the next
one. I literally copy the files into a different directory and rewrite
content. Super simple!

Latex works great with the workflow described here because it's text-based. You
are literally altering text files, so there's **nothing else** going on behind
the scenes. Word files, etc., have lot's going on under the hood, and so are
poor candidates for version control. 