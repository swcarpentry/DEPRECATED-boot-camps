[Up To Schedule](../../../README.md) - Back To [Use Version Control](../local/Readme.md) - Forward To [Collaborate](../remote/Readme.md)

# Mobility: Using Version Control at Work and Home
----

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
with the "flavor" of how such a workflow would work. 

## Setting Up the "Base" Repository

If you want to use GitHub as your repository host, you can safely skip this
step. If you want to use a university server as your repository host, you'll
have to go through these steps.

We'll start off in the home directory and create a bare repository in a new
directory.

    $ cd
    $ mkdir server
    $ cd server
    $ git init --bare myrepo

You'll see a new directory in ~/server named myrepo. If you ```cd``` into myrepo
and give an ```ls``` command, you'll see a bunch of weird stuff.

    $ cd myrepo
    $ ls
    branches  config  description  HEAD  hooks  info  objects  refs

This is git's way of storing your repository's information. You shouldn't touch
this, and you can safely ignore it.

----
### Aside: Bare Repositories

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
----

    

