Day 1 / Morning: Python
=======================

Let's start by doing two things: first, getting some files for the
morning; and second, running IPython Notebook.

You'll need to start up a shell.  On Mac OS X, this can be done
by starting up 'Terminal.app'.  On your VirtualBox, click on the
lower left, "LXTerm" icon.

Then type 'cd' just to make sure we're starting from your home directory.

To get the files, we can use 'git' to clone a repository. ::

   %% git clone https://github.com/swcarpentry/2013-02-uw-ctb.git repo

This will take the contents of that source code repository (more on that
tomorrow) and put it in the directory 'repo'.

Next, go into 'repo/notebooks'::

  %% cd repo/notebooks

and run ::

  %% ipython notebook --pylab=inline

At this point you'll be presented with a list of notebooks.  Click on
the '10 bird counting' one; the main thing you need to know is that
you can edit, etc., and the use Shift-ENTER to execute the notebook cell.

There are *two types* of notebook for the morning -- one is a completed
notebook, and the other is an empty notebook where you can type in stuff
as we go through it.  You can open both if you want -- they won't conflict
with each other.

We'll post code snippets to http://openetherpad.org/w6ffs4LN1L as we
type them in, if you want to use the empty notebook.
