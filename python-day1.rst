Introductory Python
===================

We'll be looking at the following IPython Notebooks, all of which
are under the python/ directory of the git repository:

  Scientific Python basics: `python-full.ipynb <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2012-11-scripps/master/python/python-full.ipynb>`__

  Plotting with matplotlib: `matplotlib-full.ipynb <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2012-11-scripps/master/python/matplotlib-full.ipynb>`__

  Reading and writing files: `readwrite-full.ipynb <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2012-11-scripps/master/python/readwrite-full.ipynb>`__

Running IPython Notebook
------------------------

Honestly, the hardest part is just getting things running :(.  Pick whichever
one of the solutions below works...

Running the notebook on Mac OS X using Anaconda CE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you get the Anaconda CE file (see :doc:`install` for download links!),
go run it *at the command line* -- double clicking it doesn't seem to work :(.

This should be as simple as opening up a Terminal window and typing::

   bash ~/Downloads/AnacondaCE-1.2.0-macosx.sh

and answering all the questions with the defaults ('yes' where appropriate).

Then, once it's all done installing, cd to the git directory that you
downloaded earlier, cd into the python/ subdirectory, and type ::

   ~/anaconda/bin/ipython notebook --pylab=inline

Your Web browser should pop up.  Tada!

Running the notebook on Windows using Anaconda CE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Windows, download and install Anaconda CE.  Then, open a command shell
('cmd') and change to the 'python' subdirectory of your git repository,
and type ::

   \Anaconda\Python\Scripts\ipython notebook --pylab=inline

If you're running cygwin, you might need to type ::

   /cygdrive/c/Anaconda/Python/Scripts/ipython notebook --pylab=inline

A Web browser window should open up.

Running the notebooks in the Virtual Box virtual machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start up your virtual machine (see :doc:`install` for instructions on
installing VirtualBox), and then click 'Terminal'.  Inside of terminal, run
the following commands::

   git clone https://github.com/swcarpentry/2012-11-scripps
   ls
   cd 2012-11-scripps/
   ls
   cd python
   ls
   ./run-in-vm.sh

(You can copy on your Web browser, and then paste into the Terminal in
your VM with 'ctrl-shift-V'.)

This will start up a Firefox browser pointing at IPython Notebook

IPython Notebook, a brief intro
-------------------------------

The IPython Notebook (ipynb for short) is a simple notebook interface
to Python that lets you interactively run Python code and view figures
and graphics.  You can load, save, and download notebooks as a record
of your research as well as for interaction with colleagues.

The main thing you need to know about IPython is that to execute code
in a cell, you hit Shift-ENTER.
