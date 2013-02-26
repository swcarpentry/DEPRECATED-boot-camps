Day 2 / Morning: Python scripts, and git, and github
====================================================

Etherpad: http://openetherpad.org/KTUX7IJJTx

Scripting
---------

We're going to be combining the code in two notebooks,

   `99-generate-lots-of-birds <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2013-02-uw-ctb/master/notebooks/99-generate-lots-of-birds.ipynb>`__

and

   `10-introducing-bird-counting <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2013-02-uw-ctb/master/notebooks/10-introducing-bird-counting-FULL.ipynb>`__

into scripts that we can run at the command line and place under version
control.

Notes and tips on scripts in Python
-----------------------------------

See the scripts here:

  `make-bigbirdlist.py <https://github.com/swcarpentry/2013-02-uw-ctb/blob/master/scripts/make-big-birdlist.py>`__
  `make-birdcounts.py <https://github.com/swcarpentry/2013-02-uw-ctb/blob/master/scripts/make-birdcounts.py>`__
  `plot-birdcounts.py <https://github.com/swcarpentry/2013-02-uw-ctb/blob/master/scripts/plot-birdcounts.py>`__

 - 'sys.argv' is a list of the command-line arguments.  For example, if
    you run::

         python script.py foo bar baz

    then sys.argv will contain ['script.py', 'foo', 'bar', 'baz']

 - You don't need to name Python scripts with '.py'; that's just a
   convention.

 - '#! /usr/bin/env python' is known as the she-bang line, and it tells
    the UNIX operating system to use 'python' to run the text file.

 - If you add the she-bang above, and the do 'chmod +x <filename>',
   you will be able to run the script without specifying Python.

 - 'nano' is a good, simple text editor that will get you 40% of the
   way there.  Learn emacs or vim or one of the many text editors on the
   `install page
   <http://swcarpentry.github.com/boot-camps/2013-02-25-uwash-A/>`__.
   Do not, under any circumstances, use Word.

Ways to improve these scripts:

 - put everything in a function main(args), and then put::

      if __name__ == '__main__':
         main(sys.argv)

   at the bottom.  This lets the scripts be imported and tested without
   running them.

 - write tests for them :)

 - put shebangs and chmod +x them.

 - put usage instructions in.

