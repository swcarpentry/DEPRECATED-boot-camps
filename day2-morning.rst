Day 2 / Morning: Python scripts, and git, and github
====================================================

Etherpad: http://openetherpad.org/KTUX7IJJTx

Miscellaneous notes
-------------------

To upgrade packages, you can do::

   pip install -U <package>

e.g.

   pip install -U ipython

will upgrade your IPython installation.

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


Git
------------------------------------------------------

Setting up:
-------------------------------------------------------

Possible workflow:

- Tell git who you are:

   git config --global user.name "YOUR NAME"

   git config --global user.email "YOUR EMAIL"

- Initiate a local repository in the folder you are currently working in:

   git init

- Add files to be tracked under version control. If you want to track file.txt :

   git add file.txt

- Check git status to see if you have properly added the file to version control:

   git status

- Commit file (or files). Use descriptive commit messages:

   git commit -m "Changed introductory paragraph on the manuscript to reflect new publications."

- Check git status to make sure that your commit was successful:

   git status

- To see changes between 2 files:

   git diff file1 file2

- To check the git log of commits:

   git log

- To remove a file from staging area (ignores changes):

   git reset filename

- To remove a file from the git repository (not tracking file anymore):

   git rm filename

Undoing things in git:
-----------------------------------

- If you made a commit and decide that it should not have been done then use git revert. Find the commit hash using the git log command:

   git revert 4717a5c

- If you want to take a project back in time to a certain commit and discard the history of commits that came after it use reset with that commit's hash. You can also use this to remove a commit that you have not pushed to a repository:

   git reset 4717a5c


- If you want to make a new branch based on the repository at a certain commit use:

   git checkout 4717a5c


Branching:
--------------------------------------------------

- Use branches to create parallel copies of your repository that can be changed and commited in parallel to the original

- To see which branch you are currently working on:

   git branch

- To add a new branch:

   git branch March2013Experiment

- To switch to another branch:

   git checkout March2013Experiment

- To swithc back to your original branch:

   git checkout master

- If you have been working on one branch (let's say March2013Experiment) and want to merge the contents into another branch (let's say master) use merge.

- Make sure you are on the branch you want it to be merged into:

   git branch #checks to see what branch you are on
   
   git checkout master 
   
   git merge March2013Experiment



Using a remote repository:
--------------------------------------------------

- Set up your github account see here: http://github.com

- Create new repository on github. Give it a useful name. 

- Copy the address for your repository. Should be something like: https://github.com/username/Descriptive_Repository_Name_Here.git

- In your local git repository:

   git remote add origin https://github.com/username/Descriptive_Repository_Name_Here.git

- To check that your remote has been added properly:

   git remote -v

- To send your local repository commits to your remote repository:

   git push origin master

- To get commits from the remote repository (updates local repository):

   git pull origin master

- If you want to get the changes from the remote without changing your local repository use fetch:

   git branch newbranch

   git checkout newbranch

   git fetch origin

- If you want to integrate these changes into your local copy then you need to merge (pull is like doing fetch+merge):

   git checkout master

   git merge newbranch
