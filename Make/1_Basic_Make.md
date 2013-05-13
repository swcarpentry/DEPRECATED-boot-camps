#1. Make - managing files and dependencies


For this tutorial we'll be using a set of data files (*.pdb - proteine data base) and a small program called program.awk which is written in an interpreted programming language called "awk". program.awk extracts informations about: the author of the *.pdb file, the data about the atoms and counts the total number of the atoms in the file.

Let's try running program.awk and redirect the output into a *.pdb.data file.
    $ awk -f program.awk cubane.pdb > cubane.pdb.data

If we change the author in the cubane.pdb file, the cubane.pdb.data will be out of date. In other words, cubane.pdb.data depends on cubane.pdb. Every time we change cubane.pdb, we need to update cubane.pdb.data. If we have only one file, then it'd be a short process. But what if we have many files which change as we do our research? Updating any files which depend on them will be very time consuming. And here is where Make (and other automation tools) can help a lot.

Let's create our first Makefile:

    $ nano pdbprocess.mk
    # pdbprocess.mk
    cubane.pdb.data : cubane.pdb
      awk -f program.awk cubane.pdb > cubane.pdb.data
      

The first line is a comment, which should in practice be more descriptive that the name of the file.
The second and third lines are a rule that tells Make what we want to do.
* The filename on the left of the colon in the first line is the _target_ of the rule i.e. what we want to generate.
* The target’s pre-requisites - what it depends on - are listed on the right of the colon. In our case, this is just our cubane.pdb file.
* The second line of the rule is the action. This tells Make what shell command or commands to run to bring the target up to date if it is older than any of its pre-requisites (could be any number of commands).
* 
It's important to remember that the actions in rules __must__ be indexed with a __single tab__ character, not spaces, or mixes of tabs and spaces. Make was actually written by a summer intern in the mid-70s, and its idiosyncrasies sometimes show!

Let’s run our makefile:
    $ make -f pdbprocess.mk
The output shows us it has run the command we wanted.  This happened because at least one prerequisite was newer than our cubane.pdb.data target. Make uses __last modification time__ as its age. Editing and re-saving a file, for example, will change it. 

Now, if we run the command again
    $ make -f pdbprocess.mk
This time, it doesn’t execute any commands, and tells us so. This happened since the target is newer than its prerequisites.


