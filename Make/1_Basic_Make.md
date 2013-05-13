##1. Make - managing files and dependencies


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
It's important to remember that the actions in rules __must__ be indexed with a __single tab__ character, not spaces, or mixes of tabs and spaces. Make was actually written by a summer intern in the mid-70s, and its idiosyncrasies sometimes show!

Let’s run our makefile:

    $ make -f pdbprocess.mk

The output shows us it has run the command we wanted.  This happened because at least one prerequisite was newer than our cubane.pdb.data target. Make uses __last modification time__ as its age. Editing and re-saving a file, for example, will change it. 

Now, if we run the command again:

    $ make -f pdbprocess.mk
    
This time, it doesn’t execute any commands, and tells us so. This happened since the target is newer than its prerequisites.

Now, we could extend our very simple makefile to deal with another .pdb.data target, adding the following:

    $ ethane.pdb.data : ethane.pdb
        awk -f program.awk ethane.pdb > ethane.pdb.data
And to get make to realize ethane.pdb needs to be reprocessed, we can do:

    $ touch ethane.pdb
Which updates the time stamp on the file, to the current time. Then rerun our makefile:

    $ make -f pdbprocess.mk
Why hasn’t ethane.pdb.data been updated? Because make uses first rule in makefile as its default rule.
By default, it only executes this rule. If we want make to rebuild ethane.pdb.data, we need to tell it explicitly:

    $ make -f pdbprocess.mk ethane.pdb.data


But this is only marginally better than typing individual commands. To build everything at once, we introduce a phony target. It doesn’t depend on files, so it’s never up to date. But it can depend on other things. We can edit our makefile and add this in at the top, after the comment. 

    all : cubane.pdb.data ethane.pdb.data
If we type:
    $ make -f pdbprocess.mk all
Make decides that  the ‘all’ target is out of date. It depends on the other two files, so examines them. If we now touch those two files to nudge make to rebuild, and run make again

    $ touch cubane.pdb ethane.pdb
    $ make -f pdbprocess.mk all
It rebuilds them, running both commands! Note, that order in which these commands are executed is arbitary, since there’s no dependency between the two.

###Exercise 1
Add in the last rule for methane.pdb.data, based on the rules for cubane.pdb.data and ethane.pdb.data. 
* Touch all the pdb files - touch *.pdb
* Rerun the Makefile
All *.pdb.data files should rebuild. Be sure to use a single tab when indenting the action and not spaces!

####Solution:
    # pdbprocess.mk
    all : cubane.pdb.data ethane.pdb.data methane.pdb.data
    cubane.pdb.data : cubane.pdb
        awk -f program.awk cubane.pdb > cubane.pdb.data
    ethane.pdb.data : ethane.pdb
        awk -f program.awk ethane.pdb > ethane.pdb.data
    methane.pdb.data : methane.pdb
        awk -f program.awk methane.pdb > methane.pdb.data



