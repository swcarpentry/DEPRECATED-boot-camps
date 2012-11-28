
Automation and Make Key Points
==============================

Steve Crouch, Mike Jackson, The Software Sustainability Institute.

This work is licensed under the Creative Commons Attribution License. Copyright (c) Software Carpentry and The University of Edinburgh and Southampton University 2010-2012. See http://software-carpentry.org/license.html for more information.

.. Written in reStructuredText, http://docutils.sourceforge.net/rst.html.

Prerequisites
-------------

Make

Introduction
------------

Cover AutomationMake.pptx, slides 1 to 9.

Use wget or curl e.g.
::

 wget --no-check-certificate https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/make-materials.zip
 curl --insecure https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/make-materials.zip > make-materials.zip

to get make-materials.zip from BitBucket and unzip these.

Edit cubane.pdb and change the author's name.
::

 ls -l *.pdb *.data

cubane.pdb.data is older than cubane.pdb - "out-of-date" - so needs to be recreated.

Create pdbprocess.mk.
::

 # PDB atom summarizer.
 cubane.pdb.data : cubane.pdb
         awk -f program.awk cubane.pdb > cubane.pdb.data

# denotes comment.

Makefile format is: 
 - Target - what we want to generate COLON pre-requisite files - what the target depends upon.
 - Actions - what commands to run to create or update the target if it is older than any of its prerequisites. Actions *must* be intended by TABs - a legacy of its 70s roots.

Run. -f specifies makefile. If this is omitted then make looks for a file called Makefile.
::

 make -f pdbprocess.mk

Make uses last modification time to determine if pre-requisites are newer than targets.
::

 make -f pdbprocess.mk

No change as the target is now up-to-date.

Add a rule.
::

 ethane.pdb.data : ethane.pdb
        awk -f program.awk ethane.pdb > ethane.pdb.data

touch updates a file's timestamp which makes it look as if it's been modified.
::

 touch ethane.pdb
 make -f pdbprocess.mk

Nothing happens to ethane.pdb as the first rule in the makefile, the default rule, was used.
::

 make -f pdbprocess.mk ethane.pdb.data

Introduce a phony target. 
::

 all : cubane.pdb.data ethane.pdb.data

all is not a "thing" but depends on "things" that can exist, and so can trigger their building.
::

 make -f pdbprocess.mk all
 touch cubane.pdb ethane.pdb
 make -f pdbprocess.mk all

Order of building pre-requisites is arbitrary.

Exercise 1 - add a rule 
-----------------------

Cover AutomationMake.pptx, slide 10.

Cover AutomationMake.pptx, slide 11.

Patterns
--------

Cover AutomationMake.pptx, slide 8-9 again.

Duplication and repeated code (and a makefile is code) creates maintainability issues.

Replace all target.
::

 PDBAnalysis.tar.gz : cubane.pdb.data ethane.pdb.data methane.pdb.data
        tar -czf PDBAnalysis.tar.gz cubane.pdb.data ethane.pdb.data methane.pdb.data

Run.
::

 make -f pdbprocess.mk PDBAnalysis.tar.gz

Rewrite action.
::

 tar -czf $@ cubane.pdb.data ethane.pdb.data methane.pdb.data

$@ means the target of the current rule. It is a Make automatic variable.

Rewrite action.
::

 tar -czf $@ $^

$^ means the dependencies of this rule. It is another make automatic variable.

Can use Bash wild-cards in file names.
::

 PDBAnalysis.tar.gz : *.pdb.data

Beware.
::

 rm *.pdb.data
 make -f pdbprocess.mk PDBAnalysis.tar.gz

Fails as there is no \*.pdb.data file so \*.pdb.data is used as-is.

Need to create .data files in a more manual way.
::

 make -f pdbprocess.mk cubane.pdb.data
 make -f pdbprocess.mk methane.pdb.data
 make -f pdbprocess.mk ethane.pdb.data

Output data is not just dependent upon input data but also programs that create it. 

.pdb.data files are dependent upon program.awk. 
::

 cubane.pdb.data : program.awk
 ethane.pdb.data : program.awk
 methane.pdb.data : program.awk

No dependencies for .pdb files, as these are input files - a dependency on program.awk would be a false dependency.
::

 touch program.awk
 make -f pdbprocess.mk

Pattern rules
-------------

Still duplication and repetition.

Replace .pdb.data targets and dependencies with a single target and dependency.
::

 %.pdb.data : %.pdb

% is a Make wild-card and this rule is termed a pattern rule.

Exercise 2 - simplify a rule 
----------------------------

Cover AutomationMake.pptx, slide 12.

Question: what do you notice?

Answer: program.awk is included in the processing because $^ matches all dependencies.

More on patterns
----------------

Change pattern rule action to:
::

        awk -f program.awk $< > $@

$< means use the first dependency only.

Replace the program.awk dependent rules with:
::

 %.pdb.data : program.awk

Run.
::

 touch program.awk
 make -f pdbprocess.mk

Does not rebuild. % wild-cards only matches the first rule.

Replace with a false dependency:
::

 %.pdb : program.awk
         touch $@

When program.awk is updated, timestamps of the raw data files are updated, which retriggers the build. Make is not perfect.

Exercise 3 - zip the files
--------------------------

Cover AutomationMake.pptx, slide 13.

Question: why do the .zip files need to be manually created first?

Answer: Because if \*.pdb.data.zip doesn't match any existing files, it's left as it is, it percolates down and the result is a single "\*.pdb.data.zip" file.

Macros
------

Add the program to the archive.
::

 PDBAnalysis.tar.gz : *.pdb.data program.awk

What is there are different programs? What if there are many references to this program?

Use a variable, called a Make macro.
::

 PROCESSOR=program.awk

Replace program.awk in each rule with $(PROCESSOR)

awk is a program too. Another user might only have gawk.
::

 AWKPROG=awk 

Replace awk in each rule with $(AWKPROG)

Keep macros at the top of a Makefile so they are easy to find. Or put in another file.

Move the macros to config.mk.
::

 # PDB atom summarizer configuration.
 PROCESSOR=program.awk
 AWKPROG=awk

Read the macros into the Makefile.
::

 include config.mk

A separate configuration allows for one Makefile with many configurations, no need to edit the Makefile (reduces risk of introducing a bug), separates program (Makefile) from its data.

Exercise 4 - create a macro
---------------------------

Cover AutomationMake.pptx, slide 14.

Completed makefile:
::

 # PDB atom summarizer.
 include config.mk

 $(TARFILE).tar.gz : *.pdb.data.zip $(PROCESSOR)
    tar -czf $@ $^

 %.pdb.data.zip : %.pdb.data
    zip $@ $<

 %.pdb.data : %.pdb
    $(AWKPROG) -f $(PROCESSOR) $< > $@

 %.pdb : $(PROCESSOR)
    touch $@

Completed configuration file:
::

 # PDB atom summarizer configuration.

 PROCESSOR=program.awk
 AWKPROG=awk
 TARFILE=PDBAnalysis

Conclusion
----------

Cover AutomationMake.pptx, slide 15-16.
