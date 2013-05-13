##3. Make - Macros
**Based on materials by Steve Crouch and Greg Wilson**

Let’s assume someone else wants to run our code on our data. First, we could include our awk program in the tar archive:

    PDBAnalysis.tar.gz : *.pdb.data program.awk
When we do:

    $ touch program.awk
    $ make -f pdbprocess.mk
program.awk is now included in the archive, and can be run. But let us assume we have a number of these awk programs:

    $ cp program.awk another-program.awk
And assume another-program.awk is different. Also assume our makefile is larger and contains lots of references to this program filename. Again, we have a maintenance problem. Need to change many instances of filename. How can we generalize here? We can use a make macro.

At the top of the file after comment, add:

    PROCESSOR=program.awk
And for every instance of program.awk, include instead:

    $(PROCESSOR)
This now refers to our program name as defined. Why use parenthesis? It’s yet another make idiosyncrasy.

Now, what if someone else only had gawk on their system, not awk? It’s fully compatible with awk, but is only one available. Add at top:

    AWKPROG=awk
And where ‘awk’ is mentioned:

    $(AWKPROG) -f $(PROCESSOR) $< > $@
For clarity, it’s good practice to use a macro for these kinds of things too. Put configurable stuff at the top, and it’s obvious for someone else what might need to be changed.

For clarity, we can even put these configuration options in another file. Extract the macro definitions at the top and put them in a file called config.mk:

    # config.mk
    PROCESSOR=program.awk
    AWKPROG.awk
And in our makefile, instead put:

    include config.mk

Why do this?
* Different configurations for different locations you run it
* Just need to change config file
* Avoids program duplication for different configurations
* Separation of concerns: if program is under version control, you’re only changing the aspect you need to (config or program).


###Exercise 4:
Add in a new macro to our configuration so we can change the name of our PDBAnalysis.tar.gz file if desired.


###Soluntion
pdbprocess.mk:

    # pdbprocess.mk
    include config.mk

    $(TARFILE).tar.gz : *.pdb.data.zip $(PROCESSOR)
        tar -czf $@ $^

    %.pdb.data.zip : %.pdb.data
        zip $@ $<

    %.pdb.data : %.pdb
        $(AWKPROG) -f $(PROCESSOR) $< > $@

    %.pdb : $(PROCESSOR)
        touch $@
config.mk:

    # config.mk

    PROCESSOR=program.awk
    AWKPROG=awk
    TARFILE=PDBAnalysis


Previous: [Patterns and Rules](2_Patterns_Rules.md)
