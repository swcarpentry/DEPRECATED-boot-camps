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


Previous: [Patterns and Rules](2_Patterns_Rules.md)
