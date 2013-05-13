##2. Make - patterns and rules.
**Based on materials by Steve Crouch and Greg Wilson**

Let’s look at our dependency graph again. Two big questions: 
* How are we going to handle multiple dependencies efficiently? How do we generalize?
* Can we get rid of the repeated filenames in the dependencies?

Let’s look at this now.

Now, for our tar task, we could add the following into our Makefile (replacing the all target):

    PDBAnalysis.tar.gz : cubane.pdb.data ethane.pdb.data methane.pdb.data
        tar -czf PDBAnalysis.tar.gz cubane.pdb.data ethane.pdb.data methane.pdb.data
Firstly, we can cut down the duplication a bit by changing the action:

    tar -czf $@ cubane.pdb.data ethane.pdb.data methane.pdb.data

This is better. This is a make automatic variable. It means ‘the target of the current rule’. In this case, our tar file. Now, we’d only need to change the .tar.gz filename in one place, if we wanted to. Yes, it’s cryptic! Another make idiosyncrasy!

But can we shorted the duplication of all the dependent files in the action? Yes:

    tar -czf $@ $^
...with another of make’s cryptic notations - this one means ‘the dependencies of this rule’.

But we also really want to say in our target dependencies: ‘all the files named something.pdb’. Might not know how many there will be! Can’t duplicate and modify the target every time - maintenance problem. Asking for trouble. Would like to use something like Bash’s wildcard *! The good news is we can! We can change the rule to say:

    PDBAnalysis.tar.gz : *.pdb.data
BUT when we use this, we must refer to dependencies using $^. Don’t know which filenames will match, could have many!

Also, there’s a clear dependency on our awk program. If this changes, it should properly enable a rebuild of the files, e.g. the output file format may have changed. There are some bad ways to do this. We could make our .pdb files dependent on our awk program - a false dependency, for example. But we would need to do this for all pdb files at the moment. Instead, we can add additional rules for our .pdb.data files to depend on this program, which is far more logical.
Add to end of makefile:

    cubane.pdb.data : program.awk
    ethane.pdb.data : program.awk
    methane.pdb.data : program.awk
And, pretending we’ve changed it:

    $ touch program.awk
    $ make -f pdbprocess.mk
All should rebuild!


###Using pattern rules for more generic rules
So what we have is fully functional. But it’s not easily maintainable:
* Add more pdb files to process, a lot of duplications in the rules.
* If we change the way we process pdbs, have to change more than once place.

To do in make, we use a __pattern rule__. Change our .pdb.data targets to a single new target, beginning with rule:

    %.pdb.data : %.pdb
In this rule, % is a wildcard. When expanded, it matches same on the left of the : as on the right. For example, with cubane.pdb.data, % matches on cubane on the left, and ditto on the right for cubane.pdb.

###Exercise 2:
As we did with the tar rule, simplify the awk rule.
* Instinctively, using $@ and $^ again in the right places would seem to work.
* Try this, and observe the results.

###Solution:
    # pdbprocess.mk
    PDBAnalysis.tar.gz : *.pdb.data
        tar -czf $@ $^

    %.pdb.data : %.pdb
        awk -f program.awk $^ > $@

    cubane.pdb.data : program.awk
    ethane.pdb.data : program.awk
    methane.pdb.data : program.awk

The awk command included program.awk in the processing! Why? Because $^ matches on all dependencies! Which includes program.awk in rule at the end of the file! So by using another make crypticism:

    awk -f program.awk $< > $@
We specify here we only want the first dependency in the list (only our .pdb.data file). 


Next: [Macros.](3_Macros.md)




