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
      




