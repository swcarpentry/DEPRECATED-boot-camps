Useful UNIX tools
=========

**Updated and presented by: Tracy Teal**

SSH
--------


What is ssh?
------------------

ssh stands for 'secure shell'.  It's a way of connecting to another computer
without sending your passwords in the clear.  telnet, for example, is like ssh
but your password is not encrypted when it's sent.

We're going to ssh into a Amazon instance that was set up for this class.

If you want more information on Amazon instances, Titus has a nice reference
here: http://ged.msu.edu/angus/tutorials/renting-a-computer-from-amazon.html

You can log in to this machine by typing

    ssh swc@ec2-23-21-29-27.compute-1.amazonaws.com

The password is the name of this institute.

Now you're on the Amazon instance!

Here we see the files we saw in the 1-Shell lesson and we can do all the things 
here to them that we could on our systems yesterday.  Run through a few of the 
shell commands we did yesterday.

Screen
-----------------
Now say you're on this instance and you want to run a program that's going
to take a long time, and you want to be able to log out of the instance
and still have it run.  Here is where 'screen' is your friend.

Type

    screen

You'll get the prompt back again and you can do just what you normally do.  Type 
'ls' for instance.  Now, though, type 'exit'.  You're still at a prompt.  All you 
did was exit out of screen.

Start screen again, by typing 'screen'

Do an 'ls' again.  Now we want to detach from the screen, so that the process could
keep running while we log off of Amazon or any remote machine.  To detach

   Ctrl-A  Ctrl-D

If you want to reattach to that screen, type 

   screen -r

If you have multiple screens running, it will tell you that and you have to pick 

   screen -r screen-number

That's all there is to screen.

A test example
----------------

Now let's try a test example of logging on to a remote machine, running screen, 
running a blast at the command line, detaching from the screen and then coming 
back to check on the results.

Follow the instructions above for loggin in and starting screen

Now to run command line blast

You can download command line blast from NCBI.  It's called 'blastall'

If you type 
blastall - 
You get all the options

This is a standard blastall command

blastall -e 1e-05 -p tblastx -d test_fasta -i seqs.fa  -o seqs.blast

-e is the e-value cutoff you want to use.  Any matches higher than that will not be returned
-p is the program - tblastx, blastx, blastn or tblastn
-d is the database
-i is the input file
-o is the output file
-m is the output type you want
   If you're parsing the output, then you want to use -m 8.  It outputs a tab delimited format that's easy to look through
   The default shows you all the alignments

If you do use -m 8 this is the information in each column

# Query id # Subject id # % identity # alignment length # number of mismatches # number of gap openings # position of query start # position of query end # position of subject start # position of subject end # e-value of a hit # bit score of a hit  

So, we'll run blastall on a test dataset

There are some test fasta files in the 'fasta' directory

Go in to that directory.  You can use any of these files as the input fasta file.

The blast database is: /home/swc/blastdb/nirK_ref_Rh

blastall -e 1e-05 -p tblastx -d /home/swc/blastdb/nirK_ref_Rh -i oneseq.fasta  -o seqs.blast