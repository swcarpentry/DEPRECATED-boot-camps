SSH
=========

**Updated and presented by: Tracy Teal**


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

