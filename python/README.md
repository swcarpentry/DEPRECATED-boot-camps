Introduction to Python
======================

Goals
-----

Despite what the title above might suggest, the purpose of this Software 
Carpentry lesson is specifically __not__ to teach you how to program in Python. 
While we do love Python for scientific computing, the goal of this lesson is 
actually to teach you the basic, core concepts of programming that transcend 
languages, how they fit together, and how you can use them to become a better 
scientist.

Specifically, by the end of this lesson, you will be able to:

1.	Describe and distinguish the seven core elements shared by all programming 
	languages.
2.	Use Python to write programs that use these core elements, using both the 
	core library and packages such as numpy and matplotlib.
3.	Write a program in Python that will simulate stochastic, logistic 
	population growth.

We (the instructors) are assuming that you have some background in programming, 
either in Python or another language. We will not attempt to teach you 
programming from the ground up (this would take far more than an afternoon), 
and will similarly not be able to teach you to be advanced scientific Python 
programmers in just an afternoon.

As an analogy for our goals, picture becoming a scientific programmer as 
constructing a house. First you pour the foundation, then you frame the house 
with wood beams, then you put on a roof, seal it up, and furnish it (or 
something like that). We are assuming that you arrived today with the 
foundation already poured - that is, you understand basically how to get some 
things done using small scripts or programs. This afternoon, we will construct 
the frame of the house together by putting in place the major structures and 
concepts that support all scientific programs. Once the workshop is over, 
you'll need to go out on your own and fill in this structure with drywall, 
insulation, and a couch. (We'll try to point out some useful fixtures that you 
might want to add to your house as we go along.)

As we go through this lesson, you may feel like trying to take all of this in 
is like drinking from a firehose. That's OK - the idea is to at least introduce 
you to a wide variety of topics, with the hope that you'll (a) get to reinforce 
the most important concepts during exercises, and (b) will be able to come back 
to these materials later to continue mastering the concepts.

Finally, it's worth noting that you all (unavoidably) have come with very 
different levels of background in Python programming. One-third of you are very 
experienced in basic Python and one-quarter of you are experienced with modules 
such as numpy and scipy - for those in that category, this section of the 
workshop may not be as novel as the other sections. However, we hope that the 
method of presentation in this lesson will help to solidify your existing 
knowledge. We also encourage you to take the opportunity to ask the instructors 
and volunteers about more advanced techniques that you might have heard of but 
do not know how to use well.

The Seven Core Concepts
-----------------------

As noted by Greg Wilson (the founder of Software Carpentry), every programming 
language shares [seven core elements][1]:

1.	Individual things (the number 2, the string 'hello', a matplotlib figure)
2.	Commands that operate on things (the + symbol, the `len` function)
3.	Groups of things (Python lists, tuples, and dictionaries)
4.	Ways to repeat yourself (for and while loops)
5.	Ways to make choices (if and try statements)
6.	Ways to create chunks (functions, objects/classes, and modules)
7.	Ways to combine chunks (function composition)

The lines between these are often blurry in practice - for example, a string in 
Python actually mixes some characteristics of a thing, a group, and a "chunk". 
The distinction between these categories is not particularly relevant to the 
computer - it's purely a conceptual framework that helps programmers write code 
that does what they want it to do.

Don't worry if you don't already know what all of the above examples mean - 
you'll know by the end of this lesson.

We expect that you'll find the basics of 1 and 2 fairly straightforward, and 
will go quickly through those, and will spend the most time on items 3-6. We 
won't really talk about 7, as it is not as common in scientific Python 
programming as it is in, say, shell scripting (pipes and redirection).

Starting an IPython Notebook
----------------------------

To learn about these core concepts and the Python language, we'll be using the 
IPython notebook.

To start up the notebook, open Terminal and navigate to the folder containing 
the ipynb notebook files that you wish to open (or to any directory in which 
you'd like to save a new notebook, if you're creating a new notebook from 
scratch). Once in the right directory, run the command `ipython notebook`, 
which will launch a local webserver and open your default browser. From there 
you can open an existing notebook, create a new notebook, and start working.

Note that if you are on a Windows machine, this command will probably not run 
in Git bash or Cygwin. Instead, open a Command Prompt (click on the Start menu 
and type cmd in the search box for Windows 7, or click on Run then type cmd for 
Windows XP), navigate to the appropriate directory, and run the command 
`ipython notebook`.

Asking Questions and Polls
--------------------------

As we go through this lesson, you can ask questions in two ways:

1.	If you have a question for me, just raise your hand and ask.
2.	If you have a question that you think might be restricted to just you (like 
	something on your computer isn't working), raise your hand an a volunteer 
	will come over to help you individually.



[1]: 
http://software-carpentry.org/2012/08/applying-pedagogical-principles-in-this-course/
