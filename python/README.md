Scientific Python Basics
========================

Goals
-----

Despite what the title above might suggest, the purpose of this Software 
Carpentry bootcamp is __not__ to teach you how to program in Python. While we 
do love Python for scientific computing, the goals of these modules is actually 
to teach you the basic, core concepts of programming that transcend languages, 
how they fit together, and how you can use them to become a better scientist.

By the end of these scientific Python lessons, you will be able to:

1.	Describe and distinguish the seven core elements shared by all programming 
	languages.
2.	Use Python to write simple programs that use these core elements, using 
	both the core library and scientific packages such as numpy.
3.	Make and save simple publication-quality plots using matplotlib.
4.	Read, manipulate, and save data files in csv and text formats.
5.	Write unit tests to confirm the accuracy of your Python code.
6.	Create standalone Python scripts that can be run from the command line.

We (the instructors) recognize that you all unavoidably have come with very 
different levels of background in Python programming. We expect that some of 
you might be experienced in basic Python and a few of you of you have 
experience with additional modules such as numpy and scipy - for those in that 
category, this section of the workshop may not be as novel as the other 
sections. However, we hope that the method of presentation will help to 
solidify your existing knowledge. We also encourage you to take the opportunity 
to ask the instructors and volunteers about more advanced techniques that you 
might have heard of but do not know how to use well.

For those who have no (or almost no) background in programming in any language, 
you may find that these lessons proceed quickly. We encourage you to make 
liberal use of the helpful volunteers as we proceed through these lessons. You 
may also wish to consider working together with a partner to complete the 
exercises as a team.

Regardless of your background, you will probably feel like trying to take in 
all of this material is like trying to drink from a firehose. That's OK - the 
idea is to at least introduce you to a wide variety of topics, with the hope 
that you (a) will get to reinforce the most important concepts during 
exercises, and (b) will be able to come back to these materials later to 
continue mastering the concepts.

As an analogy for our overall strategy for teaching Python, picture becoming a 
scientific programmer as constructing a house. First you pour the foundation, 
then you frame the house with wood beams, then you put on a roof, seal it up, 
and furnish it (or something like that). We are hoping that many of you have 
arrived today with the foundation already poured - that is, you understand 
basically how to get some things done using small scripts or programs in some 
language. In these lessons, we will mostly focus on constructing the frame of 
the house together by putting in place the major structures and concepts that 
support all scientific programs. Once the workshop is over, you'll need to go 
out on your own and fill in this structure with drywall, insulation, and a 
couch. (We'll try to point out some useful fixtures that you might want to add 
to your house as we go along.)

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
The distinctions between these categories is not particularly relevant to the 
computer - they are purely a conceptual framework that helps programmers write 
code that does what they want it to do.

We expect that you'll find the basics of 1 and 2 fairly straightforward. We'll 
go quickly through those and will spend the most time on items 3-6. We won't 
really talk about 7, as it is not as common in scientific Python programming as 
it is in, say, shell scripting (pipes and redirection).

Don't worry if you don't already know what all of the above examples mean - 
you'll know by the end of this lesson.

Starting an IPython Notebook
----------------------------

To learn about these core concepts and the Python language, we'll start off by 
working within the IPython notebook.

To start up a notebook, open Terminal and navigate to the folder containing the 
ipynb notebook files that you wish to open (or to any directory in which you'd 
like to save a new notebook, if you're creating a new notebook from scratch). 
Once in the directory, run the command `ipython notebook`, which will launch a 
local webserver and open your default browser. From there you can open an 
existing notebook, create a new notebook, and start working.

Note that if you are on a Windows machine, this command will probably not run 
in Git bash or Cygwin. Instead, open a Command Prompt (click on the Start menu 
and type `cmd` in the search box for Windows 7, or click on Run then type cmd 
for Windows XP), navigate to the appropriate directory, and run the command 
`ipython notebook`.

Asking Questions
----------------

As we go through this lesson, you can ask questions in two ways:

1.	If you have a question for me, just raise your hand and ask.
2.	If you have a question that you think might be restricted to just you (like 
	something on your computer isn't working), raise your hand an a volunteer 
	will come over to help you individually.


[1]: 
http://software-carpentry.org/2012/08/applying-pedagogical-principles-in-this-course/
