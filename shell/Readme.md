# The Shell

What is Shell/ Terminal?

The *shell* is a program

* presents a command line interface
* allows you to control your computer using commands entered
  with a keyboard 

There are many different flavors of Shells/Terminals

* csh  - C-shell
* tcsh - turbo C-shell
* ksh  - Korn shell 
* sh   - Bourne shell
* bash - Bourne again shell



Unix philosophy:
**Make each program do one thing well**

Most system default to the Bourne-again (bash) shell,  we will
use *bash* in these tutorials.

# Access The Shell

## Linux

Depending on the flavor of linux you use Terminal

Centos6 example

Applications Menu -> Teminal Emulator

![Linux Term](linux-term.jpg "Linux Term")
![Linux Term2](linux-term2.jpg "Linux Term2")

## Mac OS

Applications -> Utilities -> Terminal

## Install msysgit on Windows

[msysgit download link](https://code.google.com/p/msysgit/downloads/list?q=label:Featured)

eg  on a windowsXP I installed **msysGit-fullinstall-1.8.1.2-preview20130201.exe**	

double-Click to install, when finished you will see a window with a prompt

I wanted to add a shortcut on my Desktop:

    /share/msysGit/add-shortcut.tcl Desktop

# Some Data

Shell has many uses, today we are going to focus on using it to manage and setup
directories holding project data and related files.  

Keeping your projects organized is am important as
keeping an organized lab book.

[A Quick Guide to Organizinf Computational Biology Projects](http://bit.ly/AADX8F)
Noble WS (2009) A Quick Guide to Organizing Computational Biology Projects. PLoS Comput Biol 5(7): e1000424. doi:10.1371/journal.pcbi.1000424

But often we inherit badly organized and badly labelled data:

To get
the data for this test, you will need internet access and an open shell. 
Just enter the command:

    git clone -b YYYY-MM-PLACE https://github.com/USERNAME/boot-camps.git

Followed by:

    cd boot-camps
    git checkout YYYY-MM-PLACE

You should now have all the data you need for this tutorial (yea!)

** Navigation

To find out where you are in the filesystem use **pwd** (print working directory):

    pwd

msysgit prompt

/c/Documents and Settings/Administrator/My Documents/boot-camps

or (via ipython prompt)

C:\\Documents and Settings\\Administrator\\My Documents\\boot-camps

or (OSX)

/Users/cindeem/Documents/Talks/softwarecarpentry/boot-camps

or (linux)

/home/jagust/cindeem

