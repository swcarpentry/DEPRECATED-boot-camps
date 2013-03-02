* [The Shell][#The Shell]
* [Access the Shell][#access the shell]
* [pwd] [#pwd]



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

Most system default to the Bourne-again (bash) shell,  we will
use *bash* in these tutorials.

Unix philosophy:
**Make each program do one thing well**


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

## Navigation

### pwd

To find out where you are in the filesystem use **pwd** (print working directory):

    pwd

Which will give results looking like the following

using (msysgit prompt)

    /c/Documents and Settings/Administrator/My Documents/boot-camps

using (ipython prompt)

    C:\\Documents and Settings\\Administrator\\My Documents\\boot-camps

using (OSX)

    /Users/cindeem/Documents/boot-camps

using (linux)

    /home/jagust/cindeem/boot-camps

Note that the results are dependent on the OS/terminal combo you are using...
from here on out in the tutorial, I will stick to a standard OSX terminal output
(very simliar to linux output). 

### cd

**cd** allows you to change directories.

* Used above to get into boot-camps directory
* We want to move into the **shell** directory

command to move to shell directory:

    cd shell
    pwd

    /Users/cindeem/Documents/boot-camps/shell

If we wanted to go backward to the original directory we use ../ to signify going backward,
use **pwd** to get a sense of where you are:

    cd ../
    pwd
    /Users/cindeem/Documents/boot-camps/


You can be anywhere in the filesystem, if you want to go to your home directory just typing **cd** alone
will take you home:

    cd
    pwd
    /Users/cindeem

Lets get back to our shell directory, we can use a full path directive to get there:

    cd /Users/cindeem/Documents/boot-camps/shell
    pwd

    /Users/cindeem/Documents/boot-camps/shell


### touch

touch creates an empty file:

    touch myemptyfile.txt

### ls

To see the file you just created type:

    ls

    Readme.md		ex_data.txt		    linux-term.jpg
    data			generate_data.py	linux-term2.jpg
    dictionary.txt	hello			    myemptyfile.txt

You can also use a *relative* path with ls. For example if I want to see what is in data

    ls data

We could also look up a directory

    ls ../

Or put in an absolute or **full path**  eg /Users/cindeem/Documents/boot-camps/shell/data

    ls /Users/cindeem/Documents/boot-camps/shell/data

Most unix commands will also take additional arguments or flags.

For example, to get more information use the **-l** flag (stands for *long format*):

    ls -l

    total 88K
    drwxrwxr-x 9 cindeem staff  139 Mar  2 14:38 data/
    -rw-rw-r-- 1 cindeem staff 5.2K Mar  2 14:38 dictionary.txt
    -rw-rw-r-- 1 cindeem staff  236 Mar  2 14:38 ex_data.txt
    -rw-rw-r-- 1 cindeem staff 8.5K Mar  2 14:38 generate_data.py
    -rwxrwxr-x 1 cindeem staff   33 Mar  2 14:38 hello*
    -rw-rw-r-- 1 cindeem staff  41K Mar  2 14:38 linux-term2.jpg
    -rw-rw-r-- 1 cindeem staff 9.6K Mar  2 14:38 linux-term.jpg
    -rw-rw-r-- 1 cindeem staff    0 Mar  2 14:39 myemptyfile.txt
    -rw-rw-r-- 1 cindeem staff 3.4K Mar  2 14:38 Readme.md

If we look at the output we see a number of things, consider your new file *myemptyfile.txt*

**-rw-rw-r--** 

The first part gives us info about permissions on the file (we will talk about this later).
But the first character **-**, signifies that this item is a file, not a directory.

    -rw-rw-r-- 1 cindeem staff    0 Mar  2 14:39 myemptyfile.txt

In contrast, look at *data*. Data is a directory, so the first character in **drwxrwxr-x** is **d**.

    drwxrwxr-x 9 cindeem staff  139 Mar  2 14:38 data/

Both show the user, group pair **cindeem staff**.  this tells use a little about who created the file.

We also see timestamps on the file, showing it was created (or modified) on **Mar  2  14:38** (March 2 at 2:38pm).

    drwxrwxr-x 9 cindeem staff  139 Mar  2 14:38 data/

Lets create a new file, but we are going to do something odd, we are going to add a **.** to the beginning of the
filename:

    touch .hiddenfile

If we use **ls** again this file will not show up, this is beacuse of the leading **.**

However, we can see these files using the **-a** flag:

    ls -a

    ./     dictionary.txt    hello*           linux-term.jpg
    ../    ex_data.txt       .hiddenfile      myemptyfile.txt
    data/  generate_data.py  linux-term2.jpg  Readme.md

We can now see the *hidden* file (.hiddenfile). Often these files are configuration files or temporary files, 
and in general you do not edit or work with them. 
But they will be useful to know about, and sometimes you do want to access them.

Lets see one more useful **ls** flag, **-rt**.  This flag combo (-t for *"order by time"* and -r for *"reverse"*) will
order your files by the time they were last changed, but with the oldest first, and youngest last.
In the example below I will also add **-l** so we can see the timestamps:
       
    -rwxr-xr-x  1 cindeem  staff     33 Feb 27 21:28 hello
    -rw-r--r--  1 cindeem  staff   8667 Feb 27 21:28 generate_data.py
    -rw-r--r--  1 cindeem  staff    236 Feb 27 21:28 ex_data.txt
    -rw-r--r--  1 cindeem  staff   5321 Feb 27 21:28 dictionary.txt
    drwxr-xr-x  9 cindeem  staff    306 Feb 27 21:28 data
    -rw-r--r--  1 cindeem  staff  41166 Mar  2 13:32 linux-term2.jpg
    -rw-r--r--  1 cindeem  staff   9737 Mar  2 13:32 linux-term.jpg
    -rw-r--r--  1 cindeem  staff   3702 Mar  2 14:21 Readme.md
    -rw-r--r--  1 cindeem  staff      0 Mar  2 15:03 myemptyfile.txt
           
This can be a useful tool if you want to see what you have most recently changed in a large directory.

## rm

Sometimes we make things we dont want to keep. This is where **rm** (remove) comes in.  Lets get rid of our hidden file,
and then us **ls -a**  to make sure it is gone:

    rm .hiddenfile
    ls -a

**rm** can be a little scary, what if we had chosen **rm data**, would we lose all our precious data?

    rm data

    rm: data: is a directory

Luckily it will only remove local files.

rm also has a **-r** (recursive) flag. to get rid of a directory and all of its
contents you could use rm -r, but be careful, it will bite you. (with great power comes great responsibility)








        





