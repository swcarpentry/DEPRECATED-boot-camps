
Bash Shell Key Points
=====================

Steve Crouch, Mike Jackson, The Software Sustainability Institute. Steve McGough, Newcastle University.

This work is licensed under the Creative Commons Attribution License. Copyright (c) Software Carpentry and Southampton University and The University of Edinburgh and Newcastle University 2010-2012. See http://software-carpentry.org/license.html for more information.

.. Written in reStructuredText, http://docutils.sourceforge.net/rst.html.

Prerequisites
-------------

Bash

Introduction
------------

Cover BashShell.pptx, slides 1 to 4.

Bash, "Bourne again shell", is the most popular.

Download files via the command-line
-----------------------------------

wget is a simple way to download something via the command-line.
::

 wget --no-check-certificate https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/bash-materials.zip 

Alternatively, there is curl which may be needed for MacOSX.
::

 curl -L --insecure https://bitbucket.org/softwaresaved/boot-camp-edinburgh-1212/downloads/bash-materials.zip -o bash-materials.zip

Unzip.
::

 unzip bash-materials.zip

Basics
------

::

 # Bash comments start with a hash.

Use comments in bash scripts, covered below. Useful for demonstrating at the command-line too.

Who is logged on.
::

 who

If running on a standalone laptop or Cygwin this may be empty.

Who am I logged in as.
::

 whoami

Current directory, the bash name for a folder.
::

 pwd

Change into a directory. Use 'up arrow' to return to previous commands in the command 'history'.
::

 cd bash-materials
 pwd

List contents of current directory.
::

 ls 

List files ending in \*.txt.
::

 ls *.txt

\* is a wild-card. \*.txt expands to a list of matching files. This is called globbing.
::

 ls h*
 ls *_hai*

-R option recurses into sub-directories.
::

 ls -R

-F option shows directory names with / appended.
::

 ls -F

-s option shows block sizes. Block size is system dependent, often 512 bytes (0.5K).
::

 ls -s

-l option shows permissions, dates, sizes, owner, group, size in byte, creation/modification date/time, name.
::

 ls -l

-a option shows hidden files too. . is current directory and .. is parent directory
::

 ls -a 
 ls .
 ls ..

Change to parent directory.
::

 cd ..
 cd bash-materials

View file contents and parts of the file.
::

 cat haiku.txt
 head haiku.txt
 head -2 haiku.txt
 tail haiku.txt
 tail -2 haiku.txt

Copy files, make directory and copy files.
::

 cp haiku.txt another_haiku.txt
 cat another_haiku.txt
 mkdir haikus
 cp *.txt haikus

Copy directory and all its contents. -r option specifies recursion.
::

  cp -r haikus more_haikus

Remove files and directory.
::

 rm another_haiku.txt
 rm more_haikus/more_haikus.txt
 rm -rf more_haikus

-r option specifies recursion and -f option deletes files without asking.

There is no recycle bin - it is gone forever!
::

 mkdir yet_more_haikus

Empty directories can be removed with rmdir.
::

 rmdir yet_more_haikus

Word count
----------

Word count shows number lines, words, characters. It filters a file.
::

 wc haiku.txt

-l and -w options specify just lines or just words.
::

 wc -l haiku.txt
 wc -w haiku.txt

If run across multiple files, it displays a total.
::

 wc -l *.txt

Use to find out number of records in a data file if one file per line, for example.

Finding text
------------

Find text in files with grep (global/regular expression/print).

::

 grep the haiku.txt
 grep day haiku.txt
 grep is haiku.txt
 grep 'it is' haiku.txt

-w option restricts to an exact match.
::

 grep -w is haiku.txt

-n option shows lines where matches are.
::

 grep -n it haiku.txt

-i option ignores case.
::

 grep -i the haiku.txt

Can combine options.
::

 grep -wn is haiku.txt

-v option shows all non-matching lines.
::

 grep -wnv is haiku.txt

-r option recurses into sub-directories.
::

 grep  -wnr Today *

Many other options. To find out more about bash commands, check out the manual.
::

 man grep

Redirecting input and output
----------------------------

How can the matches be saved in a new file?
::

 grep -r not * > found_nots.txt
 cat found_nots.txt

> redirects output (otherwise known as standard output).
::

 ls *.txt > txt_files.txt
 cat txt_files.txt

cat by itself will echo input from what is called the standard input.
::

 cat

Exit with control-D.
::

 cat > myscript.txt
 This is a test!
 Yes it is!
 CTRL-D
 cat myscript.txt

< redirects input (the standard input).
::

 cat < haiku.txt

cat takes standard input from the file. This is not the same as "cat haiku.txt" in which cat is given a file name, even though the output/result is identical.
::

 ls idontexist.txt > output.txt
 cat output.txt

The error message is output on what is called the standard error.
::

 ls idontexist.txt 2> output.txt

Standard error is 2 and standard output is 1.
::

 ls haiku.txt 1> output.txt

To get standard output and error in the same file.
::

 ls idontexist.txt haiku.txt > output.txt 2>&1

Exercise 1 - grep 
-----------------

Cover BashShell.pptx, slide 5.

Finding files
-------------

Find all files and directories.
::

 find .

-type option finds all directories or files.
::

 find . -type d
 find . -type f

-maxdepth and -mindepth options specifies maximum and minimum depth of search.
::

 find . -maxdepth 2 -type f
 find . -mindepth 3 -type f

-perm option specifies files with specific permissions e.g. files user (u) can execute (x).
::

 find . -perm -u=x

-name option allows files with a specific name or pattern to be found.
::

 find . -name *.txt

This gives an error as the wild-card is expanded. The correct way is to use quotes.
::

 find . -name '*.txt'

-iname option ignores case.
::

 find . -iname '*.TXT'

-empty option matches empty files.
::

 find . -empty
 touch emptyfile.txt
 find . -empty

`` back-ticks allow the list of files to be passed to another command. 
::

 wc -l `find . -name '*.txt'`

Exercise 2 - find 
-----------------

Cover BashShell.pptx, slide 6.

Pipes and filters
-----------------

Count text files.
::

 find . -name '*.txt' > files.tmp
 wc -l files.tmp

All shell commands produce text output. All shell commands can take text input.

Connect the output from one command to the input of the next command by a pipe.
::

 find . -name '*.txt' | wc -l

In this context, find and wc are filters and | is a pipe.
::

 echo "Number of .txt files:" ; find . -name '*.txt' | wc -l

; separates commands. It is equivalent to running the two commands on separate lines.

It is not the same as a pipe.

Question: what does this do?
::

 ls | grep s | wc -l

Answer: counts the number of files with the letter "s" in their name.

"Little pieces loosely joined".

Exercise 3 - pipes
------------------

Cover BashShell.pptx, slide 7.

Variables
---------

Shells, like programming languages, support variables.

set shows the current variables.
::

 set

Assign a value to a variable and then see its value.
::

 MYFILE=data.txt
 echo $MYFILE
 echo "My file name is $MYFILE"

Spawn a new shell and try again.
::

 bash
 echo $MYFILE
 CTRL-D

Variables are not inherited by a new shell. Export allows a new shell to use the variables.
::

 export MYFILE
 bash
 echo $MYFILE
 CTRL-D

Bash scripts
------------

Bash supports commands similar to programming languages.

Conditional if statements
::

 NUM=1
 if [ "$NUM" -eq 1 ]; then echo "Equal"; fi

String comparisons.
::

 WORD="hello"
 if [ "$WORD" = "hello" ];  then echo "The same"; fi

Arithmetic.
::

 let NUM=$NUM+1

Save the output of a command in a variable.

::

 TEXT_FILES=`ls *.txt`
 echo TEXT_FILES

Loops.
::

 for PDB in `find . -name '*.pdb'`; do
     echo $PDB
 done

Typing in the same command sequences over and over is time-wasting, error prone, and boring. Automate.
::

 nano protein_filter.sh

Add
::

 #!/bin/bash
 DATE=`date`
 echo "Processing date: $DATE"
 for PDB in `find . -name '*.pdb'`; do
     echo $PDB
 done
 echo "Processing completed!"

Save and run.
::

 sh protein_filter.sh

chmod makes this executable.
::

 chmod +x protein_filter.sh
 ./protein_filter.sh

Exercise 4 - shell scripts
--------------------------

Cover BashShell.pptx, slide 8.

Files, directories and permissions
----------------------------------

"ls -l" shows permissions, dates, sizes, owner, group, size in byte, creation/modification date/time, name.
::

 ls -l haiku.txt

chmod can add, remove or set permissions.

Add permission to allow (+) all (a) to read (r) the file.
::

 chmod a+r haiku.txt

Remove permission for (-) all to read the file.
::

 chmod a-r haiku.txt

Add permission to allow user (u) to read the file.
::

 chmod u+r haiku.txt

Add permission to allow group (g) to write (w) the file.
::

 chmod g+w haiku.txt

Add permission to allow others (o) to execute (x) the file.
::

 chmod o+x haiku.txt

Add permission to allow group to also read and execute the file.
::

 chmod g+rx haiku.txt

Add permission to allow user and group to also read and execute the file.
::

 chmod go+rx haiku.txt

Set permission (=) explicitly to allow user, group and others to read, write and execute.
::

 chmod ugo=rwx haiku.txt

Job control
-----------

:: 

 ./counter.sh > output.txt

CTRL^Z suspends this process, or job.
::

 wc -l output.txt

jobs shows a list of jobs. -l option shows job number and process ID.
::

 jobs -l

fg resumes it in the foreground. fg can take a job number as an argument.
::

 fg 
 CTRL^Z

bg resumes it in the background. bg can take a job number as an argument.
::

 bg
 wc -l output.txt

Jobs can be started in the background by default.
::
 
./counter.sh > output.txt &

The left number is the job number and the right the process ID.

kill kills a process with a given process ID.
::

 ./counter.sh > output.txt &
 kill NNNN
 jobs

ps shows more detail on processes.
::

 ps

top shows resource consumption.
::

 top

nohup allows processes to continue even after the user logs out.
::

 bash
 nohup ./counter.sh > output.txt &
 CTRL-D
 wc -l output.txt

Secure shell
------------

Log into a remote server.
::

 ssh username@boot-camp.software-carpentry.org

Format is username AT host name.

Run a command remotely.
::

 ssh username@boot-camp.software-carpentry.org ls

Secure copy a file to a remote server.
::

 scp file.txt username@boot-camp.software-carpentry.org:

Format is username AT hostname COLON relative path.

Secure copy a file from a remote server.
::

 scp username@boot-camp.software-carpentry.org:data-files/data.txt

Script
-------

For Linux users.
::

 script
 ls -l
 CTRL-D
 cat typescript

Useful to record commands typed, commands with lots of outputs, trial-and-error when building software. 

Turn into blog or tutorial. Send exact copy of command and error message to support.

Conclusion
----------
Cover BashShell.pptx, slide 9.
