Software Carpentry Python Workshop.

**Stanford, May 6-7th**


# The Shell

## Why shell?

The shell is/has

* a program that runs other programs (usr <-> shell <-> cpu)
* text only
* cryptic commands
* many tools *only have* command line interface, esp. remote machines
* very fast and powerful commands
* allows to combine tools in powerful way
	* philosophy is to write small junks of code that do one job, and one job only
	* then combine in greater scripts
* allows to repeate reliably and reproducibly (unlike with the mouse)



## First commands

	whoami
> stickies!

whereami? No, it's
	pwd

absolute path indicated by '/' (the root directory. On Windows '\'!?)

	ls

for listing. Cryptic names.

	ls -F

This is the syntax for shell commands: 

	commandName -optionalParameters target

target can have default values. Combine two optional arguments with

	ls -F -t
	ls -Ft
	lt -tF

unless the *flag* needs its own argument

	mkdir swc

to create working directory for today. No news is good news - sometimes.
	ls, ls swc, cd swc, ls



## Tricks and shortcuts

* tab completion!
* arrow up for previous commands
* q to exit
* crtl + C to *kill*


## Data


	touch data

open files with editor

	nano data

very basic, but integrated (no need to open another program), saves in pwd

	^ = fancy for 'ctrl'

**Task. go to etherpad, add name, price, purchase date, of favorite mobile app**

> copy & paste, save file, exit nano

	cp data data.txt
	mv origin destination

also work on full directory

	rm data

Careful, no trash bin!

**Task. create backup directory in swc, copy data.txt, give it reasonable name**

relative and absolute path
	ls /Users/Berny/data
is same as
	ls data
if you are already in folder above *data*
	.. and . as special directories
	ls ~, cd ~, cd all bring you home





## Work with text files


	cat data.txt, wc data.txt
wc: lines, words, characters. One or more files.

	man wc
for manual entry (or google)

	sort data.txt
	sort data.txt > sortedData.txt

	head sortedData.txt
	rm sortedData.txt

**Task. use *tail* or optional command with *head* and *sort* to create a file that records (alphabetically) 2nd and 3rd last entry**

To do this without intermediate files, use the *pipe*

	sort data.txt | tail -3 | head -2

note that there is no argument for *tail*, *head*

**The pipe is most important concept!** Allows to combine simple commands, and create powerful scripts.

Note, output is to screen only if there is no '|' or '>' or '>>' (append).

	cut -d ',' -f 2 data.txt

**Task. combine sort, cut, uniq to get file with all prices**


## grep and find

*global regular expression print*

	grep 2013 to find all apps from 2013

**Task. Create file that only lists free apps** **Bonus: Only show names, dates**

	grep ' 0.00' to *actually* find all free apps

uses regular expressions, like end$, ^start, etc.


*find looks for files*

	find where what

How do you type the where? 'ls -a' to see '.' and '..'

	find . -type f -name "*.txt" -or -name "*Junk*"

**Task. combine find, grep to see all backup data**

	find where -exec grep -Hn pattern_in_file {};
	find . -type f | xargs grep pattern_in_file


**Task. find all occurences of "Catan"" in any text files**


## History and scripts

like lab notebook, records everything you did. Either run with '!123', or better:

Combine with tail, grep, '>' to create a script

	history | tail -20 | grep cut > listFreeAppsWithDate.sh

Use nano to clean up. To execute, run
	bash listFreeAppsWithDate.sh
	./listFreeAppsWithDate.sh
(maybe you need to run *chmod a+x listFreeAppsWithDate.sh* first)


Then use script as any other shell command.

	listFreeAppsWithData.sh | sort -k 2 -n

**Task. Use listFreeAppsWithData.sh in another script to return name of most recent free app. Try in terminal, then use history, '>' once it works to put all in a script.**

	history | tail -2 > mostRecentFreeApp.sh


## wildcard, variables

> Make copies of data. Use split -l 10 data.txt dataJunk_

	ls file
to just see file.
	ls file1 file2
to see file1 and file2. No limit on file names. Could do
	ls file*
to list all files with name fileX. 


We can do this in our scripts as well

*First*, loop

	for petName in data.txt data1.txt
	do
		echo This is $petName
		wc $petName
		cat $petName
	done

**Task. use for loop to create output with format, (recall 'echo')**

*'The free apps in fileName1 are…'*

*'The free apps in fileName2 are…'*

**Bonus: Loop over all txt files**

		
We can do this with our script as well
	replace data.txt with $1


Or wrap around for loop with
	for fileName in $*
	do
		mostRecentFreeApp.sh $fileName
	done

**Task. Create new folder *allData* and put a copy of any txt file**

$fileName is a variable. Can define our own variable

	myName=IamBatman
	echo $myName
	alias myCommand='cd /Users/Berny/swc'

Useful for commands you repeat over again, eg navigate to certain directory.