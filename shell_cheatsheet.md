# Basic Shell Commands
***

## 1. Shell Basics:

| Symbol         | Definition                                                                                                     |
|----------------|----------------------------------------------------------------------------------------------------------------|  
| `.`            | a single period refers to the current directory                                                                |  
| `..`           | a double period refers to the directory immediately above the current directory                                |  
| `~`            | refers to your home directory. _Note:_ this command does NOT work on Windows machines (Mac and Linux are okay) |  
| `cd ./dirname` | changes the current directory to the directory `dirname`                                                       |  
| `ls -F`        | tells you what files and directories are in the current directory                                              |  



## 2. Creating Things:
### a) How to create new files and directories...
* **`mkdir ./dirname`** --> makes a new directory called dirname below the current directory. _Note:_ Windows users will need to use `\` instead of `/` for the path separator
* **`nano filename`** --> if `filename` does not exist, `nano` creates it and opens the `nano` text editor. If the file exists, `nano` opens it. _Note:_ _(i)_ You can use a different text editor if you like.  In gnome Linux, `gedit` works really well too. _(ii)_ `nano` (or `gedit`) create text files. It doesn't matter what the file extension is (or if there is one)

### b) How to delete files and directories...
#### _Remember that deleting is forever. There is NO going back_
* **`rm ./filename`** --> deletes a file called `filename` from the current directory 
* **`rmdir ./dirname`** --> deletes the directory `dirname` from the current directory. _Note:_ `dirname` must be empty for `rmdir` to run.

### c) How to copy and rename files and directories...
* **`mv tmp/filename .`** --> moves the file `filename` from the directory `tmp` to the current directory. _Note:_ _(i)_ the original `filename` in `tmp` is deleted. _(ii)_ `mv` can also be used to rename files (e.g., `mv filename newname`
* **`cp tmp/filename .`** --> copies the file `filename` from the directory `tmp` to the current directory. _Note:_ _(i)_ the original file is still there



## 3. Pipes and Filters
### a) How to use wildcards to match filenames...
Wildcards are a shell feature that makes the command line much more powerful than any GUI file managers. 


** Table of commonly used wildcards **

| Wildcard               | Matches                                        |  
|------------------------|------------------------------------------------|  
| `*`                    | zero or more characters                        |  
| `?`                    | exactly one character                          |  
| `[abcde]`              | exactly one of the characters listed           |  
| `[a-e]`                | exactly one character in the given range       |  
| `[!abcde]`             | any character not listed                       |  
| `[!a-e]`               | any character that is not in the given range   |  
| `{software,carpentry}` | exactly one entire word from the options given |  



### b) That wildcards are expanded by the shell before commands are run...
### c) How to redirect a command's output to a file...
### d) How to redirect a command's input from a file...
### e) How to use the output of one command as the input to another with a pipe...
### f) That combining single-purpose filters with pipes is the most productive way to use the shell...
### g) That if a program conforms to Unix conventions, it can easily be combined with others...



## 4. Variables
### a) Assignment
* **`varname=1`** -->

### b) Indexing 
* **`varname[0]`** --> _Note:_ the shell is zero indexed.  That means you always start counting from zero

### c) Referencing
* **`${varname}` -->
* **`${varname[@]` --> 

 

## 5. Loops
NEED TO DO VARIABLE ASSIGNMENT FIRST!!!!
### a) How to repeat operations using a loop...
* **`for`** -->  
    `for filename in *.dat
    do
      mv ${filename} ${newname}
    done`
    
* **`while`** -->
    `count=0   
     while ${count} -lte 6
     do
       COMMAND HERE
     done`

### b) That the loop variable takes on a different value each time through the loop...
### c) The difference between a variable's name and its value...
### d) Why spaces and some punctuation characters shouldn't be used in files' names...
### e) How to display history and re-use commands...
* **`history`** --> displays your command history to the standard output (usually the screen)



## 6. Shell Scripts
### a) How to store shell commands in a file...
### b) How to run a shell script...
### c) How to pass filenames into a shell script...



## 7. Finding Things
### a) How to select lines matching patterns in text files...
* **`grep [options] day haiku.txt`** --> finds every instance of the string `day` in the file haiku.txt and pipes it to standard output. 
	* **`-E`** --> tells grep you will be using a regular expression. Enclose the regular expression in quotes. _Note:_ the power of `grep` comes from using regular expressions. Please see the regular expressions sheet for examples
	* **`-i`** --> makes matching case-insensitive
	* **`-n`** --> limits the number of lines that match to the first n matches
	* **`-v`** --> shows lines that do not match the pattern (inverts the match) 			
	* **`-w`** --> outputs instances where the pattern is a whole word

### b) How to find files with certain properties...
* **`find . -type d` -->
	* **`-type [df]`** --> d lists directories; f lists files
	* **`-maxdepth n`** --> `find` automatically searches subdirectories. If you don't want that, specify the number of levels below the working directory you would like to search
	* **`-mindepth n`** --> starts `find`'s search n levels below the working directory
	
### c) How to use one command's output as arguments to another command...

### d) How are text and binary files different?...

