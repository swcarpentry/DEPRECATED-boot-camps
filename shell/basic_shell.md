

# Using the shell to do more in less time
**What is the shell?**

The shell is a program that takes keyboard commands and passes them to the operating system to carry out. 

All Linux distributions supply a shell program from the GNU Project called bash. There are other kinds that can be installed in different operating systems.

Before we talk more, let's actually open the shell on your operating system.

My shell prompt looks like this:

![](shell_prompt.png)

Yours might looks different (these can be easily customized). Usually includes something like `username@machinename`, followed by the current working directory (more about that soon) and a $ sign

## Entering commands into the shell

You can just enter commands directly into the shell.

```
echo Morning Bristol
```

We just used a command called `echo` and gave it an argument called `Morning Bristol`.

If you enter a command that shell doesn't recognize, it will just report an error

```
$ gobbeltdfsf
bash: gobbeltdfsf: command not found
```

Now let's enter something useful. Let's navigate to the Desktop of your computer (more on navigation very shortly)

```
cd ~
pwd
```

What does it say?

OK, now let's type in a really useful (and for today a rather important) shell command (*PS: No need to run this if you've successfully done it already.*)

```
$ git clone -b 2013-09-bristol --single-branch https://github.com/swcarpentry/boot-camps.git
```

---

## How commands work in the shell

Commands are often followed by one or more options that modify their behavior, and further, by one or more arguments, the items upon which the command acts. So most commands look kind of like this:

```
e.g. ls -l ~/
```

```
command -option arguments
```



## Knowing where you are and seeing whats around

The first thing you want to do when you're somewhere new is get a map or figure out how to obtain directions. Since you're new to the shell, we're going to do just that.

Three really imporant commands:

* `pwd` *Acronym for print working directory*. Tell you where you are.  
* `cd` *Change directory*. Give it options for where to take you.  
* `ls` *Short for list*. List all the files and folders in your current location.  

Most operating systems have a hierarchical directory structure. The very top is called the `root directory`. It contains files, subfolders and so on. Everything else is nested below.


```
$ pwd
/Users/karthik
```

Note that I'm in my HOME directory. Yours should (hopefully) look different. 


**List all the files in this directory**

```
$ ls
Applications    Documents   Dropbox     Library     Music       Public      Desktop     Downloads         Movies      Pictures  
```
You can change the working directory at any time using the `cd` command.

```
cd /usr/bin 
pwd 
ls
```
Now change back to your home again

```
cd ~
```

Tip: `~` is a shortcut for the HOME directory for any user. My home is `/users/karthik` and I can get there two ways:

`cd /users/karthik` OR `cd ~`.

**Full versus relative paths**

In the command line you can use both full paths (much like someone's street address with post code) OR offer relative directions from one's current location. You can do the same here.

```
cd /usr
pwd
```

We're now in the `usr/` directory. Now change to `bin`/

```
cd bin
```

This is the same as doing:

```
cd /usr/bin
```

from **anywhere**.

Tip: Just plain `cd` with no options should take you back home. Try it. `cd` to some place else and type in `cd` again.


**Quick exercise**

1. Change into your HOME (or wherever you downloaded the boot-camps folder). Then into `boot-camps`. Then `shell`. List the contents of this folder. Then change back into your home again.

---
## Exploring your file system

Three really important commands

* `ls`
* `file`
* `less`

`ls` is extremely useful both for beginners and experts. `ls` can not only list the current directory contents but also contents from anywhere without changing working directories.

e.g.

```
ls /usr
```

or even multiple directories at once

```
ls ~ /usr
``

Now we can start adding more options. Recall that commands can take both options (with a `-`) followed by arguments. Let's add some to ls. 

```
cd boot-camps
ls -l
karthik:boot-camps/ (2013-09-bristol*) $ ls -l
total 13944
-rwxr-xr-x   1 karthik  staff  3772928 Sep 12 05:16 Conclusion.ppt
-rwxr-xr-x   1 karthik  staff     3827 Sep  8 15:48 HintsAndTips.md
-rwxr-xr-x   1 karthik  staff     1062 Sep  8 15:48 LICENSE.md
-rw-r--r--   1 karthik  staff     1972 Sep  8 15:48 Location.md
drwxr-xr-x  19 karthik  staff      646 Sep 11 21:07 Python
-rwxr-xr-x   1 karthik  staff     3128 Sep  8 15:48 README.md
-rw-r--r--   1 karthik  staff    11255 Sep  8 15:48 Setup.md
-rw-r--r--   1 karthik  staff  3331584 Sep  8 15:48 Welcome.ppt
drwxr-xr-x   7 karthik  staff      238 Sep  8 15:48 setup
drwxr-xr-x   8 karthik  staff      272 Sep 12 06:00 shell
drwxr-xr-x   5 karthik  staff      170 Sep  8 15:48 testing
drwxr-xr-x   7 karthik  staff      238 Sep 12 05:15 version-control
```

By adding “-l” to the command, we changed the output to the long format.

Now let's add more options

```
ls -lt
```

The `t` options now sorts by time.

Similarly you can try the following:

| Option | What it does |
| ------ | ----------- |
| -a | List all files even those that are hidden. Files starting with a `.` are considered hidden |
| -d | Only directories | 
| -F | All a trailing slash to help identify folders | 
| -h | Make file sizes human readable | 
| -l | Long format | 
| -S | Sort by file size | 
| -t | Sort by modification time | 

Try some of these. Do you see any new files that we have not discussed before? You can even combine several of these options in a single command.

What are all the extra fields in the long format?

```
$ ls -l
total 13944
-rwxr-xr-x   1 karthik  staff  3772928 Sep 12 05:16 Conclusion.ppt
-rwxr-xr-x   1 karthik  staff     3827 Sep  8 15:48 HintsAndTips.md
-rwxr-xr-x   1 karthik  staff     1062 Sep  8 15:48 LICENSE.md
-rw-r--r--   1 karthik  staff     1972 Sep  8 15:48 Location.md
drwxr-xr-x  19 karthik  staff      646 Sep 11 21:07 Python
```
* Files begin with a `-` and directories with a `d`.  
* Followed by permissions for the user, group, and everyone.   
* Permissions are in the order of read, write, and execute. If any * group is missing a permission, you'll see a `-`.  
* Ignore field 3 for now (it's the number of links to the file)  
* The owner of the file  
* What group this person belongs to  
* Size of file in bytes  (Quick question: How do you change this?)  
* Date an time the file was last modified  
* Name of file.  

**Determining file type**

```
file <filename>
```

e.g.

```
file Location.md
Location.md: ASCII English text
```

Examine files with the `less` command. Keeps the content from scrolling of the screen. You can also use the arrow keys to navigate up or down. Press enter to keep scrolling down and the `q` key to quit. 

**Quick exercise**

1. cd `HOME`
2. cd into `boot-camps` 
3. cd into a given directory
4. List directory contents with `ls -l`
5. Pick any file that looks interesting to you
6. find out what it is using `file`
7. then view it's contents using `less`

---


**Tab completion**

Bash and most other shell programs have tab completion. This means that you can begin typing in a command name or file name and just hit tab to complete entering the text. If there are multiple matches, the shell will show you all available options.

```
cd
cd bo<tab>
```

What just happened?

Try pressing `s`, then hitting tab?

When you hit the first `tab`, nothing happens. The reason is that there are multiple directories in the home directory which start with `s`. Thus, the shell does not know which one to fill in. When you hit `tab` again, the shell will list the possible choices.

Tab completion can also fill in the names of programs. For example, enter `e<tab><tab>`. You will see the name of every program that starts with an e. One of those is echo. If you enter `ech<tab>`` you will see that `tab` completion works.

**Command history**

The shell typically stores your most recent commands. View them by using the up and down arrow keys.  Typing in `history` is a great way to get a full list. Execute any (especially long commands) by referencing it with a `!`.

example: `!500` will execute line 500 from your history.


**Wildcards**

One of the biggest reasons using shell is faster than ever using a GUI file manager is that it allows for wildcards. There are special characters known as wildcards. They allow you to select files based on patterns of characters.

| Wildcard | what it means |
| -------  | -----------  | 
| * |  Matches any character | 
| ? |  Matches any single character | 
| [characters] |  Matches any character in this set | 
| ![characters] |  Matches any character NOT in this set |





### Short exercise 

Do each of the following using a single `ls` command without navigating to a different directory.

List all of the files in `/bin` that start with the letter `c`  
List all of the files in `/bin` that contain the letter `a`  
List all of the files in `/bin` that end in `o`  
BONUS: List all of the files in `/bin` that contain the letter `a`   or the letter `t`  


## Manipulating the file system

Make directories with `mkdir`

```
mkdir directory_name
```

Make as many as you like in a single call.

```
mkdir directory_name_1 directory_name_2 directory_name_3
```

Copy files with `cp`

```
cp file1 file2
```

Move files with `mv`

```
mv file1 file2
```

See the `man` command to get help on options you can use with these commands.

Remove files with `rm`

***Warning: The shell does not have a recycling bin. So any file removed with `rm` is gone forever. Use with caution.***

## Let's try out some of the commands above

First go home `cd`.  
Next create a temporary directory.


```
mkdir scratchpad
cd scratchpad
```

Make a few directories

```
mkdir dir1 dir2 dir3  
cp /etc/passwd .
```

```
ls -l
```

# Getting help

* If you're not sure where a program is located, use `which`

e.g.

```
which nano
``` 

* See the manual for any program with `man <program_name>`  
* Get help on a command with `<program_name> --help`  







