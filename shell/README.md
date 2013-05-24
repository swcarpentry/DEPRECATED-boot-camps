# Shell - crib sheet

* Even in a windows, Windows, world useful to know.
* HPC resources and cloud platform- and infrastructure-as-a-service.
* More than just file and directory management.
* Bolt together programs into powerful data processing pipelines.
* Automation.
* Bash, "Bourne again shell"

## `man` page

    $ man COMMAND

Up and down arrows to scroll.

`/` followed by search term (e.g. `/help`) then ENTER.

`q` to exit.

## `--help`

Command syntax, usage and other information.

    $ COMMAND --help

Top tip: if writing your own executables, be consistent, and provide `--help`.

## Comments

    # This is a comment. It is not executed.

## `who` and `whoami`

    $ who # who is logged on
    $ whoami # who am I logged on as

## Directories

    $ pwd       # Path to current directory (folder in Windows)
    $ ls        # List directory
    $ ls *.txt  # Wild-card
    $ ls *_hai*
    $ ls -R     # Recurse
    $ ls -F     # Append / to directories
    $ ls -l     # Permissions, date, size, owner, group...
    $ cd /      # Root directory
    $ ~         # Home directory. There's no place like ~
    $ ls -a     # Hidden files
    $ ls .      # Current directory
    $ ls ..     # Parent directory
    $ cd ..
    $ cd        # Default to home directory
    $ mkdir directory
    $ mkdir ~/directory
    $ mv directory another_directory
    $ rmdir empty_directory

# Files

    $ cat file           # View file
    $ less file          # Page through file
    $ more file          # Page through file
    $ head -2 file       # First N lines
    $ tail -3 file       # Last N lines
    $ cp file1 file2     # Copy
    $ cp *.txt directory
    $ rm file.txt        # Delete - no recycle bin.
    $ rm -r directory    # Recurse
    $ rm -rf directory   # Recurse and force - beware

## History

Up arrow browses previous commands

    $ history
    $ !NNNN   # Rerun Nth command in history

## Word count

A filter

    $ wc file 
    $ wc -l file  # Lines only
    $ wc -w file  # Words only
    $ wc -l *.txt       # Total

Use to find out number of records in a data file if one record per line.

## Regular expressions

    $ ls *.txt        # Zero or more characters
    $ ls ?o*          # Exactly one character
    $ ls a[bcde]*.txt # Exactly one of the characters listed
    $ ls a[cde]*.txt  
    $ ls *.*
    $  *.[!txt]*      # No sequence involving t, x or t

## Searching within files

Global/regular expression/print

    $ grep the haiku.txt
    $ grep day haiku.txt
    $ grep is haiku.txt
    $ grep 'it is' haiku.txt
    $ grep -w is haiku.txt   # Exact match
    $ grep -n it haiku.txt   # lines with matches
    $ grep -i the haiku.txt  # Ignore case
    $ grep -wn is haiku.txt
    $ grep -wnv is haiku.txt # Non-matching lines
    $ grep  -wnr Today *     # Recurse

## Input and output redirection

`>` redirects output (AKA standard output)

    $ grep -r not * > found_nots.txt
    $ cat found_nots.txt
    $ ls *.txt > txt_files.txt
    $ cat txt_files.txt
    $ cat                            # Echoes standard input
    Blah
    CTRL-D
    $ cat > myscript.txt
    Blah
    CTRL-D
    $ cat myscript.txt

`<` redirects input (AKA standard input).

    $ cat < haiku.txt
    $ ls idontexist.txt > output.txt
    $ cat output.txt

Error message is output on standard error.

    $ ls idontexist.txt 2> output.txt               # 2 is standard error
    $ ls haiku.txt 1> output.txt                    # 1 is standard output
    $ ls idontexist.txt haiku.txt > output.txt 2>&1

## Exercise - grep

`pdb/` contains a set of protein database files

Each `.pdb` file lists atoms in a protein

Write a single command that

* Uses `grep` to find all  hydrogen (`H`) atoms in all `.pdb` files.
* Stores these in `hydrogen.txt`.

You will need wild card, exact matches output redirection

## Searching for files

    $ find .                # Find all
    $ find . -type d        # Directories only
    $ find . -type f        # Files only
    $ find . -maxdepth 2    # Maximum depth of tree
    $ find . -mindepth 3    # Minimum depth of tree
    $ find . -name *.txt    # Fails as wild-card is expanded
    $ find . -name '*.txt'  # Name or pattern matching
    $ find . -iname '*.TXT' # Ignore case
    $ find . -empty         # Empty files only
    $ touch emptyfile.txt   # Create empty file
    $ find . -empty

`` back-ticks execute a command

    $ wc -l `find . -name '*.txt'`

## Exercise - find

Write a single command that

* Uses `find` to find all `.pdb` files.
* Uses `cat` to list their contents.
* Stores contents in `proteins.txt`.

You will need back ticks, find file name option, output redirection.

## Power of the pipe

Count text files

    $ find . -name '*.txt' > files.tmp
    $ wc -l files.tmp

`find` outputs a list of files, `wc` inputs a list of files, skip the temporary file.

    $ find . -name '*.txt' | wc -l

`|` is a pipe.

    $ echo "Number of .txt files:" ; find . -name '*.txt' | wc -l

`;` equivalent to running two commands on separate lines.

Question: what does this do?

    $ ls | grep s | wc -l

Answer: counts the number of files with `s` in their name.

    $ history | grep 'wget'

Power of well-defined modular components with well-defined interfaces,

 * Bolt together to create powerful computational and data processing workflows.
 * Good design principle applicable to programming - Python modules, C libraries, Java classes - modularity and reuse.
 * "little pieces loosely joined".

## Exercise - pipes

Write a single command that

* Uses `find` to find all `.pdb` files.
* Uses `cat` to list their contents.
* Uses `grep` to find all the hydrogen (`H`) atoms in their contents.
* Uses `wc` to count the number of hydrogen atoms found.
* Stores the count in `hydrogen_count.txt`.

You will need commands from previous exercises, back ticks, pipe.

## Variables

    $ set                            # See all variables
    $ MYFILE=data.txt
    $ echo $MYFILE
    $ echo "My file name is $MYFILE"
    $ bash                           # Spawn new shell
    $ echo $MYFILE
    CTRL-D
    $ export MYFILE                  # Export to new shells
    $ bash
    $ echo $MYFILE
    CTRL-D
    $ echo $PATH
    $ export PATH=$PATH:/path/to/bin # Common requirement

    $ let NUM=$NUM+1                 # Simple arithmetic
    $ TEXT_FILES=`ls *.txt`          # Save output in variable
    $ echo TEXT_FILES

`.bashrc` to define variables and other actions to do when logging in e.g.

    export JAVA_HOME=/opt/local/java1.5
    export PATH=$JAVA_HOME/bin:$PATH
    export ANT_HOME=/home/michelj/Software/apache-ant-1.7.0
    export PATH=$ANT_HOME/bin:$PATH

## Conditionals

    $ NUM=1
    $ if [ "$NUM" -eq 1 ]; then echo "Equal"; fi
    $ WORD="hello"
    $ if [ "$WORD" = "hello" ];  then echo "The same"; fi

## Loops.

    $ for i in `cat file`; do echo $i; done | sort | uniq

    $ for PDB in `find . -name '*.pdb'`; do
        echo $PDB
     done

## Shell scripts - automate

    $ nano protein_filter.sh
    #!/bin/bash
    DATE=`date`
    echo "Processing date: $DATE"
    for PDB in `find . -name '*.pdb'`; do
        echo $PDB
    done
    echo "Processing completed!"
    
    $ sh protein_filter.sh
    $ chmod +x protein_filter.sh # Mark as executable
    $ ./protein_filter.sh

## Exercise - shell scripts

Edit `protein_filter.sh` so that it

* Has variables `ATOM` with value `"H"` and `PDB_EXT` with value `"pdb"`.
* For each `.pdb` file it prints the file name and a count of the number of hydrogen atoms in that file.

You will need parts of your comm and from the previous exercise.

* Edit protein_filter.sh so that it takes the atom value from the command-line e.g.

    $ ./protein_filter.sh “H”
    $ ./protein_filter.sh “C”

`$1` provides access to the first command-line argument.

## Permissions

    $ ls -l # permissions, dates, sizes, owner, group, size in byte, creation/modification date/time, name.

Users, groups, others

Read, write, execute

    $ chmod a+r haiku.txt     # Add permission - all read
    $ chmod a-r haiku.txt     # Remove permission - all not read
    $ chmod u+r haiku.txt     # User read
    $ chmod g+w haiku.txt     # Group write
    $ chmod o+x haiku.txt     # Other execute
    $ chmod g+rx haiku.txt    # Group read and execute
    $ chmod go+rx haiku.txt   # Group and other read and execute
    $ chmod ugo=rwx haiku.txt # Set permission
 
## Jobs and processes

    $ ./counter.sh > output.txt
    CTRL^Z                        # Suspend
    $ wc -l output.txt
    $ jobs -l                     # Jobs, job number, process ID
    $ fg JOBNUMBER                # Resume in foreground
    CTRL^Z
    $ bg JOBNUMBER                # Resume in background
    $ wc -l output.txt            # Still working
    $ ./counter.sh > output.txt & # Start in background, job number, process ID
    $ kill PROCESSID                  
    $ jobs
    $ ps                          # Processes
    $ top                         # Resource consumption
    $ bash
    $ nohup ./counter.sh > output.txt & # Continue after log out
    $ nohup allows processes to continue even after the user logs out.
    CTRL-D
    $ wc -l output.txt

# Secure shell ***

    $ ssh username@boot-camp.software-carpentry.org
    $ ssh username@boot-camp.software-carpentry.org ls # Run remote command
    $ scp file.txt username@boot-camp.software-carpentry.org:
    $ scp username@boot-camp.software-carpentry.org:directory/file.txt . # Relative path
    $ scp -r username@boot-camp.software-carpentry.org:directory copy

# Packaging

    $ zip -r pdb.zip pdb  # Package and compress, recurse
    $ mkdir tmp
    $ cp pdb.zip tmp
    $ cd tmp
    $ unzip pdb.zip
   
    $ tar -cvf pdb.tar pdb # Create, list files, file archive
    $ mkdir tmp
    $ cp pdb.tar tmp
    $ tar -xf pdb.tar
    $ zip pdb.tar

Top tip: If preparing bundles, put your content in a directory then zip or tar up that single directory. It can be annoying if someone unzips or untars a bundle and it spews its contents all over their directory, possibly overwriting their files.

Top tip: If preparing bundles of your software put the version number or a date in the name. If someone asks for advice, you'll know exactly what version they have.

    $ ls -l
    $ gzip pdb.tar

    $ gunzip pdb.tar.fx

    $ md5sum pdb.zip

Top tip: when putting packages up for download also put up the file size and MD5 sum so people can check they've not been tampered with.

## `wget`

Download files via command-line

    $  wget --output-document=bbc.html http://news.bbc.co.uk

## script

    $ script
    $ ls -l
    $ CTRL-D
    $ cat typescript

Record commands typed, commands with lots of outputs, trial-and-error when building software. 

Send exact copy of command and error message to support.

Turn into blog or tutorial. 

## Summary

Shell and scripts

* Reduce errors by reusing tried-and-tested components.
* String together components into powerful computational and data processing pipelines.
* Avoid reinventing the wheel.
* Free up time to do research.

## Links

* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].

## Exercises

### grep

    grep -w 'H' *pdb > hydrogen.txt

### find

    cat `find . -name '*.pdb'` > proteins.txt

## Pipes

    cat `find . -name '*.pdb'` | grep -w 'H' | wc -l > hydrogen_count.txt

## Shell scripts

    DATE=`date`
    ATOM=$1
    PDB_EXT="pdb"

    echo "Processing date: $DATE"
    for PDB in `find . -name "*.$PDB_EXT"`; do
        COUNT=`grep -w $ATOM $PDB | wc -l`
        echo "$PDB $COUNT"
    done
    echo "Processing completed!"
