
introduction to scripps swc

 - where are materials
 - online commenting
 - online editing
 - schedule etc

@@ put schedule into pages

Q: do we want to do a mid-class/end-class questionnaire? make cait write one
Q: how to handle * and wildcards?
Q: wget and unzip => VM?
 - where do we post these commands?

Q: save job control for later!
Q: save ssh for 2nd day; screen.
Q: link back to exercises
Q: go back and look at questionnaire

---

does everyone have access to a command line? vm or mac os x

http://software-carpentry.org/4_0/shell/intro/

describe shell:
 - typewriter interface that dates back to beginning of computer
 - think "printer" but not laser printer :)
 - I/O of letters, numbers, punctuaton
 - ascii art
 - "clui" vs wimp
 - commands.  very hacker-esque.
 - essentially a "repl" (read evaluate print loop)
 - many languages, including matlab and python, have repl
 - so does a class of computer interfaces called "shells"
 - we tend to just call it "the shell" because whichever specific one you use (we'll be using bash) they all pretty much work the same at a basic level
 - you tell the shell to do something
 - the shell does it, potentially running programs etc.
 - waits for more input

WHY WHY WHY???

it turns out to be really easy to write shell-compatible programs.  not so
much for wimp/gui/web.  so many scientists (esp) tend to write their data
analysis programs that way.  that's why you NEED to know it.

but, more generally, the shell is its own little powerful programming language.
and you can use it to do lots of things.  which we'll now walk you through.
mostly we'll be concerned with basic files etc. as well as *automation*.
Because automation is important for science.

---

file and directory structure:

One of the most important things computers do is store data.  The
part of the computer responsible for this is the file system, which
organizes your data.  You're all familiar with this via your mac
and windows and linux gui interface.  Note that your phones don't
really have a visible file system, just as a difference.

the basic organizational structure for files/directories is like a
tree: you have directory structure with a root directory that contains
everything else, and then a bunch of directories (potentially nested).
within these directories there are files.

for example ...

let's use the shell to look around

whoami

pwd tells us where we are
you should get an answer like XXX

/ organization; / as root, then within as a separator.  files and direcotires
look the same but only a directory can have a / after it, because only a
directory can contain something within it.

ls tells us what files and directories a directory holds.

ls -F annotates the output of ls

something.something

ls -F /some/directory

absolute vs relative paths; pwd etc.

ls defaults to current directory

cd changes what directory we're currently in

special directory named ..

why doesn't it show up with ls by default?

ls -a, and secret/hidden files (config, etc.)

unix vs windows, shells etc.

----

Creating and deleting directories

mkdir will create a new directory; same absolute and relative paths apply.

ls should show new directory; but new directory is empty.

not really -- ls -a shows .files

cd newdir
nano junk

ls shows we've now created a file called junk!

ls -l shows how big things are, along iwth a lot of other info.

rm junk
ls

create file again
move up dir
try to delete dir with rm
must use rmdir

remake dir
remake file
rename with 'mv'; mv works on files and directories

mv from to

no rename

cp copies a file without moving it

can do ls with multiple paths

make copies; delete; etc.

---

shell stuff -- molecules (where can I get this data?)

@@ do with FASTA files instead!!!!!

shell expansion of *, ? -- happens AT THE COMMAND LINE

(echo)

wc

which file is shortest?
easy with 6 - how about 6000?
doing things with 6000 is just as easy as with 6!

redirect: > 

most commands output text to >

will create the file if it doesn't exist; overwrite if it does.

no screen output!
ls confirms file exists
'cat' can print out contents

sort

sort to >

head

head -1

head -20

too ... many ... intermediate ... files!!

piping!

no explicit creation or naming of an intermediate file

so unix people create simple tools that take in text and output text and
work well together. (scientists should think about doing this too...)

wc -l *.pdb > lengths

concepts of pipes and filters

tail, split, cut, and uniq

*, >, |

---

@@@ blast and organizing data
@@@ is blast and clsutall on the vm? or just do outputs?

grep and find

@@@ the tao thing 

grep finds and prints lines in files that match a pattern

grep -w
grep -n
grep -i
grep -v
grep -l

find

find finds files, kind of - not just by name, but by property.

find . -type d
find . -type f

find -empty
find -mindepth

find chains things 

find *.txt find sthe wrong thing!  put in quotes -- test with echo.

backquotes

cut

@@ mass renaming of files

what do you do with non text files?

 1. extend these commands (grep etc.) to handle those files
 2. convert files to text
 3. use a programming language like Python

---

variables and automation
