# Introduction to the Shell

## What is a shell?

Interface to the *kernel*, which is responsible for managing the
low-level functionality of your computer's operating system.

A command-line interface works through a Read-Evaluate-Print Loop.

We call these tools REPLs

## Commands to try/explain

```
# the environment
whoami

whatamidoinghere?

printenv

# pager
printenv | less

# grep
printenv | grep PWD

# creating strings of data ""
echo "Hello, $USER"

# '' and the shell
echo 'Hello, $USER'

echo 'Hello,' $USER

#delimiter
echo 'Hello,'     $USER

echo $PWD

echo $PATH

which printenv

pwd

ls .

ls ./

ls /

ls -a

ls -F

ls -l

man ls

# google ls tutorial

mkdir falafel_bin

ls falafel_bin

touch falafel_1

touch falafel_2

mv falafel_* falafel_bin

mv falafel_bin/* .

ls f*

ls -d

ls -d f*

rm falafel_bin/falafel_1

rm falafel_bin/falafel_2

rm falafel_bin

yes

Control-C

ps

top

wget/curl

find . -print

find . -type f -print
find . -type f -name "*1*"

find . -type f -exec grep Volume {} \;

find . -type f -print | xargs grep Volume

### extra stuff
top
permissions (chmod, chown, chgrp, ls -l)
```
