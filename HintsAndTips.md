# Hints and tips

## Getting help about a command - `man` pages

Type,

    $ man COMMAND

For example,

    $ man ls
    $ man grep
    $ man git

The up and down arrows on your keyboard allow you to scroll up and down the page.

To exit from the `man` page, press `q`.

## Auto-completion

Bash shells support what is called *tab completion*. If you type part of a command or file name, then press TAB, it will show the commands or files that match that prefix of the command. For example, if you type

    $ gre

then, without pressing ENTER, you press TAB you should see something like

    $ gre
    greadelf  grep           grepjar
    grefer    grep-changelog

which are all the commands that match `gre`. This also works for files. For example, if you type

    $ ls dat

then, without pressing ENTER, you press TAB you should see a list of all the files that start with `dat`. If there is only one, then the file name will be *auto-completed* e.g.

    $ ls data.txt

## Command history, or avoid having to retype commands

Bash shells support a command history - they record every command you type in. If you use up-arrow, or `CTRL-P`, at the command prompt you can scroll back to the previous command you executed. Down-arrow, or `CTRL-N`, allows you to scroll to the next command.

You can then move the cursor and edit the command. To move the cursor in the line you can do,

* Left, `CTRL-B`
* Right, `CTRL-F`
* To the start of the line, `CTRL-A`
* To the end of the line, `CTRL-E`

If you enter,

    $ history

You'll see a list of the commands in the history. Each has a number. If you enter,

    $ !NNN

where `NNN` is the number of a command in the history, that command will be rerun.

To search for a command that you've run before you can do,

    $ history | grep "COMMAND"

For example,

    $ history | grep "ls"

## A quick-start guide to the `nano` editor

nano is a simple text editor for Linux/Unix.

To start,

    $ nano file.txt

To move the cursor,

* Left, `CTRL-B`
* Right, `CTRL-F`
* Up a line, `CTRL-P`
* Down a line, `CTRL-N`
* To the start of the line, `CTRL-A`
* To the end of the line, `CTRL-E`

To delete and undelete a line,

* Delete the current line, `CTRL-K`
* Undelete the recently deleted lines, `CTRL-U`

To save a file, `CTRL-O`. You will be given the opportunity to edit the file name to save the file under a different name.

To quit, `CTRL-X`. If the file has unsaved changes you'll be given the chance to save them now.

## I don't recognise my prompt...where am I?

`GNU nano` is at the top of your window? You're in `nano`.

Your window looks like the following? You're in `vi`.

    ~
    ~
    ~    
    ~
    "somefilename" ...


Something like `somefilename  (Fundamental) ----` is at the bottom of your window? You're in `emacs` or `xemacs`.

`XEmacs: somefilename  (Fundamental) ----` is at the bottom of your window? You're in `xemacs`.

Your prompt is `>>>` ? You're in `python`.

Your prompt is `In [123]:` ? You're in `ipython`.

`--More--(...%)` is at the bottom-left of your window? You're in `more`.

`Manual page...` is at the bottom-left of your window? You're in a `man` page.

`:` is at the bottom-left of your window? You're in `less` or a `man` page.

## How do I exit from...

* `nano`, type `CTRL-X` to exit. If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` to quit without saving.
* `vi`, type `:q!` to exit without saving. If the text just appears on screen then press `ESC` then type `:q!`
* `emacs` or `xemacs`, `CTRL-X CTRL-S`. If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` then type `yes` to quit without saving.
* `python`, type `exit()` or `CTRL-D`
* `ipython`, type `exit()`, or `CTRL-D` then press `y`
* `man` page, press `q`
* `more` or `less`, press `q`
