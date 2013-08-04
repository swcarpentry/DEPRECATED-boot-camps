# Git hints and tips

## `man` page

Like many Unix/Linux commands, `git` has a `man` page,

    $ man git

You can scroll the manual page up and down using the up and down arrows.

You can search for keywords by typing `/` followed by the search term e.g. if interested in help, type `/help` and then hit enter.

To exit the manual page, type `q`.

## Command-line help

Type,

    $ git --help

and Git gives a list of commands it is able to help with, as well
as their descriptions. 

You can get more help on a specific command, by providing the command name e.g.

    $ git init --help
    $ git commit --help

## Add a repository description

You can edit the file `.git/description` and give your repository a name e.g. "My first repository".

## Ignore scratch, temporary and binary files

You can create a `.gitignore` file which lists the patterns of files you want Git to ignore. It's common practice to not add to a repository any file you can automatically create in some way e.g. C object files (`.o`), Java class (`.class`) files or temporary files e.g. XEmacs scratch files (`~`). Adding these to `.gitignore` means Git won't complain about them being untracked.

Create or edit `gitignore`,

    $ nano .gitignore

Then add patterns for the files you want to ignore, where `*` is a wildcard,

    *~
    *.o
    *.so
    *.dll
    *.exe
    *.class
    *.jar

Then, add `.gitignore` to your repository,

    $ git add .gitignore
    $ git commit -m "Added rules to ignore XEmacs scratch files and binary files"

## Add colour to `diff`

    $ git config --global color.diff auto

Previous: [Conclusions and further information](4_Conclusion.md)
