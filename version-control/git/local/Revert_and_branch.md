[Up To Schedule](../../../README.md) - Back To [Make Incremental Changes I](../../../version-control/git/local/Readme.md)  - Forward To [Plan for Mistakes](../../../python/testing/Readme.md)

----

# Make Incremental Changes II: Reverting and branching (and a bit more) in git

**Based on materials by Katy Huff, Anthony Scopatz, Joshua R. Smith, Sri
Hari Krishna Narayanan, and Matthew Gidden**


## Refresher on the basics of git

We use git to keep track of changes to the files in a particular
directory. Here are the basic commands, discussed
[previously](Readme.md).

- `git init`: Initialize a directory as a git repository.

- `git status`: Check the current status of things.

- `git add`: Add a new file to the repository, or stage the changes
  to a file.
  
- `git commit`: Commit the changes (adding, modifying, or removing
  files).
  
- `git diff`: Study the changes.

- `git log`: Summarize the history of changes.

### ![Exercise](pics/exercise.jpg) Exercise: Refresh your understanding of git
 
**Step 1**: Go back to your `~/simplestats` repository. Make a change to the
`README.md` file, or create a new file.

**Step 2**: Add and commit your changes.

**Step 3**: Study some of those differences as well as the repository log.


## Unstaging a staged file: `git reset`

There are a number of ways that you may accidentally stage a file that
you don't want to commit.  Use `git reset` to unstage a file.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git reset`

**Step 1**: Make a change to the `README.md` file and stage the change.

```
$ nano README.md
$ git add README.md
```

**Step 2**: Check the status of the repository, to see that the file
  has been staged.

```
$ git status
```

**Step 3**: Unstage the file with `git reset`; `HEAD` refers to the
  most recent commit to the repository.

```
$ git reset HEAD README.md
```

**Step 4**: Check the status again.

```
$ git status
```

## Discarding unstaged modifications: `git checkout`

If you've made changes to a file and want to just scrap those changes
and go back to the last committed version of the file, use `git
checkout`.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git checkout`

**Step 1**: Check the status of the repository, and look at your
  unstaged changes.
  
```
$ git status
$ git diff
```

**Step 2**: Discard the changes.

```
$ git checkout README.md
```

**Step 3**: Look at the status of things again.

```
$ git status
$ git diff
```

## Removing files: `git rm`

If you want to remove a file from your repository, use `git rm`.
Note that past versions of the file will remain in the repository history.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git rm`

**Step 1**: Create a file, and add and commit it to the repository.

```
$ nano READYOU.md
$ git add READYOU.md
$ git commit -m "Add READYOU.md file"
```

**Step 2**: Remove the file from the repository.

```
$ git rm READYOU.md
```

**Step 3**: Check the status.

```
$ git status
```

**Step 4**: Commit the change (removing the file).

```
$ git commit -m "Remove READYOU.md"
```

**Step 5**: Check the status again.

```
$ git status
```

What happens if you delete a file in the shell without `git rm`? Try deleting
`README.md` with `rm` rather than `git rm`.

```
$ rm README.md
```

What does `git status` say?  Oops! How can you recover this important
file?

```
$ git checkout README.md
```

Note that, just as you should use `git rm` rather than `rm` for
removing files, you should use `git mv` rather than `mv` for moving or
renaming files.


## Don't include _everything_ in the repository: the `.gitignore` file.

You probably don't want to include _every_ file in your project
directory as part of your git repository.

- Backup files automatically created by your editor (`*.bak` or `*~`)
- Very large primary data files that don't change
- Compiled code (`*.pyc`, `*.o`, `*.so`, `*.exe`)
- Files that are derived from your code (for example,
  figures/graphs/images)

If you never add them to your repository (with `git add`), then they
won't be tracked, but they'll show up in the output from `git status`,
which can be a bother.

To tell git to ignore a set of files (and so not mention them
in the status output), create a `.gitignore` file in the root of your
project directory. This should be a plain text file with file or
directory names; you can also use wildcards, like `*.bak` or `*.o`.
If you include a directory name, all files in that directory will be ignored.

### ![Exercise](pics/exercise.jpg) Exercise: Create a `.gitignore` file

**Step 1**: Create a subdirectory `Data`. Put a few data files there,
  or use `touch` to create a few files there.

**Step 2**: Use `git status`.

**Step 3**: Create a `.gitignore` file to tell git to ignore those
  data files.

**Step 4**: Use `git status` again.

**Step 5**: Add and commit the `.gitignore` file



## `git revert`: the promised "undo" button

It is possible that after many commits, you decide that you really
want to "roll back" a set of commits and start over.  It is easy to
revert your code to a previous version.

You can use `git log` and `git diff` to explore your history and
determine which version you are interested in.  Choose a version and
note the *hash* for that version. (Let's assume it's `abc456`).

**NOTE**: the version you choose will be the changes you wish to
remove, not the final point that you want to reach.  In this case, you
will remove the changes made in `abc456`, rather than "rolling back"
to `abc456`.

     git revert abc456

**Importantly**, this will not erase the intervening commits.  This
will create a new commit that is changed from the previous commit by a
change that will recreate the desired version.  This retains a
complete provenance of your software, and can be compared to the
prohibition in removing pages from a lab notebook.

### ![Exercise](pics/exercise.jpg) Exercise: Practice using `git revert`

1. Create 5 files in your directory with one line of content in each
   file.
2. Commit the files to the repository.
3. Change 2 of the 5 files and commit them.
4. Undo the changes in step 3.
5. Print out the last entry in the log.



## `git branch`: Listing, Creating, and Deleting Branches

Branches are parallel instances of a repository that can be edited and
version-controlled in parallel. They are useful for experimenting with
different ideas, for maintaining a stable core while testing
developing new features.

Without an argument, `git branch` lists the branches that
exist in your repository.

    $ git branch
    * master

The master branch is created when the repository is initialized. With an
argument, `git branch` creates a new branch with the given
name.

    $ git branch experimental
    $ git branch
    * master
      experimental

To delete a branch, use the `-d` flag.

    $ git branch -d experimental
    $ git branch
    * master

## `git checkout`: Switching Between Branches

We had previously used `git checkout` to abandon local changes. It is
also used to switch between branches.

### ![Exercise](pics/exercise.jpg) Exercise: Switch to a new branch

Create an `add_median` branch and switch to it.

    $ git branch add_median
    $ git checkout add_median
    $ git branch

How can you tell we've switched between branches? When we used the
branch command before there was an asterisk next to the master branch.
The asterisk indicates which branch you are currently in.

### ![Exercise](pics/exercise.jpg) Exercise: Copy files into your repository

Let's make sure we have a good copy of `stats.py`.

```
$ cd ~/simplestats
$ cp ~/boot-camps/python/testing/stats.py .
```

Now let's add it to our repository, but in the current branch.

```
$ git add stats.py
$ git commit -m "Adding a first version of stats.py."
```

### ![Exercise](pics/exercise.jpg) Exercise: Implement the `median()` function.

1. Use the list's `sort()` method to implement the `median()` function. (The
median is either the middle value of an odd set of numbers *or* the `mean()` of
the middle two values of an even set of numbers.)
2. Commit the changed file to your repository.

## `git merge`: Merging Branches

At some point, the `add_median` branch may be ready to become part of
the `master` branch.  In real life, we might do a lot more testing and
development.  For now, let's assume that our median function is ready
and merge this back into the master branch.  We use `git merge` to
combine the changes in two parallel branches.

```
$ git checkout master
$ git merge add_median
```

## Aside: Make your Prompt Pretty

In the next section, we'll get into the gritty details of remotes and branches
as we head toward web-based storage of your repositories. It turns out that some
folks have created a way to make this kind of navigation more convenient,
showing you what branch you're on using your bash prompt. Some super nice
properties also include color-coding when you've got changed files or when your
branch is fresh.

### ![Exercise](pics/exercise.jpg) Exercise: Update your prompt

**Step 1**: Copy the following lines into your `~/.bashrc` file, or
`~/.bash_profile` in Mac OS X (taken from a
combination of [two](http://stackoverflow.com/a/6086978)
[sources](https://gist.github.com/woods/31967)).

```
function color_my_prompt {
    local __user_and_host="\[\033[01;32m\]\u@\h"
    local __cur_location="\[\033[01;34m\]\w"
    local __git_branch='`git branch 2> /dev/null | grep -e ^* | sed -E  s/^\\\\\*\ \(.+\)$/\(\\\\\1\)\ /`'
    local __prompt_tail="\[\033[35m\]$"
    local __last_color="\[\033[00m\]"

    RED="\[\033[0;31m\]"
    YELLOW="\[\033[0;33m\]"
    GREEN="\[\033[0;32m\]"

    # Capture the output of the "git status" command.                                                                                               
    git_status="$(git status 2> /dev/null)"

    # Set color based on clean/staged/dirty.                                                                                                           
    if [[ ${git_status} =~ "working directory clean" ]]; then
        state="${GREEN}"
    elif [[ ${git_status} =~ "Changes to be committed" ]]; then
        state="${YELLOW}"
    else
        state="${RED}"
    fi

    export PS1="$__user_and_host $__cur_location ${state}$__git_branch$__prompt_tail$__last_color "
}

# Tell bash to execute this function just before displaying its prompt.                                                                              
PROMPT_COMMAND=color_my_prompt
```

**Step 2**: Source your `.bashrc` file.
(Use `~/.bash_profile` instead of `~/.bashrc` in Mac OS X)

    $ source ~/.bashrc

**Step 3**: Play around with it.

## Resources

* [git book](http://git-scm.com/book)
* [git game](http://pcottle.github.io/learnGitBranching/index.html)

----

[Up To Schedule](../../../README.md) - Back To [Don't Repeat Yourself (or Others)](../../../python/best_practice/dont_repeat_yourself.md)  - Forward To [Plan for Mistakes](../../../python/testing/Readme.md)
