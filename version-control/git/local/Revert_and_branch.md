[Up To Schedule](../../../README.md) - Back To [Planning for Mistakes](../../../python/testing) - Forward To [Collaboration](../remote)

# Make Incremental Changes II: Reverting and branching in git

----

**Based on materials by Katy Huff, Anthony Scopatz, Joshua R. Smith, Sri 
Hari Krishna Narayanan, and Matthew Gidden**

    
## Refresher on the basics of git

We use git to keep track of changes to the files in a particular
directory. Here are the basic commands, discussed
[yesterday](Readme.md).

- `git init`: Initialize a directory as a git repository.

- `git status`: Check the current status of things.

- `git add`: Add a new file to the repository, or stage the changes
  to a file.
  
- `git rm`: Remove a file from the repository.

- `git commit`: Commit the changes (adding, modifying, or removing
  files).
  
- `git diff`: Study the changes.

- `git log`: Summarize the history of changes.

### ![Exercise](pics/exercise.jpg) Exercise: Refresh your understanding of git
 
**Step 1**: Go back to your `~/simplestats` repository. Make a change to the
`README.md` file, or create a new file.

**Step 2**: Add and commit your changes.

**Step 3**: Study some of those differences as well as the repository log.


## `git revert`: the promised "undo" button

It is possible that after many commits, you decide that you really
want to "roll back" a set of commits and start over.  It is easy to
revert your code to a previous version.

You can use `git log` and `git diff` to explore your history and
determine which version you are interested in.  Choose a version and
note the *hash* for that version. (Let's assume it's `abc456`).  **NOTE**:
the version you choose will be the changes you wish to remove, not the
final point that you want to reach.  In this case, you will remove the
changes made in `abc456`, rather than "rolling back" to `abc456`.

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

Create an `add_var` branch and switch to it.

    $ git branch add_var
    $ git checkout add_var
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

### ![Exercise](pics/exercise.jpg) Exercise: Implement the `var()` function.

1. Use the list's `sort()` method to implement the `median()` function. (The
median is either the middle value of an even set of numbers *or* the `mean()` of
the middle two values in an odd set of numbers.)
2. Commit the changed file to your repository.

## `git merge`: Merging Branches

At some point, the `add_var` branch may be ready to become part of
the `master` branch.  In real life, we might do a lot more testing and
development.  For now, let's assume that our variance function is ready
and merge this back into the master branch.  We use `git merge` to
combine the changes in two parallel branches.

```
$ git checkout master
$ git merge add_var
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

[Up To Schedule](../../../README.md) - Back To [Planning for Mistakes](../../../python/testing) - Forward To [Collaboration](../remote)
