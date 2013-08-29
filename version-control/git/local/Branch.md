[Up To Schedule](../../../README.md) - Back To [Planning for Mistakes](../../../python/testing) - Forward To [Collaboration](../remote)

# Branching in Version Control
----

**Based on materials by Katy Huff, Anthony Scopatz, Joshua R. Smith, Sri 
Hari Krishna Narayanan, and Matthew Gidden**

    
## git branch : Listing, Creating, and Deleting Branches

Branches are parallel instances of a repository that can be edited and
version controlled in parallel. They are useful for pursuing various
implementations experimentally or maintaining a stable core while
developing separate sections of a code base.

Without an argument, the **branch** command lists the branches that
exist in your repository.

    $ git branch
    * master

The master branch is created when the repository is initialized. With an
argument, the **branch** command creates a new branch with the given
name.

    $ git branch experimental
    $ git branch
    * master
      experimental

To delete a branch, use the **-d** flag.

    $ git branch -d experimental
    $ git branch
    * master

## git checkout : Switching Between Branches, Abandoning Local Changes

The **git checkout** command allows context switching between branches
as well as abandoning local changes.

To switch between branches, try

    $ git branch add_stats
    $ git checkout add_stats
    $ git branch

How can you tell we've switched between branches? When we used the
branch command before there was an asterisk next to the master branch.
That's because the asterisk indicates which branch you're currently in.

### Exercise : Copy files into your repo

Let's make sure we have a good copy of `stats.py` and `test_stats.py`.

```
$ cd ~/simplestats
$ cp ~/boot-camps/python/testing/stats.py .
$ cp ~/boot-camps/python/testing/test_stats.py .
```

Now let's add them to our repo, but in the current branch.

```
$ git add *stats.py
$ git commit -m "Adding a first version of the files for mean."
```

### Exercise : Add an additional test for std() and commit the changes.

1. Write an additional test for std().  *(Ask us for a tip if necessary)*
2. Improve std() to pass this test.
3. Commit the changed files to your repo.

## git merge : Merging Branches

At some point, the `add_stats` branch may be ready to become part of
the `master` branch.  In real life, we might do a lot more testing and
development.  For now, let's assume that our mean function is ready
and merge this back to the master.  One method for combining the
changes in two parallel branches is the **merge** command.

```
$ git checkout master
$ git merge add_stats
```

## Aside: Make your Prompt Pretty

In the next section, we'll get into the gritty details of remotes and branches
as we head toward web-based storage of your repositories. It turns out that some
folks have created a way to make this kind of navigation more convenient,
showing you what branch you're on using your bash prompt. Some super nice
properties also include color-coding when you've got changed files or when your
branch is fresh.

### Exercise : Update your prompt

Step 1 : Copy the following lines into your ~/.bashrc file (taken from a
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

Step 2 : Source your bashrc (it'll change immediately)

    $ source ~/.bashrc

Step 3 : Play around with it.

## Resources

[git book](http://git-scm.com/book)

[Up To Schedule](../../../README.md) - Back To [Planning for Mistakes](../../../python/testing) - Forward To [Collaboration](../remote)
