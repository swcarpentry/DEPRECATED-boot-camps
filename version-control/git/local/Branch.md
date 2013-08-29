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

    $ git branch newbranch 
    $ git checkout newbranch 
    $ git branch

How can you tell we've switched between branches? When we used the
branch command before there was an asterisk next to the master branch.
That's because the asterisk indicates which branch you're currently in.

## git merge : Merging Branches

At some point, the experimental branch may be ready to become part of
the core or two testing branches may be ready to be combined for further
integration testing. The method for combining the changes in two
parallel branches is the **merge** command.

### Exercise : Create and Merge Branches

Step 1 : Create two new branches and list them

    $ git branch first
    $ git branch second

Step 2 : Make changes in each new branch and commit them.

    $ git checkout first
    Switched to branch 'first'
    $ touch firstnewfile
    $ git add firstnewfile
    $ git commit -am "Added firstnewfile to the first branch."
    [first 68eba44] Added firstnewfile to first branch.
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 firstnewfile
    $ git checkout second
    Switched to branch 'second'
    $ touch secondnewfile
    $ git add secondnewfile
    $ git commit -am "Added secondnewfile to the second branch."
    [second 45dd34c] Added secondnewfile to the second branch.
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 secondnewfile

Step 3 : Merge the two branches into the master branch

    $ git checkout first
    Switched to branch 'first'
    $ git merge second
    Merge made by recursive.
     0 files changed, 0 insertions(+), 0 deletions(-)
      create mode 100644 secondnewfile
    $ git checkout master
    Switched to branch 'master'
    $ git merge first
    Updating 1863aef..ce7e4b5
    Fast-forward
     0 files changed, 0 insertions(+), 0 deletions(-)
     create mode 100644 firstnewfile
     create mode 100644 secondnewfile

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
