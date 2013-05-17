# git/GitHub

The goal of this lesson is to introduce the students to [git][] via
collaboration on [GitHub][].

Yesterday, we created a collection of scripts, datafiles and output files in RSTudio. These are in a GitHub repo, but we would like to clean them up, reformat the directory structure and add README and LICENSE files. The parts about using the shell to do the cleanup are not detailed here (only the git parts).  

## Introduction

- Say some introductory stuff about version control in general, and git/GitHub
  in particular.
- Use the whiteboard to diagram the working directory + staging area + local repo + remote repo

## Setup and Signup

- Have everyone configure git:

        $ git config --global user.name "User Name"
        $ git config --global user.email "user@email.com"
        $ git config --global core.editor "nano"
        $ git config --global color.ui "true"

- Give a little tour of [GitHub][].
- Have everyone make [GitHub][] accounts.

### Make and Clone

- Have everyone fork the 2013-05-nescent repo on [GitHub][]
- via terminal / git-bash, 'git-clone' the repo to their own machines
    - Have everyone do this via the `https` URL to avoid having to set up ssh keys.
- Explain the different URL options:
    - Read/write `ssh` protocol requires [ssh keys][], which make it so you
      do not have to enter username/password.
    - `https` protocol takes username/password.
    - `git` protocol is read-only.
- `ls` to show the new directory and `cd` into it.
- Compare the local listing to what you see on [GitHub][]. Same for now, but
  not for long!

## Basics

### Local Basics

- Review the white board
- Make some local changes: we will re-organize the mess of R files from yesterday and add a README and LICENSE
- use `git status` to see changes

### Composing the Snapshot

- On the white board draw a box representing the staging area (index) and
  explain that this is where we set up the next snapshot of our project.
    - Connect the working area box and the staging box with `git add`.
- Want to see the changes you're about to add? `git diff`!
- Run `git add` for new files
- `git diff` -- now it doesn't show anything. `git diff` shows differences
- `git status` -- highlight the "Changes to be committed" section
  and git telling you how to unstage the new file.

### Taking the Snapshot

- On the white board draw a box representing the project history. Once we take
  a snapshot of the project that snapshot becomes a permanent reference point
  in the project's history that we can always go back to.
    - The history is like a photo album of changes, and each snapshot has a
      time stamp, the name of the photographer, and a description.
    - Connect the staging area to the history with `git commit`.
- Run `git commit` and explain log messages.
    - Summary message at the top, longer one below.
- `git status` -- nothing going on!

### Looking at the History

- `git log` -- Shows description of what we've done.
    - `git log --oneline` -- Abbreviated version.
- Explain the commit hash.
    - Unique identifier of that commit if we ever want to refer to it.
    - Comes from "hashing" stuff associated with the commit, like the changes
      we made.

## Sharing

- Now I want to share my changes with everyone so they can start working on
  it too.

### Remotes

- As we said back at the beginning, git uses URLs to point repositories on other
  computers, in this case [GitHub's][GitHub] servers.
- We can give these remote repositories names so that we don't have to type
  in the full URL all the time, and in fact git has already set one up for us.
- `git remote` -- there's a remote called "origin".
- `git remote -v` -- we can see that it points to the same URL we cloned from,
  git sets this up automatically.

### Branches

- On the [GitHub][] view of the repo highlight the branches pull-down -- right
  now it only has one branch called "master", this is another thing git makes
  for us by default.
- What branch are we on locally? `git branch`.
- Give a short explanation of branches and mention that we will come back to
  them later.
    - Isolated development environments.
- When git communicates with a remote repository it needs to know what branch
  is changing, in this case we'll just stick with "master".

### Pushing

- Use `push` command to send data to a remote repository, and we also have to
  specify the remote name and branch name: `git push origin master`.
- Refresh the [GitHub][] view.

### Collaborate

- Put the students into teams of two and have them initiate pull requests against each other's repos

## Developing in Branches

Often you want to leave a stable version of your code alone while you make some
potentially disruptive changes. Or you and someone else are working on the code
and you want to be able to work without worrying what others are doing.

- It's going to take a long time to get everyone's top thing to learn onto the
  list one at a time, so the room is going to break out into groups and each
  come up with their own list.
- So that they can all do this and then push their result to [GitHub][] each
  is going to work in their own, isolated branch.

### Making a New Branch

*Note: The [Learn Git Branching][] app can be a useful way to
illustrate this part of the lesson.*

- Make a new branch: `git branch matt-list`.
- `git branch` -- highlight the asterisk showing the branch we're currently on.
- Move to that branch: `git checkout matt-list`.
- `git branch` -- asterisk moved!
- Make a change and push.
    - **IMPORTANT:** have to specify new branch named when pushing, not "master".
- `git checkout master` -- show that your change is *not* in master.
- Show how to browse to the other branch on [GitHub][].
- Have each group pick a unique branch name, switch to that branch, and add
  all of their top things to learn to the list and push the result.

### Resolving the Changes

- Browse all of the new branches on [GitHub][].
- Illustrate the situation on the [Learn Git Branching][] app.
- Could just pick one of these branches as the new one "master" and move on,
  but we're more adventurous than that.
- Make sure you're in "master".
- `git fetch origin` -- without a branch name it grabs all of the new branches.
- Pick a branch and `git merge branch-name`.
    - Should be a smooth fast-forward.
    - Illustrate on the [Learn Git Branching][] app.
- Pick another branch and try to merge.
    - Resolve conflicts, add, and commit.
    - Illustrate on the [Learn Git Branching][] app.
- Repeat as necessary.
- Push the final result to [GitHub][].

[git]: http://git-scm.com/
[GitHub]: http://github.com
[ssh keys]: https://help.github.com/articles/generating-ssh-keys
[visual git]: http://marklodato.github.com/visual-git-guide/index-en.html
[Learn Git Branching]: http://pcottle.github.com/learnGitBranching/?NODEMO
