# git/GitHub

The goal of this lesson is to introduce the students to [git][] via
collaboration on [GitHub][].

## Introduction

- Say some introductory stuff about version control in general, and git/GitHub
  in particular.

*Note: If you have nowhere to draw the figures in the
[Visual Git Reference][visual git] can be a good stand-in.*

## Signup and Setup

- Have everyone make [GitHub][] accounts.
- Have everyone configure git:

        $ git config --global user.name "User Name"
        $ git config --global user.email "user@email.com"
        $ git config --global core.editor "nano"
        $ git config --global color.ui "auto"

## Basics

- Make a new local repo.
- On the white board draw a box representing the working area and explain that
  this is where you work on make changes.
- Make a new file called `top-things-to-learn.md` and put the title
  "Top Things We Want to Learn" in it.
- `git status` -- highlight the "Untracked files" section and that git tells
  you how to add a file to the next commit.
- On the white board draw a box representing the staging area (index) and
  explain that this is where we set up the next snapshot of our project.
    - Like a photographer in a studio, we're putting together a shot before
      we actually snap the picture.
    - Connect the working area box and the staging box with `git add`.
- Run `git add top-things-to-learn.md`.
- `git status` -- highlight the "Changes to be committed" section and the
  and git telling you how to unstage the new file.
- On the white board draw a box representing the project history. Once we take
  a snapshot of the project that snapshot becomes a permanent reference point
  in the project's history that we can always go back to.
    - The history is like a photo album of changes, and each snapshot has a
      time stamp, the name of the photographer, and a description.
    - Connect the staging area to the history with `git commit`.
- Run `git commit` and explain log messages.
    - Summary message at the top, longer one below.
- `git status` -- nothing going on!
- `git log` -- Shows description of what we've done.
    - `git log --oneline` -- Abbreviated version.
- Explain the commit hash.
    - Unique identifier of that commit if we ever want to refer to it.
    - Comes from "hashing" stuff associated with the commit, like the changes
      we made.
    - Can demo hashing with Python's `hashlib.sha1`.

[git]: http://git-scm.com/
[GitHub]: http://github.com
[visual git]: http://marklodato.github.com/visual-git-guide/index-en.html
