 -------------------------------

###Local Repo
Git will keep track of the _changes_ to your files, rather than keep multiple copies of the files. It saves the first version, then keeps track of subsequent changes to that version. This makes it efficient and speedy. It can recreate any version (go back in time) by adding up all the changes to get to where you want to be.

To create your own local repository, you need to first _initialize_ a repository with the infrastructure that git needs to keep track of things for you.

#####Exercise
Open the terminal window

Configure git

```bash
	git config --global user.name "User Name"
	git config --global user.email "user@email.com"
```

Selecting the default editor
* Windows
    * To use notepad++ on windows(thanks to http://stackoverflow.com/questions/1634161/how-do-i-use-notepad-or-other-with-msysgit/2486342#2486342 for that tip) `git config --global core.editor "'C:/Program Files (x86)/Notepad++/notepad++.exe' -multiInst -notabbar -nosession -noPlugin"`
    * To use notepad `git config --global core.editor "notepad.exe"`
* Mac
    * On Mac, with TextWrangler if you installed TextWrangler's command line tools then you should have an "edit" command `git config --global core.editor "edit -w"`
* Linux
    * `git config --global core.editor <your favorite editor here (nano, emacs, etc etc)>`


let's make a new repository:

`mkdir practice_git`

Navigate to the top-level directory where your code lives

and type… `git init`

It looks like nothing happened, but you have just created a repository. Now, everything inside this folder will be watched and recorded by git.

__ADVICE:__ This only has to be done once, and only for the top-level directory. DO NOT use this command multiple times within a directory or on the same folder. Things get weird.

We can see that git is watching this directory now by looking at the hidden files.

`ls -a`

look for `.git`

We don't really need to worry about what's in the `.git` directory.
But we can look inside this new directory using `cd .git` and `ls -a`, just to see what's there.

Use what you've learned. You may have noticed the file called description. 
You can describe your repository by opening the description file and replacing the text with a name for the repository. 
Mine will be called "Reproducible Science". You may call yours anything you like.
Let's not worry about the other files here. As a novice user, there's no reason to mess with them.

Navigate back up to your controlled directory ``cd ..``

On the white board draw a box representing the working area and explain that this is where you work and make changes.

#####Adding files to staging area

Create a new file `science.txt`

Now, lets add files that are inside:

On the white board draw a box representing the staging area (index) and explain that this is where we set up the next snapshot of our project.

Like a photographer in a studio, we're putting together a shot before we actually snap the picture.
Connect the working area box and the staging box with 'git add'.

`git add .` This adds __all__ the files in our repository.

But sometimes we only want to add a single file at a time.

Add your new file using `git add science.txt`

----------------------------------------
 
#####check the status of your repository

The files you've created on your machine are your local "working" copy. 
The changes your make in this local copy aren't backed up online automatically. 
Until you commit them, the changes you make are local changes. 
When you change anything, your set of files becomes different from the files in the official repository copy. 
To find out what's different about them in the terminal, try:

`git status`

This will tell you what files are currently being staged 

-- highlight the "Untracked files" section and that git tells you how to add a file to the next commit.


#####Commit changes to repository

But there's one more step! We need to commit the changes. 
Tell git _"Hey, we want you to remember the way that the files look right now"_.

On the white board draw a box representing the project history. 
Once we take a snapshot of the project that snapshot becomes a permanent reference point in the project's history that we can always go back to.
The history is like a photo album of changes, and each snapshot has a time stamp, the name of the photographer, and a description.

Connect the staging area to the history with `git commit -m "message"`.

In order to save a snapshot of the current state (revision) of the repository, we use the commit command. 
This command is always associated with a message describing the changes since the last commit and indicating their purpose. 

Git will ask you to add a commit message. 
This is just to remind you what changes you made. 

Informative commit messages will serve you well someday, so make a habit of never committing changes without at least a full sentence description.

__ADVICE: Commit often__

In the same way that it is wise to often save a document that you are working on, so too is it wise to save numerous revisions of your code. 
More frequent commits increase the granularity of your undo button.

__ADVICE: Good commit messages__

[because it's important!](http://www.commitlogsfromlastnight.com/)

There are no hard and fast rules, but good commits are atomic: they are the smallest change that remain meaningful. 
A good commit message usually contains a one-line description followed by a longer explanation if necessary.
For code, it's useful to commit changes that can be reviewed by someone in under an hour. 
Or it can be useful to commit changes that "go together" - for example, one paragraph of a manuscript, or each new function added to your script.

For example, if you work on your code all day long (add 200 lines of code, including 5 new functions and write 7 pages of your new manuscript including deleting an old paragraph), and at 3:00 you make a fatal error or deletion, but you didn't commit once, then you will have a hard time recreating the version you are looking for - because it doesn't exist!

#####Exercise

Step 1: Commit the file that you added to your repository

`git commit -m "message"`

__ADVICE:__ You must have a commit message. It's good practice and git won't let you commit without one. 

If you only want to add one file, use `git commit filename.txt -m "message"
`git commit -am "message` will add ALL tracked files.

--------------------
You can see all of the changes you have ever made to your repository by typing `git log`

#####Looking at your history
__git log__ will tell you a unique number for each commit (this will be important later), who made the changes (helpful with multiple authors), when the changes were made, and your helpful message about what you did.

git log -- Shows description of what we've done.

Explain the commit hash.

* Unique identifier of that commit if we ever want to refer to it.

Useful `git log` flags:
* -3
* --oneline
* --stat
* --since=X.minutes/hours/days/weeks/months/years or YY-MM-DD-HH:MM
* --author=<pattern>

#####differences
Git is aware of all the changes that you have made in your repository. If you change a file, you can ask git to tell you what changed (you went to lunch and forgot what you did).

Git diff will output the changes in your working directory that are not yet staged for a commit. To see how this works, make a change in your readme.rst file, but don't yet commit it.

#####Exercise
Make a small change to one of your text files (add a line or capitalize a word)

`git status` will tell us that something has changed in the local directory. Note that our changes are not staged or committed yet!

`git diff` will tell us which files changed and what was inserted and deleted. 

`git diff --stat` gives us a summary of the filename and number of insertions/deletions

`git diff -- filename` will roll back changes that are staged for the specified file, to how it looked in the most recent commit

Look at the differences. What do the messages tell you?

Stage and commit your new changes.

Look at the new changes in your log. Where are the most recent changes? The oldest?

-------------------------------------

So far, this seems like a lot of work. Why are we keeping track of all these little things??

Let's say you fatally ruin a file during an editing mistake (like when I deleted an awesome paragraph from my dissertation instead of cutting and pasting it like I meant to. Or wrote myself a mean frustrating note!) Maybe you even accidentally delete an important file (This code is old, why should I keep it?). 

If you have version control, you don't need to track down your System Administrator. You can fix your problem easily!

__NOTE:__ if you remove a file from git, you need to tell it that it will no longer be tracking that file. Use `git rm filename`

#####Exercise
Fatally change a file and delete another one

Use `ls` to see what files you have left

Use `git status` and `git diff` to see what you have done.

Stage your changes `git add filename`. Oops! We didn't actually want to do that.

Unstage your files using `git reset filename`. Now you can go back and edit them.

To discard all your most recent changes and GO BACK IN TIME, first look at your `git log` to decide what version you want to go back to. Remember the first 5-7 digits in the commit code of the version that wasn't screwed up.

use `git reset --hard versioncode`

To roll back to a specific file, use 
`git checkout version name --filename`

To roll back one version (usually I know that I messed up pretty quickly)
`git checkout master~1 path_to_file`

__NOTE:__ A hard reset is permanant. 

#####Exercise
* Create 5 files in your directory with 1-2 lines of content
* Commit files to the repository
* Change 2 of the 5 files and commit them
* Undo the changes in step 3
* Print out the last entry in the log
    * `git log -1`
* Delete a file
* Checkout that file from your last commit to get it back
    
__ADVICE__: Similar to the rationale behind frequently saving a document that you are working on, it is wise to commit revisions
to your files often. More frequent commits increase the granularity of your "undo" button.

-----------------------------
#####Branches
Branches are parallel instances of a repository that can be edited and version controlled in parallel. They are useful for pursuing various implementations experimentally or maintaining a stable core while developing separate sections of a code base.

Without an argument, the branch command lists the branches that exist in your repository.

`git branch`

So far, we just have one branch named "master"

The master branch is created when the repository is initialized. With an argument, the branch command creates a new branch with the given name.

`git branch experimental`

`git branch`

To delete a branch, use the -d flag.

`git branch -d experimental`

`git branch`

When might working on a branch be helpful?

##### Git checkout: switching between branches and abandoning local changes

The `git checkout` command allows context switching between branches as well as abandoning local changes.

To switch between branches, try:

`git branch newbranch`

`git checkout newbranch`

`git branch`

Git shows us in green with a star which branch we are on.

Working in the branch could be a safe way to make some pretty radical changes to our code without changing waht we know works. We can later merge these changes to our master. 

#####Exercise : Create and Merge Branches

Step 1 : Create two new branches and list them

`git branch first`

`git branch second`

Step 2 : Make changes in each new branch and commit them.

`git checkout first`

create a new txt file

`git add newfile`

`git commit newfile`

`git checkout second`

create second new txt file

`git add newfile2`

`git commit newfile2`

Merge the two branches into the core



Step 3 : Merge the two branches into the core

`git checkout first`

`git merge second -m "message"`

If you look at the files in first, the one you created from second is in there.

`git checkout master`

`git merge first -m "Message"`

Now look at the files and all from second and first should be in there.



