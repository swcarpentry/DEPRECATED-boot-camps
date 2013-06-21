#Teaching Git

##What is Version Control?
__any guesses?__

* keeps track of changes to a file/folder
* keeps a backup of changing files
* stores of history of the changes
* allows many people to make changes concurrently

__Do you currently use some kind of version control?__

* [The many-folder system](http://hginit.com/i/01-copies.png)
* [Final.doc](http://www.phdcomics.com/comics/archive.php?comicid=1531)
* [cup of USBs](http://www.eeweb.com/rtz/version-control)
* Dropbox
  * only keeps 30-day history
  * longer with "Packrat" but lacks concurrent change ability
* Email files to yourself
* Time Machine (or other automated backup system)


##Why should I use a version control system?
__i.e., what makes version control better than Dropbox?__

* keeping lots of copies of files clogs up your system and can be confusing
* sync copies of a project across computers
* easily share project data/code/text with collaborators
* easily pubish/share data/code/text with scientific manuscripts
* __ability to collaboratively code and merge multiple changes__
* __ability to revert to earlier copies when you (or someone else) messes up__

##What is Git?
__there are lots of version control software options...__

* mercurial (hg)
* bazaar (bzr)
* subversion (svn)
* concurent versions system (cvs)
* git (git)

we will be learning git.

Git is:

* distributed version control system (don't worry about the details)
* designed to be fast and efficient
* Has an awesome logo [octocat](http://octodex.github.com/images/original.jpg)
* FREE

Git can set up a repository on your local computer - no connection to a big server or online is needed to make it work!

###GitHub

GitHub is where code is collected online.

Show a project on GitHub (further justification for learning).

GitHub connects your local repositories to the online world and allows your project to get __social__

You can __clone__ or __fork__ other people's repositories (more on that later).

Easily see changes that have been made (and when) and who made the changes

__citable!__

__C.V./resume line!__

#####github.com?

GitHub is a site where many people store their open (and closed) source code repositories. It provides tools for browsing, collaborating on and documenting code. Your home institution may have a repository hosting system of it's own. To find out, ask your system administrator. GitHub, much like other forge hosting services (launchpad, bitbucket, googlecode, sourceforge etc.) provides :

* landing page support
* wiki support
* network graphs and time histories of commits
* code browser with syntax highlighting
* issue (ticket) tracking
* user downloads
* varying permissions for various groups of users
* commit triggered mailing lists
* other service hooks (twitter, etc.)


###Local Repo
Git will keep track of the _changes_ to your files, rather than keep multiple copies of the files. It saves the first version, then keeps track of subsequent changes to that version. This makes it efficient and speedy. It can recreate any version (go back in time) by adding up all the changes to get to where you want to be.

To create your own local repository, you need to first _initialize_ a repository with the infrastructure that git needs to keep track of things for you.

#####Exercise
Open the terminal window

Configure git

	$ git config --global user.name "User Name"
	$ git config --global user.email "user@email.com"
	$ git config --global core.editor "nano"
	$ git config --global color.ui "auto"

Navigate to the top-level directory where your code lives

and type… `git init`

It looks like nothing happened, but you have just created a repository. Now, everything inside this folder will be watched and recorded by git.

__ADVICE:__ This only has to be done once, and only for the top-level directory. DO NOT use this command multiple times within a directory or on the same folder. Things get weird.

We can see that git is watching this directory now by looking at the hidden files.

`ls -a`

look for .git

We can look inside this new directory using `cd .git` and `ls -a`

Use what you've learned. You may have noticed the file called description. You can describe your repository by opening the description file and replacing the text with a name for the repository. Mine will be called "Reproducible Science". You may call yours anything you like.
Let's not worry about the other files here. As a novice user, there's no reason to mess with them.

Navigate back up to your controlled directory ``cd ..``

On the white board draw a box representing the working area and explain that this is where you work and make changes.

#####Adding files to staging area

Now, lets add files that are inside:

On the white board draw a box representing the staging area (index) and explain that this is where we set up the next snapshot of our project.

Like a photographer in a studio, we're putting together a shot before we actually snap the picture.
Connect the working area box and the staging box with 'git add'.

`git add .` This adds __all__ the files in our repository.

But sometimes we only want to add a single file at a time.

Create a new file `science.txt`

Add it using `git add science.txt`

----------------------------------------
 
#####check the status of your repository

The files you've created on your machine are your local "working" copy. The changes your make in this local copy aren't backed up online automatically. Until you commit them, the changes you make are local changes. When you change anything, your set of files becomes different from the files in the official repository copy. To find out what's different about them in the terminal, try:

`git status`

This will tell you what files are currently being stages (draw diagram of local-->stage-->committed)

-- highlight the "Untracked files" section and that git tells you how to add a file to the next commit.

#####Commit changes to repository


But there's one more step! We need to commit the changes. Tell git _"Hey, we want you to remember the way that the files look right now"_.

On the white board draw a box representing the project history. Once we take a snapshot of the project that snapshot becomes a permanent reference point in the project's history that we can always go back to.
The history is like a photo album of changes, and each snapshot has a time stamp, the name of the photographer, and a description.

Connect the staging area to the history with `git commit`.

In order to save a snapshot of the current state (revision) of the repository, we use the commit command. This command is always associated with a message describing the changes since the last commit and indicating their purpose. 

Git will ask you to add a commit message. This is just to remind you what changes you made. 

Informative commit messages will serve you well someday, so make a habit of never committing changes without at least a full sentence description.

__ADVICE: Commit often__

In the same way that it is wise to often save a document that you are working on, so too is it wise to save numerous revisions of your code. More frequent commits increase the granularity of your undo button.

__ADVICE: Good commit messages__

[because it's important!](http://www.commitlogsfromlastnight.com/)

There are no hard and fast rules, but good commits are atomic: they are the smallest change that remain meaningful. A good commit message usually contains a one-line description followed by a longer explanation if necessary.

For example, if you work on your code all day long (add 200 lines of code, including 5 new functions and write 7 pages of your new manuscript including deleting an old paragraph), and at 3:00 you make a fatal error or deletion, but you didn't commit once, then you will have a hard time recreating the version you are looking for - because it doesn't exist!

#####Exercise

Step 1: Commit the file that you added to your repository

`git commit -m "message"

__ADVICE:__ You must have a commit message. It's good practice and git won't let you commit without one. 

If you only want to add one file, use `git commit filename.txt -m "message"

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

Look at the differences. What do the messages tell you?

Stage and commit your new changes.

Look at the new changes in your log. Where are the most recent changes? (top) The oldest? (bottom)

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

#####git clone : Copying a Repository

Yesterday, you checked out a git type repository at <https://github.com/swcarpentry/boot-camps/2013-06-wise-beginners>

`git clone -b 2013-06-wise-beginners git@github.com:swcarpentry/boot-camps.git`

When you clone the Original repository, the one that is created on your local machine is a copy, and will behave as a fully fledged local repository locally. However, with the right configuration, it will be able to pull changes from collaborators to your local machine and push your changes to the Original repository. We'll get to that soon, but for now, let's __fork__ the repository from GitHub.

#####Exercise : Cloning a Repository from GitHub

Step 1 : Pick any repository you like. There are many cool projects hosted on github. Take a few minutes here, and pick a piece of code.

I like: https://github.com/gemgon/Octocat.git

Step 2 : Clone it. If you didn't find anything cool, you can chose the "instructional" Spoon-Knife repository:

`git clone git@github.com/octocat/Spoon-Knife.git`

Step 3 : You should see many files download themselves onto your machine. Let's make sure it worked. Change directories to the source code and list the contents.

`cd Spoon-Knife`

`ls` 

`git pull` : Pulling updates from the Original Repository

Updating your repository is like voting. You should update early and often especially if you intend to contribute back to the upstream repository and particularly before you make or commit any changes. This will ensure you're working with the most up-to-date version of the repository. Updating won't overwrite any changes you've made locally without asking, so don't get nervous. When in doubt, update.

`git pull`

Already up-to-date.
Since we just pulled the repository down, we will be up to date unless there has been a commit by someone else to the Original repository in the meantime.

=========================

###Remote Repo

Advantages to using and participating on remote repositories:

####github password

Setting up github requires a github user name and password. Please take a moment to create a free github account (if you want to start paying, you can add that to your account some other day).

####git remote : Steps for Forking a Repository

A key step to interacting with an online repository that you have forked is adding the original as a remote repository. By adding the remote repository, you inform git of a new option for fetching updates and pushing commits.

The git remote command allows you to add, name, rename, list, and delete repositories such as the original one upstream from your fork, others that may be parallel to your fork, and so on.

Exercise : Fork Our GitHub Repository

While you probably already have a copy of the SWC-bootcamp repository, GitHub doesn't know about it yet. You'll need to tell github you want to have an official fork of this repository.

Step 1 : Go to our repository () from your browser, and click on the Fork button. Choose to fork it to your username rather than any organizations.

Step 2 : Clone it. From your terminal :

`git clone https://github.com/YOU/boot-camps.git`

`cd boot-camps`

Step 3 :

`git remote add upstream https://github.com/USERNAME/boot-camps.git`

`git remote -v`

All repositories that are clones begin with a remote called origin.

#####git fetch : Fetching the contents of a remote

Now that you have alerted your repository to the presence of others, it is able to pull in updates from those repositories. In this case, if you want your master branch to track updates in the original SWC-bootcamp repository, you simply `git fetch` that repository into the master branch of your current repository.

The fetch command alone merely pulls down information recent changes from the original master (upstream) repository. By itself, the fetch command does not change your local working copy. To update your local working copy to include recent changes in the original (upstream) repository, it is necessary to also merge.

#####git merge : Merging the contents of a remote

To incorporate upstream changes from the original master repository (in this case USERNAME/boot-camps) into your local working copy, you must both fetch and merge. The process of merging may result in conflicts, so pay attention. This is where version control is both at its most powerful and its most complicated.

#####Exercise : Fetch and Merge the Contents of Our GitHub Repository

Step 1 : Fetch the recent remote repository history

`git fetch upstream`

Step 2 : Make certain you are in the YYYY-MM-PLACE branch and merge the upstream YYYY-MM-PLACE branch into your YYYY-MM-PLACE branch

`git checkout 2013-06-wise-beginners`

`git branch` will show you what branch you are on, to confirm

`git merge upstream\2013-06-wise-beginners`

Step 3 : Check out what happened by browsing the directory.

#####git pull : Pull = Fetch + Merge

The command `git pull` is the same as executing `git fetch` followed by `git merge`. Though it is not recommened for cases in which there are many branches to consider, the pull command is shorter and simpler than fetching and merging as it automates the branch matching. Specificially, to perform the same task as we did in the previous exercise, the pull command would be :

`git pull upstream`

Already up-to-date.

When there have been remote changes, the pull will apply those changes to your local branch, unless there are conflicts with your local changes.

#####git push : Sending Your Commits to Remote Repositories

The `git push` command pushes commits in a local working copy to a remote repository. The syntax is `git push [remote] [local branch]`. Before pushing, a developer should always pull (or fetch + merge), so that there is an opportunity to resolve conflicts before pushing to the remote.

#####Exercise : Push a change to github

We'll talk about conflicts later, but first, since we have no conflicts and are up to date, we can make a minor change and send our changes to your fork, the "origin."

`git push origin 2013-06-wise-beginners`

If you have permission to push to the upstream repository, sending commits to that remote is exactly analagous.

`git push upstream 2013-06-wise-beginners`

In the case of the YYYY-MM-PLACE code, new developer accounts will not allow this push to succeed. You're welcome to try it though.

#####git merge : Conflicts

This is the trickiest part of version control, so let's take it very carefully.

In the YYYY-MM-PLACE code, you'll find a file called `Readme.md`. This is a standard documentation file that appears rendered on the landing page for the repository in github. To see the rendered version, visit your fork on github, (https://github.com/YOU/boot-camps/tree/YYYY-MM-PLACE/README.md).

For illustration, let's imagine that, suddenly, each of the developers on the YYYY-MM-PLACE code would like to welcome visitors in a language other than English. Since we're all from so many different places and speak so many languages, there will certainly be disagreements about what to say instead of "Welcome."

I, for example, speak fluent LOLcat speak, so I'll push (to the upstream repository) my own version of Welcome on line 5 of Readme.md. 

You may speak another language, perhaps even English, however, and may want to replace the Tamil word 'O Hai' with an equivalent word that you prefer (welcome, willkommen, bienvenido, benvenuti, etc.).

You'll want to start a new branch for development. It's a good convention to think of your master branch (in this case your YYYY-MM-PLACE branch) as the "production branch," typically by keeping that branch clean of your local edits until they are ready for release. Developers typically use the master branch of their local fork to track other developers changes in the remote repository until their own local development branch changes are ready for production.

#####Exercise : Experience a Conflict

__Step 1 :__ Make a new branch, edit the readme file in that branch, and commit your changes.

`git branch development`

`git checkout development`

Switched to branch 'development'

`nano Readme.md`

`git commit -m "Changed the welcome message to ... `

__Step 2 :__ Mirror the remote upstream repository in your master branch (in this case your YYYY-MM-PLACE branch) by pulling down my changes

`git checkout 2013-06-wise-beginners`

Switched to branch '2013-06-wise-beginners'

`git fetch upstream`

`git merge upstream/2013-06-wise-beginners`

	Updating 43844ea..3b36a87
	Fast-forward
	README.rst |   2 +-
	1 files changed, 1 insertions(+), 1 deletions(-)

__Step 3 :__ You want to push it to the internet eventually, so you pull updates from the upstream repository, but will experience a conflict.

`git merge development`

	Auto-merging Readme.md
	CONFLICT (content): Merge conflict in Readme.md
	Automatic merge failed; fix conflicts and then commit 	the result.

#####git resolve : Resolving Conflicts

Now what?

Git has paused the merge. You can see this with the git status command.

The only thing that has changed is the `Readme.md` file. Opening it, you'll see something like this at the beginning of the file.

	<<<<<<< HEAD
	O Hai Dere!!!1!
	================
	=======
	Bula!
	========
	>>>>>>> devo

The intent is for you to edit the file, knowing now that I wanted the Welcome to say O Hai. If you want it to say Willkommen, you should delete the other lines. However, if you want to be inclusive, you may want to change it to read O Hai and Willkommen. Decisions such as this one must be made by a human, and why conflict resolution is not handled more automatically by the version control system.

""O Hai and Bula!!1!""

This results in a status To alert git that you have made appropriate alterations,

`git add Readme.md`

`git commit -am "message"`

`git push origin 2013-06-wise-beginners`


=========================================

###Repos with multiple authors
#####Make and Clone


* Explain that much like a browser navigates to websites using a URL, git talks to remote repositories using a URL.
* Explain the different URL options:
  * Read/write ssh protocol requires ssh keys, which make it so you don't have to enter username/password.
  * https protocol takes username/password.
  * git protocol is read-only.
* We can `Fork` repositories by clicking on the button
* We can `Watch` repositories
  * options to get or ignore messages when the repo changes or is updated
  * We can use `Settings` to 
    * delete a repository
    * rename a repository
    * transfer ownership, change permissions, make public/private
* We can look at our [network](https://github.com/blog/39-say-hello-to-the-network-graph-visualizer)

#####Make a new repository on github
* Make a new demo repo on GitHub explaining the process as you go (probably on your personal account).
  * Click on green `new repository` button on your account
  * Have GitHub put in a `README` so it can be cloned.

* Now we want to get a copy of this down on all of our computers

#####Exercise
 `git clone`
 
* Have everyone do this via the https URL.
* `ls` to show the new directory and cd into it.
* Compare the local listing to what you see on GitHub. Same for now, but not for long!

#####Local Basics

__IMPORTANT:__ Make sure you tell people not to make their own local changes, that will make things really complicated later when people pull. Alternatively, you can go ahead and let them do whatever they like and use it as a teaching moment on `git reset --hard` in a few minutes when it's time to start the collaboration.


Make a new file called `to-learn.md` and put the title "Top Things We Want to Learn" in it.

`git add to-learn.md`

`git status`

`git commit to-learn.md`

`git log`


#####Exercise - Previewing Changes
The file we're making is going to be a list of the top things everyone wants to learn in the boot camp.
Add your item (e.g. everyone's names) and save.

`git status` -- point out that now we have a modified file instead of an untracked file, but the procedure for adding it to the next snapshot is the same.

Want to see the changes you're about to add? `git diff`!

* `git add`
* `git diff` -- now it doesn't show anything. git diff shows differences between the working area and the staging area.
* `git commit -m`

#####Remotes

As we said back at the beginning, git uses URLs to point repositories on other computers, in this case GitHub's servers.

We can give these remote repositories names so that we don't have to type in the full URL all the time, and in fact git has already set one up for us.

* `git remote` -- there's a remote called "origin".
* `git remote -v` -- we can see that it points to the same URL we cloned from, git sets this up automatically.

#####Branches

On the GitHub view of the repo highlight the branches pull-down -- right now it only has one branch called "master", this is another thing git makes for us by default.

What branch are we on locally? git branch.

Give a short explanation of branches and mention that we will come back to them later.

* Isolated development environments.

When git communicates with a remote repository it needs to know what branch is changing, in this case we'll just stick with "master".

#####Pushing

Use push command to send data to a remote repository, and we also have to specify the remote name and branch name: 

`git push origin master`

Refresh the GitHub view.

#####Pulling

__IMPORTANT:__ If students have been making local commits, this is the time at which they will need to use `git reset --hard` to get back in sync with you.

`pull` is the reciprocal command, must specify remote and branch.

Have everyone `git pull origin master`

#####Collaborate

Pick a student to add their top thing to learn to the list:

* Add them to the collaborator list on the demo repo.
  * Settings --> Collaborators
edit, save, `add`, `commit`, `push`
Have everyone `pull` changes from the origin.

#####Rebase
__No Conflict__

* Have another student add their thing and `push`.
* Make a change to the README file before pulling.
* Try to push.
* On the white board draw the situation: my repo and the remote repo have different development histories and git doesn't know how to pull things together.
* It would be nice if I could move my local change after the remote change. (Draw picture.) There's a command for that!
* `git fetch origin` -- This gets information from the remote repository but doesn't integrate it with your local repo like pull does.
* `git rebase origin/master` -- `origin/master` is how we specify the fetched data for the remote named "origin" and it's branch named "master".
  * This replays my local changes on top of the state of the remote repo.
* `git log --oneline` -- Now my change is after the remote change.
* `git push origin master`
* Have everyone `pull`.

__With Conflicts__

* Again, have a student add their thing and `push`.
* Before pulling make a change in the same place in the same file.
* Try to `rebase` as above.
* Explain the conflict message git prints out.
* Show the conflict messages in the file and how to clean it up.
* Continue the `rebase` and `push` the result.
* Have everyone `pull`.

#####Developing in Branches

Often you want to leave a stable version of your code alone while you make some potentially disruptive changes. 
Or you and someone else are working on the code and you want to be able to work without worrying what others are doing.

It's going to take a long time to get everyone's top thing to learn onto the list one at a time, so the room is going to break out into groups and each come up with their own list.

So that they can all do this and then push their result to GitHub each is going to work in their own, isolated branch.

####Working on branches

* Make a new branch: `git branch sarah-list`.
* `git branch` -- highlight the asterisk showing the branch we're currently on.
* Move to that branch: `git checkout sarah-list`.
* `git branch` -- asterisk moved!
* Make a change (`add`, `commit`) and `git push origin sarah-list`.
  * __IMPORTANT:__ have to specify new branch named when pushing, not "master".

`git checkout master` -- show that your change is not in master.

Show how to browse to the other branch on GitHub.

#####Exercise
Have each group/table pick a unique branch name, switch to that branch, and add all of their top things to learn to the list and push the result.

#####Resolving the Changes

* Browse all of the new branches on GitHub.
* Illustrate the situation on the board
* Could just pick one of these branches as the new one "master" and move on, but we're more adventurous than that.
* Make sure you're in "master".
* `git fetch origin` -- without a branch name it grabs all of the new branches.
* Pick a branch and `git merge branch-name -m "message"`.
  * Should be a smooth fast-forward.
* Pick another branch and try to merge.
  * Resolve conflicts, add, and commit.
* Repeat as necessary.
* Push the final result to GitHub.
  * `git push origin master`

###Show some scientific repositories for work flow

####Review workflow on whiteboard again

### Add to concept map of skills learned

