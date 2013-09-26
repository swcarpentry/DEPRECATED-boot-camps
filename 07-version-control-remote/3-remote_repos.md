#Teaching Git - 3

###Remote Repo

GitHub allows us to set up personal and public repositories. Repos are public by default. 

####github password

Setting up github requires a github user name and password. Please take a moment to create a free github account (if you want to start paying, you can add that to your account some other day).

####git remote : Steps for Forking a Repository

A key step to interacting with an online repository that you have forked is adding the original as a remote repository. By adding the remote repository, you inform git of a new option for fetching updates and pushing commits.

The git remote command allows you to add, name, rename, list, and delete repositories such as the original one upstream from your fork, others that may be parallel to your fork, and so on.

#####git clone : Copying a Repository

Yesterday, you checked out a git type repository at <https://github.com/swcarpentry/boot-camps/2013-09-msu>

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


Exercise : Fork Our GitHub Repository

While you probably already have a copy of the SWC-bootcamp repository, GitHub doesn't know about it yet. You'll need to tell github you want to have an official fork of this repository.

Step 1 : Go to our repository (https://github.com/swcarpentry/boot-camps/) from your browser, and click on the Fork button. Choose to fork it to your username rather than any organizations.

Step 2 : Clone it. From your terminal :

`git clone https://github.com/YOU/boot-camps.git`

`cd boot-camps`

Step 3 :

`git remote add upstream https://github.com/swcarpentry/boot-camps.git`

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

`git checkout 2013-09-msu`

`git branch` will show you what branch you are on, to confirm

`git merge upstream\2013-09-msu`

Step 3 : Check out what happened by browsing the directory.

#####git pull : Pull = Fetch + Merge

The command `git pull` is the same as executing `git fetch` followed by `git merge`. Though it is not recommened for cases in which there are many branches to consider, the pull command is shorter and simpler than fetching and merging as it automates the branch matching. Specifically, to perform the same task as we did in the previous exercise, the pull command would be :

`git pull upstream`

Already up-to-date.

When there have been remote changes, the pull will apply those changes to your local branch, unless there are conflicts with your local changes.

#####git push : Sending Your Commits to Remote Repositories

The `git push` command pushes commits in a local working copy to a remote repository. The syntax is `git push [remote] [local branch]`. Before pushing, a developer should always pull (or fetch + merge), so that there is an opportunity to resolve conflicts before pushing to the remote.

#####Exercise : Push a change to github

We'll talk about conflicts later, but first, since we have no conflicts and are up to date, we can make a minor change and send our changes to your fork, the "origin."

`git push origin 2013-09-msu`

In the case of the YYYY-MM-PLACE code, new developer accounts will not allow this push to succeed. You're welcome to try it though.

#####git merge : Conflicts

This is the trickiest part of version control, so let's take it very carefully.

In the YYYY-MM-PLACE code, you'll find a file called `Readme.md`. This is a standard documentation file that appears rendered on the landing page for the repository in github. To see the rendered version, visit your fork on github, (https://github.com/YOU/boot-camps/tree/2013-09-msu/README.md).

For illustration, let's imagine that, suddenly, each of the developers on the 2013-09-msu code would like to welcome visitors in a language other than English. Since we're all from so many different places and speak so many languages, there will certainly be disagreements about what to say instead of "Welcome."

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

`git checkout 2013-09-msu`

Switched to branch '2013-09-msu'

`git fetch upstream`

`git merge upstream/2013-09-msu`

	Updating #######..#######
	Fast-forward
	README.rst |   2 +-
	1 files changed, 1 insertions(+), 1 deletions(-)

__Step 3 :__ You want to push it to the internet eventually, so you pull updates from the upstream repository, but will experience a conflict.

`git merge development`

	Auto-merging Readme.md
	CONFLICT (content): Merge conflict in Readme.md
	Automatic merge failed; fix conflicts and then commit 	the result.

# Exercise

With a friend (or a small group) take turns:

* Select one person to be the 'upstream' person, the others are the 'forkers'
* Have the upstream make a repository on github
* Have the forkers fork the upstream's repository
* Everyone clones their own version of the repository
* Have forkers add upstream's repository as a remote
* Have upstream add a few lines to the README, commit and push the changes
* Have Forkers pull the changes, make changes of their own and commit them
* Have upstream make another change, commit and push
* Forkers pull this new change creating a conflict and resolve the conflict
* Trade off roles 