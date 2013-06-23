#Teaching Git - 3

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

####Show some scientific repositories for work flow

####Review workflow on whiteboard again

#### Add to concept map of skills learned

