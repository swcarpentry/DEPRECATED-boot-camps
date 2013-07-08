# Version Control and Git - crib sheet

[Git](http://git-scm.com/)

## Git command-line

    $ git --version
    $ man git
    $ git help
    $ git help --all
    $ git help checkout

## Working with a local repository

Repository for directories and files.

    $ mkdir papers
    $ cd papers
    $ git init

Working directory.

Git configuration files in `.git` directory.

    $ ls -A

### Git configuration

Who we are.

    $ git config --global user.name "Your Name"
    $ git config --global user.email "yourname@yourplace.org"

Editor.

    $ git config --global core.editor nano
    $ git config --global core.editor vi
    $ git config --global core.editor xemacs
    
Global configuration.

    $ cat ~/.gitconfig
    [user]
	name = Your Name
	email = yourname@yourplace.org
    [core]
	editor = nano
    $ git config -l

### Working with files in the repository

[Markdown](http://daringfireball.net/projects/markdown/syntax) file.

    $ nano rfc-template.md
    # Title

    by Authors

    ## Overview

    ## Key points

    * One
    * Two
    * Three...

    ## References
    
    * Reference
    * ...

    $ git status rfc-template.md

Untracked - working directory but unknown to Git.

Add to staging area / index / cache / "loading dock".

    $ git add rfc-template.md
    $ git status rfc-template.md

Changes to be committed - in staging area.

    $ git commit

Commit message - provide the "why" about a change.

Top tip: Write useful commit messages, not "made a change".

Commit shows number of files changed and the number of lines inserted or deleted.

    $ git status rfc-template.md

`nothing to commit` means file is in the repository, working directory up-to-date, no uncommitted changes in staging area.

    $ git log

Commit identifier (AKA revision number) uniquely identifies the changes made in this commit, author, date, and message.

    $ git log --relative-date
    $ cp rfc-template.md data-rfc.md
    $ git add data-rfc.md

Provide commit message via command-line.

    $ git commit -m "Create data RFC from rfc-template.md."

Edit.

    # Super distributed data concept RFC

    by Mike Jackson, Mario Antonioletti, Amrey Krause

    ## Overview

    A distributed data management platform.

    ## Key points

    * Revolutionary.
    * Novel.
    * Modest.

    ## References
    
    * WS-DAIR specification

    $ git status data-rfc.md

Modified - changed but not staged or commited.

    $ git add data-rfc.md
    $ git commit -m "Added title, authors, overview, points, references" data-rfc.md 

Top tip: good commits are atomic. Code should be reviewable in an hour.

What we know about software development - code reviews work. Fagan (1976) discovered that a rigorous inspection can remove 60-90% of errors before the first test is run. M.E., Fagan (1976). [Design and Code inspections to reduce errors in program development](http://www.mfagan.com/pdfs/ibmfagan.pdf). IBM Systems Journal 15 (3): pp. 182-211.

What we know about software development - code reviews should be about 60 minutes long. Cohen (2006) discovered that all the value of a code review comes within the first hour, after which reviewers can become exhausted and the issues they find become ever more trivial. J. Cohen (2006). [Best Kept Secrets of Peer Code Review](http://smartbear.com/SmartBear/media/pdfs/best-kept-secrets-of-peer-code-review.pdf). SmartBear, 2006. ISBN-10: 1599160676. ISBN-13: 978-1599160672.

    $ mkdir images
    $ cd images

Use wget to pull your photo from the web or find a suitable alternative (call it me.jpg or use whatever extension).

    $ wget http://software-carpentry.org/img/bootcamps/2012-12-uta.png
    $ mv 2012-12-uta.png > authors.png
    $ cd ..

Add hyper-link.   

    ![Group photo](images/authors.png "This is us!")

    $ git add images
    $ git commit -m "Added an images directory and image of authors." images
    $ git commit -m "Added link to image of authors." data-rfc.md

Top tip: Commit anything that cannot be automatically recreated e.g. `.tex`, `.java`, `.py`, `.c`, not `.pdf`, `.dvi`, `.dll`, `.jar`, `.exe`. Reduce risk of source-binary divergence.

### Discarding changes

Make more changes.

    $ git diff data-rfc.md

`-` line removed, `+` line added, `-` and `+` line edited.

Throw away changes or revert.

    $ git checkout -- data-rfc.md
    $ git status data-rfc.md

### History

    $ git log
    $ git log data-rfc.md

Globally-unique commit identifier.

    $ git diff COMMITID
    $ git diff OLDER_COMMITID NEWER_COMMITID
    $ git log
    $ git checkout COMMITID
    $ ls
    $ git checkout master
    $ ls

Undo and redo for directories and files.

Top tip: Commit often increases the granularity of "undo". 

DropBox and GoogleDrive also preserve every version, they delete old versions after 30 days, or, for GoogleDrive, 100 revisions. 

DropBox allows for old versions to be stored for longer but you have to pay for this. 

Using revision control the only bound is space available!

Commit changes to sets of files and rollback to exact state.

### Tags

Nicknames for commit identifiers

    $ git tag BOOT_CAMP
    $ git tag

Make more changes.

    $ git add data-rfc.md
    $ git commit -m "..." data-rfc.md
    $ git checkout BOOT_CAMP
    $ git checkout master

Top tip: tag significant "events" e.g. submitted papers, released versions.

### Branches

    $ git status data-rfc.md

`master` is a branch name.

Multiple sets of changes to files and directories - parallel instances.

Useful for bug-fixing releases while working on product, saves user waiting ror next full release.

Create branch for release.

    -o---o---o                               master
              \
               o                             release

Continue developing on default, master, branch.

    -o---o---o---o---o---o                   master
              \
               o                             release

Fix bug in release branch and release bug-fixed version.

    -o---o---o---o---o                       master
              \
               o---o---o                     release

Merge changes into master.

    -o---o---o---o---o---M                   master
              \         /
               o---o---o                     release

Continue developing.

    -o---o---o---o---o---M---o---o           master
              \         /
               o---o---o                     release

Popular model of release branch, master (up-to-date stable) branch, feature-specific and/or developer specific branches.

             0.1      0.2        0.3
              o---------o------o------    release
             /         /      /
     o---o---o--o--o--o--o---o---o---    master
     |            /     /
     o---o---o---o     /                 fred
     |                /
     o---o---o---o---o                   kate

## Working from multiple locations with a remote repository

Keep track of changes, a lab notebook for code and documents.

Roll back changes.

Delete repository means lose all files and all changes.

How to access repository from multiple locations e.g. desktop, laptop?

### GitHub and BitBucket

* [GitHub](http://github.com) 
* [BitBucket](http://bitbucket.com)
* [Launchpad](https://launchpad.net)
* [GoogleCode](http://code.google.com)
* [SourceForge](http://sourceforge.net)

Repositories and project infrastructure e.g. wiki, issue (ticket) and bug trackers, release management, commit-triggered e-mails etc.

All rights are reserved on unlicensed repositories - you should licence it to let people know what they can and cannot do with it.

BitBucket offers free, private repositories to researchers. 

GitHub and others offer pricing plans to host private repositories.

Git =/= GitHub

[Sign-up for free GitHub account](https://github.com/signup/free)

[Sign-up for free BitBucket account](https://bitbucket.org/account/signup/)

### Create new repository

GitHub:
* [Log in](https://github.com).
* Click on Create a new repo icon on top right, next to user name.
* Enter Repository name: `bootcamp`.
* Check Public option is selected.
* Check Initialize this repository with a README is unselected.
* Click Create Repository.

BitBucket:
* [Log in](https://bitbucket.com).
* Click on Create  icon on top left, next to Bitbucket logo.
* Enter Repository name: `bootcamp`.
* Check private repository option is ticked.
* Check repository type is `Git`.
* Check Initialize this repository with a README is unselected.
* Click Create Repository.

Push `master` branch to GitHub:

    $ git remote add origin https://github.com/USERNAME/bootcamp.git
    $ git push -u origin master

Push `master` branch to BitBucket:

    $ git remote add origin https://USERNAME@bitbucket.org/USERNAME/bootcamp.git
    $ git push -u origin --all

`origin` is an alias for repository URL. `-u` sets local repository to track remote repository.

GitHub, click Code tab and click Network tab.

BitBucket, click Source tab and click Commits tab.

### Cloning a remote repository

    $ cd ..
    $ rm -rf papers
    $ git clone https://github.com/USERNAME/bootcamp.git
    $ git clone https://USERNAME@bitbucket.org/USERNAME/bootcamp.git
    $ cd bootcamp
    $ git log
    $ ls -A

### Push changes to a remote repository

    # Make changes, add, commit.
    $ git push

GitHub, click Code tab and click Network tab.

BitBucket, click Source tab and click Commits tab.

Always pull before push.

### Pull changes from a remote repository

    $ cd ..
    $ ls
    $ git clone https://github.com/USERNAME/bootcamp.git anotherbootcamp
    $ git clone https://USERNAME@bitbucket.org/USERNAME/bootcamp.git anotherbootcamp
    $ ls

Pretend clones are on separate machines. 3 repositories - one remote, 2 local on our 'separate machines'.

    $ cd bootcamp
    $ nano data-rfc.md

Make changes, add, commit.

    $ git push
    $ cd ../anotherbootcamp
    $ git fetch
    $ git diff origin/master

Compare `master` branch with `origin/master` branch, alias for branch on remote repository.

    $ git merge origin/master
    $ cat data-rfc.md
    $ git log

pull = fetch + merge

    $ nano data-rfc.md

Make changes, add, commit.

    $ git push
    $ cd ../bootcamp
    $ git pull
    $ cat data-rfc.md
    $ git log

### Conflicts and resolution

    $ nano data-rfc.md

Make changes, add, commit.

    $ git push
    $ cd ../anotherbootcamp

Make changes to same lines, add, commit.

    $ git push

Push fails. Fetch changes made to remote repository.

    $ git fetch
    $ git diff origin/master
    $ git merge origin/master

Merge done file-by-file, line-by-line. Conflict - file has two changes affecting the same line.

    $ git status

Unmerged.

    $ cat data-rfc.md

`<<<<<<< HEAD` lines from local version.

`=======` delimiter

`>>>>>>> 71d34decd32124ea809e50cfbb7da8e3e354ac26` lines from remote version.

Keep local, or keep remote, or combine both. Remove all the mark-up.

    $ git add data-rfc.md
    $ git commit -m "Resolved conflict in data-rfc.md by ...."
    $ git push

GitHub, click Code tab and click Network tab.

BitBucket, click Source tab and click Commits tab.

DropBox and GoogleDrive don't do this. No work is lost.

## Collaborating with colleagues

Form into pairs and swap GitHub / BitBucket user names.

One of you (Owner) share your repository with your partner.

Owner, on GitHub, click on the Settings tab, click on Collaborators, add partner's GitHub name.

Owner, on BitBucket, click Share link, add partner's BitBucket name.
Both,

    $ git clone https://github.com/OWNERUSERNAME/bootcamp.git
    $ git clone https://USERNAME@bitbucket.org/USERNAME/bootcamp.git 

Now:

* Both edit same file, add and commit.
* Owner push.
* Partner push, fetch, merge, resolve conflicts (if any), add, commit, push.
* Owner pull.
* Both edit same file, add and commit.
* Owner push, fetch, merge, resolve conflicts (if any), add, commit, push.
* Partner pull.
* Repeat! Try adding and editing new files too.

## More on branching

    $ git branch

`*` signifies current branch.

    $ git branch new_template
    $ git branch
    $ git checkout new_template
    $ git branch
    $ nano template-rfc.md

Change, add, commit.

    $ git checkout master
    $ ls
    $ git checkout new_template
    $ ls
    $ git checkout master
    $ ls
    $ git merge new_template
    $ git branch -D new_template

## Push the Software Carpentry repository to BitBucket

Create "swc" repository.

    cd ~
    cd boot-camps
    git push -u http://bitbucket.org/USERNAME/swc.git 2013-07-bath

## Summary

* Keep track of changes, a lab notebook for code and documents.
* Roll back changes to any point in the history of changes, "undo" and "redo" for files.
* Back up entire history of changes in various locations.
* Work on files from multiple locations.
* Identify and resolve conflicts when the same file is edited within two repositories without losing any work.
* Collaboratively work on code or documents or any other files.

## Links

* K. Ram  (2013) "git can facilitate greater reproducibility and increased transparency in science", Source Code for Biology and Medicine 2013, 8:7 doi:[10.1186/1751-0473-8-7](http://dx.doi.org/10.1186/1751-0473-8-7) - survey of the range of ways in which version control can help research.
* [Visual Git Reference](http://marklodato.github.com/visual-git-guide/index-en.html) - pictorial representations of what Git commands do.
* [Pro Git](http://git-scm.com/book) - the "official" online Git book.
* [Version control by example](http://www.ericsink.com/vcbe/) - an acclaimed online book on version control by Eric Sink.
* [Git commit policies](http://osteele.com/posts/2008/05/commit-policies) - images on what Git commands to with reference to the working directory, staging area, local and remote repositories.
* [Gitolite](https://github.com/sitaramc/gitolite) - a way for you to host
your own multi-user Git repositories. Your collaborators send you their public SSH keys then they can pull and push from/to the repositories.
* [A successful Git branching model](http://nvie.com/posts/a-successful-git-branching-model/)
* G. Wilson, D. A. Aruliah, C. T. Brown, N. P. Chue Hong, M. Davis, R. T. Guy, S. H. D. Haddock, K. Huff, I. M. Mitchell, M. Plumbley, B. Waugh, E. P. White, P. Wilson (2012) "[Best Practices for Scientific Computing](http://arxiv.org/abs/1210.0530)", arXiv:1210.0530 [cs.MS].

## Git hints and tips

### `man` page

Like many Unix/Linux commands, `git` has a `man` page,

    $ man git

You can scroll the manual page up and down using the up and down arrows.

You can search for keywords by typing `/` followed by the search term e.g. if interested in help, type `/help` and then hit enter.

To exit the manual page, type `q`.

### Add a repository description

You can edit the file `.git/description` and give your repository a name e.g. "My first repository".

### Ignore scratch, temporary and binary files

You can create a `.gitignore` file which lists the patterns of files you want Git to ignore. It's common practice to not add to a repository any file you can automatically create in some way e.g. C object files (`.o`), Java class (`.class`) files or temporary files e.g. XEmacs scratch files (`~`). Adding these to `.gitignore` means Git won't complain about them being untracked.

Create or edit `gitignore`,

    $ nano .gitignore

Then add patterns for the files you want to ignore, where `*` is a wildcard,

    *~
    *.o
    *.so
    *.dll
    *.exe
    *.class
    *.jar

Then, add `.gitignore` to your repository,

    $ git add .gitignore
    $ git commit -m "Added rules to ignore XEmacs scratch files and binary files"

### Add colour to `diff`

    $ git config --global color.diff auto
