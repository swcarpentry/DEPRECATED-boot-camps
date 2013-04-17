## Working from multiple locations with a remote repository

We're going to set up a remote repository that we can use from multiple locations. The remote repository can also be shared with colleagues, if we want to.

### Bitbucket

[Bitbucket](http://bitbucket.org) is a web service which allows users to set up  their private and public source code Git repositories. It provides tools for browsing, collaborating on and documenting code. Your organisation may also offer support for hosting Git repositories - ask your local system administrator. Bitbucket, like other services such as [Launchpad](https://launchpad.net), [GitHub](https://github.com),
[GoogleCode](http://code.google.com), and [SourceForge](http://sourceforge.net) provides a wealth of resources to support projects including:

* Time histories changes to repositories
* Commit-triggered e-mails
* Browsing code from within a web browser, with syntax highlighting
* Software release management
* Issue (ticket) and bug tracking
* Download
* Varying permissions for various groups of users
* Other service hooks e.g. to Twitter.

**Note** Bitbucket allows you to set up private repositories for free. However, GitHub's free repositories have public licences **by default**. If you don't want to share (in the most liberal sense) your stuff with the world and you want to use GitHub (instead of BitBucket), you will need to pay for the private GitHub repositories. GitHub and others offer pricing plans different pricing plans, for more details, check out their websites.

### Get an account

Setting up Bitbucket at first requires a user name and password. Let's now  [create a free one](https://bitbucket.org/account/signup/). 

### Create a new repository

Now, we can create a repository on Bitbucket,

* Log in to [Bitbucket](https://github.com)
* Click on the Create  icon on the top left, next to the Bitbucket logo
* Enter Repository name: `bootcamp`
* Make sure that the private repository option is ticked
* Make sure that you select Git in the `repository type` section
* Make sure the Initialize this repository with a README is *unselected*
* Click Create Repository

You'll get a page with new information about your repository. We already have our local repository and we will be copying it to Bitbucket. On the page which just opened (`Add some code`) click `I have code I want to import`. Copy the following line (from Bitbucket) and paste it into your command line

    $ git remote add origin https://USERNAME@bitbucket.org/USERNAME/bootcamp.git

This sets up an alias, `origin`, to correspond to the URL of our new repository on Bitbucket.

Now copy and paste the second line,

    $ git push -u origin --all   # to push changes for the first time
    Password for 'https://USERNAME@bitbucket.org': 
    Counting objects: 3, done.
    Writing objects: 100% (3/3), 220 bytes, done.
    Total 3 (delta 0), reused 0 (delta 0)
    remote: bb/acl: USERNAME is allowed. accepted payload.
    To https://USERNAME@bitbucket.org/USERNAME/bootcamp.git
    * [new branch]      master -> master
    Branch master set up to track remote branch master from origin.


This *pushes* our `master` branch to the remote repository, named via the alias `origin` and creates a new `master` branch in the remote repository.

Now, on Bitbucket, click the Source tab and we should see our code and click the Commits tab we should see our complete history of commits.  

Our local repository is now available on Bitbucket. So, anywhere we can access Bitbucket, we can access our repository.

### Cloning a remote repository

Now, let's do something drastic!

    $ cd ..
    $ rm -rf papers

Gulp! We've just wiped our local repository! But, as we've a copy on Bitbucket we can just copy, or *clone* that,

    $ git clone https://USERNAME@bitbucket.org/USERNAME/bootcamp.git
    Cloning into 'bootcamp'...
    Password for 'https://USERNAME@bitbucket.org':
    remote: Counting objects: 12, done.
    remote: Compressing objects: 100% (4/4), done.
    remote: Total 12 (delta 0), reused 0 (delta 0)
    Unpacking objects: 100% (12/12), done.


In Bitbucket, there is also an option for called *forking* which essentially allows to do the same thing but forking is done in the web service (that is not using the command line and creating a repository on your local file system).

Now, if we change into `bootcamp` we can see that we have our repository,

    $ cd bootcamp
    $ git log

and we can see our Git configuration files too,

    $ ls -A
    common  .git  journal.txt

But where is the `papers` directory, you might ask? `papers` was the directory that held our local repository but was not a part of it.

### Push changes to a remote repository

We can use our cloned repository just as if it was a local repository so let's make some changes to our files and commit these.

Having done that, how do we send our changes back to the remote repository? We can do this by *pushing* our changes,

    $ git push

If we now check our Bitbucket page we should be able to see our new changes under the Commit tab.

Before you push to a remote repository you should always pull so you have the most up-to-date copy of the remote repository. So, on that note...

### Pull changes from a remote repository

So, we can now access our repository from multiple locations and make changes to it. But how do we get the latest changes. One way is simply to clone the repository every time but this is inefficient, especially if our repository is very large. So, Git allows us to get the latest changes down from a repository. 

So, first, let us leave our current local repository,

    $ cd ..
    $ ls
    bootcamp

And let us clone our repository again, but this time specify the local directory name,

    $ git clone https://USERNAME@bitbucket.org/USERNAME/bootcamp.git bootcamp2
    Cloning into 'bootcamp2'...


So we now have two clones of our repository,

    $ ls
    $ bootcamp anotherbootcamp

Let's pretend these clones are on two separate machines! So we have 3 versions of our repository - our two local versions, on our separate machines (we're still pretending!) and one on Bitbucket. So let's go into one of our clones, make some changes, commit these and push these to Bitbucket:

    $ cd bootcamp
    $ nano journal.txt
    $ git add journal.txt
    $ git commit -m "Added some notes to the sections" journal.txt
    $ git push

Now let's change to our other repository and *fetch* the changes from our remote repository,

    $ cd ../anotherbootcamp
    $ git fetch

We can now see what the differences are by doing,

    $ git diff origin/master

which compares our current, `master` branch, with an `origin/master` branch which is the name of the `master` branch in `origin` which is the alias for our cloned repository, the one on Bitbucket.

We can then *merge* these changes into our current repository, which merges the branches together,

    $ git merge origin/master

And then we can check that we have our changes,

    $ cat journal.txt
    $ git log

As a short-hand, we can do a Git *pull* which does a *fetch* then a *merge*,

    $ nano journal.txt
    $ git add journal.txt
    $ git commit -m "Added some references" journal.txt
    $ git push
    $ cd ../bootcamp
    $ git pull

And then check that we have our changes,

    $ cat journal.txt
    $ git log

### Conflicts and how to resolve them

Let's continue to pretend that our two local, cloned, repositories are hosted on two different machines, and make some changes to our file, and push these to Bitbucket,

    $ nano journal.txt
    $ git add journal.txt
    $ git commit -m "Rewrote the title" journal.txt
    $ git push

Now let us suppose, at a later, date, we use our other repository (for example, we may have been working on a local repository on our laptop, and now are using the one on our workstation) and we come up with a better idea for a title,

    $ cd ../anotherbootcamp
    $ git add journal.txt
    $ git commit -m "Changed the first author" journal.txt
    $ git push

Our push fails, as we've not yet pulled down our changes from our remote repository. Before pushing we should always pull, so let's do that...

    $ git pull

and we get

    Auto-merging journal.txt
    CONFLICT (content): Merge conflict in journal.txt
    Automatic merge failed; fix conflicts and then commit the result.

As we saw earlier, with the fetch and merge, a pull pulls down changes from the repository and tries to merge these. It does this on a file-by-file basis, merging files line by line. We get a *conflict* when if a file has changes that affect the same lines and those changes can't be seamlessly merged. If we look at the status,

    $ git status

we can see that our file is listed as `Unmerged` and if we look at `journal.txt`, we may see something like,

    <<<<<<< HEAD 
    Title: A paper about proteines
    =======
    Title: A paper everything but proteines
    >>>>>>> 71d34decd32124ea809e50cfbb7da8e3e354ac26 

The mark-up shows us the parts of the file causing the conflict and the versions they come from. We now need to manually edit the file to *resolve* the conflict. This means removing the mark-up and doing one of

* Keep the local version, which, here, is the one marked-up by `HEAD` i.e. `Title: A paper about proteines`
* Keep the remote version, which, here, is the one marked-up by the commit identifier i.e. `Title: A paper everything but proteines`
* Or keep a combination of the two e.g. `Title: A paper about proteines and eveything else`

We edit the file. Then commit our changes e.g.

    $ git add journal.txt
    $ git commit -m "Resolved conflict in journal.txt by rewriting title to combine best of both originals"

Now if we push,

    $ git push

All goes well. If we now go to Bitbucket and click on the Overview tab we can see where our repository diverged and came together again.

This is where version control proves itself better than DropBox or GoogleDrive, this ability to merge text files line-by-line and highlight the conflicts between them, so no work is ever lost.

## The story so far...

So far, we've now seen how we can,

* Host our private repository on Bitbucket
* Copy, or clone, our remote repository onto a local machine
* Make changes in a local repository and push these to a remote repository
* Fetch and merge, or pull, changes from a remote repository into our local repository
* Identify and resolve conflicts when the same file is edited within two repositories

We now have everything we need to collaborate with our colleagues!

Previous: [Tracking our changes with a local repository](1_Local.md) Next: [Collaborating with our colleagues](3_Collaboration.md)
