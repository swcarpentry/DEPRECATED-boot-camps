#Teaching Git - 1

##What is Version Control?
__any guesses?__

__Some questions__

* How easy would it be to grab the code you used to develop a project/manuscript from last week? Last year? Last 5 years?
* How much work would you lose if your laptop was stolen/damaged?
* If you deleted part of your manuscript or analysis, then later realized that you wanted to keep it in, would you have to rewrite it?
* If you were collaborating with someone on a manuscript or some code and someone rewrote something, how easy would it be to track down who changed it and why?

__Version control...__

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

