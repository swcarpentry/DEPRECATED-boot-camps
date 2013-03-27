---
root: ..
layout: base
title: University of Oxford Doctoral Training Centres Software Carpentry Boot Camp
---

<div>
<a href="http://software-carpentry.org">
<img src="{{page.root}}/logos/software-carpentry-with-hammer.svg" style="float: right; height: 60px;">
</a>
<h1>Oxford DTCs</h1>
</div>

 * Dates: 9-10 May 2013
 * Times: 9am - 5pm

## Schedule

Timings are to be confirmed, but the sessions will be:

* Registration and introduction
* What we actually know about developing software, and why we believe it's true
* Version control (with [git][])
* The basics of [Python][], Functions and Testing
* Scientific programming in Python: [NumPy][] and [SciPy][]
* Using relational databases
* Applying software carpentry in Matlab: coding style, code organisation, testing, databases
* Hopefully some flexible time to discuss attendee interests
* Putting it all together

[git]: http://git-scm.com/
[Python]: http://python.org/
[NumPy]: http://www.numpy.org/
[SciPy]: http://www.scipy.org/

## Installation

The venue will provide Windows desktops with the necessary software installed.  However, if you have your own laptop we recommend that you bring it with you to work on during the boot camp.  In this case you should ensure the following software is installed before you arrive, so that we can get started promptly with the tutorials.

 * Bash (the particular Unix shell we'll be using)
 * Git (for version control)
 * [Python][] (2.6 or 2.7) with [NumPy][], [SciPy][], the [iPython notebook][] and [matplotlib][]
 * [Matlab][]
 * Firefox with SQLite Manager Add-on (for relational databases)

The subsections below explain how to install these tools on different operating systems.
You can also save time on the day by creating an account beforehand at [GitHub][].

[iPython notebook]: http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html
[matplotlib]: http://matplotlib.org/
[Matlab]: http://people.maths.ox.ac.uk/gilesm/matlab.html
[GitHub]: https://github.com/

### All platforms (Linux, Mac OS X, Windows)

#### Python

Installing everything you need on your own can be a bit difficult so we recommend just installing the [Enthought Python Distribution][EPD], which has an [academic][EPD Acad] version for Mac, Windows, and Linux if you sign up with your Oxford University email address.  There is also a cut-down [free][EPD Free] version for anyone.

For other options check the Python4Astronomers page on [installing scientific Python][astpy].

[EPD]: http://www.enthought.com/products/epd.php
[EPD Free]: http://www.enthought.com/products/epd_free.php
[EPD Acad]: http://www.enthought.com/products/edudownload.php
[astpy]: http://python4astronomers.github.com/installation/python_install.html

#### Firefox

You can download Firefox from the [Mozilla web site][mozilla] (if it's not available via your normal system installation process), and then the SQLite Manager add-on can be installed through the "Add-ons Manager" in Firefox.

[mozilla]: http://www.mozilla.org/

### Linux

The default shell is usually bash but if not you can get to bash by opening a terminal and typing `bash`.

If git is not already available on your machine you should be able to install it via your distro's package manager (e.g. `apt-get` or `yum`).

Many different text editors suitable for programming are available.  If you don't already have a favourite, you could look at [Kate].

[Kate]: http://kate-editor.org/

### Max OS X

The default shell in Mac OS X is bash.

For git, either install [Xcode][] and the command line tools (from the Download preferences pane) or [install just git][Mac git].

Many different text editors suitable for programming are available.  We recommend [Text Wrangler][] or [Sublime Text][].

[Xcode]: https://developer.apple.com/xcode/
[Mac git]: http://code.google.com/p/git-osx-installer/downloads/list?can=3
[Text Wrangler]: http://www.barebones.com/products/textwrangler/
[Sublime Text]: http://www.sublimetext.com/

### Windows

For a bash shell and git, the simplest route is to install '[Git Bash][]' following the instructions [here][gitbash].  Alternatively, you can use [Cygwin][] and install its git package.

The Enthought Python Distribution includes the [SciTE] editor.  [Notepad++] is a popular free code editor for Windows.

[Git Bash]: http://msysgit.github.com/
[gitbash]: https://openhatch.org/missions/windows-setup/install-git-bash
[Cygwin]: http://www.cygwin.com/
[SciTE]: http://www.scintilla.org/SciTE.html
[Notepad++]: http://notepad-plus-plus.org/

## Contact

This workshop is being organised by [Jonathan Cooper][].

[Jonathan Cooper]: http://www.cs.ox.ac.uk/people/jonathan.cooper

