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
 * Location: Computer Room C21, Department of Zoology, South Parks Road, Oxford, OX1 3PS

## Schedule

The sessions will run (roughly) as follows.
Drinks and edibles will be provided for the coffee breaks, but you will need to make your own lunch arrangements.
Links to the materials for each session are also given, although these will not be finalised until a day or two before the workshop.

* Day 1, Thursday 9th May
    * 09:00 Registration and introduction
    * 09:30 Version control (with [git][]) \[[materials][m vcs]]
    * _11:00 Coffee break_
    * 11:30 Version control continued
    * 12:00 The basics of [Python][] \[[materials][m python]]
    * _13:00 Lunch_
    * 14:00 Organising code: Python functions \[[materials][m func]]
    * _15:00 Coffee break_
    * 15:30 Testing your code \[[materials][m test]]
* Day 2, Friday 10th May
    * 09:00 Scientific programming in Python: [NumPy][] and [SciPy][] \[[materials][m numpy]]
    * _11:00 Coffee break_
    * 11:30 Using relational databases \[[materials][m sql]]
    * _13:00 Lunch_
    * 14:00 Applying software carpentry in MATLAB \[[materials][m matlab]]
    * _15:00 Coffee break_
    * 15:30 Some flexible time to discuss attendee interests (hopefully!)
    * 16:15 What we actually know about developing software, and why we believe it's true
    * 16:30 Putting it all together
    * 16:45 Feedback

[git]: http://git-scm.com/
[Python]: http://python.org/
[NumPy]: http://www.numpy.org/
[SciPy]: http://www.scipy.org/

[m vcs]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/version-control/git-oxford
[m python]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/python/2a-PythonVariables
[m func]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/python/2c-PythonFlowControl
[m test]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/python/testing-oxford
[m numpy]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/python/numpy
[m sql]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/sql/sql.md
[m matlab]: https://github.com/swcarpentry/boot-camps/tree/2013-05-oxford-dtc/matlab

## Before you arrive

The venue will provide Windows desktops with the necessary software installed.
However, if you have your own laptop we recommend that you bring it with you to work on during the boot camp.
In this case you will need to ensure you have installed some software on your laptop, following the instructions below.
Please do this well in advance of turning up (not the night before and certainly not on the first morning!).
Jonathan will be available in the DTC from 9.30am-12.30pm on Wednesday to troubleshoot any issues.

You can also save time on the day by creating a free account beforehand at [Bitbucket][].

A pre-camp [questionnaire][q] is available at [http://goo.gl/FFr6R][q].
If you could take a minute to fill this in it will help us get an idea as to what level we should pitch our material.

[q]: http://goo.gl/FFr6R
[GitHub]: https://github.com/
[Bitbucket]: https://bitbucket.org/

## Installation

If bringing your own machine, please ensure the following software is installed before you arrive, so that we can get started promptly with the tutorials.

 * Bash (the particular Unix shell we'll be using)
 * Git (for version control)
 * For convenience when using git, you will also want a friendly command-line-based editor, such as [nano][].
 * [Python][] (2.6 or 2.7) with [NumPy][], [SciPy][], [iPython][] and [matplotlib][]
 * [Matlab][]
 * Firefox with SQLite Manager Add-on (for relational databases)

The subsections below explain how to install these tools on different operating systems.

[iPython]: http://ipython.org/
[iPython notebook]: http://ipython.org/ipython-doc/dev/interactive/htmlnotebook.html
[matplotlib]: http://matplotlib.org/
[Matlab]: http://people.maths.ox.ac.uk/gilesm/matlab.html
[nano]: http://www.nano-editor.org/

### All platforms (Linux, Mac OS X, Windows)

#### Python

Installing everything you need on your own can be a bit difficult so we recommend just installing the [Enthought Python Distribution][EPD], which has an [academic][EPD Acad] version for Mac, Windows, and Linux if you sign up with your Oxford University email address.  There is also a cut-down [free][EPD Free] version for anyone.

For other options check the Python4Astronomers page on [installing scientific Python][astpy].

[EPD]: http://www.enthought.com/products/epd.php
[EPD Free]: http://www.enthought.com/products/epd_free.php
[EPD Acad]: http://www.enthought.com/products/edudownload.php
[astpy]: http://python4astronomers.github.com/installation/python_install.html

#### Firefox

You can download Firefox from the [Mozilla web site][mozilla] (if it's not available via your normal system installation process), and then the SQLite Manager add-on can be installed through the "Add-ons Manager" in Firefox.  Make sure you get [version 0.8.0][] as earlier versions do not work properly with the latest Firefox.

[mozilla]: http://www.mozilla.org/
[version 0.8.0]: https://addons.mozilla.org/en-US/firefox/addon/sqlite-manager/versions/

### Linux

The default shell is usually bash but if not you can get to bash by opening a terminal and typing `bash`.

If git is not already available on your machine you should be able to install it via your distro's package manager (e.g. `apt-get` or `yum`).  Similarly for nano.

Many different text editors suitable for programming are available.  If you don't already have a favourite, you could look at [Kate].

[Kate]: http://kate-editor.org/

### Max OS X

The default shell in Mac OS X is bash.

For git, either install [Xcode][] and the command line tools (from the Download preferences pane) or [install just git][Mac git].  I believe Mac OS X comes with pico, which is very similar to nano.

Many different text editors suitable for programming are available.  We recommend [Text Wrangler][] or [Sublime Text][].

[Xcode]: https://developer.apple.com/xcode/
[Mac git]: http://code.google.com/p/git-osx-installer/downloads/list?can=3
[Text Wrangler]: http://www.barebones.com/products/textwrangler/
[Sublime Text]: http://www.sublimetext.com/

### Windows

For a bash shell and git, the simplest route is to install '[msysgit][]' following the instructions [here][gitbash].  We provide a script [swc-windows-installer.py][] to install nano for you; download it and double-click to run, assuming you're using the Enthought Python Distribution.

Alternatively, you can use [Cygwin][] and install its git and nano packages.  This is probably a better route if you are likely to use the command-line much long term.

The full Enthought Python Distribution includes the [SciTE] editor (the free version only includes [IDLE][]).  [Notepad++] is a popular free code editor for Windows.


[msysgit]: http://msysgit.github.com/
[gitbash]: https://openhatch.org/missions/windows-setup/install-git-bash
[Cygwin]: http://www.cygwin.com/
[SciTE]: http://www.scintilla.org/SciTE.html
[Notepad++]: http://notepad-plus-plus.org/
[IDLE]: http://docs.python.org/2/library/idle.html
[swc-windows-installer.py]: https://github.com/jonc125/boot-camps/blob/2013-05-oxford-dtc/setup/swc-windows-installer.py

## Venue directions

1. Enter at main reception on South Parks Road.  You may need to show your University card.
2. Go straight ahead into Zoology.
3. Go up the stairs to the left, following signs for Darwins Cafe & Teaching Labs.
4. Do a u-turn right at the top through Darwins gallery.
5. Take the first right and the room (C21) is in front of you.

Coffee breaks & toilets are just down the stairs by the door to C21.

## Contact

This workshop is being organised by [Jonathan Cooper][].

[Jonathan Cooper]: http://www.cs.ox.ac.uk/people/jonathan.cooper

