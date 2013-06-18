# Manual set-up

You should ensure you have the following software and tools available: 

* Web browser (preferably Firefox or Google Chrome)
* Bash shell
* A text editor e.g. nano, vi or emacs
* [Git](http://git-scm.com/)
* Python 2.7 or higher
* [numpy](http://numpy.scipy.org/), [scipy](http://scipy.org) and [matplotlib](ttp://matplotlib.org/) for Scientific Python
* [nosetests](https://nose.readthedocs.org/en/latest/) framework for Testing

If you experience any problems, please arrive at 8:30 on Day 1 and the instructors and helpers will help you deal with the issue. 

## To install:

### Bash:  
**Mac:**  
The default shell in Mac OS X is bash.
**Windows:**  
In order to practice working with shell, scripting and Git on Windows, we will use [Cygwin](http://www.cygwin.com/.
Once you download and run the installation package, make sure that you install:
* A text editor (nano - probably the easiest choice if you're new to *nix environments, emacs, vi etc.)
* Git
* python-nose  
To install the above packages simply type in their names in the "Select packages" window during the installation process (note that Git will be listed in the "Devel" section).

**Linux:**  
The default shell is usually bash but if not you can get to bash by opening a terminal and typing `bash`.

### Git
**Mac:**   
Install Xcode and the command line tools (from the Download preferences pane) or install just Git.

**Windows:**  
Install Git using cygwin packages (see above).

**Linux:**  
If git is not already available on your machine you can try to install it via your distro's package manager (e.g. apt-get).


### Python and Scientific Python packages
We recommend the all-in-one scientific Python installer [Anaconda CE](http://continuum.io/anacondace.html). Installation on Mac and Linux requires using the shell and if you aren't comfortable doing the installation yourself just download the installer and we'll help you at the boot camp.

### Python Nose

**Mac and Linux:**  
`easy_install nose` or pip `install nose`

To check the installation run `nosetests` and the output should be similar to:

    ..
    ----------------------------------------------------------------------
    Ran 2 tests in 0.015s

    OK

**Windows:**  
Python Nose can be installed via cygwin (please see above).

   
