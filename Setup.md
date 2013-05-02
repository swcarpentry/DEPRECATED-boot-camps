# Manual set-up

You should ensure you have the following software and tools available. 

* Web browser (preferably Firefox or Google Chrome)
* Bash shell
* A text editor e.g. nano, vi or emacs
* Make
* [Git](http://git-scm.com/)
* Python 2.7 or higher

Below there are setup instructions for Windows (XP, Vista) RedHat/Scientific Linux 6 and Ubuntu and Mac OSX.
If you have any problems then please send an email to: admin [_at_] software-carpentry.org.

## To install under Windows 

In order to work shell, scripting, make and git on Windows, you 'll need [Cygwin](http://www.cygwin.com/]).
Once you download and run the installation package make sure that you install:
* A text editor (during the tutorials we will be using "nano")
* Make
* Git
* Python 2.7 or higher

To install the above modules simply type in their names in the "Select packages" window during the installation process (note that Make and Git will be listed in the 'Devel' section).

## To install under RedHat/Scientific Linux 6

Scientific Linux 6 already comes with shell and vi text editor. To install the other packages run,

    $ sudo su -
    # yum install nano
    # yum install git
    # yum install python
    
## To install under Ubuntu

Ubuntu 11.04 and above already comes with shell, vi and nano text editors. To install the other packages run,

    $ sudo su -
    # apt-get install git
    # apt-get install python
   
