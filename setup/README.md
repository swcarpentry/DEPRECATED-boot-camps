## Setup

### Overview

#### Editor

When you're writing code, it's nice to have a text
editor that is optimized for writing code, with features
like automatic color-coding of key words.
The default text editor on Mac OS X and Linux is usually set to Vim,
which is not famous for being intuitive.
If you accidentally find yourself stuck in it,
try typing the escape key,
followed by `:q!` (colon, lower-case q, exclamation mark),
then hitting Return
to return to the shell.

#### The Bash Shell

Bash is a commonly-used shell. Using a shell gives you
more power to do more tasks more quickly with your
computer.

#### Git

Git is a state-of-the-art version control system. It
lets you track who made changes to what when and has
options for easily updating a shared or public version of
your code on [github.com](https://github.com/).

#### Python

Python is becoming very popular in scientific computing,
and it's a great language for teaching general programming concepts due to its easy-to-read syntax.
We teach with Python version 2.7,
since it is still the most widely used.
Installing all the scientific packages for Python individually can be a bit difficult,
so we recommend an all-in-one installer.

---

### Windows

#### Python

Download and install  [Anaconda CE](http://continuum.io/anacondace.html).

Use all of the defaults for installation
_except_ make sure to check **Make Anaconda the default Python**.

#### Git Bash

Install Git for Windows by downloading and running
[the installer](http://msysgit.github.io/).
This will provide you with both Git and Bash in the Git Bash program.

#### Software Carpentry Installer

After installing Python and Git Bash:

- Download the [installer](http://files.software-carpentry.org/SWCarpentryInstaller.exe).
  (This installer requires an active internet connection.)
- If the file opens directly in the browser select **File &rarr; Save Page As**
  to download it to your computer.
- Double click on the file to run it.

#### Editor

`nano` is the editor installed by the Software Carpentry Installer.
It is a basic editor integrated into the lesson material.

[Notepad++](http://notepad-plus-plus.org/) is a
popular free code editor for Windows.
Be aware that you must add its installation directory to your system path
in order to launch it from the command line
(or have other tools like Git launch it for you).
Please ask your instructor to help you do this.

---

### Mac OS X

#### Bash

The default shell in all versions of Mac OS X is bash,
so no need to install anything.  You access bash from
the Terminal (found
in `/Applications/Utilities`).  You may want
to keep Terminal in your dock for this workshop.

#### Editor

We recommend
[Text Wrangler](http://www.barebones.com/products/textwrangler/) or
[Sublime Text](http://www.sublimetext.com/).
In a pinch, you can use `nano`,
which should be pre-installed.

#### Git

Install Git for Mac 10.8 (Mavericks or newer) by downloading and running
[the installer](http://git-scm.com/downloads).  For older
versions of OS X (10.6-10.8, before Mavericks) use the most recent available
installer for Snow Leopard [available
here](http://sourceforge.net/projects/git-osx-installer/files/).

#### Python

Download and install [Anaconda CE](http://continuum.io/anacondace.html).

Use all of the defaults for installation _except_ make sure to check
**Make Anaconda the default Python**.

---

### Linux

#### Bash

The default shell is usually `bash`,
but if your machine is set up differently
you can run it by opening a terminal and typing `bash`.
There is no need to install anything.

#### Git

If Git is not already available on your machine you can try
to install it via your distribution's package manager
(e.g. `apt-get` or `yum`).

#### Editor

[Kate](http://kate-editor.org/) is one option for Linux users.
In a pinch, you can use `nano`,
which should be pre-installed.

#### Python

We recommend the all-in-one scientific Python installer
[Anaconda](http://continuum.io/downloads.html).
(Installation requires using the shell and if you aren't
comfortable doing the installation yourself just
download the installer and we'll help you at the boot
camp.)

- Download the installer that matches your operating
  system and save it in your home folder.
- Open a terminal window.
- Type `bash Anaconda-` and then press
  tab. The name of the file you just downloaded should
  appear.
- Press enter. You will follow the text-only
  prompts. When there is a colon at the bottom of the
  screen press the down arrow to move down through the
  text. Type `yes` and press enter to approve
  the license. Press enter to approve the default
  location for the files. Type `yes` and
  press enter to prepend Anaconda to
  your `PATH` (this makes the Anaconda
  distribution the default Python).
