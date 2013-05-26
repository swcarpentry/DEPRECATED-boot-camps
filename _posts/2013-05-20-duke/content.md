<a id="location"></a> Location
==============================

CRTP Classroom, 2nd floor, Hock Plaza

2424 Erwin Rd, Durham, NC

<a id="schedule"></a> Schedule
==============================

### Day 1

* 9:00-10:30 R Basics, RStudio basics, workspace, working directory, projects in RStudio
* 10:30-11:00 *coffee break*
* 11:00-12:30 Basic care and feeding of data, special emphasis on data.frames
* 12:30-1:30 *lunch*
* 1:30-3:00 Data aggregation, reshaping
* 3:00-3:30 *coffee break*
* 3:30-4:30 putting it all together: a bit of file-fu with R, organizing your workflow, strategies for 
    reproducibility and share ability

optional content (will be covered by example as we go but might be tackled explicitly briefly): 
* lattice graphics

### Day 2

* 9:00-10:30 using the shell, part 1
* 10:30-11:00 *coffee break*
* 11:00-12:30 version control with Git
* 12:30-1:30 *lunch*
* 1:30-3:00 Make and special topics
* 3:00-3:30 *coffee break*
* 3:30-4:30 special topics


<a id="install"></a> Installation
=================================

You need to install and test R, RStudio, Git, and a Bash shell before the workshop. Download some files 
we will use, just in case we have any network problems. It is recommended to also install a few 
add-on packages.

### R and RStudio

* Download and install [R, a free software environment for statistical computing and graphics](http://www.r-project.org) from [CRAN](http://cran.rstudio.com), the Comprehensive R Archive Network. It is _highly recommended_ to install a precompiled binary distribution for your operating system -- use the links up at the top of the CRAN page linked to above!

* Install RStudio, a powerful user interface for R: <http://www.rstudio.com/ide/download/>

### Testing testing

* Do whatever is appropriate for your OS to launch RStudio. You should get a window similar to the screenshot you see [here](http://www.rstudio.com/ide/), but yours will be more boring because you haven't written any code or made any figures yet!

* Put your cursor in the pane labelled Console, which is where you interact with the live R process. Create a simple object with code like `x <- 2 * 4` (followed by enter or return). Then inspect the `x` object by typing `x` followed by enter or return. Obviously you should see the value 8 print to screen. If yes, you are good to go.

### Gapminder stuff

* We will work with some of the data from the [Gapminder project](http://www.gapminder.org). Here is an excerpt I have prepared for our use. Please save this file on your computer prior to the workshop and note the location!
  - <http://www.stat.ubc.ca/~jenny/notOcto/STAT545A/examples/gapminder/data/gapminderDataFiveYear.txt>

### Add-on packages

* Installing add-on packages. R is an extensible system and many people share useful code they have developed as a _package_ via CRAN and github. To install a package, for example the [`plyr` package](http://plyr.had.co.nz) for data aggregation, here is one way to do it in the R console (there are others).

    ```
    install.packages("plyr", dependencies = TRUE)
    ```
    Another package you may wish to play around with soon is [`knitr`](http://yihui.name/knitr/), which facilitates the creation of dynamic reports. You could install it in the same way.
    ```
    install.packages("knitr", dependencies = TRUE)
    ```

### Unix Bash Shell

* Please follow the [Software Carpentry instructions](http://software-carpentry.org/bootcamps/setup.html) 
    in the section titled "The Unix Bash Shell" to get a Bash shell on your machine.

* If you're using **Windows**, you'll want to install `nano`, an easy-to-use
    command-line text editor.
    [Download it here](http://www.nano-editor.org/dist/v2.2/NT/nano-2.2.6.zip),
    extract to a directory (like C:/Program Files/nano), and add that directory
    to your Windows path:
    
    * Right-click on **My Computer** and select **Properties**
    
    * Select the **Environment** page

    * In the System Variables area, find the PATH variable

    * Add the directory you extracted nano in to the end of the path, separating
        it from the previous directory with a semicolon.


### Git

* Follow the [Software Carpentry instructions](http://software-carpentry.org/bootcamps/setup.html) 
in the section titled "Git" to install Git.

* If you don't have a GitHub account, please go to <https://github.com/> and register for a free account.

* You'll have the choice of working with Git through either the command line, a GUI, or both, whatever 
you're comfortable with (we'll demonstrate both.) If you want to use the command line, no additional 
setup is required. The GUI we'll be using is SourceTree, which is free and supports Windows and OS X. The 
installer can be downloaded here: <http://www.sourcetreeapp.com/>. 

    When you first run SourceTree, enter the username and e-mail address you used for your GitHub account 
and 
click Next. For the remaining settings in the setup wizard, the default options will work fine.


### Make

* **Linux**: you should already have Make installed.

* **OS X**: follow the instructions under GNU Make on the 
    [Software Carpentry instruction page](http://software-carpentry.org/bootcamps/setup.html)
    to install the Apple developer tools.

* **Windows**: you can download make from this link:
    <https://github.com/msysgit/msysgit/blob/master/bin/make.exe?raw=true>.
    Download this file and place it in the "bin" directory where you installed
    Git bash.

* To test: open a terminal, type "make", and press Enter. You should see 
    the following message:

        make: *** No targets specified and no makefile found.  Stop.

    This means that Make was successfully installed. Otherwise, you'll see this
    error message:

        bash: make: command not found


### Further resources

The above is enough preparation but here are some links if you are interested in reading a bit further.

* How to Use RStudio:
    - <http://www.rstudio.com/ide/docs/>
* RStudio Public Discussion & Troubleshooting Guide:
  - <http://support.rstudio.org/help/>
* How to Install R:
    - <http://cran.r-project.org/doc/manuals/R-admin.html>
    - <http://cran.stat.sfu.ca/doc/FAQ/R-FAQ.html#How-can-R-be-installed_003f>
* R FAQ:
    - <http://cran.r-project.org/doc/FAQ/R-FAQ.html>
* More about add-on packages in the R Installation and Administration Manual
     - <http://cran.r-project.org/doc/manuals/R-admin.html#Add_002don-packages>

     

<!-- Notes from an October 2012 workshop
  ["R carpentry - Finding Help"](../modules/r-carpentry-finding-help.html)
  (quite rough at this point) -->
  
<!-- we don't need this if we remove Q4c, right? -->
