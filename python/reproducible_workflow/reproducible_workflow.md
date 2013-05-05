Creating a Reproducible Workflow
================================

Introduction
------------
We're now in the home stretch of the workshop - congratulations! Up to this
point, we've talked quite a bit about how to make your code efficient (good
programming practice), accurate (testing), and maintainable (modularization +
version control).

Now we're going to talk about a final and very important concept
known as reproducibility. This idea has a long tradition here at Stanford,
starting with the Earth Scientist Jon Claerbout and Prof. David Donoho at our
stats department, who wrote in a paper with Buckheit the following:

	 An article about a computational result is advertising, not
	 scholarship. The actual scholarship is the full software environment,
	 code and data, that produced the result.

One of Donoho's former student's Victoria Stodden has written extensively about
the idea of reproducibility in scientific software - you may want to look up
[some of her papers][1] for reference.

For our purposes, we can summarize the goal of reproducibility in two related 
ways, one technical and one colloquial.

In a technical sense, your goal is to __have a complete chain of custody (ie, 
provenance) from your raw data to your finished results and figures__. That is, 
you should _always_ be able to figure out precisely what data and what code 
were used to generate what result - there should be no "missing links".

If you have ever had the experience of coming across a great figure that you
made months ago and having no idea how in the world you made it, then you
understand why provenance is important. Or, worse, if you've ever been unable
to recreate the results that you once showed on a poster or (gasp) published in
a paper...

In a colloquial sense, I should be able to sneak into your lab late at night, 
delete everything except for your raw data and your code, and __you should be 
able to run a single command to regenerate EVERYTHING, including all of your 
results, tables, and figures in their final, polished form__. Think of this as 
the "push button" workflow. This is your ultimate organizational goal as a 
computational scientist. Importantly, note that this rules out the idea of 
manual intervention at any step along the way - no tinkering with figure axes 
in a pop-up window, no deleting columns from tables, no copying data from one 
folder to another, etc. All of that needs to be fully automated.

As an added bonus, if you couple this with a version control system that tracks 
changes over time to your raw data and your code, you will be able to instantly 
recreate your results from any stage in your research (the lab presentation 
version, the dissertation version, the manuscript version, the Nobel Prize 
committee version, etc.). Wouldn't that be nice?

To illustrate these ideas, we'll set up a small but realistic research project 
that follows a reproducible workflow. Just like you would in your own research 
projects, we'll go through the following key steps:

1.	Create a clear and useful directory structure for our project.
2.	Set up (and use) Git to track our changes.
3.	Add the raw data to our project.
4.	Write our code to perform the analysis, including tests.
5.	Push the button and watch the magic.

One other aspect of this is that someone else should be able to look at your
analysis and understand what you did and to some degree also why (informative
git commit messages can be helpful in that regard...). In the most common case,
the 'someone else' is you 6 or 12 or 18 months down the line. 

One final note - the workflow that we're following here is just a suggestion. 
Organizing code and data is an art, and a room of 100 scientists will give you 
101 opinions about how to do it best. Consider the below a useful place to get 
started, and don't hesitate to tinker and branch out as you get a better feel 
for this process. You also might want to review [William Noble's paper][2] on 
this topic for more ideas.

1.	Setting up the project directory
------------------------------------

Let's create a project (a reasonably self-contained set of code, data, and 
results to answer a discrete scientific question) that will analyze responses
in auditory receptors of grasshopper neurons.

We begin by creating a directory called `grasshopper` in a convenient place on
our hard drive. You might want to create a main directory called `Projects` or
`Research` in your home folder or in your Documents folder to hold the
directories for all of your individual research projects.

Now, within the `grasshopper` directory, create four subdirectories:

        .
	├── data
	├── doc
	├── results
	├── src

The `data` directory will hold all of the raw data associated with the project,
which in this case will be files containing spike-time data and a file
containing auditory stimuli. 

The `man` folder, short for manuscript, will (someday) contain the manuscript
that we'll write describing the results of our analysis (you were planning on
using version control for your manuscript too, weren't you?). The `results`
folder will contain the results of our analysis, including both tables and
figures, and the `src` directory will contain all of our code.

In a more complex project, each of these directories may have additional 
subdirectories to help keep things organized.

For bonus points, do this all from the command line.

2.	Initialize a Git repository
-------------------------------

Since we want to use version control to track the development of our project, 
we'll start off right away by initializing an empty Git repository within this 
directory. To do this, open a Terminal window, navigate to the main 
`grasshoppers` directory, and run the command `git init`.

As you add things to the project directory, and modify old things, you'll want 
to frequently commit your changes as we discussed in the Git tutorial.

3.	Add raw data
----------------

Often, we start a project with a particular data file, or set of data files.
In this case, we will look at a data set that is available on the [CRCNS
data-sharing website][crcns.org]. For the purposes of this tutorial, we will
look at a small subset of the data available for download
[here][http://arokem.org/data/grasshoppers.zip].

Download the data and move the files into the `data` subdirectory.

Now we reach an interesting question - should your `data` directory be placed 
under version control (ie, should you `git add` and `git commit` these files)? 
Although you might automatically think that this is necessary, in principle our 
raw data should never change - that is, there's only one version (the original 
version!), and it will never be updated. As a result, it's not necessarily 
useful to place this file under version control for the purpose of tracking 
changes to it.

A reasonable rule of thumb for getting started is that if the file is
realatively small, go ahead and commit it to the Git repository, as you won't
be wasting much hard disk space. Additionally, the file will then travel with
your code, so if you push your repository to Github (for example) and one of
your collaborators clones a copy, they'll have everything they need to generate
your results.

However, if your file is relatively large (the stimuli files here are about 3.5
MB) AND is backed up elsewhere, you might want to avoid making a duplicate copy
in the `.git` directory.

In either case, you'll want to ensure that every one of your data files has
some sort of metadata associated with it to describe where it came from, how it
got to you, the meaning of the columns, etc. There are many formats for
metadata that vary from simple to very complex. At a minimum, make sure to create a
`README.txt` file that describes your data as best you can. This is also
helpful because empty directories do not get 

Copy and paste the text below into a `README.txt` file and place it in the data 
subdirectory. Remember that this is a bare-bones description - in your own 
work, you'll want to include as much information as you have.

	Data from grasshopper auditory receptor recordings, donwloaded from the
	CRCNS ia-1 data-set. Includes files with spike times and files with
	auditory stimuli. For full methods, refer to Rokem et al. J Neurophys (2006). 

At this point, your project directory should look like this:

	.
	├── data
	│   ├── README.txt
	│   ├── grasshopper_spike_times1.txt
	│   ├── grasshopper_spike_times2.txt
	│   ├── grasshopper_stimulus1.txt
	│   ├── grasshopper_stimulus2.txt
	├── doc
	├── results
	├── src

Add the readme file to your git repository.

What about the case in which your raw data is hosted elsewhere, on a SQL 
server, for example, or a shared hard drive with your lab? Now your data is 
living somewhere else, and you don't necessarily have direct control over its 
provenance (what if someone changes it while you weren't looking?). In this 
situation, you should try to make your `runall.py` script (see below) make a 
copy of the metadata associated with the dataset (it does have metadata, 
doesn't it?), which hopefully will include something like a version number and 
a last-updated date, and store this along with your results. That way you'll at 
least have some information on the version of the data that was used. If 
there's no metadata, try to shame your collaborators into creating some. If all 
else fails, at least record the date on which your analysis was run so that, in 
principle, you could later try to find out what state the raw data was in on 
that date. If you're really nervous about the data changing, though, you might 
want to look into making yourself a local copy.

4. Write code to perform analysis
---------------------------------

Now for the real work - writing the code that will perform our analysis. We'd
like to generate two kinds of outputs. First, we want to compute the
spike-triggered average stimulus (or STA) for each stimulus and plot it. For
each STA, we would like to save a figure, with an informative file-name. 

#### Modules and tests

We have just finished writing code that will read the data from files and
compute the STA.

We still need to test that code. Let's use the fact that we have one naive (but
slow) way of computing the STA and another more sophisticated (but less
obvious) way of computing the STA and compare them to each other. Though not
completely PhD-proof, it stands to argue that if you get the same result in
these two distinct ways, you might be getting it right.

Let's start with the test file. We'll call it `test_analysis.py` and it will contain the
following code, which we have copied from our data_analysis_and_exploration
notebook:

    #!/usr/bin/env python

    import numpy as np
    import numpy.testing as npt # We'll need this for testing
    import matplotlib.mlab as mlab # We'll use this to read the stimulus file

    from files import get_data
    from analysis import sta1, sta2, volt2dB

    def test_sta():
       """
       Test that the STA comes out to be the same when computed two different ways
       """

       spikes, header = get_data('../data/grasshopper_spike_times1.txt')

       rec_arr = mlab.csv2rec('../data/grasshopper_stimulus1.txt', delimiter=' ',
                              names = ['time', 'volt'])

       est_sta1 = sta1(spikes, rec_arr['time'],
                        volt2dB(rec_arr['volt'], header['intensity (dB)']))

       est_sta2 = sta1(spikes, rec_arr['time'],
                        volt2dB(rec_arr['volt'], header['intensity (dB)']))

       # Note - we are using numpy's testing module here: 
       npt.assert_array_equal(est_sta1, est_sta2)

This will, of course fail immediately if we run it at this point, because none
of the other code has been filled in yet. But in a sense, this provides us with
a "contract" of what we expect to be there.

Next, let's write the two files: `analysis.py` and `files.py`. In 'analysis.py'
I put the following:

    import numpy as np

    def volt2dB(volt, max_dB):
        """ 
        Convert from voltage to dB SPL

        Parameters 
        ----------
        volt : 1d array
            The stimulus in volts 

        max_dB : int or float 
            The maximal dB value presented in the experiment

        Notes
        -----
        `max_dB` can be gleaned from the header information, see `get_data`
        """
        stim = 20  * (np.log10(volt / 2.0e-5))
        return max_dB - stim.max() + stim


    def sta1(spike_times, time_arr, stim_arr, num_bins=200): 
        """ 
        Calculate the STA by looping over the spikes

        Parameters
        ----------
        spike_times : 1d array
            Times of occurence of action potentials (micro-seconds)

        time_arr : 1d array 
            The time-bins in the stimulus presentation

        stim_arr : 1d array
            The stimulus values (in dB) at each time bin

        num_bins : int
            The number of time-bins to calculate the STA for. 

        Return 
        ------
        STA : 1d array
            The spike-triggered average stimulus 

        """

    def sta1(spike_times, time_arr, stim_arr, num_bins=200): 
        """ 
        Calculate the STA by looping over the spikes

        Parameters
        ----------
        spike_times : 1d array
            Times of occurence of action potentials (micro-seconds)

        time_arr : 1d array 
            The time-bins in the stimulus presentation

        stim_arr : 1d array
            The stimulus values (in dB) at each time bin

        num_bins : int
            The number of time-bins to calculate the STA for. 

        Return 
        ------
        STA : 1d array
            The spike-triggered average stimulus 

        """
        result = []
        # We enumerate the spike times:
        for spike_idx, time in enumerate(spike_times):
            idx = np.where(time_arr == time)[0]
            # Only use spikes with sufficient 'history':
            if (idx - num_bins) > 0:
                result.append(stim_arr[idx-num_bins:idx])
        return np.mean(result, 0)

    def sta2(spike_times, time_arr, stim_arr, num_bins=200):
        """ 
        Calculate the STA by allocating a binary array

        Parameters
        ----------
        spike_times : 1d array
            Times of occurence of action potentials (micro-seconds)

        time_arr : 1d array 
            The time-bins in the stimulus presentation

        stim_arr : 1d array
            The stimulus values (in dB) at each time bin

        num_bins : int
            The number of time-bins to calculate the STA for. 

        Return 
        ------
        STA : 1d array
            The spike-triggered average stimulus 

        """
        raster = np.zeros(time_arr.shape)
        spike_idx = spike_times/time_arr[1]
        # Indice arrays need  to be integers:
        spike_idx = spike_idx.astype(int)
        raster[spike_idx] = 1
        idx = np.where(raster==1)
        idx_w_len = np.array([idx[0] - count for count in range(num_bins, 0, -1)])
        return np.mean(stim_arr[idx_w_len], 1)

And in `files.py` let's put the following:

    import numpy as np

    def get_data(fname): 
        """ 
        Read spike time data and header information from file

        Parameters 
        ----------
        fname : str
            Path to the file 

        Returns
        -------
        spikes : 1-d array
            Contains the spike times in the order they were 
            recorded in the file

        header : dict
            Header information from the data file
        """
        f = open(fname, 'r')
        header = {}
        l = f.readline()
        while l.startswith('#'):  # Same as l[0]=='#' 
            line_parts = l.split(':')
            # We split the key again, so that we don't get the '#'
            try:
                header[line_parts[0].split('# ')[-1]] = float(line_parts[1])
            except:
                header[line_parts[0].split('# ')[-1]] = line_parts[1]
            l = f.readline()

        spikes = np.loadtxt(fname)
        f.close()
        return spikes, header

This should make the test now pass. Just to be sure we did everything right, go
ahead and run `nosetests` from the Terminal and make sure that the functions
pass. 

Hmm. This doesn't work. Maybe we have a mistake? 

It turns out that we have neglected to ignore the spikes that occur in bins
before num_bins. Good things we have tests in place!

When we add the line:

     raster[spike_idx[np.where(spike_idx<num_bins)]] = 0

The tests now pass. Convince yourself that you understand why this does what
it's supposed to be doing

Now that everything passes, we are happy and quite confident that we can read files
and calculate STA without error.

Our project should now look like this:

	.
	├── data
	│   ├── README.txt
	│   ├── grasshopper_spike_times1.txt
	│   ├── grasshopper_spike_times2.txt
	│   ├── grasshopper_stimulus1.txt
	│   ├── grasshopper_stimulus2.txt
	├── doc
	├── results
	├── src
	│   ├── analysis.py 
	│   ├── files.py 
	│   ├── test_analysis.py 

This would be a good time to add all these new source files to our version control and
make a commit. 

Note, of course, that this is not the normal workflow for these steps. Normally, 
you'd spend days/weeks/months working in the `src` directory, writing code and 
tests, generating results, looking at the results, writing new code and tests, 
generating new results, etc. This iterative cycle isn't unlike writing a 
paper - you spew out a draft that's not too bad, then go back and revise it, 
then spew out some new material, revise that, etc.

Different people have different favorite approaches and tools for this 
iterative cycle.

One strategy is to simultaneously work on three files at once - a module like 
`analysis.py`, a file to test the functions in your module like 
`test_analysis.py`, and a third script that calls various functions from 
your module and runs them so that you can see whether your code is doing what 
you want it to do. I sometimes call this `scratch.py` or something like that, 
and fill up my Terminal window with hundreds of lines of `python scratch.py` as 
I modify my code and look at the results that are saved, results that are 
printed to the Terminal window, and errors that pop up.

Another strategy is to take advantage of the IPython notebook, where you can 
write your code in individual cells and easily execute each one sequentially, 
make sure each cell executes properly, and review the values of variables after 
each step. This can be great and efficient in the event that you really don't 
have any idea what your final code will look like. The downside here is that, 
at the end, you'll be left with one enormous notebook file (probably without 
unit tests), and you'll need to go back at the end to properly modularize your 
code into functions and separate files, similar to the structure that we're 
using in this exercise, so that you can have a fully reproducible workflow. 
Plus you may end up writing your unit tests at the end (you are going to write 
them, aren't you?) rather than iteratively with your code as you go. All that 
said, though, this is a great strategy if you think you need to feel your way 
around for a while.


#### The runall script

Now that we have our core functions and tests in place, it's time to create the 
"button" for our push-button workflow - the `runall.py` script. The idea is 
that you will be able to start with an empty results directory, execute the 
line `python runall.py` in Terminal, and have our table and figures saved in the 
`results` directory.

Create a new text file called `runall.py` and copy and paste the following code 
into it.

    #!/usr/bin/env python
    import os

    import numpy as np
    import matplotlib.mlab as mlab
    import matplotlib.pyplot as plt

    from files import get_data
    from analysis import sta1, sta2, volt2dB

    if __name__=="__main__":
        data_path = '../data/'
        results_path = '../results/'
        data_list = os.listdir(data_path)
        for file_num in range(1, len(data_list)/2+1):    
            spikes, header = get_data(
                '../data/grasshopper_spike_times%s.txt'%file_num)

            rec_arr = mlab.csv2rec('../data/grasshopper_stimulus%s.txt'%file_num,
                                   delimiter=' ',
                                   names = ['time', 'volt'])

            sta = sta2(spikes, rec_arr['time'],
                        volt2dB(rec_arr['volt'], header['intensity (dB)']))


            fig, ax = plt.subplots(1)
            ax.plot(sta)
            ax.set_xlabel('Time (before spike, msec)')
            ax.set_ylabel('Amplitude (dB)')
            xticks = np.array([0,50,100,150,200])
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(s*50./1000.) for s in xticks[::-1]])
            fig.savefig(os.path.join(results_path, 'figure%s.png'%file_num))


We won't go over this code in too much detail. Note just that the entire code
is inside the `if __name__ == "__main__"` conditional. This is useful in cases
in which a file might be run directly by the python interpreter (as we did
here), or imported from other modules. Only when the file is called directly
from the python interpreter (`python runall.py`) is the `__name__` variable set
to be `'__main__'`, so that block of code only gets executed in that case. 

In this script we have set up variables that define the locations of the `data`
and `results` directories (relative to the `src` directory where our
`runall.py` file is located) as well as the names of our input data file and
the table and figure that we will create. Note that the code discovers how many
data files are in the directory, so when you add more data files, they will
automatically be incorporated into the analysis, without any need to edit this
code. 

Now, go back to your Terminal window, navigate to the `src` directory, and run 
the command `python runall.py`.  Figures should appear in the `results` directory.

Don't forget to add `runall.py` to your git repo.

5. Run the push button analysis
-------------------------------

Now with everything in place, we're ready for the magic. Just for good measure, 
delete any files that are hanging around in your `results` directory. Then, 
execute `python runall.py` from your `src` subdirectory and marvel at your 
fully reproducible workflow! (If you've been making a lot of changes to your 
code, and aren't quite sure what's in your `results` directory, you may want to 
periodically clear out this folder and re-run everything to make sure that 
everything is regenerating properly.)

At this point, your directory should look like the below.

	.
	├── data
	│   ├── README.txt
	│   ├── grasshopper_spike_times1.txt
	│   ├── grasshopper_spike_times2.txt
	│   ├── grasshopper_stimulus1.txt
	│   ├── grasshopper_stimulus2.txt
	├── doc
	├── results
	│   ├── figure1.png
	│   ├── figure2.png
	├── src
	│   ├── analysis.py 
	│   ├── files.py 
	│   ├── test_analysis.py
	│   ├── runall.py 

	
At this point, a natural question to ask is whether you need to add the 
contents of your `results` directory to your git repository. The answer should 
be obvious - you do not need to do this, since the files in your `results` 
directory contain no unique information on their own. Everything you need to 
create them is contained in the `data` and `src` directories. One exception to 
this, though, might be if your analysis takes a very long time to run and the 
outputs are fairly small in size, in which case you may want to periodically 
commit (so that you can easily recover) the results associated with 
"intermediate" versions of your code.

While some of your projects might be nearly this simple, most will probably be
more complex, sometimes significantly so. You will eventually come across the
need to deal with modules that are shared across multiple projects, running the
same analysis on multiple sets of parameters simultaneously, running analyses
on multiple computers, etc. While we don't have time to go into these extra
bits in detail, feel free to ask the instructors about any specific issues that
you expect to encounter in the near future.

Good luck with your reproducible research!


[1]:
http://www.stanford.edu/~vcs/Papers.html

[2]: 
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424

[3]: http://knb.ecoinformatics.org/software/eml/

[4]: http://knb.ecoinformatics.org/morphoportal.jsp
