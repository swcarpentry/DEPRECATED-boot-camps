Installing Python packages; useful Python packages
--------------------------------------------------

(You will want to do a 'git pull origin master' in your 2012-11-scripps
directory to get the latest data.)

----

The easiest way to install Python packages is to use pip, the Python
package installer.

pip and Anaconda on OS X
~~~~~~~~~~~~~~~~~~~~~~~~

If you're using anaconda on OS X, you have pip installed; but you will
need to refer to the python package installer by a full path name::

   ~/anaconda/bin/pip

More generally, you'll need to prefix many of the commands below with
'~/anaconda/bin/'.  You can set this as a default for your current
shell by doing::

   export PATH=~/anaconda/bin:$PATH

(or you can add that command to the file ~/.bashrc using nano.  Ask a TA
for help!)

pip on VirtualBox
~~~~~~~~~~~~~~~~~

If you're running on VirtualBox, you may need to install pip first
by doing::

   sudo apt-get install python-pip

at the terminal.

(If you're running any other way, you should already have pip installed.)

----

OPTIONAL: Using virtualenv
~~~~~~~~~~~~~~~~~~~~~~~~~~

Down the road, if you're running on a machine where you don't have
sysadmin access, you can use a python package called 'virtualenv' to
set up your own installation of Python into which you can install your
own packages.  Once virtualenv is installed (by a sysadmin,
presumably) It's as simple as ::

   python -m virtualenv NAME

where NAME is the name of your workspace, e.g. ::

   python -m virtualenv env

followed by ::

   . env/bin/activate

From that point on, you will be able to use pip to install things
within this workspace, and Python (again from within that workspace)
will be able to access and use those installed packages.

Python packages
~~~~~~~~~~~~~~~

There are, literally, *thousands* of Python packages.  The basic deal
is this: Python comes with "batteries included", which means that you
can do amazing numbers of things with just a basic Python install.
The anaconda install and VirtualBox virtual machine come with tons
*more* stuff.  But there's always the need to use an updated version
of something, or a little package that someone wrote that addresses
*just your concern*... so you'll always need to install stuff.

Here's how to install and use some potentially useful packages from 
my lab, but there's a whole world of Python packages out there.
See http://docs.python.org/2/library for packages that come included
with Python, and http://pypi.python.org/pypi for the Python package
index for third-party packages.

screed
~~~~~~

Screed is a little Python package from Titus's lab that reads in
DNA sequences -- more explicitly, it's a FASTA and FASTQ parser.
You can see some documentation here:

   http://screed.readthedocs.org/en/latest/

But how do *you* use it?

To install screed directly from github, do::

   pip install git+https://github.com/ged-lab/screed.git

*Using screed:*

screed can read FASTA and FASTQ files, as well as gzip or bzip2 versions
of those files.  For example, in the python directory there is a file
called '25k.fq.gz'; check it out::

   gunzip -c 2012-11-scripps/python/25k.fq.gz | less

(type 'q' to get out of less, and space bar to scroll through the file.)

**Note:**

All of the below screed commands are in the `using-screed.ipynb notebook <http://nbviewer.ipython.org/urls/raw.github.com/swcarpentry/2012-11-scripps/master/python/using-screed.ipynb>`__.

screed, in a nutshell, lets you read in all that data and access it
in Python. Try::

   import screed
   for record in screed.open('/path/to/2012-11-scripps/python/25k.fq.gz'):
      print record.name
      print record.sequence
      print record.accuracy
      break

A couple of points here.

First, there are 25,000 sequences in this file.  You might want to avoid
printing them all out (hence the 'break' command at the end of the loop!)

Second, you can use this for short read data or genomic sequences or
whatever.  We've mostly designed it for short-read data but it works
fine for genome-scale data (which is, after all, rather smaller than
most short-read data...)

Third, you can open any kind of sequence file with this command.

This can be a simple and handy way to extract a particular sequence
from a large file -- ::

   for record in screed.open('/path/to/2012-11-scripps/python/25k.fq.gz'):
      if record.name == '@895:1:4:1596:8538/2':
         break

   # do stuff with record

You can even pull out a list::

   list_of_names = ['@895:1:4:1596:8538/2', '@895:1:4:1596:6003/2']
   list_of_records = []

   for record in screed.open('/path/to/2012-11-scripps/python/25k.fq.gz'):
      if record.name in list_of_names:
         list_of_records.append(record)

   # do stuff with list_of_records

(You might want to use a 'set' here, note.)

So how is this stuff useful!?

Well, here's one simple example -- ::

   n = 0.
   m = 0.
   for record in screed.open('/path/to/2012-11-scripps/python/25k.fq.gz'):
      n += len(record.sequence)
      m += record.sequence.count('G') + record.sequence.count('C')

   print '%.3f G/C content' % (m / n,)

You can also do your quality trimming, or analysis of the first bases,
or... whatever.

Another example -- ::

   outfp = open('out.fa', 'w')
   for record in screed.open('/path/to/2012-11-scripps/python/25k.fq.gz'):
      outfp.write('>%s\n%s\n' % (record.name, record.sequence))

This converts FASTQ to FASTA.

(Does anyone want to see random access?)

blastparser
~~~~~~~~~~~

blastparser is another little Python package from Titus's lab
that reads in BLAST output and makes it accessible to Python.
This is really the only documentation :).

To install blastparser directly from github, do::

   pip install git+https://github.com/ged-lab/blastparser.git

blastparser is both less mature and more complicated to use than
screed, because BLAST files are more complicated than FASTA files.

Before we move forward, let's look at a BLAST output file -- check out
2012-11-scripps/python/sample-blast.txt::

   less python/sample-blast.txt

Each query is a record; each record has a bunch of hits; each hit has
a bunch of matches!

Here's how blastparser does it::

   import blastparser
   fp = open('python/sample-blast.txt')
   for record in blastparser.parse_fp(fp):
       for hit in record.hits:
           for match in hit.matches:
               print record.query_name, hit.subject_name
	       print match.subject_start, match.query_start
	       print match.subject_end, match.query_end
       break

A few things to cover --

 * figuring out what is stored in each object
 * print out to csv
