#!/usr/bin/env python

# The line above tells the shell to run this script with Python.
import sys

from meananimals import mean_animals

# sys.argv contains the arguments the user entered at the command line when
# calling this script. See more at http://docs.python.org/2/library/sys.html.
# Another great way of putting together command line interfaces for scripts
# is the argparse module: http://docs.python.org/2/library/argparse.html

# Try running this script by typing "./mean_sighted.py big_animals.txt Elk"
# at the command line.

if len(sys.argv) != 3:
    print 'Usage: mean_sighted.py <filename> <species>'
else:
    filename = sys.argv[1]
    species = sys.argv[2]
    print mean_animals(filename, species)
