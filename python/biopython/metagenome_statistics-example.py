#!/usr/bin/env python
'''This script retrieves a metagenome_statistics data structure from the MG-RAST API and
plots a graph using data from the web interface'''

import urllib, json, sys
import numpy as np

# retrieve the data by sending at HTTP GET request to the MG-RAST API
ACCESSIONNUMBER = "mgm4440613.3"   # this is a public job
some_url = "http://api.metagenomics.anl.gov/api2.cgi/metagenome_statistics/%s?verbosity=full" % ACCESSIONNUMBER
jsonobject = urllib.urlopen(some_url).read()

# convert the data from a JSON structure to a python data type, a dict of dicts.
jsonstructure = json.loads(jsonobject)

# get the elements of the data that we want out of the dict of dicts..
spectrum = np.array( jsonstructure["qc"]["kmer"]["15_mer"]["data"], dtype="float")
lengthdistribution = np.array( jsonstructure["length_histogram"]["upload"], dtype="int")

# display the data table
np.savetxt(sys.stderr, spectrum, fmt="%d", delimiter="\t")

# plot the length distribution graph
import matplotlib.pyplot as plt
plt.plot(lengthdistribution[:, 0], lengthdistribution[:, 1])
plt.xlabel("length")
plt.ylabel("number of reads")
plt.show()
