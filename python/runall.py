#!/usr/bin/env python

'''
Script to create all results for camera_analysis project.
'''

import numpy as np
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
from mean_sightings import get_sightings


# ---------------------------------------------------------------------------
# Declare variables
# ---------------------------------------------------------------------------

# Set paths to data and results directories. Note that this method of relative 
# paths only works on *nix - for Windows, see os.path module.
data_dir = '../data/'
results_dir = '../results/'

# Set name of data file, table, and figure
data_name = 'sightings_tab_lg.csv'
table_name = 'spp_table.csv'
fig_name = 'spp_fig.png'

# Set names of species to count
spp_names = ['Fox', 'Wolf', 'Grizzly', 'Wolverine']


# ---------------------------------------------------------------------------
# Perform analysis 
# ---------------------------------------------------------------------------

# Declare empty list to hold counts of records
spp_recs = []

# Get total number of records for each species
for spp in spp_names:
    totalrecs, meancount = get_sightings(data_dir + data_name, spp)
    spp_recs.append(totalrecs)

print spp_names
print spp_recs

# ---------------------------------------------------------------------------
# Save results as table 
# ---------------------------------------------------------------------------

# Put two lists into a recarray table
table = np.array(zip(spp_names, spp_recs),
                 dtype=[('species', 'S12'), ('recs', int)])

# Save recarray as csv
ml.rec2csv(table, results_dir + table_name)


# ---------------------------------------------------------------------------
# Save results as figure 
# ---------------------------------------------------------------------------

# Set up figure with one axis
fig, ax = plt.subplots(1, 1)

# Create bar chart: args are location of left edge, height, and width of bar
ax.bar([0,1,2,3], spp_recs, 0.8)

# Place tick marks in center of each bar
ax.set_xticks([0.4, 1.4, 2.4, 3.4])

# Set limits to give some white space on either side of first/last bar 
ax.set_xlim([-0.2, 4])

# Add species names to x axis
ax.set_xticklabels(spp_names)

# Save figure
fig.savefig(results_dir + fig_name)
