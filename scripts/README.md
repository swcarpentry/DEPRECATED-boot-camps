To run our simple pipeline:

    python make-big-birdlist.py long-birds.csv

    python make-birdcounts.py long-birds.csv long-birds.counts

    python plot-birdcounts.py long-birds.counts plot.pdf

These will make a random list of birds, states, and days; count the
number of birds seen each day; and plot the resulting distribution by
day of year.  The result will be placed in 'plot.pdf'.
