from matplotlib.pyplot import *
import numpy
from cPickle import load
from datetime import datetime

birdcount_by_day = load(open(sys.argv[1]))

plotme = []

for day in birdcount_by_day:  # note, iterating over dictionaries gives you keys
    date = datetime.strptime(day + ' 2013', '%B %d %Y')
    day_of_year = date.strftime('%j')
    
    # trick: we need to convert day_of_year into an integer
    day_of_year = int(day_of_year)
    
    # retrieve birdcount
    count = birdcount_by_day[day]
    
    # now add day_of_year and birdcount
    plotme.append((day_of_year, count))

plotme.sort()

plotme = numpy.array(plotme)

plot(plotme[:,0], plotme[:,1])
savefig('test.pdf')
