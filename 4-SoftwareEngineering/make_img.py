import numpy as np
import pylab

#read in the image
img = np.genfromtxt('spec_example.dat')
#Display the image
pylab.imshow(img, origin = 'lower', interpolation = 'nearest')
pylab.xlim(0, 150)

#Collapse the spectrum along x axis
img_collapse = np.sum(img, axis = 1)
#create and array to plot against
y = np.arange(np.shape(img_collapse)[0])

#Make a different figure
#pylab.figure()
#pylab.plot(img_collapse, y)

#--------OR---------
#Plot on the same figure
#Create a twin y  axis
pylab.twiny()
#Plot to new axis
pylab.plot(img_collapse, y, 'k', lw = 2)
pylab.xlim(0, 60000)
