#! /Users/t/anaconda/bin/python
import numpy
from pylab import *

num1 = numpy.loadtxt('numbers.dat', delimiter=',')
num2 = numpy.loadtxt('numbers2.dat', delimiter=',')
num3 = numpy.loadtxt('numbers3.dat', delimiter=',')
num4 = num3*2

plot(num1[:,0], num1[:,1], label='num1')
plot(num2[:,0], num2[:,1], label='num2')
plot(num3[:,0], num3[:,2], label='num3')
plot(num4[:,0], num4[:,2], label='num4')
legend()
savefig('make_figure.pdf')
print 'I have now saved make_figure.pdf'

