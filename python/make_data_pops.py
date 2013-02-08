#!/usr/bin/python

'''
Make and save logistic growth population data
'''

from __future__ import division
import numpy as np

def logistic_growth(r, K, n0, p=0.1):
    # Set up population array
    n = np.zeros((npops, nt))
    n[:,0] = n0

    # Loop through all time steps
    steps = np.arange(1, nt)
    for i in steps:
        cat = (np.random.rand(npops) < p)  # Random catastrophe 10% of time
        ni = np.round(n[:,i-1] + r*n[:,i-1]*(1 - n[:,i-1]/K))
        ni[cat] = n0
        n[:,i] = ni
    
    return n

if __name__ == '__main__':
    nt = 15
    npops = 50
    result = logistic_growth(0.6, 100, 10)
    np.savetxt('data_pops.csv', result, '%i', ',')
