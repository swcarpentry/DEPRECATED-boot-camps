#!/usr/bin/python

# This script counts the number of Elk that were sighted
# The input to this script is a file name and an animal name

import sys

file = sys.argv[1]
animal = sys.argv[2]


newfile = open(file, 'r') 
filelines = newfile.readlines()
    
for line in filelines:
    d, t, a, n = line.split()
    if a == animal:
        print a, n


