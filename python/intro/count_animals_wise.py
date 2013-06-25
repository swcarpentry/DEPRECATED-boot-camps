#!/usr/bin/python

# This script counts the number of Elk that were sighted

import sys

file = sys.argv[1]

newfile = open(file, 'r') 
filelines = newfile.readlines()

elk_list = []
    
for line in filelines:
    d, t, a, n = line.split()
    if a == 'Elk':
        num = int(n)
        if num > 10:
            elk_list.append(line)

print elk_list

