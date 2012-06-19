#!/usr/bin/python

import animals

import sys

filename = sys.argv[1]
animal = sys.argv[2]

try:
  mean_count = animals.main(filename, animal)
except:
  print "There were no", animal

print "The mean count of", animal, "is", mean_count