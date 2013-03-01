import sys

import animals

if len(sys.argv) != 3:
    print 'Usage: python mean_animals.py filename kind'
    sys.exit()

filename = sys.argv[1]
kind = sys.argv[2]

mean = animals.mean_animals(filename, kind)

print mean
