import csv
from cPickle import dump
import sys

def read_birdlist2(filename):
    birdlist = []

    fp = file(filename, 'rb')
    reader = csv.reader(fp)

    n = 0
    for bird, state, day in reader:
        birdlist.append((bird, state, day))

    return birdlist

def make_birdcount_by_day(birdlist2):
    birdcount_by_day = {}

    for bird, state, day in birdlist2:
        birdcount_by_day[day] = 0

    for bird, state, day in birdlist2:
        birdcount_by_day[day] = birdcount_by_day[day] + 1

    return birdcount_by_day

birdlist = read_birdlist2(sys.argv[1])
d = make_birdcount_by_day(birdlist)

dump(d, open(sys.argv[2], 'w'))
