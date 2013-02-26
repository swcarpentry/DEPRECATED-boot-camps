#! /usr/bin/env python
import sys
import random, datetime

birdnames = """\
Common Teal
Long-tailed Duck
Little Grebe
Albatross
Flamingo
Heron
Stork
Shoebill
Vulture
Bald Eagle""".splitlines()

states = """Kentucky
Missouri
Kansas
Maine
Michigan
New York""".splitlines()

def generate_random_day(year):
    doy = int(random.gauss(200, 50))
    doy = min(doy, 365)
    doy = max(doy, 1)
    datestr = "%s %s" % (doy, year)
    date = datetime.datetime.strptime(datestr, "%j %Y")
    new_date = date.strftime("%B %d")
    return new_date.lower()

def generate_random_observation():
    bird = random.choice(birdnames).lower()
    state = random.choice(states).lower()
    day = generate_random_day(2013)
    return bird, state, day

import csv
fp = file(sys.argv[1], 'wb')
w = csv.writer(fp)

for i in range(10000):
    row = generate_random_observation()
    w.writerow(row)
fp.close()
