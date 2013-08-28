#!/usr/bin/python

import sys

x = []
if len(sys.argv) == 2:
	stream = open(sys.argv[1])
else:
	stream = sys.stdin
for line in stream:
	if line.strip() == "":
		continue

	try:
		num = float(line)
		x.append(num)
	except ValueError:
		pass

print "average is", sum(x) /float(len(x))
