from random import random, getrandbits
from types import CodeType as code

print "I'll exhibit random runtime failures..."
n = 1
while 0 < n:
    print "\r" + str(n)
    if random() < 0.001:
	# Generate some non-functioning code.
        exec code(0, 0, 0, 0, "hello thw", (), (), (), "", "", 0, "") 
        n = -1000
    n += 1

