#!/usr/bin/python

import sys

# The # indicates that this is a comment
# This is our toy function, it is suppose to return the integer 2
# no matter what value v is passed as an argument
def return_v(v):
    return v

# The simplest type of test is to print out the value and compare it by eye to the expected value
for line in sys.stdin:
    v = int(line)
    print "Return value of return_v(", v, "):", return_v(v)
