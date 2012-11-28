import sys

import mean
from mean import read_file, get_mean

filename = sys.argv[1]
grades = read_file(filename)
mean = get_mean(grades)
print 'mean of ', filename, 'is', mean
