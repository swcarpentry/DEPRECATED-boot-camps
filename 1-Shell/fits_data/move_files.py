#!/usr/env/bin python

import shutil

with open('issues_to_file.txt') as f:
    for line in f:
        dirname, filename = line.split()
        shutil.move(filename, dirname)
