#!/usr/bin/env python

import sys, os
from optparse import OptionParser

if __name__ == '__main__':
    usage  = "usage: skeleton.py [<options>] [<arguments>]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input",  dest="infile", default=None, help="Input filename")
    parser.add_option("-o", "--output", dest="outfile", default=None, help="Output filename")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=True, help="Verbose")
    (opts, args) = parser.parse_args()
    if len(args) == 0 :
        sys.exit("No arguments supplied!")
    argument = args[0]

    if opts.verbose: 
        sys.stdout.write("Argument: %s\n" % argument ) 
#     Maybe do something here

    if opts.verbose: 
        sys.stdout.write("Done. \n")
