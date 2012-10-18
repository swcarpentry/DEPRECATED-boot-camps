#!/usr/bin/env python

import sys

from argparse import ArgumentParser


def subset_by_name(match, infile, outfile):
    with open(outfile, 'w') as outf:
        with open(infile, 'r') as inf:
            # transfer the column headings
            outf.write(inf.readline())

            for line in inf:
                name = line.split(',')[0]
                if match in name:
                    outf.write(line)


def parse_args(args=None):
    d = 'Output a CSV subset of the exoplanet data filtering on name.'
    parser = ArgumentParser(description=d)
    h = ('Only lines which contain this string in the name field will be'
         'included in the subset.')
    parser.add_argument('match', type=str, help=h)
    parser.add_argument('infile', type=str, help='Name of input CSV file.')
    parser.add_argument('outfile', type=str, help='Name of output CSV file.')
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
    subset_by_name(args.match, args.infile, args.outfile)


if __name__ == '__main__':
    sys.exit(main())
