#!/usr/bin/env python

import sys
import argparse
import csv

import matplotlib.pyplot as plt


def read_mass_period(df):
    mass = []
    period = []

    with open(df, 'r') as f:
        csvr = csv.DictReader(f)

        for line in csvr:
            if line['mass']:
                mass.append(float(line['mass']))

            if line['period']:
                period.append(float(line['period']))

    return mass, period


def plot_data(df):
    figname = df[:-3] + 'png'

    mass, period = read_mass_period(df)

    fig, (m_ax, p_ax) = plt.subplots(2, 1)

    # plot mass histogram
    m_ax.hist(mass, bins=50)
    m_ax.set_xlabel('Mass ($M_{Earth}$)')

    # plot period histogram
    p_ax.hist(period, bins=50)
    p_ax.set_xlabel('Period (Days)')

    fig.savefig(figname)


def parse_args(args=None):
    d = ('Plot histograms of planet mass and planet period for each'
         'given CSV data file.')
    parser = argparse.ArgumentParser(d)
    parser.add_argument('datafiles', type=str, nargs='+',
                        help='CSV files with data to plot.')
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)

    for df in args.datafiles:
        plot_data(df)


if __name__ == '__main__':
    sys.exit(main())
