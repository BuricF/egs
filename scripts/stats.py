#! /usr/bin/python -tt


# This file is part of the EGS - Extended Gray-Scott Model simulation package
#
# Copyright (C) 2014 Filip Buric
#
# EGS is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# EGS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 


""" 
Perform various statistics (measures) on the time evolution 
of the reaction-diffusion system.
Data is read from the simulation output (.dat files).

Usage:
    
stats.py - Compute statistics on data series and plot results in image files.
stats.py DIR - Additionally save statistics in DIR for easier subsequent use.
"""

__author__ = "Filip Buric"
__date__ = "19 June 2014"


import matplotlib.pyplot as plt
import numpy as np
import math
import sys


class StatsCollector:
    """ 
    The StatsCollector class encasulates statistic data series computed 
    for simulation data.

    Its usage is to invoke data collection on each system generation, after
    selecting the desired statistics function.
    These functions are defined outside the class and bound to the statistic
    member reference at object initialisation.
    
    Plotting is also done by this class since it is statistics-dependent,
    though title, axes labels, and saving is done externally.
    """

    def __init__(self, stats_function, no_series=1, output_file = ''):
        """Initialize with relevant statistics function,
           number of data series and file to save data.
           If the filename is invalid, data series is output to stdout."""
        self.statistic = stats_function
        self.s = no_series
        self.outfile = None
        self.stat_data = []  # each line is concentrations vector c = (c_1, .. ,c_s)
                             # columns are the time evolution of each component of these vectors

        if output_file != '':
            try:
                self.outfile = open(output_file, 'w')
            except IOError:
                print('[Error] invalid output file specifed. using stdout..')
                self.outfile = sys.stdout


    def collect(self, data):
        """ Collect statistics for given generation."""
        values = self.statistic(data)
        self.stat_data.append(values)

        if self.outfile:
            for v in values:
                self.outfile.write(str(v) + '\t')

            self.outfile.write('\n')


    def plot(self):

        if self.outfile:
            self.outfile.close()

        plt.figure()

        size = len(self.stat_data)
        stats = np.array(self.stat_data)

        colors = ['k', 'g', 'r', 'b']

        for x in xrange(self.s):
            plt.plot(range(1, size+1), stats[:, x], color=colors[x+1], label='V'+str(x+1), )

        plt.legend()



def crest_spread(data):
    """ Measure the amount of "crest"/"peak" values. (See def. in thesis report)
        for one generation

        data is a list of species conentration matrices

        returns a list of total peak count for each species
    """

    S = len(data)     # number of species
    N = len(data[0])  # size of matrices

    crest = [0.0 for x in xrange(S)]

    # for each species
    for sp in xrange(S):

        # compute mean concentration
        mean = sum([sum(row) for row in data[sp]]) / float(N*N)

        # compute deviation
        dev = sum([sum([abs(x-mean) for x in row]) for row in data[sp]]) / float(N*N)

        # count cells in desired range
        count = 0
        for row in data[sp]:
            for c in row:
                if c > mean + dev:
                    count +=1

        # record percentage of found cell
        crest[sp] = count / float(N*N)

    return crest


def total_concentration(data):
    """ Measure the total concentrations of chemical species in one generation,
    normalized by grid size N^2."""

    S = len(data)                               # number of species
    N = len(data[0])                            # size of matrices

    conc = [0.0 for x in xrange(S)]            # total concentrations

    for sp in xrange(S):                       # for each species
        for i in xrange(N):                    # sum up matrix values line by line
            conc[sp] += sum([float(x) for x in data[sp][i]])

    return [c/float(N*N) for c in conc]



def main():

    print "\n[ Statistics ]"

    datafile = [
                open('simdata/catalyst1.dat', 'r'),
                open('simdata/catalyst2.dat', 'r'),
                open('simdata/catalyst3.dat', 'r')
                ]

    S = len(datafile)   # get no. of species

    if (len(sys.argv)>=2):
        # detected saving path command line argument
        crests = StatsCollector(crest_spread, S, sys.argv[1]+'/crests.dat')
        totconc = StatsCollector(total_concentration, S, sys.argv[1]+'/totconc.dat')
    else:
        crests = StatsCollector(crest_spread, S)
        totconc = StatsCollector(total_concentration, S)

    # get size and pass line
    N = int(datafile[0].readline())
    for i in xrange(1, S):
        datafile[i].readline()

    done = False
    iter = 0
    while True:

        # feedback
        if iter % 100 == 0:
            sys.stdout.write('\rgen %d' % iter)
            sys.stdout.flush()

        matrices = []     # concentration matrices for one generation

        for species in xrange(S):

            m = []
            try:
                for i in xrange(N):

                    vals = (datafile[species].readline()).split()

                    if not vals:        # missing line in one file: terminate
                        raise Exception('missing line')

                    m.append([float(x) for x in vals])

                vals = datafile[species].readline()      # pass gap

                if not vals:
                    raise Exception('missing line')

                matrices.append(m)

            except Exception:
                done = True
                break

        if done:
            break

        # Call collectors on matrices
        crests.collect(matrices)
        totconc.collect(matrices)

        iter += 1


    # Plot statistics gathered by collectors
    crests.plot()
    plt.title('Total wave crest length (spread)')
    plt.xlabel('generation')
    plt.ylabel('% cells')
    plt.savefig('simdata/crest_spread.png')

    totconc.plot()
    plt.title('Total concentration')
    plt.xlabel('generation')
    plt.ylabel('concentration')
    plt.savefig('simdata/total_conc.png')

    plt.show()

    for dfile in datafile:
        dfile.close()


if __name__ == '__main__':
    main()
