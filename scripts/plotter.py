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
    Plot the superimposed species concentration matrices at each generation
    and save this as a movie file.
    Optionally, plot the matrices at a single generation as an image file.
    
    While configured, the plotting of the U concentrations is disabled.
    See comments in code to enable it.
    The code is meant to be easily extendable to more species.
    
    To indicate current progress, the plotter will output when it reaches
    every 100th generation to not delay script with constant printing.
    
    Note that the plotter overwrites previous output files (with the same name).
    
    Usage:
        plotter.py - Execution without parameters plots the entire data series.
        plotter.py final - Plots only the final concentration matrices (state),
                           from the .final files.
        plotter.py GEN   - Plots only the specified generation GEN.
        plotter.py GEN1 GEN2 - Plots the matrices from generation GEN1 up to and
                               including GEN2 as a single files.
"""

__author__ = "Filip Buric"
__date__ = "19 June 2014"


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys


def main():
    
    print "\n[ Plotting system evolution ]"

    resolution = 100        # DPI
    skip_info = False       # to skip plotting colorbars, and title

    start_gen = 0           # at which generation to start plotting

    plotfinal = False       # flag if only final state is to be plotted: only outputs image, not movie file
    plotgen = []            # generation interval to optionally plot as single files

    if len(sys.argv) > 1:
        if sys.argv[1] == "final":
            plotfinal = True

        else:
            try:
                if len(sys.argv) > 2:
                    plotgen = [int(sys.argv[1]), int(sys.argv[2])]
                else:
                    plotgen = [int(sys.argv[1]), int(sys.argv[1])]
                
                print 'gen. interval: [{0}, {1}]'.format(plotgen[0], plotgen[-1])

            except:
                print '[error:] invalid numbers as generations to plot. Plotting only first generation.'
                plotgen = [0]

    if not plotfinal:
        datafile = [open('simdata/fuel.dat', 'r'),
                    open('simdata/catalyst1.dat', 'r'),
                    open('simdata/catalyst2.dat', 'r'),
                    open('simdata/catalyst3.dat', 'r')
                    ]
    else:
        datafile = [open('simdata/fuel.final', 'r'),
                    open('simdata/catalyst1.final', 'r'),
                    open('simdata/catalyst2.final', 'r'),
                    open('simdata/catalyst3.final', 'r')
                    ]

    # Define the colormaps for each species
    colmap = []
    colmap.append(mpl.colors.LinearSegmentedColormap.from_list('my_cmap1',['white','black'],256)) # U
    colmap.append(mpl.colors.LinearSegmentedColormap.from_list('my_cmap2',['white','green'],256)) # V1
    colmap.append(mpl.colors.LinearSegmentedColormap.from_list('my_cmap3',['white','red'],256))   # V2
    colmap.append(mpl.colors.LinearSegmentedColormap.from_list('my_cmap4',['white','blue'],256))  # V3


    # Set up plotting
    fig = plt.figure()
    plt.hold(True)

    # Default behaviour: plot everything as movie
    if not plotfinal and not plotgen:
        writer = animation.FFMpegWriter(fps=10, codec ='libx264')
        writer.setup(fig, 'simdata/composite.mp4', dpi=100)

    # Get matrix size and pass that first line in all files
    for dfile in datafile[1:]:
        N = int(dfile.readline())


    # Go through files and plot found matrices (generations) until EOF (unknown no. of entries)        
    gen = 0                 # generation 0 is the initial state
    halt = False

    while not halt:

        # feedback
        if gen % 100 == 0:
            sys.stdout.write('\rgen %d' % gen)
            sys.stdout.flush()

        # For each chemical species, plot data.
        # Skip U data. Append 0 to the list below to enable it.            
        # Note: need reverse order in species list for in-order V_i colorbars        
        for fno in [3,2,1]:

            # Read and plot (imshow()) one matrix (current generation) from each file
            C = []                                     # concentration matrix
            for i in xrange(N):                        # read N matrix lines
                line = datafile[fno].readline()
                if not line:                          
                    halt = True               # stop ploting if no more data
                    break

                # Check for less characters than required entries and stop if data malformed
                if len(line) < N:
                    print "Malformed data: row " + str(gen) + ":" + str(rowno) + " is too short. Aborting."
                    exit(1)
                
                # Skip current generation if unwanted (need to do it here, after traversing file lines)
                if (gen >= start_gen) and ((not plotgen) or ((gen >= plotgen[0]) and (gen <= plotgen[-1]))):

                    row = [float(x) for x in (line.rstrip()).split()]       # get numeric list
                    C.append(row)

            # Plot (if wanted generation and if data of correct size)
            if (gen >= start_gen) and ((not plotgen) or ((gen >= plotgen[0]) and (gen <= plotgen[-1]))):

                # if not plotfinal:
                if (not skip_info):
                    plt.title('gen: '+str(gen))

                if (len(C) == N):

                    C = np.array(C)
                    # Note!: It's important to explicitly specify vmin and vmax, because of weird default behaviour in matplotlib
                    #plt.imshow(C, origin='lower', cmap = colmap[fno], alpha=1.0/len(datafile), vmin=0.0, vmax=C.max(), interpolation = 'none')
                    plt.imshow(C, origin='lower', cmap = colmap[fno], alpha=0.33, vmin=0.0, vmax=C.max(), interpolation = 'none')               

                    if (not skip_info):
                        # To reduce colorbar stacking, plot U colorbar horizontally beneath the matrix plot
                        if fno == 0:
                            cbar = plt.colorbar(orientation='horizontal', pad=0.01, shrink=0.5)
                            cbar.set_label('[U]')
                        else:
                            cbar = plt.colorbar(orientation='vertical', shrink=0.5)
                            cbar.set_label('[V'+str(fno)+']')                                
                    else:
                        plt.axis('off')

            # Pass gap btw matrices if not at EOF already
            if not halt:
                line = datafile[fno].readline()
                if not line:
                    halt = True
                    break

        if (gen >= start_gen) and ((not plotgen) or ((gen >= plotgen[0]) and (gen <= plotgen[-1]))):

            if not plotfinal and not plotgen:
                # default behaviour: save movie frame
                writer.grab_frame()
                plt.clf()

            elif not plotgen:
                plt.title('final generation')
                plt.savefig('simdata/composite_final.png', dpi=resolution)

            else:
                plt.savefig('simdata/composite_'+str(gen)+'.png', dpi=resolution)
                plt.clf()


        if (plotgen and gen == plotgen[-1]):
            # reached end of given plotting interval
            break
        
        else:
            gen += 1


    # Clean up
    # Make sure to close the output movie file 
    if not plotfinal  and (not plotgen):
        writer.finish()
    
    for dfile in datafile:
        dfile.close()

    print


if __name__ == '__main__':
    main()
