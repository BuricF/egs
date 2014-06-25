 #!/bin/sh

##################################################################################
#
# This script controls the simulation
#
# Because simulation paramters are hardcoded, the program is recompiled 
# and the current main source file and dynamics header are saved with a timestap.
#
# ! Note ! By default, the script removes old simulation data and movies 
# so these should be saved. Disable this behaviour by setting:
#
#       cleanup_old = false
#
# Executions of the simulation and statistics scripts are timed with 'time'.
#
##################################################################################

echo ""


cleanup_old=true

if [ "$cleanup_old" = true ] ; then
    # Clean up old simulation data
    rm simdata/fuel.dat
    rm simdata/catalyst1.dat
    rm simdata/catalyst2.dat
    rm simdata/catalyst3.dat

    # Clean up olf movie files
    rm simdata/composite.mp4
fi


# Recompile simulator (abort on error)
set -e
g++ -O3 src/egs_cycle.cpp -o bin/egs_cycle -lgsl -lgslcblas -lm -Wall -Wextra -funsafe-loop-optimizations -xT -fipa-pta -ftracer -march=native -fopenmp


# Save a copy of the main source file and dynamics header in order to record simulation parameters
# Rathen inelegant hack - in lieu of a proper versioning system
timestamp="$(date +"%s")"
cp src/egs_cycle.cpp simdata/egs_cycle_"$timestamp".cpp
cp src/include/egs_cycle_dynamics.hpp simdata/egs_cycle_dynamics_"$timestamp".hpp


# Run simulation
echo ""
date
time bin/egs_cycle
set +e


###############################
# Perform various analysis on the simulation data and plot animations of the time evolution (continue on error)
###############################

# Statistics
#time python scripts/stats.py simdata/stats

# Email statistics
#echo "Simulation results for 3-hypercycle simulation" | mutt -a simdata/total_conc.png -s "sim results" -- <email>

# Plot animation
time python scripts/plotter.py

# View animations/plots
#vlc simdata/composite.mp4
#gwenview simdata/composite_final.png
