# EGS - Extended Gray-Scott Model simulation package

This package simulates and analyzes the evolution of a more general version
of the Gray-Scott reaction-diffusion model,
namely one with 3 autocatalytic chemical species
(i.e. a catalytic hypercycle as described in [ref](http://studentarbeten.chalmers.se/publication/203896-pattern-formation-and-chemical-evolution-in-extended-gray-scott-models)).


> **Note**
>
> This is a fork of the simulator created for my master's thesis project.
> The original repository is frozen to reflect the thesis and any future work will
> be done here.

## Requirements

The package was developed on Ubuntu Linux and should be relatively easy to use
on other UNIX (i.e. POSIX) systems.
The simulation control scripts are written in Bash but they are not necessary
for compiling or running the simulator. They simply automate certain tasks.

The source code is quite portable so it shouldn't be too difficult to
use on Windows with minor adjustments to things like hard-coded file paths.

* C++
  * [GCC](https://gcc.gnu.org) >= 4.6.3
  * [GSL](https://www.gnu.org/software/gsl) >= 1.15

* Python 2.7
  * [matplotlib](https://matplotlib.org) >= 1.2.1

## Using the package

There are two simulators provided:

* `egs_cycle` - Basic simulator in which the chemical species are defined only through
                the initial conditions. The simulation is not perturbed once started.
* `egs_cycle_evol` - Evolutionary simulation.
                     The chemical species are introduced at different times to study how
                     such perturbations affect the system at long times.

These may be compiled from `src/` and run without any arguments since the
parameters are hard-coded, except for the initial conditions.

The **PDE model** is defined in `src/include/egs_cycle_dynamics.hpp`.

To compile and run the simulators, run:

~~~
egs_run.sh
~~~

and, respectively:

~~~
egs_evol_run.sh
~~~

These scripts will also run the plotting and statistics scripts (see below). 


### Initial conditions 

The area and initial concentration of the chemical species is specified through a
`.cond` file. 
It is assumed that outside these specified
areas the grid holds constant [U] = 1 and [V_i] = 0 concentrations.

These files are in the `cond/` subdirectory.
The `.cond` file is hard-coded in the simulator source code.

Simulation data is saved to the `simdata/` directory and this is where the scripts
assume the data is located and where their own output is saved.  On POSIX
systems, this directory may be replaced with a symbolic or hard
link should the data need saving on a larger drive.

### Visualization and statistics

Two scripts are included to visualize the simulation outcomes:

* `scripts/plotter.py` -  Plots the grid at each output generation and save these 
                          as a movie file, or, optionally, as separate frames.
                          The script provides documentation on options.

* `scripts/stats.py` - Plots statistics (e.g. total concentrations over simulation time) 
                       collected from the simulation output `.dat` files (see script for options).

These may be called without any arguments and the default behavior is plotting
data (movie, statistics) for the entire length of the simulation.

Minimal Doxygen documentation has been generated for the simulators (start at
`doc/html/index.html`) though the code is meant to be self-documenting and
to be inspected/modified when in use, at least for setting parameters.


## Notes on implementation

The various design decisions (like hard-coding parameters and recompiling for each run)
were made to prioritize speed of execution, while keeping the architecture modular and
extensible.

### GSL

The program requires the GNU Scientific Library (GSL) for random number
generation.  GSL version 1.15 was used during development but the few functions
employed should exist in much earlier version. The compiler flags to link with
the library are:

~~~
-lgsl -lgslcblas -lm
~~~

### OpenMP

Optionally, portions of the code are parallelizable with the
[OpenMP](http://www.openmp.org/) library, which should be included with your compiler.
The program has been compiled with GCC versions 4.6.3 in the latest iteration. 
Earlier version should be fine, but the status of OpenMP implementation should be checked,
as well as the various compiler flags. The compiler flags used were:

~~~
-O3 -Wall -Wextra -funsafe-loop-optimizations -xT -fipa-pta -ftracer -march=native -fopenmp
~~~

Only mathematically safe optimizations were enabled.

The statistics and plotting are done with Python 2.7 scripts,
using matplotlib library version >= 1.2.1
