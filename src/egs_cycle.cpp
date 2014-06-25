/*
This file is part of the EGS - Extended Gray-Scott Model simulation package

Copyright (C) 2014 Filip Buric

EGS is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

EGS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/


/** @file egs_cycle.cpp

   Extended Gray-Scott Model Simulator
   
   The system features one 2D reaction chamber modelled as a grid 
   with periodic boundaries (torus topology).
   
   The differential equations are integrated using the forward Euler method
   with the optional half-step (midpoint method).
   
   Initial conditions are loaded from a file in the cond/ directory 
   which specifies the concentration and spread of each species present initially
   (this markup format is described in the .cond files).
   
   If species are to appear through mutation, the convention is 
   that these are to be \f$V_2\f$ and \f$V_3\f$ (\f$V_1\f$ is always the "original" species).
   
   The system and reaction-diffusion paramters are hardcoded since it proved
   easier to manage when reparameterising and extending the models, 
   or when defining experiment parameters.
   
   Degradation rates \f$k_i\f$ and the fuel feed rate F are cosntant across the grid.
   
   The location of the generated simulation data directory is hardcoded at simdata/
   Status messages and events are printed to stdout.
   The format of the concentration matrix evolution data for each species is a
   .dat file containing the matrix at each generation, separated by an empty line.
   The numerical values are stored as text, not binary, for human inspection.
   The first line is the number of rows M (height) of the matrices.
   
   Errors during execution are meant to be loud and terminating so no hidden errors
   affect the results.
   
   File paths are relative to current working directory, which is set by the controlling script.
   
   As a note on execution speed, the bottleneck was profiled to be 
   the discrete Laplacian operator.
*/


#include <cstdio>       // stdio was chosen over iostream since it should have less overhead
#include <cstdlib>
#include <vector>
#include <set>
#include <ctime>        // use system time to set random seed
#include <cmath>
#include <cfloat>

#include <omp.h>        // OpenMP control function

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// effective null concentration (based on Avogadro constant): def'd here since needed in includes
#define LOWLIM 1.0e-24

// EGS header files
#include "include/operators.hpp"
#include "include/egs_cycle_dynamics.hpp"       // includes concentration operators used in experiments
#include "include/system.hpp"                   // handles saving and loading of data to the hard drive
#include "include/memory.hpp"                   // handles management of concentration matrices (allocation, swapping,..)


using namespace std;



int main() {

    // Parameters ===========================

    // Reaction ==============
    const float Dc = 1.4e-9;     // V diffusion coefficient
    const float Du = 2 * Dc;     // U diffusion coefficient

    // fuel feed rate
    const float F = 0.04;

    // catalyst degradation rates
    const float k1 = 0.06;
    const float k2 = 0.065;
    const float k3 = 0.067;
    
    // reaction rate matrix
    float r[][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}
    };
    
    // 2sqrt(r_ii*r_jj) abbreviations for convenience
    #define     sqs12       (2.0*sqrt(r[0][0]*r[1][1]))
    #define     sqs13       (2.0*sqrt(r[0][0]*r[2][2]))
    #define     sqs23       (2.0*sqrt(r[1][1]*r[2][2]))

    // Set received/given catalytic rates and impose C2 constraint
    r[0][1] = 0.1;                              // r12
    if (r[0][1] > sqs12) r[0][1] = sqs12;       // check for too large rate
    r[1][0] = sqs12 - r[0][1];                  // r21 (C2)
    
    r[1][2] = 0.1;                              // r23
    if (r[1][2] > sqs23) r[1][2] = sqs23;       // check for too large rate
    r[2][1] = sqs23 - r[1][2];                  // r32 (C2)  
    
    r[2][0] = 0.9;                              // r31
    if (r[2][0] > sqs13) r[2][0] = sqs13;       // check for too large rate
    r[0][2] = sqs13 - r[2][0];                  // r13 (C2)


    char icondfile[] = "cond/mix.cond";         // Initial conditions specification file
    const bool initial_randomization = true;    // If the initial conditions are to be perturbed
    bool noisy_dynamics = true;                 // Enable noise term in reaction-diffusion dyanmics

    const float noise_stddev = sqrt(1e-20 / Du);  // use higher diffusion value as most significant 
    const float noise_level = 1.0e6;        // noise amplification \lambda
    // ========================


    // Simulation =============
    const int N = 150;                // grid length (no. of columns)
    const int M = 1*N;                // grid height (no. of rows)
    const float L = 0.01;             // system size
    float dt = 0.05;                  // integration time step
    const bool midpoint = false;                    // if to use the midpoint method instead of forward Eurler (extra half step)
    
    // square root of grid cell particle capacity: discretisation used in computing noise amount
    // The noise functions uses this value so no need to square to the 2D cell capacity
    const float cell_cap_1D = (L/N) * 1.0e6;
    
    const int GENLEN = (int)ceil((L*N)/(dt*dt));    // length of a generation heuristic
    const int NITER = 100 * GENLEN;                  // simulation length: no. of iterations

    const bool resume = false;                      // Load final state of previous simulation and continue from there
                                                    // (Note that should ensure the same parameters between simulations)
    const bool skip_save = false;                   // save only initial and final system states
    // ========================


    // Mutant =================    
    const char mscheme = 'M';         // mutation scheme: mut to be placed at {m=min, M=max, r=rand} w.t. concentration
    const bool m_internal = false;    /* true: The mutant quantity is to be a fraction of the wild type (mconc = v*mfract).
                                               This will also cause it to be subtracted from the w.t. quantity at that position.
                                         false: Constant mutant concentration (mfract), not to be subtracted from w.t. quantity. */

    const float mfract = 0.2;         // fraction or concentration of mutated catalyst (see m_internal above)
    const int MSIZE = 2;              // radius of mutant drop i.e. range of Moore neighbourhood (0: singe cell)

    // whether V2 and/or V3 appear in the system as mutants
    const bool mutation_event1 = false;     // V2
    const bool mutation_event2 = false;     // V3
    
    const int MTIME1 = 150 * GENLEN;        // iteration for first mutation event 
    const int MTIME2 = 155 * GENLEN;        // iteration for second mutation event
    // ========================


    // Experiment: perturbation with noise =====
    const int PERTURB_DURATION = 20;        // how long is the noise perturbation to be, counted in generations
    set<int> perturbation_times;            // keeps the generation numbers when noise perturbations are to begin
    //perturbation_times.insert(500);
    // ===========================================

   
    
    // == Initialization =====================================
    const float hsq = (L/N)*(L/N);   // precompute h^2 for speed

    // set up GSL random number generator =============
    const gsl_rng_type * T;
    T = gsl_rng_taus;
    gsl_rng* rgen  = gsl_rng_alloc(T);
    if (!rgen) {
        perror("main(): Error allocating random number generator.");
        return 1;
    }
    const unsigned int rand_seed = (unsigned)time(0); 
    gsl_rng_set(rgen, rand_seed);
    
    // save random no. generator seed for this simulation (determines all pseudo-random numbers)
    FILE *ssave = fopen("simdata/seed", "w");
    if (!ssave) {
        perror("Could not open a file to save random seed. Printing on screen");
        printf("%d", rand_seed);
    }
    else {
        fprintf(ssave, "%d", rand_seed);
        fclose(ssave);                      // close early in case simulation crashes
    }
    // ==============================================
    
    
    // Show numeric details ==================
    printf("\n[ Extended Gray-Scott model simulation ]\n\n");
    printf("dx = %e\n", (L/N));
    printf("dt = %e\n", dt);
    printf("Du = %e\n", Du);
    printf("cell cap = %e particles\n", cell_cap_1D*cell_cap_1D);
        
    printf("CFL number for u: %f\n", (Du*dt) / hsq );
    printf("CFL number for v_i: %f\n", (Dc*dt) / hsq );
    
    if (noisy_dynamics) {
        printf("noise std. dev. = %e\n", noise_stddev);
        printf("noise amplif. = %e\n", noise_level);
    }
    
    printf("\nN = %d for %d generations (genlen = %d)..\n", N, NITER/GENLEN, GENLEN);
    
    printf("r12 = %f\tr21 = %f\n", r[0][1], r[1][0]);
    printf("r23 = %f\tr32 = %f\n", r[1][2], r[2][1]);
    printf("r31 = %f\tr13 = %f\n", r[2][0], r[0][2]);
    // ======================================


    // OpenMP settings =======
    int MAX_THREADS = 2;        // maximum number of allowed parallel threads to compute dynamics
    
    #ifdef _OPENMP
    omp_set_num_threads(MAX_THREADS);
    printf("[OpenMP: max %d threads]\n", MAX_THREADS);
    #endif
    // =======================

    /* 
     * Allocate data matrices: dynamic since expected to be of large size
     * 
     * Note: Raw arrays were preferred over std::vector or similar for less overhead
     * and the possibility of adaopting them to other paralellising schemes other than OpenMP.
     * 
     * {..}_new matrices hold computed values for next time step but these roles swap every iteration.
     * 
     */
    double **u, **v1, **v2, **v3;
    double **u_new, **v1_new, **v2_new, **v3_new;       
    
    try {
        gridalloc(&u, M, N);
        gridalloc(&v1, M, N);
        gridalloc(&v2, M, N);
        gridalloc(&v3, M, N);
        gridalloc(&u_new, M, N);
        gridalloc(&v1_new, M, N);
        gridalloc(&v2_new, M, N);
        gridalloc(&v3_new, M, N);
    }
    catch(int e) {
        perror("main(): Cannot allocate memory. Aborting.");
        return 1;
    }

    // Open files for saving fuel and catalyst values
    datafile ffuel("simdata/fuel.dat", "fuel.dat", u, M, N, resume);
    datafile fcat1("simdata/catalyst1.dat", "catalyst1.dat", v1, M, N, resume);
    datafile fcat2("simdata/catalyst2.dat", "catalyst2.dat", v2, M, N, resume);
    datafile fcat3("simdata/catalyst3.dat", "catalyst3.dat", v3, M, N, resume);
    // ======================================

    
    // Simulation ===========================   
    if (resume) {
        // load final state of old data
        load_state(u, M, N, "simdata/fuel.final");   
        load_state(v1, M, N, "simdata/catalyst1.final");
        load_state(v2, M, N, "simdata/catalyst2.final");
        load_state(v3, M, N, "simdata/catalyst3.final");
    }

    // Read initial conditions (concentrations) from file
    else 
        if (load_initial_conditions(icondfile, u, v1, v2, v3, M, N, initial_randomization, cell_cap_1D, noise_stddev, noise_level, rgen)) {
            perror("main(): Could not get initial conditions! Aborting.");
            return 1;
        }

    // write initial state
    ffuel.dump();
    fcat1.dump();
    fcat2.dump();
    fcat3.dump();

    /* Record species parameters to file
     * The first four arguments give:  iteration number, 
     * and three 1/0 if species i exists or not in the system at that time; irrelevant here but format kept consistent
     * with the evolution simulation for plotting the cycle diagram
     */
    FILE *paramrec = fopen("simdata/cycleparams.dat", "w");    
    if (!paramrec) paramrec = stdout;
    fprintf(paramrec, "%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                      0, 1, 1, 0, 
                      k1, k2, k3, 
                      r[0][0], r[0][1], r[0][2],
                      r[1][0], r[1][1], r[1][2],
                      r[2][0], r[2][1], r[2][2]);
    fclose(paramrec);

    // generation at which a noise perturbation started
    int pert_start = 0;       
    
    // Step the evolution =========================
    for (int iter=1; iter<=NITER; ++iter) {

        if (iter % (100*GENLEN) == 0) {
            printf("gen %d \r", iter / GENLEN);
            fflush(stdout);
        }
        
        int int_steps = 1;  // number of integration steps
        if (midpoint) {
            ++int_steps;
            dt /= 2.0;      // midpoint method half-step
        }

        
        // Experiment : noise perturbation ===============================
        if (!perturbation_times.empty()) {
            printf("\n noisy dynamics\n");
            if (iter % GENLEN == 0) {               // only check once per generation
                
                // if current generation marked as start of event, start noise
                if (perturbation_times.find(iter / GENLEN) != perturbation_times.end()) {   
                    noisy_dynamics = true;
                    pert_start = iter / GENLEN;
                }
                
                // if curr. gen. marked as end of event, stop noise
                if (iter / GENLEN == pert_start + PERTURB_DURATION) {
                    noisy_dynamics = false;
                }
            }
        }
        // ================================================================


        for (int istep=0; istep<int_steps; ++istep) {

            if ((midpoint) && (istep == 1))     
                dt *= 2.0;                              // readjust time step after half step if using midpoint method

            #pragma omp parallel for collapse(2)
            for (int i=0; i<M; ++i)
                for (int j=0; j<N; ++j) {

                    // deterministic dynamics
                    if (!noisy_dynamics) {
                        u_new[i][j] =   u[i][j] + dt *  react_u(u, v1, v2, v3, i, j, hsq, M, N, r, F, Du);
                        v1_new[i][j] = v1[i][j] + dt * react_v1(u, v1, v2, v3, i, j, hsq, M, N, r, F, k1, Dc);
                        v2_new[i][j] = v2[i][j] + dt * react_v2(u, v1, v2, v3, i, j, hsq, M, N, r, F, k2, Dc);
                        v3_new[i][j] = v3[i][j] + dt * react_v3(u, v1, v2, v3, i, j, hsq, M, N, r, F, k3, Dc);
                    }
                    // noisy dynamics
                    else {
                        u_new[i][j] =   u[i][j] + dt *  react_u(u, v1, v2, v3, i, j, hsq, M, N, r, F, Du, cell_cap_1D, noise_stddev, noise_level, rgen);
                        v1_new[i][j] = v1[i][j] + dt * react_v1(u, v1, v2, v3, i, j, hsq, M, N, r, F, k1, Dc, cell_cap_1D, noise_stddev, noise_level, rgen);
                        v2_new[i][j] = v2[i][j] + dt * react_v2(u, v1, v2, v3, i, j, hsq, M, N, r, F, k2, Dc, cell_cap_1D, noise_stddev, noise_level, rgen);
                        v3_new[i][j] = v3[i][j] + dt * react_v3(u, v1, v2, v3, i, j, hsq, M, N, r, F, k3, Dc, cell_cap_1D, noise_stddev, noise_level, rgen);
                    }
                    
                    // hard checks that computed values stay within defined range (accumulation may creep in and have possible noise terms)
                    if (u_new[i][j] > 1.0) u_new[i][j] = 1.0;
                    if (u_new[i][j] < LOWLIM) u_new[i][j] = 0.0;

                    if (v1_new[i][j] > 1.0) v1_new[i][j] = 1.0;
                    if (v1_new[i][j] < LOWLIM) v1_new[i][j] = 0.0;

                    if (v2_new[i][j] > 1.0) v2_new[i][j] = 1.0;
                    if (v2_new[i][j] < LOWLIM) v2_new[i][j] = 0.0;
                    
                    if (v3_new[i][j] > 1.0) v3_new[i][j] = 1.0;
                    if (v3_new[i][j] < LOWLIM) v3_new[i][j] = 0.0;
                }

            swap_mat(&u, &u_new);
            swap_mat(&v1, &v1_new);
            swap_mat(&v2, &v2_new);
            swap_mat(&v3, &v3_new);
        }

        // Mutantion events:  generate new species and place perturbation thereof
        if (mutation_event1 && (iter == MTIME1) && (!resume)) {
            mutation_event(v1, v2, M, N, mfract, MSIZE, mscheme, m_internal);
            printf("mutation 1: k2 = %f\ts2 = %f\tr12 = %f\tr21 = %f\n", k2, r[1][1], r[0][1], r[1][0]);
        }
        else if (mutation_event2 && (iter == MTIME2) && (!resume)) {
            mutation_event(v2, v3, M, N, mfract, MSIZE, mscheme, m_internal);        
            printf("mutation 2: k3 = %f\ts3 = %f\tr13 = %f\tr31 = %f\tr23 = %f\tr32 = %f\n", k3, r[2][2], r[0][2], r[2][0], r[1][2], r[2][1]);
        }
            
        // Write current generation to file, unless set to skip
        if (((!skip_save) && (iter % GENLEN == 0)) || (iter == NITER)) {
            ffuel.dump();
            fcat1.dump();
            fcat2.dump();
            fcat3.dump();
        }

    }
    // ========================================
    

    // Write final state ===
    ffuel.dump_final();
    fcat1.dump_final();
    fcat2.dump_final();
    fcat3.dump_final();
    // =====================

    
    // Clean up ================
    gsl_rng_free(rgen);
    gridfree(&u, M);
    gridfree(&v1, M);
    gridfree(&v2, M);
    gridfree(&v3, M);
    gridfree(&u_new, M);
    gridfree(&v1_new, M);
    gridfree(&v2_new, M);
    gridfree(&v3_new, M);
    // ========================

    return 0;
}
