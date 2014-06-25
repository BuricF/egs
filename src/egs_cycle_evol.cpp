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



/** @file egs_cycle_evol.cpp

   Extended Gray-Scott Model Evlution Simulator

   This version is built upon the simple simulator,
   running the dyanmic for a specified length of time.
   
   The simulation features subsequent mutation events until all species
   go extinct or the alloted times runs out. The roles of wild type and mutant
   alternate between \f$V_1\f$, \f$V_2\f$, \f$V_3\f$ in a deterministic way 
   (see comments and formal algorithm outlined in the thesis report).
   
   This version only features deterministic dynamics beyond the initial conditions randomisation
   but can easily be adapted following the code in the simple simulator.
*/

#include <cstdio>        // stdio was chosen over iostream since it should have less overhead
#include <cstdlib>
#include <vector>
#include <ctime>        // use system time to set random seed
#include <cmath>
#include <cfloat>

#include <omp.h>        // OpenMP control function

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// EGS header files
#include "include/operators.hpp"                // various "operators" that manipulate the state of the concentration matrices
#include "include/egs_cycle_dynamics.hpp"       // includes concentration operators used in experiments
#include "include/system.hpp"                   // handles saving and loading of data to the hard drive
#include "include/memory.hpp"                   // handles management of concentration matrices (allocation, swapping,..)

using namespace std;

// effective null concentration (based on Avogadro constant)
#define LOWLIM 1.0e-24


int main() {

    // set up random number generator =============
    const gsl_rng_type * T;
    T = gsl_rng_taus;
    gsl_rng* rgen  = gsl_rng_alloc(T);
    if (!rgen)
        fprintf(stderr, "main(): Error allocating random number generator.\n");
    const unsigned int rand_seed = 1391087569; //(unsigned)time(0); 
    gsl_rng_set(rgen, rand_seed);
    
    // save rng seed for this simulation (determines all pseudo-random numbers)
    FILE *ssave = fopen("simdata/seed", "w");
    if (!ssave) ssave = stdout;
    fprintf(ssave, "%d", rand_seed);
    fclose(ssave);                      // close early in case simulation crashes
    // ============================================
    

    
    // Parameters ===========================

    // Reaction ==============
    const float Dc = 1.4e-9;     // V diffusion coefficient
    const float Du = 2 * Dc;     // U diffusion coefficient

    // fuel feed rate
    const float F = 0.04;

    // parameter ranges
    const float k_low = 0.05;
    const float k_high = 0.07;
    const float r_low = 0.0;
    const float r_high = 2.0;
    
    const float k_stddev = 0.01;    // offset of new k drawn from N(0, k_stddev^2) 
    const float r_stddev = 0.01;    // offset of new r drawn from N(0, k_stddev^2)
    
    // initial catalyst degradation rates
    float k1 = 0.065;
    float k2 = 0.066;
    float k3 = 0.065;         // irrelevant, will be set by mutation
    
    // reaction rate matrix  (s_i == r[i][i])    
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
    r[0][1] = 1 + 0.3;          //sqs12 - 0.1;  // r12
    if (r[0][1] > sqs12) r[0][1] = sqs12;       // check for too large rate
    r[1][0] = sqs12 - r[0][1];                  // r21 (C2)
    
    // These rate will be reassigned when V3 is generated randomly 
    // but kept in case different mix given as init. cond.
     r[0][2] = 0.1;                              // r13
     if (r[0][2] > sqs13) r[0][2] = sqs13;       // check for too large rate
     r[2][0] = sqs13 - r[0][2];                  // r31 (C2)
 
     r[1][2] = sqs23 - 0.1;                      // r23
     if (r[1][2] > sqs23) r[1][2] = sqs23;       // check for too large rate
     r[2][1] = sqs23 - r[1][2];                  // r32 (C2)  


    char* icondfile = "cond/v1v2_mix.cond";      // Initial conditions specification file
    const bool initial_randomization = true;     // If the initial conditions are to be perturbed
    const bool noisy_dynamics = false;           // Enable noise term in reaction-diffusion dyanmics

    const float noise_stddev = 1e-20 / Du;  // use higher diffusion value as most significant
    const float noise_level = 1.0e3;        // noise amplification \lambda
    // ========================


    // Simulation =============
    const int N = 150;             // grid length (no. of columns)
    const int M = 1*N;             // grid height (no. of rows)
    const float L = 0.01;          // system size
    float dt = 0.05;               // integration time step
    const bool midpoint = false;                    // if to use the midpoint method instead of forward Eurler (extra half step)
    
    // square root of grid cell particle capacity: discretisation used in computing noise amount
    // The noise functions uses this value so no need to square to the 2D cell capacity
    const float cell_cap_1D = (L/N) * 1.0e6;
    
    const int GENLEN = (int)ceil((L*N)/(dt*dt));     // length of a generation heuristic
    const int NITER = 20000 * GENLEN;                // simulation length: no. of iterations
   
    const bool resume = false;                       // Load final state of previous simulation and continue from there
                                                     // (Note that should ensure the same parameters between simulations)
    const bool skip_save = false;                    // save only initial and final system states
    // ========================


    // Mutant =================    
    const char mscheme = 'M';         // mutation scheme: mut to be placed at {m=min, M=max, r=rand} w.t. concentration
    const bool m_internal = false;    /* true: The mutant quantity is to be a fraction of the wild type (mconc = v*mfract).
                                               This will also cause it to be subtracted from the w.t. quantity at that position.
                                         false: Constant mutant concentration (mfract), not to be subtracted from w.t. quantity. */

    const float mfract = 0.2;        // fraction or concentration of mutated catalyst (see m_internal above)
    const int MSIZE = 2;             // radius of mutant drop i.e. range of Moore neighbourhood (0: singe cell)
    
    const int EXTINCT_CHECK = 5 * GENLEN;      // need to check (often) if one of the species disappeared
    const int CHECK_DELAY = 150 * GENLEN;      // start checking only after a certain time has passed since simuation started
    // ========================


    // == Initialization =====================================
    const float hsq = (L/N)*(L/N);   // precompute h^2 for speed
    
    // Show numeric details ==================
    printf("\n[ Extended Gray-Scott model evolution simulation ]\n\n");
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
        perror("main(): Cannot allocate memory. Aborting.\n");
        return 1;
    }

    // Open files for saving fuel and catalyst values
    datafile ffuel("simdata/fuel.dat", "fuel.dat", u, M, N, resume);
    datafile fcat1("simdata/catalyst1.dat", "catalyst1.dat", v1, M, N, resume);
    datafile fcat2("simdata/catalyst2.dat", "catalyst2.dat", v2, M, N, resume);
    datafile fcat3("simdata/catalyst3.dat", "catalyst3.dat", v3, M, N, resume);
    // ======================================

    
    // Simulation ==============================================================================

    if (resume) {
        // load final state of old data
        load_state(u, M, N, "simdata/fuel.final");         
        load_state(v1, M, N, "simdata/catalyst1.final");
        load_state(v2, M, N, "simdata/catalyst2.final");
        load_state(v3, M, N, "simdata/catalyst3.final");
    }

    // Read initial perturbations (concentrations) from file
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

    
    // Species bookkeeping ==========
    bool x[3] = {true,true,false};             // marks which species are present in the system
    int eventno = 0;                           // event number
    bool cont = true;                          // flag to stop sim in certain conditions
    
    FILE *evollog = fopen("simdata/evol.log", "w");           // record parameters of species during mutation events
    if (!evollog) evollog = stdout;
    
    FILE *paramrec = fopen("simdata/evolparams.dat", "w");    // record parameters of species during mutation events
    if (!paramrec) paramrec = stdout;
    // ==============================
    
    
    // Step the evolution =========================
    for (int iter=1; iter<=NITER; ++iter) if (cont) {

        if (iter % GENLEN == 0) {
            printf("gen %d \r", iter / GENLEN);
            fflush(stdout);
        }
        
        int int_steps = 1;  // number of integration steps
        if (midpoint) {
            ++int_steps;
            dt /= 2.0;      // midpoint method half-step
        }

        for (int istep=0; istep<int_steps; ++istep) {

            if ((midpoint) && (istep == 1))  
                dt *= 2.0;                              // readjust time step after half step if using midpoint method

            #pragma omp parallel for collapse(2)
            for (int i=0; i<M; ++i)
                for (int j=0; j<N; ++j) {

                    // deterministic dynamics
                    u_new[i][j] =   u[i][j] + dt *  react_u(u, v1, v2, v3, i, j, hsq, M, N, r, F, Du);
                    v1_new[i][j] = v1[i][j] + dt * react_v1(u, v1, v2, v3, i, j, hsq, M, N, r, F, k1, Dc);
                    v2_new[i][j] = v2[i][j] + dt * react_v2(u, v1, v2, v3, i, j, hsq, M, N, r, F, k2, Dc);
                    v3_new[i][j] = v3[i][j] + dt * react_v3(u, v1, v2, v3, i, j, hsq, M, N, r, F, k3, Dc);

                    // hard checks
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
        

        
        // Regular check if a mutant should be introduced ( <=> when a species has gone extinct)
        if ((iter > CHECK_DELAY) && (iter % EXTINCT_CHECK == 0)) {
        
            // Get total concentrations (local variables as optimization)
            double v1_tot = 0.0;
            double v2_tot = 0.0;
            double v3_tot = 0.0;
            
            #pragma omp parallel for collapse(2) schedule(static) reduction(+:v1_tot,v2_tot,v3_tot)
            for (int i=0; i<M; ++i)
                for (int j=0; j<N; ++j) {
                    if (v1[i][j] > LOWLIM) 
                        v1_tot += v1[i][j];
        
                    if (v2[i][j] > LOWLIM)
                        v2_tot += v2[i][j];
                    
                    if (v3[i][j] > LOWLIM)
                        v3_tot += v3[i][j];
                }
            
            
            /* Mark which species present in the system:
             * can compare absolutely with zero since quantities <= LOWLIM were truncated at exact zero
             */
            x[0] = (v1_tot != 0.0);     
            x[1] = (v2_tot != 0.0);     
            x[2] = (v3_tot != 0.0);
                       
            // If not all species present
            if (!(x[0] && x[1] && x[2])) {
            
                // In case of complete extinction
                if (!(x[0] || x[1] || x[2])) {
                    printf("All 3 species extinct: stopping simulation.\n");
                    cont = false;
                    break;
                }
                
                ++eventno;                         // Another mutation event occurs

                // Assign roles (wild type, mutant based on w.t., and the third "original" species)
                double ** v[3] = {v1, v2, v3};     // indexes for easier role assignment
                float * k[3] = {&k1, &k2, &k3};    // (not used from beginning given dereferencing overhead)
                int o = -1;                        // role indexes
                int w = -1;
                int m = -1;
                
                for (int i=0; i<3; ++i) {
                    if (!(x[i]))
                        m = i;          // "mutant" is the last not-present species
                    else {
                        o = w;          // "original" is the first present species,
                                        // can be unset (-1) if only one species present
                        w = i;          // "wild type" is the last present species, always set to an existing species
                                        // always the species under mutation
                    }
                }
                
                /* Generate mutant values
                 * 
                 * Convention:  Only _received_ catalytic rates are set by mutation.
                 *              The back-reactions (offered catalytic support) is set according to Constraint 2
                 */
                
                *k[m] = *k[w] + gsl_ran_gaussian(rgen, k_stddev);
                r[m][m] = r[w][w] + gsl_ran_gaussian(rgen, r_stddev);
                r[w][m] = r[w][w] + gsl_ran_gaussian(rgen, r_stddev);
                    
                // impose limits
                if (*k[m] < k_low) *k[m] = k_low;
                if (*k[m] > k_high) *k[m] = k_high;
                if (r[m][m] < r_low) r[m][m] = r_low;
                if (r[m][m] > r_high) r[m][m] = r_high;
                if (r[w][m] < r_low) r[w][m] = r_low;
                if (r[w][m] > r_high) r[w][m] = r_high;
                    
                // impose Constraint 2
                r[m][w] = 2*sqrt(r[w][w]*r[m][m]) - r[w][m];
                
                // There are 2 species in the system: additional catalytic link with "original"
                if (o != 1) {

                    r[o][m] = r[o][w] + gsl_ran_gaussian(rgen, r_stddev);
                    
                    // impose limits
                    if (r[o][m] < r_low) r[o][m] = r_low;
                    if (r[o][m] > r_high) r[o][m] = r_high;
                    
                    // impose C2
                    r[m][o] = 2*sqrt(r[o][o]*r[m][m]) - r[o][m];
                }
                
                // Record and output mutation parameters
                if (o != 1) {   // 3 species
                    
                    fprintf(evollog, "mutation %d: o = %d\tw = %d\tm = %d;\tk_m = %f\tr_mm = %f\tr_wm = %f\tr_mw = %f\tr_om = %f\tr_mo = %f\n", 
                            eventno, o+1, w+1, m+1, *k[m], r[m][m], r[w][m], r[m][w], r[o][m], r[m][o]);
                    printf("mutation %d: o = %d\tw = %d\tm = %d;\tk_m = %f\tr_mm = %f\tr_wm = %f\tr_mw = %f\tr_om = %f\tr_mo = %f\n", 
                            eventno, o+1, w+1, m+1, *k[m], r[m][m], r[w][m], r[m][w], r[o][m], r[m][o]);
                } else {        // 2 species
                    
                    fprintf(evollog, "mutation %d: w = %d\tm = %d;\tk_m = %f\tr_mm = %f\tr_wm = %f\tr_mw = %f\n", 
                            eventno, w+1, m+1, *k[m], r[m][m], r[w][m], r[m][w]);
                    printf("mutation %d: w = %d\tm = %d;\tk_m = %f\tr_mm = %f\tr_wm = %f\tr_mw = %f\n", 
                            eventno, w+1, m+1, *k[m], r[m][m], r[w][m], r[m][w]);
                }
            
                // Place mutant perturbation
                mutation_event(v[w], v[m], M, N, mfract, MSIZE, mscheme, m_internal);
                
                
                // Record all cycle parameters (including existence vector to distinguish 1->2 and 2->3 species cases)
                fprintf(paramrec, "%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
                               iter/GENLEN, x[0], x[1], x[2], 
                               k1, k2, k3, 
                               r[0][0], r[0][1], r[0][2],
                               r[1][0], r[1][1], r[1][2],
                               r[2][0], r[2][1], r[2][2]);
            }
        }
        
        // write current generation to file
        if (((!skip_save) && (iter % GENLEN == 0)) || (iter == NITER)) {
            ffuel.dump();
            fcat1.dump();
            fcat2.dump();
            fcat3.dump();
        }

    }
    // ========================================
    

    // Clean up ================
    fclose(evollog);
    fclose(paramrec);
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
