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


/** @file egs_cycle_dynamics.hpp

    Dynamics for the 3-component hypercyle extension of the Gray-Scott model
    (reaction-diffusion equations)

        U + 2 V1  -r[0][0]->  3V1
        U + 2 V2  -r[1][1]->  3V2
        U + 2 V3  -r[2][2]->  3V3

        U + V1 + V2 -r21-> 2V1 + V2
        U + V1 + V2 -r12-> 2V2 + V1
        
        U + V1 + V3 -r31-> 2V1 + V3
        U + V1 + V3 -r13-> 2V3 + V1
        
        U + V2 + V3 -r32-> 2V2 + V3
        U + V2 + V3 -r23-> 2V3 + V2
        
        V1 -k1-> P
        V2 -k2-> P
        V3 -k3-> P

    There are two separate versions of the reaction-diffusion functions,
    with and witout noise, respectively, as an optimisation.
    Uses 5-point stencil Laplacian discretisation (for speed)
    
    While technically a concentration operator, the initial conditions loading procedure
    was included here so that this file would describe the initial value problem.
    
 */

// #include "operators.hpp"

/*
 * Noiseless (deterministic) versions
 */

/**
 *  Computes the change in U (deterministic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] Du        Diffusion coefficient of U
 */
inline double react_u(double** u, double** v1, double** v2, double** v3,
                      int i, int j, float hsq, int M, int N,
                      float r[][3], float F, float Du) {

    double u_new;

    // diffusion
    u_new = Du * laplacian5_per_grid(u, i, j, M, N, hsq);

    // reaction
    u_new += -u[i][j] * (r[0][0]*v1[i][j]*v1[i][j] + r[1][1]*v2[i][j]*v2[i][j] + r[2][2]*v3[i][j]*v3[i][j] +
                         (r[0][1] + r[1][0])*v1[i][j]*v2[i][j] +
                         (r[1][2] + r[2][1])*v2[i][j]*v3[i][j] +
                         (r[0][2] + r[2][0])*v1[i][j]*v3[i][j]
                        );

    // feed
    u_new += F * (1.0 - u[i][j]);

    return u_new;
}


/**
 *  Computes the change in V1 (deterministic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k1        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V1
 */
inline double react_v1(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k1, float Dc) {

    double v1_new;

    // diffusion
    v1_new = Dc * laplacian5_per_grid(v1, i, j, M, N, hsq);

    // reaction
    v1_new += u[i][j] * (r[0][0]*v1[i][j]*v1[i][j] + r[1][0]*v1[i][j]*v2[i][j] + r[2][0]*v1[i][j]*v3[i][j]);

    // feed
    v1_new += -(F + k1) * v1[i][j];

    return v1_new;
}


/**
 *  Computes the change in V2 (deterministic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k2        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V2
 */
inline double react_v2(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k2, float Dc) {

    double v2_new;

    // diffusion
    v2_new = Dc * laplacian5_per_grid(v2, i, j, M, N, hsq);

    // reaction
    v2_new += u[i][j] * (r[1][1]*v2[i][j]*v2[i][j] + r[0][1]*v1[i][j]*v2[i][j] + r[2][1]*v2[i][j]*v3[i][j]);

    // feed
    v2_new += -(F + k2) * v2[i][j];

    return v2_new;
}


/**
 *  Computes the change in V3 (deterministic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k3        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V3
 */
inline double react_v3(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k3, float Dc) {

    double v3_new;

    // diffusion
    v3_new = Dc * laplacian5_per_grid(v3, i, j, M, N, hsq);

    // reaction
    v3_new += u[i][j] * (r[2][2]*v3[i][j]*v3[i][j] + r[1][2]*v2[i][j]*v3[i][j] + r[0][2]*v1[i][j]*v3[i][j]);

    // feed
    v3_new += -(F + k3) * v3[i][j];

    return v3_new;
}





/*
 *  Versions with noise (stochastic dynamics)
 */

/**
 *  Computes the change in U (stochasti version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] Du        Diffusion coefficient of U
 */
inline double react_u(double** u, double** v1, double** v2, double** v3,
                      int i, int j, float hsq, int M, int N,
                      float r[][3], float F, float Du,
                      float cell_cap_1D, float noise_stddev, float noise_level, gsl_rng* rgen) {

    double u_new;

    // diffusion
    u_new = Du * laplacian5_per_grid(u, i, j, M, N, hsq);

    // reaction
    u_new += -u[i][j] * (r[0][0]*v1[i][j]*v1[i][j] + r[1][1]*v2[i][j]*v2[i][j] + r[2][2]*v3[i][j]*v3[i][j] +
                         (r[0][1] + r[1][0])*v1[i][j]*v2[i][j] +
                         (r[1][2] + r[2][1])*v2[i][j]*v3[i][j] +
                         (r[0][2] + r[2][0])*v1[i][j]*v3[i][j]
                        );

    // feed
    u_new += F * (1.0 - u[i][j]);

    // noise
    u_new += gauss_noise(u[i][j], cell_cap_1D, noise_stddev, noise_level, rgen);

    return u_new;
}


/**
 *  Computes the change in V1 (stochastic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k1        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V1
 */
inline double react_v1(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k1, float Dc,
                       float cell_cap_1D, float noise_stddev, float noise_level, gsl_rng* rgen) {

    double v1_new;

    // diffusion
    v1_new = Dc * laplacian5_per_grid(v1, i, j, M, N, hsq);

    // reaction
    v1_new += u[i][j] * (r[0][0]*v1[i][j]*v1[i][j] + r[1][0]*v1[i][j]*v2[i][j] + r[2][0]*v1[i][j]*v3[i][j]);

    // feed
    v1_new += -(F + k1) * v1[i][j];

    // noise
    v1_new += gauss_noise(v1[i][j], cell_cap_1D, noise_stddev, noise_level, rgen);

    return v1_new;
}


/**
 *  Computes the change in V2 (stochastic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k2        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V2
 */
inline double react_v2(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k2, float Dc,
                       float cell_cap_1D, float noise_stddev, float noise_level, gsl_rng* rgen) {

    double v2_new;

    // diffusion
    v2_new = Dc * laplacian5_per_grid(v2, i, j, M, N, hsq);

    // reaction
    v2_new += u[i][j] * (r[1][1]*v2[i][j]*v2[i][j] + r[0][1]*v1[i][j]*v2[i][j] + r[2][1]*v2[i][j]*v3[i][j]);

    // feed
    v2_new += -(F + k2) * v2[i][j];

    // noise
    v2_new += gauss_noise(v2[i][j], cell_cap_1D, noise_stddev, noise_level, rgen);

    return v2_new;
}



/**
 *  Computes the change in V3 (stochastic version)
 * 
 * @param[in] u,v1,v2,v3        Concentration matries
 * @param[in] i,j       Position in grid
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$ used by the Laplacian
 * @param[in] M,N       Size (height = rows, length = columns) of the grid
 * @param[in] r         Reaction rates matrix
 * @param[in] F         Fuel feed rates
 * @param[in] k3        V1 degradation rate
 * @param[in] Dc        Diffusion coefficient of V3
 */
inline double react_v3(double** u, double** v1, double** v2, double** v3,
                       int i, int j, float hsq, int M, int N,
                       float r[][3], float F, float k3, float Dc,
                       float cell_cap_1D, float noise_stddev, float noise_level, gsl_rng* rgen) {

    double v3_new;

    // diffusion
    v3_new = Dc * laplacian5_per_grid(v3, i, j, M, N, hsq);

    // reaction
    v3_new += u[i][j] * (r[2][2]*v3[i][j]*v3[i][j] + r[1][2]*v2[i][j]*v3[i][j] + r[0][2]*v1[i][j]*v3[i][j]);

    // feed
    v3_new += -(F + k3) * v3[i][j];

    // noise
    v3_new += gauss_noise(v3[i][j], cell_cap_1D, noise_stddev, noise_level, rgen);
    
    return v3_new;
}


/**
 * Loads initial conditions or fail with a non-null return value.
 * 
 * Notes:
 *    1. Replaces all previous data.
 *    2. Any unkown catalyst type numbers (>=4) are ignored.
 *
 * @param[in] cond_file The file describing the inital conditions
 * @param[in] u,v1,v2,v3 Concentration matrices
 * @param[in] nrows,ncols Size of matrices
 * @return 0 = success, 1 = failure to open file
 */
int load_initial_conditions(char* cond_file, double** u, double** v1, double** v2, double** v3, int nrows, int ncols,
                         bool initial_randomization, float cell_cap_1D, float noise_stddev, float noise_level, gsl_rng* rgen) {

    #pragma omp parallel for collapse(2) schedule(static)
    for (int i=0; i<nrows; ++i)
        for (int j=0; j<ncols; ++j) {
            u[i][j] = 1.0;
            v1[i][j] = 0.0;
            v2[i][j] = 0.0;
            v3[i][j] = 0.0;
        }

    FILE *conds = fopen(cond_file, "r");
    if (!conds)
        return 1;

    char c;
    char buff[256];
    int read;

    double conc;    // concentration
    int cid;        // catalyst id
    int xlo, ylo, xhi, yhi;
    double** species;

    while (!feof(conds)) {

        // set defaults to account for malformed inputs
        species = NULL;
        cid = 1;
        conc = 0.0;
        read = 0;

        c = fgetc(conds);
        switch(c) {

        // catalyst
        case 'C':
            if (!fscanf(conds, "%d", &cid))
                fprintf(stderr, "load_initial_conditions(): No catalyst type found. Defaulting to type %d\n", cid);

            switch(cid) {
            case 1: species = v1; break;
            case 2: species = v2; break;
            case 3: species = v3; break;
            default: perror("load_initial_conditions(): Unknown chemical species found! Ignoring.\n"); break;
            }
            break;

        // fuel
        case 'F':
            species = u;
            break;

        // comment line
        case '#':
            char* foo = fgets(buff, 256, conds);  // stops at 256 chars or newline
            if (foo[0] == '#')                    // getting the compiler to shut up
                break;
            break;

        }

        // read coordinates and concentration values
        if (species) {      // allow for ignoring unkown species

            read += fscanf(conds, "%lf", &conc);
            read += fscanf(conds, "%d %d", &xlo, &ylo);
            read += fscanf(conds, "%d %d", &xhi, &yhi);
            if (read < 5)
                perror("load_initial_conditions(): Not all values correctly read. Used default concentration.\n");

            // silently limit out of bounds coordinates to system size
            if (xhi >= nrows) xhi = nrows-1;
            if (yhi >= ncols) yhi = ncols-1;

            // assign area with chemical species + noise
            for (int i=xlo; i<=xhi; ++i)
                for (int j=ylo; j<=yhi; ++j) {

                    species[i][j] = conc;
                    if (initial_randomization)
                        species[i][j] += gauss_noise(1, cell_cap_1D, noise_stddev, noise_level, rgen);

                    if (species[i][j] < 0.0)
                        species[i][j] = 0.0;
                    if (species[i][j] > 1.0)
                        species[i][j] = 1.0;
                }
        }
    }

    return 0;
}