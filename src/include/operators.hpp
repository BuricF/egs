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


/** @file operators.hpp
 * 
 * This header contains the implementations of "operators" 
 * (term very loosely defined),
 * which are used in the evolution rules of the RD system
 * or to modify the system state in some way.
 */


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Discrete Laplace operators  ========================

/**
 * Discrete Laplacian with 5-point method and periodic boundaries.
 * 
 * Note on implementation: 
 *   Boundary tests done explicitly (not with modulos) for speed
 *   since this is the simulation bottleneck (profiled).
 *   The code is cumbersome but faster.
 * 
 *
 * @param[in] u         Concentration matrix
 * @param[in] i,j       Grid coordinates, i: M, j:N
 * @param[in] M,N       Matrix size
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$
 * @return      Laplacian approximation
 */
inline double laplacian5_per_grid(double** u, int i, int j, int M, int N, float hsq) {

    double lap = -4*u[i][j];

    //u[i-1][]
    if (i == 0)
        lap += u[M-1][j];
    else
        lap += u[i-1][j];

    //u[i+1][]
    if (i == M-1)
        lap += u[0][j];
    else
        lap += u[i+1][j];

    // u[][j-1]
    if (j == 0)
        lap += u[i][N-1];
    else
        lap += u[i][j-1];

    // u[][j+1]
    if (j == N-1)
        lap += u[i][0];
    else
        lap += u[i][j+1];

    return lap / hsq;
}


/**
 * Discrete Laplacian with 5-point method and periodic boundaries:
 * square matrix version
 * 
 * The square matrix case was defined as a separate function as an optimization
 * to reduce the overhead of a redundant argument.
 * 
 * Note on implementation: 
 *   Boundary tests done explicitly (not with modulos) for speed
 *   since this is the simulation bottleneck (profiled).
 *   The code is cumbersome but faster.
 * 
 *
 * @param[in] u         Concentration matrix
 * @param[in] i,j       Grid coordinates
 * @param[in] N       Matrix size
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$
 * @return      Laplacian approximation
 */
inline double laplacian5_per_grid(double** u, int i, int j, int N, float hsq) {

    double lap = -4*u[i][j];

    //u[i-1][]
    if (i == 0)
        lap += u[N-1][j];
    else
        lap += u[i-1][j];

    //u[i+1][]
    if (i == N-1)
        lap += u[0][j];
    else
        lap += u[i+1][j];

    // u[][j-1]
    if (j == 0)
        lap += u[i][N-1];
    else
        lap += u[i][j-1];

    // u[][j+1]
    if (j == N-1)
        lap += u[i][0];
    else
        lap += u[i][j+1];

    return lap / hsq;
}



/**
 * Discrete Laplacian with 9-point method and periodic boundaries.
 * 
 * Note on implementation: 
 *   Boundary tests done explicitly (not with modulos) for speed
 *   since this is the simulation bottleneck (profiled).
 *   The code is cumbersome but faster.
 * 
 *
 * @param[in] u         Concentration matrix
 * @param[in] i,j       Grid coordinates
 * @param[in] M,N       Matrix size
 * @param[in] hsq       Spatial discretisation step squared \f$ h^2 = {\Delta x}^2 = \frac{L}{N}^2 \f$
 * @return      Laplacian approximation
 */
inline double laplacian9_per(double** u, int i, int j, int M, int N, float hsq) {

    if (i != 0) {

        if (i != M-1) {  // i: interior

            if (j != 0) {

                if (j != N-1)       // j: interior
                    return (4*u[i-1][j] + u[i-1][j-1] + 4*u[i][j-1] + u[i+1][j-1] + 4*u[i+1][j] + u[i+1][j+1] + 4*u[i][j+1] + u[i-1][j+1] - 20*u[i][j]) / (6*hsq*hsq);

                else                    // j: last colunmn
                    return (4*u[i-1][N-1] + u[i-1][N-2] + 4*u[i][N-2] + u[i+1][N-2] + 4*u[i+1][N-1] + u[i+1][0] + 4*u[i][0] + u[i-1][0] - 20*u[i][N-1]) / (6*hsq*hsq);

            }
            else {        // j: first colunmn
                return (4*u[i-1][0] + u[i-1][N-1] + 4*u[i][N-1] + u[i+1][N-1] + 4*u[i+1][0] + u[i+1][1] + 4*u[i][1] + u[i-1][1] - 20*u[i][0]) / (6*hsq*hsq);
            }

        }
        else {  // i: last line
            if (j != 0) {

                if (j != N-1)       // j: interior
                    return (4*u[M-2][j] + u[M-2][j-1] + 4*u[M-1][j-1] + u[0][j-1] + 4*u[0][j] + u[0][j+1] + 4*u[M-1][j+1] + u[M-2][j+1] - 20*u[M-1][j]) / (6*hsq*hsq);

                else                     // j: last colunmn
                    return (4*u[M-2][N-1] + u[M-2][N-2] + 4*u[M-1][N-2] + u[0][N-2] + 4*u[0][N-1] + u[0][0] + 4*u[M-1][0] + u[M-2][0] - 20*u[M-1][N-1]) / (6*hsq*hsq);

            }
            else {        // j: first colunmn
                return (4*u[M-2][0] + u[M-2][N-1] + 4*u[M-1][N-1] + u[0][N-1] + 4*u[0][0] + u[0][1] + 4*u[M-1][1] + u[M-2][1] - 20*u[M-1][0]) / (6*hsq*hsq);
            }
        }

    }

    else {      // i: first line
        if (j != 0) {

                if (j != N-1)       // j: interior
                    return (4*u[M-1][j] + u[M-1][j-1] + 4*u[0][j-1] + u[1][j-1] + 4*u[1][j] + u[1][j+1] + 4*u[0][j+1] + u[M-1][j+1] - 20*u[0][j]) / (6*hsq*hsq);

                else                    // j: last colunmn
                    return (4*u[M-1][N-1] + u[M-1][N-2] + 4*u[0][N-2] + u[1][N-2] + 4*u[1][N-1] + u[1][0] + 4*u[0][0] + u[M-1][0] - 20*u[0][N-1]) / (6*hsq*hsq);

            }
        else {        // j: first colunmn
                return (4*u[M-1][0] + u[M-1][N-1] + 4*u[0][N-1] + u[1][N-1] + 4*u[1][0] + u[1][0] + 4*u[0][0] + u[M-1][0] - 20*u[0][0]) / (6*hsq*hsq);
        }
    }


}



// Noise operators ===================

/**
 * Gaussian noise with sigma = noise_level * sqrt(target_value) and mean = 0
 *
 * @param[in] target    Base value of the concentration to perturb with noise
 * @param[in] cell_cap  Square root of grid cell particle capacity
 * @param[in] noise_var Variance of gaussian noise
 * @param[in] level     Noise level present in the system
 * @param[in] rgen      GSL random number generator
 * @return      Amount of noise
 */
inline double gauss_noise(double target, float cell_cap_1D, float noise_stddev, float level, gsl_rng* rgen) {

    return level * (sqrt(target)/cell_cap_1D) * gsl_ran_gaussian_ziggurat(rgen, noise_stddev);
}


/**
 * Uniform linear noise in the range (-target*level, +target*level)
 * (uniform: range as function of target)
 * 
 * Note: used for some tests but rand() is not to be trusted.
 *
 * @param[in] target >=0 Base value of the concentration to perturb with noise
 * @param[in] level >=0 Noise level present in the system
 * @return Amount of noise
 *
 */
inline double unif_lin_noise(double target, double level) {

    if (rand() % 2)
        return - target * (double)rand()/((double)RAND_MAX/level);
    else
        return target * (double)rand()/((double)RAND_MAX/level);

}



// Concentration operators =============

/**
 * Sets the concentration of a given area to a given value.
 * Any invalid parameters leave the matrix unaltered.
 * 
 * Notes:
 * 1. The limits are virtual coordinates that go beyond array bounds (e.g. (-32, ncols+100))
 *    and which wrap around while looping over area.
 *    This approach is for long wrap-arounds that cause ambiguity in resulting area to be set.
 * 2. As long as type(a) == type(conc), any data type could be passed.
 *
 * @param[out] a        The concentration matrix
 * @param[in] conc      The concentration to set
 * @param[in] xlo,ylo   Lower corner of area to set (inclusive); x: matix lines, y: columns
 * @param[in] xhi,yhi   Higher corner of area to set (inclusive); x: matix lines, y: columns
 * @param[in] nrows, ncols      Size of matrix (rows, columns)
 *
 */
template <typename NumericType>
void set_concentration(NumericType** a, NumericType conc, int xlo, int ylo, int xhi, int yhi, int nrows, int ncols)
{
    int xi, yi;

    for (int i=xlo; i<=xhi; ++i)
        for (int j=ylo; j<=yhi; ++j) {
                        
            if (i>=0)
                xi = i % nrows;
            else
                xi = (i+nrows-1) % nrows; 
            
            if (j>=0)
                yi = j % ncols;
            else
                yi = (j+ncols-1) % ncols;
            
            a[xi][yi] = conc;            
        }            
}


/**
 * Copy concentrations from one matrix to the other (replacing previous values)
 * The caller should check destination matrix dimensions to be no greater than the source.
 *
 * @param[in] src       Source matrix
 * @param[in] nrows,ncols       Source matrix size
 * @param[out] dest     Destination matrix, which should be smaller or equal to the source.
 */
void clone_values(double** src, double** dest, int nrows, int ncols)
{
    for (int i=0; i<nrows; ++i)
        for (int j=0; j<ncols; ++j)
            dest[i][j] = src[i][j];
}


/**
 * Returns total concentration. Written for convenience.
 * 
 * @param[in] mat       Cocentration matrix
 * @param[in] nrows,ncols      Matrix size
 * @return      Total concentration (simple sum)
 */
inline double total_conc(double **mat, int nrows, int ncols)
{
    double tot = 0.0;
    
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (int i=0; i<nrows; ++i)
        for (int j=0; j<ncols; ++j)
            tot += mat[i][j];
    
    return tot;
}




// Mutation operators ===================================

/**
 * Returns the coordinates at which a mutation event should take place,
 * given scheme and bounds. Not an operator in the sense of this file,
 * but seems like the best place to put it.
 *
 * @param[in] mscheme   The mutation scheme
 * @param[in] v         The data matrix of the species underogoing mutation
 * @param[in] xlo,xhi,ylo,yhi   Inclusive bounds in which to place event
 * @param[in] LOWLIM    The effective null concentration (passed since def'd in main file)
 * @return Coordinates as a vector<int> of 2 elements
 */
std::vector<int> get_mutation_coords(char mscheme, double** v, int xlo, int xhi, int ylo, int yhi) 
{
    std::vector<int> coords(2,0);

    switch(mscheme) {

    // mutant at a maximum of the catalyst
    case 'M': {
        // look for a peak
        double cmax = LOWLIM;
        for (int x=xlo; x<=xhi; ++x)
            for (int y=ylo; y<=yhi; ++y)
                if (v[x][y] > cmax) {
                    cmax = v[x][y];
                    coords[0] = x;
                    coords[1] = y;
                }
    }
        break;

    // mutant at a minimum of the catalyst
    case 'm': {
        // look for a peak
        double cmin = 1.0;
        for (int x=xlo; x<xhi; ++x)
            for (int y=ylo; y<yhi; ++y)
                if ((v[x][y] < cmin) && (v[x][y] > LOWLIM)) {
                    cmin = v[x][y];
                    coords[0] = x;
                    coords[1] = y;
                }
    }
        break;

    // mutant at random position
    case 'r': {
        // set up random number generator
        const gsl_rng_type * T;
        T = gsl_rng_taus;
        gsl_rng * rgen  = gsl_rng_alloc(T);
        if (!rgen)
            fprintf(stderr, "Error allocating random number generator.\n");
        gsl_rng_set(rgen, (unsigned)time(0));

        coords[0] = xlo + gsl_rng_uniform_int(rgen, xhi);
        coords[1] = gsl_rng_uniform_int(rgen, yhi);
    }
        break;
    }

    return coords;
}



/**
 * Introduces a mutant "drop" of given parameters into the system at certain coordinates, according to given scheme.
 * 
 * @param[in] wt        Wild type species concentration matrix, upon which the mutant parameters are based
 * @param[in] mut       Mutant species concentration matrix, into which the drop is placed
 * @param[in] M,N       Size of the concentration matrices
 * @param[in] mfract    Concentration of mutant drop or fraction of converted wild type, based on @param m_internal
 * @param[in] msize     The size (Moore radius) of the mutant drop to place 
 * @param[in] mscheme   Way in which to place the drop (M = at maximum peak of wild type,  m = at a minimum of wt, r = random)
 * @param[in] m_internal        "Internal" mutation: convert @param mfract of the wild type to mutant
 */
void mutation_event(double** wt, double** mut, int M, int N, float mfract, int msize, char mscheme, bool m_internal)
{
    float mconc = 0.0;   // mutant concentration to place
    int i,j;             // position of mutant event

    // Get event coordinates by different schemes
    std::vector<int> coords = get_mutation_coords(mscheme, wt, 0, M-1, 0, N-1);
    i = coords[0];
    j = coords[1];

    // Place mutant drop
    if (m_internal)
        mconc = mfract * wt[i][j];
    else
        mconc = mfract;

    mut[i][j] = mconc;
    if (m_internal) {
        // converted catalyst
        wt[i][j] -= mconc;
        if (wt[i][j] < LOWLIM)
            wt[i][j] = 0.0;
    }

    // periodic boundary conditions
    if (msize)
        for (int x= -msize; x<= msize; ++x)
            for (int y= -msize; y<= msize; ++y) {
                
                mut[(i+M+x)%M][(j+N+y)%N] = mconc;
                
                if (m_internal) {
                    // convert catalyst
                    wt[(i+M+x)%M][(j+N+y)%N] -= mconc;
                    if (wt[(i+M+x)%M][(j+N+y)%N] < LOWLIM)
                        wt[(i+M+x)%M][(j+N+y)%N] = 0.0;
                }
            }
}