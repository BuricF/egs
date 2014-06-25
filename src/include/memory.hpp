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



/** @file memory.hpp
 * 
 * This header contains memory management (allocation/freein) functions.
 * Ugly mix of C/C++ because raw data structures were thought to have less overhead,
 * and because earlier on, some pure C mathematical libraries were used.
 */

#include <cstdlib>
#include <vector>


/**
 * Allocates a grid of generic numeric type.
 * 
 * @param[in] a         Pointer to NumericType** matrix
 * @param[in] nrows,ncols       Size of the matrix to allocate
 * @return 0    on success
 * @throw 1     on allocation failure
 */
template <typename NumericType>
int gridalloc(NumericType*** a, int nrows, int ncols) throw(int){
    
    *a = new NumericType*[nrows];

    if (!(*a))
        //return 1;
        throw 1;

    for (int i=0; i<nrows; ++i) {
        (*a)[i] = new NumericType[ncols];

        if (!(*a)) {
            // clean up pointer array
            delete[] *a;
            //return 1;
            throw 1;
        }
    }

    return 0;
}


/**
 * Frees a grid of generic numeric type
 */
template <typename NumericType>
void gridfree(NumericType*** a, int nrows){
    for (int i=0; i<nrows; ++i)
        delete[] (*a)[i];
    delete[] (*a);
}



/**
   Standard pointer swapping function.
   After call, the two supplied pointers will reference the other's matrix
 */
inline void swap_mat(double*** u, double*** v) {
    double **aux = *u;
    *u = *v;
    *v = aux;
}



/**
 *  Allocate a square matrix of side N
*/
int matalloc(double*** a, int N) {
    *a = (double **)malloc(N * sizeof(double *));
    if (!(*a))
        return 0;
    for (int i=0; i<N; ++i) {
        (*a)[i] = (double *)calloc(N, sizeof(double));

        if (!(*a)[i])
            return 0;
    }

    return 1;
}


/**
 *  Free square matrix of side N
*/
void matfree(double*** a, int N) {
    for (int i=0; i<N; ++i)
        free((*a)[i]);
    free(*a);
}
