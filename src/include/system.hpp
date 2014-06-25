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



/** @file system.hpp
 *
 * This header contains classes and functions used to interact with
 * the system (e.g. reading/writing data from files)
 *
 */

// These will be included through the main file before this header
// #include <cstdio>
// #include <cstdlib>
#include <string>

/**
 * Write data matrix to _open_ file, since it needs
 * to be called during the simulation to save data.
 *
 * On error, the matrix content is dumped to stdout (terminal screen on most systems)
 *
 * @param[in] a         The concentration matrix for a species
 * @param[in] n         The size of the matrix
 * @param[out] fout     The output file
 */
void fdump(double** a, int n, FILE *fout) {
    if (!fout) {
        perror("Output file not specified. Dumping to stdout.\n");
        fout = stdout;
    }

    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j)
            fprintf(fout, "%.15le ", a[i][j]);
        fprintf(fout, "\n");
    }

    fprintf(fout, "\n");  // gap between matrices
}


/**
 * Write nrows x ncols data matrix to _open_ file, since it needs
 * to be called during the simulation to save data.
 *
 * The format of the concentration matrix evolution data for each species is a
 * .dat file containing the matrix at each generation, separated by an empty line.
 * The numerical values are stored as text, not binary, for human inspection.
 * 
 * On error, the matrix content is dumped to stdout (terminal screen on most systems)
 *
 * @param[in] a         The concentration matrix for a species
 * @param[in] nrows,ncols The size of the matrix
 * @param[out] fout     The output file
*/
void griddump(double** a, int nrows, int ncols, FILE *fout) {
    if (!fout) {
        perror("Output file not specified. Dumping to stdout.\n");
        fout = stdout;
    }

    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j)
            fprintf(fout, "%.15le ", a[i][j]);
        fprintf(fout, "\n");
    }

    fprintf(fout, "\n"); // gap between matrices
}



/**
 * This class encapsulates a stdio FILE (not std::fstream for less overhead),
 * along with associated operations, for easier use in main program.
 *
 * Notes:
 *  1. no copies are permitted (copy constructor disabled) since this is
 *     meant as one data source per file.
 *
 *  2. the output file is opened for writing, i.e. the file is replaced.
 *  TODO: implement resume function
 */
class datafile {

private:
    FILE *file;
    FILE *endfile;          ///< will keep the last system state for easy resuming

    int nrows, ncols;       ///< size of the matrix
    double** data;
    std::string filepath;

    datafile(const datafile& that);  // hide copy constructor to prevent use

public:

    datafile(const char* filename, double** datamatrix, int rows, int cols, bool append) :
        nrows(rows), ncols(cols), data(datamatrix) {

        if (append)
            file = fopen(filename, "a");
        else
            file = fopen(filename, "w");

        filepath = filename;

        std::string name(filename);
        name = name.substr(0, name.length() - 4) + ".final";
        endfile = fopen(name.c_str(), "w");

        // write data matrix size as first entry in both datafiles
        fprintf(file,  "%d\n", nrows);
        fprintf(endfile,  "%d\n", nrows);
    }


    /// specifies a fallback filename to use in case the first location is not writeable
    datafile(const char* filename, const char* fallback_name, double** datamatrix, int rows, int cols, bool append) :
        nrows(rows), ncols(cols), data(datamatrix) {

            if (append)
                file = fopen(filename, "a");
            else
                file = fopen(filename, "w");

            if (!file) {
                if (append)
                    file = fopen(fallback_name, "a");
                else
                    file = fopen(fallback_name, "w");

                filepath = fallback_name;
            }
            else
                filepath = filename;

            std::string name(filename);
            name = name.substr(0, name.length() - 4) + ".final";
            endfile = fopen(name.c_str(), "w");

            // write data matrix size as first entry in both datafiles
            fprintf(file,  "%d\n", nrows);
            fprintf(endfile,  "%d\n", nrows);
    }

    ~datafile() {
        fclose(file);
        fclose(endfile);
    }
    
    /// Save matrix to .final file - meant to save last system state
    void dump_final() {
        griddump(data, nrows, ncols, this->endfile);
        fflush(this->endfile);
    }

    void dump() {
        griddump(data, nrows, ncols, this->file);
        fflush(this->file);
    }

    // test for successful access to file
    bool isopen() { return (file != NULL); }


    std::string get_path() {
        return filepath;
    }

};



/**
 * Load the state (concentration distribution) of a chemical species from file.
 *
 * @param[out] a The square concentration matrix for a species
 * @param[in]  nrows,ncols The size of the matrix
 * @param[in] filepath
 * //@param[in] from_end The procedure reads the concentration matrix starting at the end of the file.
 *                     This is used when resuming from outputs of previous simulations.
 *
 * @return true if no errors occured during read, else false.
 *
 * Notes:
 *   If the size of the data (nrows or ncols) in filepath is less than that of the matrix a,
 *   the remaining entries are set to zero, and error signalled with return parameter.
 */
bool load_state(double** a, int nrows, int ncols, const char* filepath) {

    bool no_errs = true;

    FILE *sfile = fopen(filepath, "r");
    if (!sfile) {
        perror("Could not open state file. Data left unaltered.\n");
        return false;
    }

    // read first matrix of file
    for (int i=0; i<nrows; ++i) {
        for (int j=0; j<ncols; ++j) {
            if (fscanf(sfile, "%le", &a[i][j]) < 1)
                a[i][j] = 0.0;  // set to null on error
                no_errs = false;
            }
        }

    fclose(sfile);

    // DEPRECATED: read last matrix of file
    /*
    else {
        int clen = sizeof(char);       // assuming file uses char of same size
        char c;
        long cpos = 0;

        // get behind the final newline: undefined where SEEK_END is
        fseek(sfile, -clen, SEEK_END);
        cpos = ftell(sfile);           // remember attempt

        if (fgetc(sfile) != '\n') {
            fseek(sfile, -1, cpos);    // try one byte back from prev attempt
        }


        for (int no_nl=0; no_nl<nrows; ++no_nl) {   // need to read nrows

            c = fgetc(sfile);

            while (c != '\n') {
                fseek(sfile, -2 * clen, SEEK_END);
                cpos = ftell(sfile);

            }
        }
    }
    */

    return no_errs;
}

