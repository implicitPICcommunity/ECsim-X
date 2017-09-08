/*******************************************************************************************
  Probes.h    
  -------------------
  developers: Diego Gonzalez
 ********************************************************************************************/

#ifndef PROBES_H
#define PROBES_H

/**
 * 
 * Class for probes. The values of the probes are computed and stored every cycle
 * @data vie sep  8 09:52:20 CEST 2017
 * @author Diego Gonzalez
 * @version 1.0
 *
 */

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Alloc.h"


class Probes {
public:
    Probes(Collective* col, Grid* grid);
    ~Probes();

    int getNP() { return np;};
    double getX(int i) { return x[i];};
    double getY(int i) { return y[i];};
    double getZ(int i) { return z[i];};
    void setValues(EMfields3D* EMf);
    void writeFiles();

    void print();

private:
    void init(Grid* grid);
    void writeHeader();
    double interpolate(int p, double*** field);

private:
    int rank;
    bool active;  
    const char* probesfile;
    char savedir[256];
    int np;     // Number of probes
    int nval;   // Number of values in each probe
    int npt;    // Number of time steps of each probe
    int nptmax; // Maximum number of time steps of each probe
    int npmax;  // Maximum number of probes
    int cycle;  // Current cycle
    double dt;  // time step
    int*   id;  // id of the probes
    double* x;  // x-coordinates of the probes
    double* y;  // y-coordinates of the probes
    double* z;  // z-coordinates of the probes
    double**** W; //Wegiths of each probe for interpolation
    int* indx;  // x-index of the cell where the probes are
    int* indy;  // y-index of the cell where the probes are
    int* indz;  // z-index of the cell where the probes are

    double*** data; // Historical data in each probe
};



#endif
