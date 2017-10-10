/*******************************************************************************************
  FilePartID.h    
  -------------------
  developers: Diego Gonzalez
 ********************************************************************************************/

#ifndef FILEPARTID_H
#define FILEPARTID_H

/**
 * 
 * Class for reading a file with the IDs of several particles of several species
 * @data Mon Oct  9 18:12:07 CEST 2017
 * @author Diego Gonzalez
 * @version 1.0
 *
 */

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "Collective.h"

#define NS_MAX 10
#define NP_MAX 5000

class FilePartID {
public:
    FilePartID(Collective* col);
    ~FilePartID();

    int getNS() { return ns;};
    int*  getSpec() { return spec;};
    long long getNP(int s) { return np[s];};
    long long* getTagR(int s) { return tagR[s];};
    long long* getTagP(int s) { return tagP[s];};

    void print();

private:
    int   ns;  // Number of species
    int        spec[NS_MAX];  // Indices of different species
    long long    np[NS_MAX];  // Number of particles of each species
    long long* tagR[NS_MAX];  // TagR for each species and each particle
    long long* tagP[NS_MAX];  // TagP for each species and each particle
 
};



#endif
