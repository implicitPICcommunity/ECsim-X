/***************************************************************************
  PetscSolver.h  -  Class to solve the Maxwell equations with PETSc
  -------------------
begin                : 2016
copyright            : 2016 KU Luveun
developers           : Diego Gonzalez
 ***************************************************************************/
#ifdef __WITH_PETSC__

#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H


#include <petscdmda.h>
#include <petscmat.h>
#include <petscksp.h>
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "VirtualTopology3D.h"

typedef struct {
    Grid *grid;
    VirtualTopology3D *vct;
    EMfields3D *EMf;
    Collective *col;
  } CtxSolver;

class PetscSolver {
public:
  PetscSolver(EMfields3D *EMf, Grid* grid, VirtualTopology3D *vct, Collective *col);
  ~PetscSolver();

  void solveE();

private:
  // Dimension of the system
  int dim;
  // Global vectors
  Vec globX, globB;
  // Coefficient Matrix 
  Mat A;
  // Solver
  KSP ksp;
  // Preconditioner
  PC pc;
  // Context for the A*x product
  CtxSolver ctx;
};

#endif

extern PetscErrorCode MyMatMult(Mat A,Vec x,Vec y);
#endif
