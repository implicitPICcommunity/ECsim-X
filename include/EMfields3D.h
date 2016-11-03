/*!************************************************************************* 
 * EMfields3D.h - ElectroMagnetic fields definition 
 * ------------------- 
 *  begin : May 2008 
 *  copyright : KU Leuven 
 *  developers : Stefano Markidis, Giovanni Lapenta, Diego Gonzalez
 *  ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H


#include <iostream>
#include <sstream>

#include <math.h>
#include <mpi.h>

#include "hdf5.h"

#include "Alloc.h"
#include "Basic.h"
#include "Grid3DCU.h"
#include "TransArraySpace3D.h"
#include "GMRES.h"
#include "Collective.h"
#include "ComNodes3D.h"
#include "ComInterpNodes3D.h"
#include "asserts.h"
#include "BCStructure.h"
#include "NeighboringNodes.h"

using std::cout;
using std::cerr;
using std::endl;

/*! Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver @date May 2008 @par Copyright: (C) 2008 KUL @author Stefano Markidis, Giovanni Lapenta. @version 3.0 */

class EMfields3D                // :public Field
{
  public:
    /*! constructor */
    EMfields3D(Collective * col, Grid * grid);
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    void init(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoubleGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! init beam */

    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Perpare the memory to calculate E with PETSc */
    void startEcalc(VirtualTopology3D * vct, Collective *col);
    /*! Frees the memory allocated in startEcalc and apply BC to the field */
    void endEcalc(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid);
    /*! Image of Poisson Solver (for SOLVER) */
    void MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Calculate Magnetic and Electric field */
    void calculateEB(Grid*, VirtualTopology3D*, Collective*);
    void productMassE(double*, double*, double*, double***, double***, double***, int, int, int);
    void smoothEth(double value, VirtualTopology3D * vct, Collective *col);
    void smoothJh(double value, VirtualTopology3D * vct, Collective *col);
    void centers2nodesB(VirtualTopology3D*, Grid*, Collective*);

    void UpdateRHOcs(Grid * grid);

    /*! communicate ghost for densities and interp rho from node to center */
    void interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid);
    void communicateGhostP2G(VirtualTopology3D * vct);
    void communicateGhostP2G_Test3(VirtualTopology3D * vct);
    /*! set to 0 all the densities fields */
    void setZeroDensities();
    /*! Sum rhon over species */
    void sumOverSpecies(VirtualTopology3D * vct);

    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction X to current density field at node X,Y,Z */

    void addMxx(double, int, int, int, int);
    void addMxy(double, int, int, int, int);
    void addMxz(double, int, int, int, int);
    void addMyx(double, int, int, int, int);
    void addMyy(double, int, int, int, int);
    void addMyz(double, int, int, int, int);
    void addMzx(double, int, int, int, int);
    void addMzy(double, int, int, int, int);
    void addMzz(double, int, int, int, int);
    void addJxh(double weight[][2][2], int X, int Y, int Z);
    void addJyh(double weight[][2][2], int X, int Y, int Z);
    void addJzh(double weight[][2][2], int X, int Y, int Z);


    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorRightX(double **susxx, double **susyx, double **suszx);
    void sustensorLeftX (double **susxx, double **susyx, double **suszx);
    void sustensorRightY(double **susxy, double **susyy, double **suszy);
    void sustensorLeftY (double **susxy, double **susyy, double **suszy);
    void sustensorRightZ(double **susxz, double **susyz, double **suszz);
    void sustensorLeftZ (double **susxz, double **susyz, double **suszz);


    /*! get Electric Field component X defined on node(indexX,indexY,indexZ) */
    double &getEx(int indexX, int indexY, int indexZ) const;
    /*! get Electric field X component array */
    double ***getEx();
    double ***getExth();
    double ***getEFxth();
    /*! get Electric field X component cell array without the ghost cells */
    double ***getExc();
    /*! get Electric Field component Y defined on node(indexX,indexY,indexZ) */
    double &getEy(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Y component array */
    double ***getEy();
    double ***getEyth();
    double ***getEFyth();
    /*! get Electric field Y component cell array without the ghost cells */
    double ***getEyc();
    /*! get Electric Field component Z defined on node(indexX,indexY,indexZ) */
    double &getEz(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Z component array */
    double ***getEz();
    double ***getEzth();
    double ***getEFzth();
    /*! get Electric field Z component cell array without the ghost cells */
    double ***getEzc();
    /*! get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field X component array */
    double ***getBx();
    /*! get Magnetic field X component cell array without the ghost cells */
    double getBxc(int i, int j, int k);
    /*! get Magnetic field X component cell array without the ghost cells */
    double ***getBxc();
    /*! get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Y component array */
    double ***getBy();
    /*! get Magnetic field Y component cell array without the ghost cells */
    double getByc(int i, int j, int k);
    /*! get Magnetic field Y component cell array without the ghost cells */
    double ***getByc();
    /*! get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Z component array */
    double ***getBz();
    /*! get Magnetic field Z component cell array without the ghost cells */
    double getBzc(int i, int j, int k);
    /*! get Magnetic field Z component cell array without the ghost cells */
    double ***getBzc();
    /*! get density on cell(indexX,indexY,indexZ) */
    double &getRHOc(int indexX, int indexY, int indexZ) const;
    /*! get density array on center cell */
    double ***getRHOc();
    /*! get density on nodes(indexX,indexY,indexZ) */
    double &getRHOn(int indexX, int indexY, int indexZ) const;
    /*! get density array on nodes */
    double ***getRHOn();
    /*! SPECIES: get density on nodes(indexX,indexY,indexZ) */
    double &getRHOns(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get density on center cell(indexX,indexY,indexZ) */
    double &getRHOcs(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get density on center cell(indexX,indexY,indexZ) */
    double ****& getRHOcs();
    /*! SPECIES: get density array on nodes */
    double ****getRHOns();
    /*! SPECIES: get density array on nodes */
    double ***& getRHOns(int is);
    /*! SPECIES: get density array on cells without the ghost cells */
    double ***getRHOcs(int is);
    double ***& getRHOcs(int is, int dummy);

    /** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx_ext(int indexX, int indexY, int indexZ) const;
    /** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy_ext(int indexX, int indexY, int indexZ) const;
    /** get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz_ext(int indexX, int indexY, int indexZ) const;

    /** get Magnetic Field component X */
    double ***getBx_ext();
    /** get Magnetic Field component Y */
    double ***getBy_ext();
    /** get Magnetic Field component Z */
    double ***getBz_ext();

    double ***&getBxTot();
    double ***&getByTot();
    double ***&getBzTot();

    void UpdateFext(int cycle);
    double getFext();

    /*! get Jx(X,Y,Z) */
    double &getJxh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction X */
    double ***getJxh();
    /*! get Jy(X,Y,Z) */
    double &getJyh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Y */
    double ***getJyh();
    /*! get Jz(X,Y,Z) */
    double &getJzh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Z */
    double ***getJzh();
    /*! get Jzs(X,Y,Z,is) */

    /*! get the electric field energy */
    double getEenergy();
    /*! get the magnetic field energy */
    double getBenergy();


    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col);

    /* ********************************* // VARIABLES ********************************* */
  private:
    /*! light speed */
    double c;
    /* 4*PI for normalization */
    double FourPI;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    /*! delt = c*th*dt */
    double delt;
    /*! number of particles species */
    int ns;
    /*! GEM challenge parameters */
    double B0x, B0y, B0z, delta;
    /** Earth Model parameters */
    double B1x, B1y, B1z;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! Boundary electron speed */
    double ue0, ve0, we0;

    // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
    /*! number of cells - X direction, including + 2 (guard cells) */
    int nxc;
    /*! number of nodes - X direction, including + 2 extra nodes for guard cells */
    int nxn;
    /*! number of cell - Y direction, including + 2 (guard cells) */
    int nyc;
    /*! number of nodes - Y direction, including + 2 extra nodes for guard cells */
    int nyn;
    /*! number of cell - Z direction, including + 2 (guard cells) */
    int nzc;
    /*! number of nodes - Z direction, including + 2 extra nodes for guard cells */
    int nzn;
    /*! local grid boundaries coordinate */
    double xStart, xEnd, yStart, yEnd, zStart, zEnd;
    /*! grid spacing */
    double dx, dy, dz, invVOL;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;

    /*! Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ex;
    /*! Exth: implicit electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Exth, ***EFxth;
    /*! Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ey;
    /*! Eyth: implicit electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Eyth, ***EFyth;
    /*! Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
    double ***Ez;
    /*! Ezth: implicit electric field Z-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ezth, ***EFzth;
    /*! Bxc: magnetic field X-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Bxc;
    /*! Byc: magnetic field Y-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Byc;
    /*! Bzc: magnetic field Z-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Bzc;
    /*! Bxn: magnetic field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Bxn;
    /*! Byn: magnetic field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Byn;
    /*! Bzn: magnetic field Z-component (indexX, indexY, indexZ), defined on nodes */
    double ***Bzn;

    // *************************************
    // TEMPORARY ARRAY
    // ************************************
    /*!some temporary arrays (for calculate hat functions) */
    double ***tempXC;
    double ***tempYC;
    double ***tempZC;
    double ***tempXN;
    double ***tempYN;
    double ***tempZN;
    /*! other temporary arrays (in MaxwellSource) */
    double ***tempC;
    double ***tempX;
    double ***tempY;
    double ***tempZ;
    double ***temp2X;
    double ***temp2Y;
    double ***temp2Z;
    /*! and some for MaxwellImage */
    double ***imageEX;
    double ***imageEY;
    double ***imageEZ;
    double ***imageBX;
    double ***imageBY;
    double ***imageBZ;
    double ***vectX;
    double ***vectY;
    double ***vectZ;
    double ***arr;


    // *******************************************************************************
    // *********** SOURCES **
    // *******************************************************************************

    /*! Charge density, defined on central points of the cell */
    double ***rhoc;
    /*! Charge density, defined on nodes */
    double ***rhon;
    /*! Implicit charge density, defined on central points of the cell */
    double ***rhoh;
    /*! SPECIES: charge density for each species, defined on nodes */
    double ****rhons;
    /*! SPECIES: charge density for each species, defined on central points of the cell */
    double ****rhocs;
    /*! Implicit current density X-component, defined on nodes */
    double ***Jxh;
    /*! Implicit current density Y-component, defined on nodes */
    double ***Jyh;
    /*! Implicit current density Z-component, defined on nodes */
    double ***Jzh;
   
    double ***J_partx, ***J_party, ***J_partz;

    /*! External magnetic field component-X, defined on nodes */
    double***  Bx_ext;
    /*! External magnetic field component-Y, defined on nodes */
    double***  By_ext;
    /*! External magnetic field component-Z, defined on nodes */
    double***  Bz_ext;
    /*! Mass matrix */
    double ****Mxx;
    double ****Mxy;
    double ****Mxz;
    double ****Myx;
    double ****Myy;
    double ****Myz;
    double ****Mzx;
    double ****Mzy;
    double ****Mzz;

    NeighboringNodes NeNo; 


    /*! Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic condition */

    /*! Boundary condition for electric field 0 = perfect conductor 1 = magnetic mirror */
    /*! Boundary Condition EM Field: FaceXright */
    int bcEMfaceXright;
    /*! Boundary Condition EM Field: FaceXleft */
    int bcEMfaceXleft;
    /*! Boundary Condition EM Field: FaceYright */
    int bcEMfaceYright;
    /*! Boundary Condition EM Field: FaceYleft */
    int bcEMfaceYleft;
    /*! Boundary Condition EM Field: FaceZright */
    int bcEMfaceZright;
    /*! Boundary Condition EM Field: FaceZleft */
    int bcEMfaceZleft;


    /*! GEM Challenge background ion */
    double *rhoINIT;
    double *rhoINJECT;
    /*! Drift of the species */
    bool *DriftSpecies;

    /*! RESTART BOOLEAN */
    int restart1;
    /*! String with the directory for the restart file */
    string RestartDirName;
    /*! Case */
    string Case;
   
    double Fext;

    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;

};

typedef EMfields3D Field;

#endif
