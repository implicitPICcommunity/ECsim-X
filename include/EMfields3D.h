/*!************************************************************************* 
 * EMfields3D.h - ElectroMagnetic fields definition 
 * ------------------- 
 *  begin : May 2008 
 *  copyright : KU Leuven 
 *  developers : Stefano Markidis, Giovanni Lapenta, Diego Gonzalez
 *  ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H

#define ZERO_J0 2.40482555769577
#define ZERO_J1 3.83170597020751

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

    /* Set all the arrays to zero */
    void setAllzero();

    /*! initialize the electromagnetic fields with constant values */
    void init(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoubleGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initForceFree(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void tangential_discontinuity(VirtualTopology3D * vct, Grid * grid, Collective *col); 
    void mag_KHI(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void ion_cyclotron(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initWB8(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initMirror(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initWB8_6(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initBessel1D(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initBessel2D(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTEmode(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTMmode(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initHalfPlanet(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initPulse(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initFluxRope(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void slabs_instability(VirtualTopology3D * vct, Grid * grid, Collective *col);

    /* Fix boundary conditions in cylindrical coordinate */
    void fixGhostB();
    void fixBC_Cyl_Image(Grid* grid, double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz);
    void fixBC_Cyl_Source(double ***SoX, double ***SoY, double ***SoZ);
    void fixBC_Image(VirtualTopology3D*, double***, double***, double***, double***, double***, double***, double***, double***, double***);
    void fixBC_ImageC(VirtualTopology3D*, double***, double***, double***, double***, double***, double***, double***, double***, double***);
    void fixBC_Source(VirtualTopology3D*, Collective*, double***, double***, double***);
    void fixBC_SourceC(VirtualTopology3D *vct, Collective *col, double ***SoX, double ***SoY, double ***SoZ) ;
    void fixBC_PlanetImage  (VirtualTopology3D *vct, Collective *col, Grid *grid, double***, double***, double***, double***, double***, double***, double***, double***, double***);
    void fixBC_PlanetImageC (VirtualTopology3D *vct, Collective *col, Grid *grid, double***, double***, double***, double***, double***, double***, double***, double***, double***);
    void fixBC_PlanetSource (VirtualTopology3D *vct, Collective *col, Grid *grid, double***, double***, double***);
    void fixBC_PlanetSourceC(VirtualTopology3D *vct, Collective *col, Grid *grid, double ***SoX, double ***SoY, double ***SoZ) ;

    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Perpare the memory to calculate E with PETSc */
    void startEcalc(VirtualTopology3D * vct, Collective *col);
    /*! Frees the memory allocated in startEcalc and apply BC to the field */
    void endEcalc(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid);
    /*! Image of Poisson Solver (for SOLVER) */
    void MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col);
    void RhoImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /*! Calculate Magnetic and Electric field */
    void calculateEB(Grid*, VirtualTopology3D*, Collective*);
    void productMassE(double*, double*, double*, double***, double***, double***, int, int, int);
    double productMassEV(double ***vX, int i, int j, int k, Grid* grid);
    void smoothEth(double value, VirtualTopology3D * vct, Collective *col);
    void smoothE(double value, VirtualTopology3D * vct, Collective *col);
    void smooth(double value, double ***Vx,double ***Vy,double ***Vz,VirtualTopology3D * vct, Collective *col);
    
    void smoothBth(double value, double ***Bxth, double ***Byth, double ***Bzth, VirtualTopology3D * vct, Collective *col);
    void smoothJh(double value, VirtualTopology3D * vct, Collective *col);
    void centers2nodesB(VirtualTopology3D*, Grid*, Collective*);

    void UpdateRHOcs(Grid * grid);

    /*! communicate ghost for densities and interp rho from node to center */
    void interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid);
    void communicateGhostP2G(VirtualTopology3D * vct);
    void communicateGhostP2G_JP(VirtualTopology3D * vct);
    void communicateGhostP2G_Rho(VirtualTopology3D * vct);
    /*! set to 0 all the densities fields */
    void setZeroDensities();
    void setZeroRho();
    void setZeroMass();
    /*! Sum rhon over species */
    void sumOverSpecies(VirtualTopology3D * vct);

    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
  #pragma acc routine seq
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of particles at node X,Y,Z */
  #pragma acc routine seq
    void addN(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of charge at cell X,Y,Z */
    void addQrem(double val, int X, int Y, int Z, int ns);

#pragma acc routine seq
    void addMass(double val[3][3], int, int, int, int);
    void addMassXX(double val, int, int, int, int);
#pragma acc routine seq
    void addJxh(double weight[][2][2], int X, int Y, int Z, int is);
#pragma acc routine seq
    void addJyh(double weight[][2][2], int X, int Y, int Z, int is);
#pragma acc routine seq
    void addJzh(double weight[][2][2], int X, int Y, int Z, int is);
#pragma acc routine seq
  void addJxh_1Dacc(double weight000, double weight001, double weight010, double weight011, double weight100, double weight101, double weight110, double weight111, int X, int Y, int Z, int is);
#pragma acc routine seq
  void addJyh_1Dacc(double weight000, double weight001, double weight010, double weight011, double weight100, double weight101, double weight110, double weight111, int X, int Y, int Z, int is);
#pragma acc routine seq
  void addJzh_1Dacc(double weight000, double weight001, double weight010, double weight011, double weight100, double weight101, double weight110, double weight111, int X, int Y, int Z, int is); 

  void preFetchMoments();
  void PrefetchAfterGatherMoments_extra();

    #pragma acc routine seq
    void addJx(double weight[][2][2], int X, int Y, int Z, int ns);
    #pragma acc routine seq
    void addJy(double weight[][2][2], int X, int Y, int Z, int ns);
    #pragma acc routine seq
    void addJz(double weight[][2][2], int X, int Y, int Z, int ns);

    /*! add an amount of EF - direction X  at node X,Y,Z */
    #pragma acc routine seq
    void addEFx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of EF - direction Y at node X,Y,Z */
    #pragma acc routine seq
    void addEFy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of EF - direction Z at node X,Y,Z */
    #pragma acc routine seq
    void addEFz(double weight[][2][2], int X, int Y, int Z, int is);
    
    /*! add an amount of pressure density - direction XX to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPxx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XY to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPxy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPxz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YY to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPyy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPyz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
    #pragma acc routine seq
    void addPzz(double weight[][2][2], int X, int Y, int Z, int is);

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
    /*! get Electric field X component cell array without the ghost cells */
    double ***getExc();
    /*! get Electric Field component Y defined on node(indexX,indexY,indexZ) */
    double &getEy(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Y component array */
    double ***getEy();
    double ***getEyth();
    /*! get Electric field Y component cell array without the ghost cells */
    double ***getEyc();
    /*! get Electric Field component Z defined on node(indexX,indexY,indexZ) */
    double &getEz(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Z component array */
    double ***getEz();
    double ***getEzth();
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
    /*! SPECIES: get number of particles array on nodes */
    double ***& getNns(int is);
    /*! SPECIES: get density array on cells without the ghost cells */
    double ***getRHOcs(int is);
    double ***& getRHOcs(int is, int dummy);
    double ***getRhotilde(int is);

    double ***getDivE();
    double ***getQrem(int is);
    double ***getTmpX();

    /** get external magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx_ext(int indexX, int indexY, int indexZ) const;
    /** get external magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy_ext(int indexX, int indexY, int indexZ) const;
    /** get external magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz_ext(int indexX, int indexY, int indexZ) const;

    /** get external magnetic Field component X */
    double ***getBx_ext();
    /** get external magnetic Field component Y */
    double ***getBy_ext();
    /** get external magnetic Field component Z */
    double ***getBz_ext();
    /** get external magnetic Field component X. Centers */
    double ***getBxc_ext();
    /** get external magnetic Field component Y. Centers */
    double ***getByc_ext();
    /** get external magnetic Field component Z. Centers */
    double ***getBzc_ext();
    /** get external electric Field component X */
    double ***getEx_ext();
    /** get external electric Field component Y */
    double ***getEy_ext();
    /** get external electric Field component Z */
    double ***getEz_ext();

    double ***&getBxTot();
    double ***&getByTot();
    double ***&getBzTot();

    void UpdateFext(int cycle);
    double getFext();


    /*! get pressure tensor XX for species */
    double ****getpXXsn();
    /*! get pressure tensor XY for species */
    double ****getpXYsn();
    /*! get pressure tensor XZ for species */
    double ****getpXZsn();
    /*! get pressure tensor YY for species */
    double ****getpYYsn();
    /*! get pressure tensor YZ for species */
    double ****getpYZsn();
    /*! get pressure tensor ZZ for species */
    double ****getpZZsn();

    /* get the pressure tensor */
    double ***getPxx(int is);
    double ***getPxy(int is);
    double ***getPxz(int is);
    double ***getPyy(int is);
    double ***getPyz(int is);
    double ***getPzz(int is);

    /*! get Jxh(X,Y,Z) */
    double &getJxh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction X */
    double ***getJxh();
    /*! get Jyh(X,Y,Z) */
    double &getJyh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Y */
    double ***getJyh();
    /*! get Jzh(X,Y,Z) */
    double &getJzh(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Z */
    double ***getJzh();

    /*! get Jxs(X,Y,Z) */
    double &getJxs(int indexX, int indexY, int indexZ, int is) const;
    /*! get current -Direction X */
    double ****getJxs();
    double ***getJxs(int is);
    double ***getJx_ext();
    /*! get Jy(X,Y,Z) */
    double &getJys(int indexX, int indexY, int indexZ, int is) const;
    /*! get current -Direction Y */
    double ****getJys();
    double ***getJys(int is);
    double ***getJy_ext();
    /*! get Jz(X,Y,Z) */
    double &getJzs(int indexX, int indexY, int indexZ, int is) const;
    /*! get current -Direction Z */
    double ****getJzs();
    double ***getJzs(int is);
    double ***getJz_ext();

    /*! SPECIES: get Energy Fluxes  */
    double ***&getEFxs(int is);
    double ***&getEFys(int is);
    double ***&getEFzs(int is);

   
    double getJxcs (int i,int j,int k, int is);
    double getJycs (int i,int j,int k, int is);
    double getJzcs (int i,int j,int k, int is);

    /*! get the electric field energy */
    double getEenergy();
    double getEenergyRemoved(bool save);
    /*! get the magnetic field energy */
    double getBtotenergy();
    double getBextenergy();
    double getBintenergy();


    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col);

    void testCurl(Grid* grid);
    void testDivergence(Grid* grid);
    void computeRhotilde(VirtualTopology3D * vct, Collective *col, Grid* grid, int is);
    void testJdotE(Grid* grid, double deltaK);
    void PoyntingFluxTest();

    void updateLagrangianFrame(VirtualTopology3D * vct, Collective *col);
    void setLagrangianFrame(VirtualTopology3D* vct, Collective* col, bool mode);
    void setZeroCurrent(VirtualTopology3D * vct);

    void chargeFromCurrent(Grid* grid, VirtualTopology3D * vct);

    void SetLambda(Grid *grid, Collective *col);
    double ***GetLambda();

    void setRhoInit(Grid* grid,Collective* col);

    void setFext(int v) { Fext = v; };

    void divergenceOfE(Grid* grid, VirtualTopology3D * vct);

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
    double dx, dy, dz;
    double ***invVOLn, ***invVOLc;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;
    /** source center - X direction   */
    double x_center;
    /** source center - Y direction   */
    double y_center;
    /** source center - Z direction   */
    double z_center;
    /** Characteristic length */
    double L_square;
    /** Coil Parameter: magnetic coil diameter    */
    double coilD;
    /** Coil Parameter: magnetic coil spacing   */
    double coilSpacing;

    /*! Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ex;
    /*! Exth: implicit electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Exth;
    /*! Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ey;
    /*! Eyth: implicit electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Eyth;
    /*! Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
    double ***Ez;
    /*! Ezth: implicit electric field Z-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ezth;
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
    /*! SPECIES: charge density for each species, defined on nodes */
    double ****rhons;
    /*! SPECIES: charge density for each species, defined on central points of the cell */
    double ****rhocs;
    /* Rhotilde: density to compute the charge of each particles to achieve rho when interpolate from particles to nodes */
    double ****rhotilde;
    /*! Implicit current density X-component, defined on nodes */
    double ***Jxh;
    /*! Implicit current density Y-component, defined on nodes */
    double ***Jyh;
    /*! Implicit current density Z-component, defined on nodes */
    double ***Jzh;
    /*! SPECIES: Implicit current density X-component, defined on nodes */
    double ****Jxhs;
    /*! SPECIES: Implicit current density Y-component, defined on nodes */
    double ****Jyhs;
    /*! SPECIES: Implicit current density Z-component, defined on nodes */
    double ****Jzhs;
   
    /*! SPECIES: current density X-component, defined on nodes */
    double ****Jxs;
    /*! SPECIES: current density Y-component, defined on nodes */
    double ****Jys;
    /*! SPECIES: current density Z-component, defined on nodes */
    double ****Jzs;
    /*! SPECIES: number of particles, defined on nodes */
    double ****Nns;
    /*! SPECIES: Energy Flux density component-X for species, defined on nodes */
    double ****EFxs;
    /*! SPECIES: Energy Flux density component-Y for species, defined on nodes */
    double ****EFys;
    /*! SPECIES: Energy Flux component-Z for species, defined on nodes */
    double ****EFzs;


    /** Divergence of E */
    double ***divE;

    /** Charge removed in each cell */
    double ****Qrem;
   
    /*! External magnetic field component-X, defined on nodes */
    double***  Bx_ext;
    /*! External magnetic field component-Y, defined on nodes */
    double***  By_ext;
    /*! External magnetic field component-Z, defined on nodes */
    double***  Bz_ext;
    /*! External magnetic field component-X, defined on centers */
    double***  Bxc_ext;
    /*! External magnetic field component-Y, defined on centers */
    double***  Byc_ext;
    /*! External magnetic field component-Z, defined on centers */
    double***  Bzc_ext;
    /*! External electric field component-X, defined on nodes */
    double***  Ex_ext;
    /*! External electric field component-Y, defined on nodes */
    double***  Ey_ext;
    /*! External electric field component-Z, defined on nodes */
    double***  Ez_ext;

    /*! Initial zero current */
    int zeroCurrent;
    /*! External current: used to asure that Jt = 0 in the initialization */
    double ***Jx_ext;
    double ***Jy_ext;
    double ***Jz_ext;

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

    /*! SPECIES: pressure tensor component-XX, defined on nodes */
    double ****pXXsn;
    /*! SPECIES: pressure tensor component-XY, defined on nodes */
    double ****pXYsn;
    /*! SPECIES: pressure tensor component-XZ, defined on nodes */
    double ****pXZsn;
    /*! SPECIES: pressure tensor component-XZ, defined on nodes */
    double ****pYYsn;
    /*! SPECIES: pressure tensor component-YZ, defined on nodes */
    double ****pYZsn;
    /*! SPECIES: pressure tensor component-ZZ, defined on nodes */
    double ****pZZsn;

    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;

    /** Cylindrical coordinates */
    int Cylindrical;
    /** Bulk velocity (for the Lagrangian reference frame */
    double u_bulk, v_bulk, w_bulk;

    double ***Lambda;
    int damping;

    double Eenergy;
    double Benergy;
    double EenergyRemoved;

};

typedef EMfields3D Field;
extern double besselJ(int alpha, double x);

#endif
