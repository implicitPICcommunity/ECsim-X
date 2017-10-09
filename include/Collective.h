/***************************************************************************
  Collective.h  -  Stefano Markidis, Giovanni Lapenta
  -------------------------------------------------------------------------- */


/*! Collective properties. Input physical parameters for the simulation.  Use ConfigFile to parse the input file @date Wed Jun 8 2011 @par Copyright: (C) 2011 K.U. LEUVEN @author Pierre Henri, Stefano Markidis @version 1.0 */

#ifndef Collective_H
#define Collective_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ConfigFile.h"
#include "input_array.h"
#include "hdf5.h"
#include "mpi.h"

#ifdef BATSRUS
#include "InterfaceFluid.h"
#endif

using std::cout;
using std::endl;
using std::ofstream;
using namespace std;

class Collective
#ifdef BATSRUS
: public InterfaceFluid
#endif
  {
  public:
    /*! constructor: initialize physical parameters with values */
    Collective(int argc, char **argv);
    /*! destructor */
    ~Collective();
    /*! read input file */
    void ReadInput(string inputfile);
    /*! read the restart input file from HDF5 */
    int ReadRestart(string inputfile);
    /*! Print physical parameters */
    void Print();
    /*! save setting in a file */
    void save();
    /*! get the physical space dimensions */
    int getDim();
    /*! Get length of the system - direction X */
    double getLx();
    /*! Get length of the system - direction Y */
    double getLy();
    /*! Get length of the system - direction Z */
    double getLz();
    /*! Get object center - direction X */
    double getx_center();
    /*! Get object center - direction Y */
    double gety_center();
    /*! Get object center - direction Z */
    double getz_center();
    /*! Get object size - cubic box */
    double getL_square();
    /*! Get object center for injection - direction X */
    double getx_center_inj(int ind);
    /*! Get object center for injection - direction Y */
    double gety_center_inj(int ind);
    /*! Get object center for injection - direction Z */
    double getz_center_inj(int ind);
    /*! Get object size for injection - cubic box */
    double getL_square_inj(int ind);
    /** getgeometry coild - Coil Diameter */
    double getcoilD();
    /** getgeometry coild - Coil Spacing */
    double getcoilSpacing();
    /** getgeometry coild - Coil thickness */
    double getcoilThickness();
    /** Get cylindrical coordinates */
    int getCylindrical();
    /** Get lagrangian reference frame */
    int getLagrangian();
    /** Get the bulk velocity */
    double getUbulk();
    double getVbulk();
    double getWbulk();

    /*! Get the number of cells - direction X */
    int getNxc();
    /*! Get the number of cells - direction Y */
    int getNyc();
    /*! Get the number of cells - direction Z */
    int getNzc();

    int  getXLEN()      {return (XLEN);};
    int  getYLEN()      {return (YLEN);};
    int  getZLEN()      {return (ZLEN);};

    bool getPERIODICX() {return (PERIODICX);};
    bool getPERIODICY() {return (PERIODICY);};
    bool getPERIODICZ() {return (PERIODICZ);};

    /*! Get the grid spacing - direction X */
    double getDx();
    /*! Get the grid spacing - direction Y */
    double getDy();
    /*! Get the grid spacing - direction Z */
    double getDz();
    /*! get the light speed */
    double getC();
    /*! get the time step */
    double getDt();
    /*! get the decentering parameter */
    double getTh();
    /*! get the Smoothing value */
    double getSmooth();
    /*! get the number of time cycles */
    int getNcycles();
    /*! get the number of species */
    int getNs();
    /*! get the number of particles for different species */
    long getNp(int nspecies);
    /*! get the number of particles per cell */
    int getNpcel(int nspecies);
    /*! get the number of particles per cell - direction X */
    int getNpcelx(int nspecies);
    /*! get the number of particles per cell - direction Y */
    int getNpcely(int nspecies);
    /*! get the number of particles per cell - direction Z */
    int getNpcelz(int nspecies);
    /*! get maximum number of particles for different species */
    long getNpMax(int nspecies);
    /*! NpMax/Np is the ratio between the maximum number of particles allowed on a processor and the number of particles */
    double getNpMaxNpRatio();
    /*! get charge to mass ratio for different species */
    double getQOM(int nspecies);
    /*! get background charge for GEM challenge */
    double getRHOinit(int nspecies);
    /*! get rho injection */
    double getRHOinject(int nspecies);
    /*! get thermal velocity - X direction */
    double getUth(int nspecies);
    /*! get thermal velocity - Y direction */
    double getVth(int nspecies);
    /*! get thermal velocity - Z direction */
    double getWth(int nspecies);
    /*! get Drift velocity - Direction X */
    double getU0(int nspecies);
    /*! get Drift velocity - Direction Y */
    double getV0(int nspecies);
    /*! get Drift velocity - Direction Z */
    double getW0(int nspecies);

    int getTP_Ns();
    /*! Test particles: get the number of particles for different species */
    long getTP_Np(int nspecies);
    /*! Test particles: get the number of particles per cell */
    int getTP_Npcel(int nspecies);
    /*! Test particles: get the number of particles per cell - direction X */
    int getTP_Npcelx(int nspecies);
    /*! Test particles: get the number of particles per cell - direction Y */
    int getTP_Npcely(int nspecies);
    /*! Test particles: get the number of particles per cell - direction Z */
    int getTP_Npcelz(int nspecies);
    /*! Test particles: get maximum number of particles for different species */
    long getTP_NpMax(int nspecies);
    /*! Test particles: get charge to mass ratio for different species */
    double getTP_QOM(int nspecies);
    /*! Test particles: get background charge for GEM challenge */
    double getTP_RHOinit(int nspecies);
    /*! Test particles: get thermal velocity - X direction */
    double getTP_Uth(int nspecies);
    /*! Test particles: get thermal velocity - Y direction */
    double getTP_Vth(int nspecies);
    /*! Test particles: get thermal velocity - Z direction */
    double getTP_Wth(int nspecies);
    /*! Test particles: get Drift velocity - Direction X */
    double getTP_U0(int nspecies);
    /*! Test particles: get Drift velocity - Direction Y */
    double getTP_V0(int nspecies);
    /*! Test particles: get Drift velocity - Direction Z */
    double getTP_W0(int nspecies);
    /*! Test particles: Get object center - direction X */
    double getTP_x_center(int is);
    /*! Test particles: Get object center - direction Y */
    double getTP_y_center(int is);
    /*! Test particles: Get object center - direction Z */
    double getTP_z_center(int is);
    /*! Test particles: Get object size - cubic box */
    double getTP_L_square(int is);

    /*! get the boolean value for TrackParticleID */
    bool getTrackParticleID(int nspecies);
    /*! get SaveDirName */
    string getSaveDirName();
    /*! get last_cycle */
    int getLast_cycle();
    /*! get RestartDirName */
    string getRestartDirName();

    /*! get Case type */
    string getCase();
    /*! get particle initialization type */
    string getPartInit();
    /*! get particle injection type */
    string getPartInj();
    /*! get output writing method */
    string getWriteMethod();
    /*! get simulation name */
    string getSimName();
    /*! get probes file name */
    string getProbesFile();
    /*! get poisson correction flag */
    string getPoissonCorrection();

    /*! get initial solution flag */
    bool getSolInit();
    /*! get static fields flag */
    bool getStaticField();

    /*! get Boundary Condition Particles: FaceXright */
    int getBcPfaceXright();
    /*! get Boundary Condition Particles: FaceXleft */
    int getBcPfaceXleft();
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceYright();
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceYleft();
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceZright();
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceZleft();

    /*! get Boundary Condition Electrostatic Potential: FaceXright */
    int getBcPHIfaceXright();
    /*! get Boundary Condition Electrostatic Potential:FaceXleft */
    int getBcPHIfaceXleft();
    /*! get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceYright();
    /*! get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceYleft();
    /*! get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceZright();
    /*! get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceZleft();

    /*! get Boundary ConditionElectric Field: FaceXright */
    int getBcEMfaceXright();
    /*! get Boundary Condition Electric Field: FaceXleft */
    int getBcEMfaceXleft();
    /*! get Boundary Condition Electric Field: FaceYright */
    int getBcEMfaceYright();
    /*! get Boundary Condition Electric Field: FaceYleft */
    int getBcEMfaceYleft();
    /*! get Boundary Condition Electric Field: FaceZright */
    int getBcEMfaceZright();
    /*! get Boundary Condition Electric Field: FaceZleft */
    int getBcEMfaceZleft();

    /*! get RESTART */
    int getRestart_status();

    /*! get the sheet thickness */
    double getDelta();
    /*! get the amplitude of the magnetic field along x */
    double getB0x();
    /*! get the amplitude of the magnetic field along y */
    double getB0y();
    /*! get the amplitude of the magnetic field along z */
    double getB0z();
    /*! get the amplitude of the magnetic field 1 along x */
    double getB1x();
    /*! get the amplitude of the magnetic field 1 along y */
    double getB1y();
    /*! get the amplitude of the magnetic field 1 along z */
    double getB1z();

    /*! get the boolean value for verbose results */
    bool getVerbose();

    /*! get the velocity of injection of the plasma from the wall */
    double getVinj();

    /*! get the converging tolerance for CG solver */
    double getCGtol();
    /*! get the converging tolerance for GMRES solver */
    double getGMREStol();
    /*! get the numbers of iteration for the PC mover */
    int getNiterMover();

    /*! output of fields */
    int getFieldOutputCycle();
    /*! output of particles */
    int getParticlesOutputCycle();
    /*! output of test particles */
    int getTestParticlesOutputCycle();
    /*! output of restart */
    int getRestartOutputCycle();
    /*! output of diagnostics */
    int getDiagnosticsOutputCycle();

    /*! Boundary condition selection for BCFace for the electric field components */
    int bcEx[6], bcEy[6], bcEz[6];
    /*! Boundary condition selection for BCFace for the magnetic field components */
    int bcBx[6], bcBy[6], bcBz[6];

    /*! get initfile */
    string getinitfile();

    int getNdim() { return ndim;}

    int getZeroCurrent();

	  double getMemPartGb();

    string getLambdaDamping();
  private:
    /*! inputfile */
    string inputfile;
    /*! Restart HDF5 file */
    string initfile;
    /*! light speed */
    double c;
    /*! 4 pi */
    double fourpi;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    /*! number of time cycles */
    int ncycles;
    /*! physical space dimensions */
    int dim;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;
    /*! object center - X direction */
    double x_center;
    /*! object center - Y direction */
    double y_center;
    /*! object center - Z direction */
    double z_center;
    /*! object size - assuming a cubic box */
    double L_square;
    /*! object size - assuming a cubic box */
    double *L_square_inj;
    /*! object center - X direction */
    double *x_center_inj;
    /*! object center - Y direction */
    double *y_center_inj;
    /*! object center - Z direction */
    double *z_center_inj;
    /** geometry coil - Coil Diameter */
    double coilD;
    /** geometry coil - Coil Spacing */
    double coilSpacing;
    /** geometry coil - Coil Thickness */
    double coilThickness;
    /*! number of cells - X direction */
    int nxc;
    /*! number of cells - Y direction */
    int nyc;
    /*! number of cells - Z direction */
    int nzc;
    /*! grid spacing - X direction */
    double dx;
    /*! grid spacing - Y direction */
    double dy;
    /*! grid spacing - Z direction */
    double dz;

    int XLEN;          /*! Number of subdomains in the X direction */
    int YLEN;          /*! Number of subdomains in the Y direction */
    int ZLEN;          /*! Number of subdomains in the Z direction */
    bool PERIODICX;    /*! Periodicity in the X direction */
    bool PERIODICY;    /*! Periodicity in the Y direction */
    bool PERIODICZ;    /*! Periodicity in the Z direction */

    /*! number of species */
    int ns;
    /*! number of particles per cell */
    int *npcel;
    /*! number of particles per cell - X direction */
    int *npcelx;
    /*! number of particles per cell - Y direction */
    int *npcely;
    /*! number of particles per cell - Z direction */
    int *npcelz;
    /*! number of particles array for different species */
    long *np;
    /*! maximum number of particles array for different species */
    long *npMax;
    /*! max number of particles */
    double NpMaxNpRatio;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! charge to mass ratio array for different species */
    double *rhoINIT;
    /*! density of injection */
    double *rhoINJECT;
    /*! thermal velocity - Direction X */
    double *uth;
    /*! thermal velocity - Direction Y */
    double *vth;
    /*! thermal velocity - Direction Z */
    double *wth;
    /*! Drift velocity - Direction X */
    double *u0;
    /*! Drift velocity - Direction Y */
    double *v0;
    /*! Drift velocity - Direction Z */
    double *w0;

    /*! Test particles: number of species */
    int TP_ns;
    /*! Test particles: number of particles per cell */
    int *TP_npcel;
    /*! Test particles: number of particles per cell - X direction */
    int *TP_npcelx;
    /*! Test particles: number of particles per cell - Y direction */
    int *TP_npcely;
    /*! Test particles: number of particles per cell - Z direction */
    int *TP_npcelz;
    /*! Test particles: number of particles array for different species */
    long *TP_np;
    /*! Test particles: maximum number of particles array for different species */
    long *TP_npMax;
    /*! Test particles: charge to mass ratio array for different species */
    double *TP_qom;
    /*! Test particles: charge to mass ratio array for different species */
    double *TP_rhoINIT;
    /*! Test particles: thermal velocity - Direction X */
    double *TP_uth;
    /*! Test particles: thermal velocity - Direction Y */
    double *TP_vth;
    /*! Test particles: thermal velocity - Direction Z */
    double *TP_wth;
    /*! Test particles: Drift velocity - Direction X */
    double *TP_u0;
    /*! Test particles: Drift velocity - Direction Y */
    double *TP_v0;
    /*! Test particles: Drift velocity - Direction Z */
    double *TP_w0;
    /*! object center - X direction */
    double *TP_x_center;
    /*! object center - Y direction */
    double *TP_y_center;
    /*! object center - Z direction */
    double *TP_z_center;
    /*! object size - assuming a cubic box */
    double *TP_L_square;

    /*! Case type */
    string Case;
    /*! Fields initialization type */
    string FieldsInit;
    /*! Particle initialization type */
    string PartInit;
    /*! Particle injection type */
    string PartInj;
    /*! Output writing method */
    string wmethod;
    /*! Simulation name */
    string SimName;
    /*! Poisson correction flag */
    string PoissonCorrection;

    /*! HDF5 initial solution flag */
    bool SolInit;
    /*! Static field flag */
    bool StaticField;

    /*! TrackParticleID */
    bool *TrackParticleID;
    /*! SaveDirName */
    string SaveDirName;
    /*! RestartDirName */
    string RestartDirName;
    /*! get inputfile */
    string getinputfile();
    /*! restart_status 0 --> no restart; 1--> restart, create new; 2--> restart, append; */
    int restart_status;
    /*! last cycle */
    int last_cycle;

    /*! Boundary condition on particles 0 = exit 1 = perfect mirror 2 = riemission */
    /*! Boundary Condition Particles: FaceXright */
    int bcPfaceXright;
    /*! Boundary Condition Particles: FaceXleft */
    int bcPfaceXleft;
    /*! Boundary Condition Particles: FaceYright */
    int bcPfaceYright;
    /*! Boundary Condition Particles: FaceYleft */
    int bcPfaceYleft;
    /*! Boundary Condition Particles: FaceYright */
    int bcPfaceZright;
    /*! Boundary Condition Particles: FaceYleft */
    int bcPfaceZleft;


    /*! Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic Condition */
    /*! Boundary Condition Electrostatic Potential: FaceXright */
    int bcPHIfaceXright;
    /*! Boundary Condition Electrostatic Potential:FaceXleft */
    int bcPHIfaceXleft;
    /*! Boundary Condition Electrostatic Potential:FaceYright */
    int bcPHIfaceYright;
    /*! Boundary Condition Electrostatic Potential:FaceYleft */
    int bcPHIfaceYleft;
    /*! Boundary Condition Electrostatic Potential:FaceZright */
    int bcPHIfaceZright;
    /*! Boundary Condition Electrostatic Potential:FaceZleft */
    int bcPHIfaceZleft;

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


    /*! GEM Challenge parameters */
    /*! current sheet thickness */
    double delta;
    /* Amplitude of the field */
    double B0x;
    double B0y;
    double B0z;
    double B1x;
    double B1y;
    double B1z;


    /*! boolean value for verbose results */
    bool verbose;
    /*! RESTART */
    bool RESTART1;
    /*! SOLINIT */
    bool SOLINIT1;

    /*! velocity of the injection from the wall */
    double Vinj;

    /*! CG solver stopping criterium tolerance */
    double CGtol;
    /*! GMRES solver stopping criterium tolerance */
    double GMREStol;
    /*! mover predictor correcto iteration */
    int NiterMover;

    /*! Output for field */
    int FieldOutputCycle;
    /*! Output for particles */
    int ParticlesOutputCycle;
    /*! Output for test particles */
    int TestParticlesOutputCycle;
    /*! restart cycle */
    int RestartOutputCycle;
    /*! Output for diagnostics */
    int DiagnosticsOutputCycle;

    /** Cylindrical coordinates */
    int Cylindrical;
    /** Lagrangian frame of reference */
    int Lagrangian;
    /** Number of dimension of the problem */
    int ndim;
    /** Bulk velocity (for the Lagrangian reference frame */
    double u_bulk, v_bulk, w_bulk;

    /** Initial current equal zero */
    int zeroCurrent;

    /** Memory per process for the particles (Gb) . If -1, it is not used */
    double MemPartGb;

    /** Usage of probes */
    string probesFile;
    
    /** Lambda damping */
    string LambdaDamping;
};

#endif
