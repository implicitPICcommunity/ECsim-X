/*******************************************************************************************
  Particles3Dcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
  -------------------
developers: Stefano Markidis, Giovanni Lapenta, Diego Gonzalez
 ********************************************************************************************/

#ifndef Part3DCOMM_H
#define Part3DCOMM_H

#include "Particles.h"
#include "NeighboringNodes.h"

/**
 * 
 * Abstract class for particles of the same species, in a 2D space and 3component velocity with communications methods
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3Dcomm:public Particles {
public:
  /** constructor */
  Particles3Dcomm();
  /** destructor */
  ~Particles3Dcomm();
  /** allocate particles */
  void  allocate(bool test, int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid);

  /** calculate the weights given the position of particles */
  void calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid);
  /** Compute the implicit current and the mass matrix */
  void computeMoments(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** Compute the mass matrices for the charge initialization */
  void computeMassMatrices(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** Compute the current (for writting the results) */
  void computeCurrent(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  void computeCurrentBar(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  void computeCharge(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes */
  int communicate(VirtualTopology3D * ptVCT);
  /** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferXright(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferYright(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZleft(double *b_, long long np, VirtualTopology3D * vct);
  /** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
  void bufferZright(double *b_, long long np, VirtualTopology3D * vct);
  /** Delete the a particle from a list(array) and pack the list(array) */
  void del_pack(long long np, long long *nplast);

  /** method to debuild the buffer received */
  int unbuffer(double *b_);

  /** resize the receiving buffer */
  void resize_buffers(int new_buffer_size);
  /** a method to compute how many particles are not in the right domain */
  int isMessagingDone(VirtualTopology3D * ptVCT);
  /** calculate the maximum number exiting from this domain */
  int maxNpExiting();
  /** calculate the weights given the position of particles */
  // void calculateWeights(double*** weight, double xp, double yp, double zp,int ix, int iy, int iz, Grid* grid);
  /** get X-position array for all the particles */
  double *getXall() const;
  /** get Y-position array for all the particles */
  double *getYall() const;
  /** get Z-position array for all the particles */
  double *getZall() const;
  /** get u (X-velocity) array for all the particles */
  double *getUall() const;
  /** get v (Y-velocity) array for all the particles */
  double *getVall() const;
  /** get w (Z-velocity) array for all the particles */
  double *getWall() const;
  /** get X-position array for all the particles by reference */
  double *& getXref();
  /** get Y-position array for all the particles by reference */
  double *& getYref();
  /** get Z-position array for all the particles by reference */
  double *& getZref();
  /** get u (X-velocity) array for all the particles by reference */
  double *& getUref();
  /** get v (Y-velocity) array for all the particles by reference */
  double *& getVref();
  /** get w (Z-velocity) array for all the particles by reference */
  double *& getWref();
  /** get q array for all the particles by reference */
  double *& getQref();
  /** get the ID array   */
  long long *getParticleIDall() const;
  long long *getParticleRankall() const;
  /** get X-position of particle with label indexPart */
  double getX(long long indexPart) const;
  /** get Y-position of particle with label indexPart */
  double getY(long long indexPart) const;
  /** get Z-position of particle with label indexPart */
  double getZ(long long indexPart) const;
  /** get u (X-velocity) of particle with label indexPart */
  double getU(long long indexPart) const;
  /** get v (Y-velocity) of particle with label indexPart */
  double getV(long long indexPart) const;
  /** get w (Z-velocity) of particle with label indexPart */
  double getW(long long indexPart) const;
  /** get ID of particle with label indexPart */
  long long getParticleID(long long indexPart) const;
  long long getParticleRank(long long indexPart) const;
  /**get charge of particle with label indexPart */
  double getQ(long long indexPart) const;
  /** get charge of array for ID particles */
  double *getQall() const;
  /** get the number of particles of this subdomain */
  long long getNOP() const;
  /** return the Kinetic energy */
  double getKe();
  /** return the Kinetic energy removed from the system*/
  double getKremoved();
  /** return the maximum kinetic energy */
  double getMaxVelocity();
  /** return the total charge */
  double getTotalCharge();
  /** return the total number of particles */
  long long getTotalNumPart();
  /** return energy distribution */
  unsigned long *getVelocityDistribution(int nBins, double maxVel);
  /** return the momentum */
  double getP();
  /** Print particles info: positions, velocities */
  void Print(VirtualTopology3D * ptVCT) const;
  /** Print the number of particles of this subdomain */
  void PrintNp(VirtualTopology3D * ptVCT) const;
  /** Add distributions in this iteration to the total */
  void Add_vDist3D();
  void Write_vDist3D(string SaveDirName);

  double *ubar, *vbar, *wbar;
protected:
  /** number of species */
  int ns;
  /** maximum number of particles of this species on this domain. used for memory allocation */
  long long npmax;
  /** number of particles of this species on this domain */
  long long nop;
  /** total number of particles */
  long long np_tot;
  /** number of particles per cell */
  int npcel;
  /** number of particles per cell - X direction */
  int npcelx;
  /** number of particles per cell - Y direction */
  int npcely;
  /** number of particles per cell - Z direction */
  int npcelz;
  /** charge to mass ratio */
  double qom;
  /** recon thick */
  double delta;
  /** thermal velocity  - Direction X*/
  double uth;
  /** thermal velocity  - Direction Y*/
  double vth;
  /** thermal velocity  - Direction Z*/
  double wth;
  /** u0 Drift velocity - Direction X */
  double u0;
  /** v0 Drift velocity - Direction Y */
  double v0;
  /** w0 Drift velocity - Direction Z */
  double w0;
  /** Positions arra - X component */
  double *x;
  /** Positions array - Y component */
  double *y;
  /** Positions array - Z component */
  double *z;
  /** Velocities array - X component */
  double *u;
  /** Velocities array - Y component */
  double *v;
  /** Velocities array - Z component */
  double *w;

  /** Rank */
  int rank;

  /** TrackParticleID */
  bool TrackParticleID;
  /** Local particle ID, number of process where it was created and number of paricles created by this process*/
  long long *partID;
  long long *partRank;
  long long npTracked;
  /** number of variables to be stored in buffer for communication for each particle  */
  int nVar;
  /** Charge array */
  double *q;

  /** Initial charge density */
  double rhoINIT;
  /** Injection charge density */
  double rhoINJECT;

  /** Simulation domain lengths */
  double xstart, xend, ystart, yend, zstart, zend, invVOL;
  /** time step */
  double dt;
  /** Lx = simulation box length - x direction   */
  double Lx;
  /** Ly = simulation box length - y direction   */
  double Ly;
  /** Lz = simulation box length - z direction   */
  double Lz;
  /** center of object immersed - x direction */
  double x_center;
  /** center of object immersed - y direction */
  double y_center;
  /** center of object immersed - z direction */
  double z_center;
  /** size of immersed box - cube */
  double L_square;
  /** grid spacings */
  double dx, dy, dz;
  /** number of grid 
          nodes */
  int nxn, nyn, nzn;
  /** buffers for communication */
  /** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
  int buffer_size;
  /** smaller buffer size */
  int buffer_size_small;
  /** buffer with particles going to the right processor - Direction X */
  double *b_X_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_X_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction X */
  double *b_X_LEFT;
  /** pointer to the buffer for resizing */
  double *b_X_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Y */
  double *b_Y_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Y_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Y */
  double *b_Y_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Y_LEFT_ptr;
  /** buffer with particles going to the right processor - Direction Z */
  double *b_Z_RIGHT;
  /** pointer to the buffer for resizing */
  double *b_Z_RIGHT_ptr;
  /** buffer with particles going to the left processor - Direction Z */
  double *b_Z_LEFT;
  /** pointer to the buffer for resizing */
  double *b_Z_LEFT_ptr;

  /** number of particles exiting per cycle*/
  int npExitXright;
  /** number of particles exiting to X-LEFT per cycle*/
  int npExitXleft;
  /** number of particles exiting to Y-RIGHT per cycle*/
  int npExitYright;
  /** number of particles exiting to Y-LEFT per cycle*/
  int npExitYleft;
  /** number of particles exiting to Z-RIGHT per cycle*/
  int npExitZright;
  /** number of particles exiting to Z-LEFT per cycle*/
  int npExitZleft;
  /** total number of particles exiting per cycle */
  int npExit;
  /** number of particles not in the right domain   */
  int rightDomain;


  /** bool for communication verbose */
  bool cVERBOSE;
  /** Boundary condition on particles:
          <ul>
          <li>0 = exit</li>
          <li>1 = perfect mirror</li>
          <li>2 = riemission</li>
          <li>3 = periodic condition </li>
          </ul>
          */
  /** Boundary Condition Particles: FaceXright */
  int bcPfaceXright;
  /** Boundary Condition Particles: FaceXleft */
  int bcPfaceXleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceYright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceYleft;
  /** Boundary Condition Particles: FaceYright */
  int bcPfaceZright;
  /** Boundary Condition Particles: FaceYleft */
  int bcPfaceZleft;
  /** speed of light in vacuum */
  double c;
  /** restart variable for loading particles from restart file */
  int restart;
  /** Number of iteration of the mover*/
  int NiterMover;
  /** velocity of the injection of the particles */
  double Vinj;
  /** removed charge from species */
  double Q_removed;
  /** density of the injection of the particles */
  double Ninj;
  /** kinetic energy removed */
  double K_removed;

  NeighboringNodes NeNo;

  /** Seed for random numbers generator */
  long long seed;
  /** Cylindrical coordinates */
  int Cylindrical;
  /** Number of dimension of the problem */
  int ndim;
  /** Lagrangian reference frame */
  int Lagrangian;
  /** Bulk velocity (for the Lagrangian reference frame */
  double u_bulk, v_bulk, w_bulk;
};

#endif
