/*******************************************************************************************
  Particles3Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta.
 ********************************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "ComParticles3D.h"
#include "Alloc.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "MPIdata.h"

#include "Particles3Dcomm.h"
//#include "Memory.h"

#include "hdf5.h"
#include <vector>
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-32

#include "MyClock.h"
extern MyClock *clocks;

/**
 * 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles3Dcomm::Particles3Dcomm() {
  // see allocate(int species, Collective* col, VirtualTopology3D* vct, Grid* grid)
  
}

/** deallocate particles */
Particles3Dcomm::~Particles3Dcomm() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
  // deallocate buffers
  delete[]b_X_RIGHT;
  delete[]b_X_LEFT;
  delete[]b_Y_RIGHT;
  delete[]b_Y_LEFT;
  delete[]b_Z_RIGHT;
  delete[]b_Z_LEFT;
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(bool test, int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid) {

  /* initialize random generator with different seed on different processor */
  seed = vct->getCartesian_rank()+ 2;

  // info from collectiveIO
  ns = species;
  if (test) {
    npcel = col->getTP_Npcel(species);
    npcelx = col->getTP_Npcelx(species);
    npcely = col->getTP_Npcely(species);
    npcelz = col->getTP_Npcelz(species);
    qom = col->getTP_QOM(species);
    uth = col->getTP_Uth(species);
    vth = col->getTP_Vth(species);
    wth = col->getTP_Wth(species);
    u0 = col->getTP_U0(species);
    v0 = col->getTP_V0(species);
    w0 = col->getTP_W0(species);
    TrackParticleID = 1;
  } else {
    npcel = col->getNpcel(species);
    npcelx = col->getNpcelx(species);
    npcely = col->getNpcely(species);
    npcelz = col->getNpcelz(species);
    qom = col->getQOM(species);
    uth = col->getUth(species);
    vth = col->getVth(species);
    wth = col->getWth(species);
    u0 = col->getU0(species);
    v0 = col->getV0(species);
    w0 = col->getW0(species);
    TrackParticleID = col->getTrackParticleID(species);
  }

  npTracked = 0,
  dt = col->getDt();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();
  invVOL = 1.0/(dx*dy*dz);
  delta = col->getDelta();

  c = col->getC();

  if (TrackParticleID)
    nVar = 9;
  else
    nVar = 7;

  // This if is necessary to restart with H5hut-io
  if (initnpmax==0){
    nop   = 0;
  } else {
    nop   = initnpmax;
  }
  long ncproc = int(col->getNxc()/col->getXLEN()) *
                int(col->getNyc()/col->getYLEN()) *
                int(col->getNzc()/col->getZLEN());
  npmax = ncproc * npcel * col->getNpMaxNpRatio();


  if (col->getMemPartGb() != -1) {
		double memGb = col->getMemPartGb();
		double mem = memGb*(1000*1024*1024)/col->getNs();
		double szpp = sizeof(double)*(nVar + 0.05*nVar*6/col->getNpMaxNpRatio());
		double npmax_0 = mem/szpp;
		//printf ("With %g Gb I can use %.0f particles\n", memGb, npmax_0);
		if (vct->getCartesian_rank() == 0) {
            printf("\nWARNING: particle memory set in the input file (%.3f Gb). The total memory allocated might not be big enough for holding all the particles"
                   "(npcelx * npcely * npcelz * number of cells)\n\n",memGb);
    }
		npmax = (long long) npmax_0;
  }

  rhoINIT   = col->getRHOinit(species);
  rhoINJECT = col->getRHOinject(species);

  // info for mover
  NiterMover = col->getNiterMover();
  // velocity of the injection from the wall
  Vinj = col->getVinj();
  Ninj = col->getRHOinject(species);
  // info from Grid
  xstart = grid->getXstart();
  xend = grid->getXend();
  ystart = grid->getYstart();
  yend = grid->getYend();
  zstart = grid->getZstart();
  zend = grid->getZend();

  // info from VirtualTopology3D
  cVERBOSE = vct->getcVERBOSE();
  Cylindrical = col->getCylindrical();
  ndim = col->getNdim();

  Lagrangian = col->getLagrangian();
  u_bulk = col->getUbulk();
  v_bulk = col->getVbulk();
  w_bulk = col->getWbulk();

  // boundary condition for particles
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceZright = col->getBcPfaceZright();
  bcPfaceZleft = col->getBcPfaceZleft();
  // //////////////////////////////////////////////////////////////
  // ////////////// ALLOCATE ARRAYS /////////////////////////
  // //////////////////////////////////////////////////////////////
  // BUFFERS
  // the buffer size should be decided depending on number of particles
  // the buffer size should be decided depending on number of particles
  buffer_size = (int) (.05 * npmax * nVar + 1); // max: 5% of the particles in the processors is going out
  buffer_size_small = (int) (.01 * npmax * nVar + 1); // max 1% not resizable 

  if (vct->getCartesian_rank() == 0) {
    long long memreq = (npmax*nVar + buffer_size*6);
    memreq *= sizeof(double);
    printf("---------------------------------------------------\n"
					 " Species %d:\n"
					 "   Maximum number of particles per process: %ld\n"
           "   Memory needed for the particles: %.3f Mb/process\n"
           "---------------------------------------------------\n",species, npmax, ((double) memreq)/(1024*1024));
  }

  //double vm1, res1, vm0, res0;
  //process_mem_usage(vm0, res0);

  // positions
  x = new double[npmax];
  y = new double[npmax];
  z = new double[npmax];
  // velocities
  u = new double[npmax];
  v = new double[npmax];
  w = new double[npmax];
  #ifdef __JDOTE_TEST__
  ubar = new double[npmax];
  vbar = new double[npmax];
  wbar = new double[npmax];
  #endif

  // charge
  q = new double[npmax];
  // ID
  if (TrackParticleID) {
    partID   = new long long[npmax];
    partRank = new long long[npmax];
  }

  b_X_RIGHT = new double[buffer_size];
  b_X_RIGHT_ptr = b_X_RIGHT;    // alias to make the resize
  b_X_LEFT = new double[buffer_size];
  b_X_LEFT_ptr = b_X_LEFT;      // alias to make the resize
  b_Y_RIGHT = new double[buffer_size];
  b_Y_RIGHT_ptr = b_Y_RIGHT;    // alias to make the resize
  b_Y_LEFT = new double[buffer_size];
  b_Y_LEFT_ptr = b_Y_LEFT;      // alias to make the resize
  b_Z_RIGHT = new double[buffer_size];
  b_Z_RIGHT_ptr = b_Z_RIGHT;    // alias to make the resize
  b_Z_LEFT = new double[buffer_size];
  b_Z_LEFT_ptr = b_Z_LEFT;      // alias to make the resize

  // if RESTART is true initialize the particle in allocate method
  restart = col->getRestart_status();
  if (restart != 0) {
    cout << " Restart not implemented." << endl; 
  }

  K_removed = 0;
  //process_mem_usage(vm1, res1);
  //printf("Memmory actually allocated: %g %g\n", vm1-vm0, res1-res0);

}

/** calculate the weights given the position of particles 0,0,0 is the left,left, left node */
void Particles3Dcomm::calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid) {
  double xi[2], eta[2], zeta[2];
  xi[0] = xp - grid->getXN(ix - 1, iy, iz);
  eta[0] = yp - grid->getYN(ix, iy - 1, iz);
  zeta[0] = zp - grid->getZN(ix, iy, iz - 1);
  xi[1] = grid->getXN(ix, iy, iz) - xp;
  eta[1] = grid->getYN(ix, iy, iz) - yp;
  zeta[1] = grid->getZN(ix, iy, iz) - zp;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        weight[i][j][k] = xi[i] * eta[j] * zeta[k] * invVOL;
}



/** Compute the mass matrix and the implicit current */
void Particles3Dcomm::computeMoments(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  double alphaU, alphaV, alphaW;
  double Bxl, Byl, Bzl;
  double beta = 0.5*qom*dt;

  double ***Bx     = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By     = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz     = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());
  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());
  double Fext      = EMf->getFext();

 //#pragma omp parallel
  {
    //#pragma omp for
    for (register long long ip = 0; ip < nop; ip++)
    {
      #ifdef __PROFILING__
      clocks->start(6);
      #endif
      const int ix = 2 + int (floor((x[ip] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[ip] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[ip] - zstart) * inv_dz));
      double weight[2][2][2];
      double temp[2][2][2];
      calculateWeights( weight, x[ip], y[ip], z[ip], ix, iy, iz, grid);

      // B field: Grid -> particle
      Bxl = 0.0;
      Byl = 0.0;
      Bzl = 0.0;
      for (int i=0; i<=1; i++)
        for (int j=0; j<=1; j++)
          for (int k=0; k<=1; k++) {
            Bxl += weight[i][j][k] * (Bx[ix-i][iy-j][iz-k] + Fext*Bx_ext[ix-i][iy-j][iz-k]);
            Byl += weight[i][j][k] * (By[ix-i][iy-j][iz-k] + Fext*By_ext[ix-i][iy-j][iz-k]);
            Bzl += weight[i][j][k] * (Bz[ix-i][iy-j][iz-k] + Fext*Bz_ext[ix-i][iy-j][iz-k]);
          }

      // Rotation matrix alpha:
      double omcx = beta * Bxl;
      double omcy = beta * Byl;
      double omcz = beta * Bzl;
      double oxx = omcx*omcx;
      double oyy = omcy*omcy;
      double ozz = omcz*omcz;
      double oxy = omcy*omcx;
      double oxz = omcx*omcz;
      double oyz = omcy*omcz;
      double denom = 1.0 / (1.0 + oxx + oyy + ozz);
    
      double alpha[3][3];
      alpha[0][0] = ( 1    + oxx)*denom;
      alpha[0][1] = ( omcz + oxy)*denom;
      alpha[0][2] = (-omcy + oxz)*denom;

      alpha[1][0] = (-omcz + oxy)*denom;
      alpha[1][1] = ( 1    + oyy)*denom;
      alpha[1][2] = ( omcx + oyz)*denom;

      alpha[2][0] = ( omcy + oxz)*denom;
      alpha[2][1] = (-omcx + oyz)*denom;
      alpha[2][2] = ( 1    + ozz)*denom;

      double qau = q[ip]*(alpha[0][0]*u[ip] + alpha[0][1]*v[ip] + alpha[0][2]*w[ip]);
      double qav = q[ip]*(alpha[1][0]*u[ip] + alpha[1][1]*v[ip] + alpha[1][2]*w[ip]);
      double qaw = q[ip]*(alpha[2][0]*u[ip] + alpha[2][1]*v[ip] + alpha[2][2]*w[ip]);

      //double vbEx = u[ip];
      //double vbEy = v[ip];
      //double vbEz = w[ip];
      //double edotb = vbEx * omcx + vbEy * omcy + vbEz * omcz;
      //qau = q[ip]*(vbEx + (vbEy * omcz - vbEz * omcy + edotb * omcx)) * denom;
      //qav = q[ip]*(vbEy + (vbEz * omcx - vbEx * omcz + edotb * omcy)) * denom;
      //qaw = q[ip]*(vbEz + (vbEx * omcy - vbEy * omcx + edotb * omcz)) * denom;


      #ifdef __PROFILING__
      clocks->stop(6);
      #endif

//     #ifdef __PROFILING__
//     clocks->start(7);
//     #endif
//     // add charge density               
//     for (int ii = 0; ii < 2; ii++)
//       for (int jj = 0; jj < 2; jj++)
//         for (int kk = 0; kk < 2; kk++) {
//           temp[ii][jj][kk] = q[ip]* weight[ii][jj][kk];
//         }
//     EMf->addRho(temp, ix, iy, iz, ns);
//     #ifdef __PROFILING__
//     clocks->stop(7);
//     #endif

      #ifdef __PROFILING__
      clocks->start(8);
      #endif
      // add implicit current density - X
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = qau* weight[ii][jj][kk];
      EMf->addJxh(temp, ix, iy, iz, ns);

      // add implicit current density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = qav* weight[ii][jj][kk];
      EMf->addJyh(temp, ix, iy, iz, ns);

      // add implicit current density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = qaw* weight[ii][jj][kk];
      EMf->addJzh(temp, ix, iy, iz, ns);
      #ifdef __PROFILING__
      clocks->stop(8);
      #endif

      #ifdef __PROFILING__
      clocks->start(9);
      #endif
      // Mass Matrix:
      //  Each node of the cell in which the particle is:
      for (int i=0; i<2; i++) 
        for (int j=0; j<2; j++) 
          for (int k=0; k<2; k++) {
            int ni = ix-i;
            int nj = iy-j;
            int nk = iz-k;

            // Only half of the possible nodes  (M is symetric)
            for (int c=0; c<14; c++) {
            //for (int c=0; c<27; c++) {
              int n2i = ni + NeNo.getX(c);
              int n2j = nj + NeNo.getY(c);
              int n2k = nk + NeNo.getZ(c);

              int i2 = ix - n2i;
              int j2 = iy - n2j;
              int k2 = iz - n2k;

              // Check if this node is one of the cell in which the particle is:
              if (i2>=0 && i2<2 && j2>=0 && j2<2 && k2>=0 && k2<2) {
                //Mass[i][j][j] 
                double qww = q[ip]*beta*weight[i][j][k]*weight[i2][j2][k2];
                double val[3][3];
                for (int ind1=0; ind1<3; ind1++)
                   for (int ind2=0; ind2<3; ind2++) val[ind1][ind2] = alpha[ind2][ind1]*qww;

                EMf->addMass(val, ni, nj, nk, c);
              }

            }
          }
        #ifdef __PROFILING__
        clocks->stop(9);
        #endif
    }
  }
}

/** Compute the mass matrix for the charge initialization */
void Particles3Dcomm::computeMassMatrices(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;

  EMf->setZeroMass();

 //#pragma omp parallel
  {
    //#pragma omp for
    for (register long long ip = 0; ip < nop; ip++)
    {
      const int ix = 2 + int (floor((x[ip] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[ip] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[ip] - zstart) * inv_dz));
      double weight[2][2][2];
      calculateWeights( weight, x[ip], y[ip], z[ip], ix, iy, iz, grid);

      //  Each node of the cell in which the particle is:
      for (int i=0; i<2; i++) 
        for (int j=0; j<2; j++) 
          for (int k=0; k<2; k++) {
            int ni = ix-i;
            int nj = iy-j;
            int nk = iz-k;

            // Only half of the possible nodes  (M is symetric)
            for (int c=0; c<14; c++) {
            //for (int c=0; c<27; c++) {
              int n2i = ni + NeNo.getX(c);
              int n2j = nj + NeNo.getY(c);
              int n2k = nk + NeNo.getZ(c);

              int i2 = ix - n2i;
              int j2 = iy - n2j;
              int k2 = iz - n2k;

              // Check if this node is one of the cell in which the particle is:
              if (i2>=0 && i2<2 && j2>=0 && j2<2 && k2>=0 && k2<2) {
                //Mass[i][j][j] 
                double ww = weight[i][j][k]*weight[i2][j2][k2];
                EMf->addMassXX(ww, ni, nj, nk, c);
              }

            }
          }
    }
  }
}
/** Compute only the charge density **/
void Particles3Dcomm::computeCharge(Field * EMf, Grid * grid, VirtualTopology3D * vct) {

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
 //#pragma omp parallel
  {
    //#pragma omp for
    for (register long long ip = 0; ip < nop; ip++)
    {
      const int ix = 2 + int (floor((x[ip] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[ip] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[ip] - zstart) * inv_dz));
      double weight[2][2][2];
      double temp[2][2][2];
      calculateWeights( weight, x[ip], y[ip], z[ip], ix, iy, iz, grid);

      // add charge density               
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
            temp[ii][jj][kk] = q[ip]*weight[ii][jj][kk];
          }
      EMf->addRho(temp, ix, iy, iz, ns);

      // add number of particles          
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
          }
      EMf->addN(weight, ix, iy, iz, ns);
    }
  }
}


/** Compute the currents (for writting the results) */
void Particles3Dcomm::computeCurrent(Field * EMf, Grid * grid, VirtualTopology3D * vct) {

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
 //#pragma omp parallel
  {
    //#pragma omp for
    for (register long long ip = 0; ip < nop; ip++)
    {
      const int ix = 2 + int (floor((x[ip] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[ip] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[ip] - zstart) * inv_dz));
      double weight[2][2][2];
      double temp[2][2][2];
      calculateWeights( weight, x[ip], y[ip], z[ip], ix, iy, iz, grid);

      // add current density - X
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*u[ip]*weight[ii][jj][kk];
      EMf->addJx(temp, ix, iy, iz, ns);

      // add current density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*v[ip]*weight[ii][jj][kk];
      EMf->addJy(temp, ix, iy, iz, ns);

      // add current density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*w[ip]*weight[ii][jj][kk];
      EMf->addJz(temp, ix, iy, iz, ns);

      // Pxx - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * u[ip] * u[ip] * weight[ii][jj][kk];
      EMf->addPxx(temp, ix, iy, iz, ns);
      // Pxy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * u[ip] * v[ip] * weight[ii][jj][kk];
      EMf->addPxy(temp, ix, iy, iz, ns);
      // Pxz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * u[ip] * w[ip] * weight[ii][jj][kk];
      EMf->addPxz(temp, ix, iy, iz, ns);
      // Pyy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * v[ip] * v[ip] * weight[ii][jj][kk];
      EMf->addPyy(temp, ix, iy, iz, ns);
      // Pyz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * v[ip] * w[ip] * weight[ii][jj][kk];
      EMf->addPyz(temp, ix, iy, iz, ns);
      // Pzz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip] * w[ip] * w[ip] * weight[ii][jj][kk];
      EMf->addPzz(temp, ix, iy, iz, ns);

    }
  }
}
/** Compute the currents (for writting the results) */
void Particles3Dcomm::computeCurrentBar(Field * EMf, Grid * grid, VirtualTopology3D * vct) {

  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
 //#pragma omp parallel
  {
    //#pragma omp for
    for (register long long ip = 0; ip < nop; ip++)
    {
      const int ix = 2 + int (floor((x[ip] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[ip] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[ip] - zstart) * inv_dz));
      double weight[2][2][2];
      double temp[2][2][2];
      calculateWeights( weight, x[ip], y[ip], z[ip], ix, iy, iz, grid);

      // add charge density               
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
            temp[ii][jj][kk] = q[ip]*weight[ii][jj][kk];
          }
      EMf->addRho(temp, ix, iy, iz, ns);

      // add current density - X
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*ubar[ip]*weight[ii][jj][kk];
      EMf->addJx(temp, ix, iy, iz, ns);

      // add current density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*vbar[ip]*weight[ii][jj][kk];
      EMf->addJy(temp, ix, iy, iz, ns);

      // add current density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = q[ip]*wbar[ip]*weight[ii][jj][kk];
      EMf->addJz(temp, ix, iy, iz, ns);

    }
  }
}


/** communicate buffers */
int Particles3Dcomm::communicate(VirtualTopology3D * ptVCT) {
  // allocate buffers
  MPI_Status status;
  int new_buffer_size;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4, avail5, avail6;
  for (int i = 0; i < buffer_size; i++) {
    b_X_RIGHT[i] = MIN_VAL;
    b_X_LEFT[i] = MIN_VAL;
    b_Y_RIGHT[i] = MIN_VAL;
    b_Y_LEFT[i] = MIN_VAL;
    b_Z_RIGHT[i] = MIN_VAL;
    b_Z_LEFT[i] = MIN_VAL;
  }
  npExitXright = 0, npExitXleft = 0, npExitYright = 0, npExitYleft = 0, npExitZright = 0, npExitZleft = 0, npExit = 0, rightDomain = 0;
  long long np_current = 0, nplast = nop - 1;

  while (np_current < nplast+1){

    // BC on particles
    if (x[np_current] < 0 && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
    else if (x[np_current] > Lx && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft); 
    if (y[np_current] < 0 && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
    else if (y[np_current] > Ly && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft); 
    if (z[np_current] < 0 && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);
    else if (z[np_current] > Lz && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);

    // if the particle exits, apply the boundary conditions add the particle to communication buffer
    if (x[np_current] < xstart || x[np_current] >xend){
      // communicate if they don't belong to the domain
      if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitXleft+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferXleft(b_X_LEFT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitXleft++;
      } 
      else if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitXleft++;
      } 
      else if (x[np_current] > xend && ptVCT->getXright_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitXright+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferXright(b_X_RIGHT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitXright++;
      }
      else if (x[np_current] > xend && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitXright++;
      }

    } else  if (y[np_current] < ystart || y[np_current] >yend){
      // communicate if they don't belong to the domain
      if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitYleft+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferYleft(b_Y_LEFT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitYleft++;
      }
      else if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL){
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitYleft++;
      }
      else if (y[np_current] > yend && ptVCT->getYright_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitYright+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferYright(b_Y_RIGHT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitYright++;
      }
      else if (y[np_current] > yend && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL){
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitYright++;
      }
    } else  if (z[np_current] < zstart || z[np_current] >zend){
      // communicate if they don't belong to the domain
      if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitZleft+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferZleft(b_Z_LEFT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);

        npExitZleft++;
      } 
      else if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitZleft++;
      }
      else if (z[np_current] > zend && ptVCT->getZright_neighbor_P() != MPI_PROC_NULL){
        // check if there is enough space in the buffer before putting in the particle
        while(((npExitZright+1)*nVar)>=buffer_size){
          resize_buffers((int) (buffer_size*2)); 
        }
        // put it in the communication buffer
        bufferZright(b_Z_RIGHT,np_current,ptVCT);
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);

        npExitZright++;
      }
      else if (z[np_current] > zend && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitZright++;
      }
    }  else {
      // particle is still in the domain, procede with the next particle
      np_current++;
    }

  }

  nop = nplast + 1;
  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  npExitingMax = reduceMaxNpExiting(npExitingMax);

  /*****************************************************/
  /* SEND AND RECEIVE MESSAGES */
  /*****************************************************/

  new_buffer_size = npExitingMax * nVar + 1;

  if (new_buffer_size > buffer_size) {
    cout << "resizing the receiving buffer" << endl;
    resize_buffers(new_buffer_size);
  }

  if (npExitingMax > 0) {
    communicateParticles(new_buffer_size, b_X_LEFT, b_X_RIGHT, b_Y_LEFT, b_Y_RIGHT, b_Z_LEFT, b_Z_RIGHT, ptVCT);

    // UNBUFFERING
    // message from XLEFT
    avail1 = unbuffer(b_X_RIGHT);
    avail2 = unbuffer(b_X_LEFT);
    avail3 = unbuffer(b_Y_RIGHT);
    avail4 = unbuffer(b_Y_LEFT);
    avail5 = unbuffer(b_Z_RIGHT);
    avail6 = unbuffer(b_Z_LEFT);
    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4 + avail5 + avail6;
    availALL = reduceNumberParticles(avail);
    if (availALL < 0)
      return (-1);              // too many particles coming, save data nad stop simulation
  }

  return (0);                   // everything was fine


}
/** resize the buffers */
void Particles3Dcomm::resize_buffers(int new_buffer_size) {
  cout << "RESIZING FROM " << buffer_size << " TO " << new_buffer_size << endl;
  // resize b_X_LEFT
  double *temp = new double[buffer_size];

  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_X_LEFT_ptr[i];
  // delete[] b_X_LEFT_ptr;
  delete[]b_X_LEFT;
  b_X_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_X_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_X_LEFT[i] = MIN_VAL;

  // resize b_X_RIGHT 
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_X_RIGHT_ptr[i];
  // delete[] b_X_RIGHT_ptr;
  delete[]b_X_RIGHT;
  b_X_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_X_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_X_RIGHT[i] = MIN_VAL;

  // resize b_Y_RIGHT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Y_RIGHT_ptr[i];
  // delete[] b_Y_RIGHT_ptr;
  delete[]b_Y_RIGHT;
  b_Y_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Y_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Y_RIGHT[i] = MIN_VAL;

  // resize b_Y_LEFT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Y_LEFT_ptr[i];
  // delete[] b_Y_LEFT_ptr;
  delete[]b_Y_LEFT;
  b_Y_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Y_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Y_LEFT[i] = MIN_VAL;

  // resize b_Z_RIGHT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Z_RIGHT_ptr[i];
  // delete[] b_Z_RIGHT_ptr;
  delete[]b_Z_RIGHT;
  b_Z_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Z_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Z_RIGHT[i] = MIN_VAL;

  // resize b_Z_LEFT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Z_LEFT_ptr[i];
  // delete[] b_Z_LEFT_ptr;
  delete[]b_Z_LEFT;
  b_Z_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Z_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Z_LEFT[i] = MIN_VAL;

  delete[]temp;

  b_X_RIGHT_ptr = b_X_RIGHT;
  b_Y_RIGHT_ptr = b_Y_RIGHT;
  b_Z_RIGHT_ptr = b_Z_RIGHT;
  b_Y_LEFT_ptr = b_Y_LEFT;
  b_X_LEFT_ptr = b_X_LEFT;
  b_Z_LEFT_ptr = b_Z_LEFT;

  buffer_size = new_buffer_size;
}
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  if (x[np_current] < 0)
    b_[npExitXleft * nVar] = x[np_current] + Lx;  // this applies to the the leftmost processor
  else
    b_[npExitXleft * nVar] = x[np_current];
  b_[npExitXleft * nVar + 1] = y[np_current];
  b_[npExitXleft * nVar + 2] = z[np_current];
  b_[npExitXleft * nVar + 3] = u[np_current];
  b_[npExitXleft * nVar + 4] = v[np_current];
  b_[npExitXleft * nVar + 5] = w[np_current];
  b_[npExitXleft * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitXleft * nVar + 7] = partID[np_current];
    b_[npExitXleft * nVar + 8] = partRank[np_current];
  }

}
/** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXright(double *b_, long long np_current, VirtualTopology3D * vct) {
  if (x[np_current] > Lx)
    b_[npExitXright * nVar] = x[np_current] - Lx; // this applies to the right most processor
  else
    b_[npExitXright * nVar] = x[np_current];
  b_[npExitXright * nVar + 1] = y[np_current];
  b_[npExitXright * nVar + 2] = z[np_current];
  b_[npExitXright * nVar + 3] = u[np_current];
  b_[npExitXright * nVar + 4] = v[np_current];
  b_[npExitXright * nVar + 5] = w[np_current];
  b_[npExitXright * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitXright * nVar + 7] = partID[np_current];
    b_[npExitXright * nVar + 8] = partRank[np_current];
  }
}
/** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitYleft * nVar] = x[np_current];
  if (y[np_current] < 0)
    b_[npExitYleft * nVar + 1] = y[np_current] + Ly;
  else
    b_[npExitYleft * nVar + 1] = y[np_current];
  b_[npExitYleft * nVar + 2] = z[np_current];
  b_[npExitYleft * nVar + 3] = u[np_current];
  b_[npExitYleft * nVar + 4] = v[np_current];
  b_[npExitYleft * nVar + 5] = w[np_current];
  b_[npExitYleft * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitYleft * nVar + 7] = partID[np_current];
    b_[npExitYleft * nVar + 8] = partRank[np_current];
  }
}
/** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYright(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitYright * nVar] = x[np_current];
  if (y[np_current] > Ly)
    b_[npExitYright * nVar + 1] = y[np_current] - Ly;
  else
    b_[npExitYright * nVar + 1] = y[np_current];
  b_[npExitYright * nVar + 2] = z[np_current];
  b_[npExitYright * nVar + 3] = u[np_current];
  b_[npExitYright * nVar + 4] = v[np_current];
  b_[npExitYright * nVar + 5] = w[np_current];
  b_[npExitYright * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitYright * nVar + 7] = partID[np_current];
    b_[npExitYright * nVar + 8] = partRank[np_current];
  }
}
/** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitZleft * nVar] = x[np_current];
  b_[npExitZleft * nVar + 1] = y[np_current];
  if (z[np_current] < 0)
    b_[npExitZleft * nVar + 2] = z[np_current] + Lz;
  else
    b_[npExitZleft * nVar + 2] = z[np_current];
  b_[npExitZleft * nVar + 3] = u[np_current];
  b_[npExitZleft * nVar + 4] = v[np_current];
  b_[npExitZleft * nVar + 5] = w[np_current];
  b_[npExitZleft * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitZleft * nVar + 7] = partID[np_current];
    b_[npExitZleft * nVar + 8] = partRank[np_current];
  }
}
/** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZright(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitZright * nVar] = x[np_current];
  b_[npExitZright * nVar + 1] = y[np_current];
  if (z[np_current] > Lz)
    b_[npExitZright * nVar + 2] = z[np_current] - Lz;
  else
    b_[npExitZright * nVar + 2] = z[np_current];
  b_[npExitZright * nVar + 3] = u[np_current];
  b_[npExitZright * nVar + 4] = v[np_current];
  b_[npExitZright * nVar + 5] = w[np_current];
  b_[npExitZright * nVar + 6] = q[np_current];
  if (TrackParticleID) {
    b_[npExitZright * nVar + 7] = partID[np_current];
    b_[npExitZright * nVar + 8] = partRank[np_current];
  }
}
/** This unbuffer the last communication */
int Particles3Dcomm::unbuffer(double *b_) {
  long long np_current = 0;
  // put the new particles at the end of the array, and update the number of particles
  while (b_[np_current * nVar] != MIN_VAL) {
    x[nop] = b_[nVar * np_current];
    y[nop] = b_[nVar * np_current + 1];
    z[nop] = b_[nVar * np_current + 2];
    u[nop] = b_[nVar * np_current + 3];
    v[nop] = b_[nVar * np_current + 4];
    w[nop] = b_[nVar * np_current + 5];
    q[nop] = b_[nVar * np_current + 6];
    if (TrackParticleID) {
      partID[nop]   = (long) b_[nVar * np_current + 7];
      partRank[nop] = (long) b_[nVar * np_current + 8];
    }
    np_current++;
    // these particles need further communication
    if (x[nop] < xstart || x[nop] > xend || y[nop] < ystart || y[nop] > yend || z[nop] < zstart || z[nop] > zend)
      rightDomain++;            // the particle is not in the domain
    nop++;
    if (nop > npmax) {
      cout << "Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
      return (-1);              // end the simulation because you dont have enough space on the array
    }
  }
  return (0);                   // everything was fine
}
/** Delete the a particle from the array and pack the the array, update the number of 
 * particles that are exiting
 * For deleting the particle from the array take the last particle and put it
 * in the position of the particle you want to delete
 * @param np = the index of the particle that must be deleted
 * @param nplast = the index of the last particle in the array
 */
void Particles3Dcomm::del_pack(long long np_current, long long *nplast) {
  x[np_current] = x[*nplast];
  y[np_current] = y[*nplast];
  z[np_current] = z[*nplast];
  u[np_current] = u[*nplast];
  v[np_current] = v[*nplast];
  w[np_current] = w[*nplast];
  q[np_current] = q[*nplast];
  if (TrackParticleID) {
    partID[np_current]   = partID[*nplast];
    partRank[np_current] = partRank[*nplast];
  }
  npExit++;
  (*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles3Dcomm::isMessagingDone(VirtualTopology3D * ptVCT) {
  int result = 0;
  result = reduceNumberParticles(rightDomain);
  if (result > 0 && cVERBOSE && ptVCT->getCartesian_rank() == 0)
    cout << "Further Comunication: " << result << " particles not in the right domain" << endl;
  return (result);

}
/** calculate the maximum number exiting from this domain */
int Particles3Dcomm::maxNpExiting() {
  int maxNp = 0;
  if (npExitXright > maxNp)
    maxNp = npExitXright;
  if (npExitXleft > maxNp)
    maxNp = npExitXleft;
  if (npExitYright > maxNp)
    maxNp = npExitYright;
  if (npExitYleft > maxNp)
    maxNp = npExitYleft;
  if (npExitZright > maxNp)
    maxNp = npExitZright;
  if (npExitZleft > maxNp)
    maxNp = npExitZleft;
  return (maxNp);
}
/** return X-coordinate of particle array */
double *Particles3Dcomm::getXall()  const {
  return (x);
}
/** return Y-coordinate  of particle array */
double *Particles3Dcomm::getYall()  const {
  return (y);
}
/** return Z-coordinate  of particle array*/
double *Particles3Dcomm::getZall()  const {
  return (z);
}
/** get X-velocity of particle with label indexPart */
double *Particles3Dcomm::getUall()  const {
  return (u);
}
/** get Y-velocity of particle with label indexPart */
double *Particles3Dcomm::getVall()  const {
  return (v);
}
/**get Z-velocity of particle with label indexPart */
double *Particles3Dcomm::getWall()  const {
  return (w);
}
/**get ID of particle with label indexPart */
long long *Particles3Dcomm::getParticleIDall()  const {
  return (partID);
}
/**get proc of particle with label indexPart */
long long *Particles3Dcomm::getParticleRankall()  const {
  return (partRank);
}
/**get charge of particle with label indexPart */
double *Particles3Dcomm::getQall()  const {
  return (q);
}
/** return X-coordinate of particle array as reference */
double *& Particles3Dcomm::getXref() {
  return (x);
}
/** return Y-coordinate  of particle array as reference */
double *& Particles3Dcomm::getYref() {
  return (y);
}
/** return Z-coordinate  of particle array as reference */
double *& Particles3Dcomm::getZref() {
  return (z);
}
/** get X-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getUref() {
  return (u);
}
/** get Y-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getVref() {
  return (v);
}
/**get Z-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getWref() {
  return (w);
}
/**get charge of particle with label indexPart as reference */
double *& Particles3Dcomm::getQref() {
  return (q);
}
/** return X-coordinate of particle with index indexPart */
double Particles3Dcomm::getX(long long indexPart)  const {
  return (x[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getY(long long indexPart)  const {
  return (y[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getZ(long long indexPart)  const {
  return (z[indexPart]);
}
/** get u (X-velocity) of particle with label indexPart */
double Particles3Dcomm::getU(long long indexPart)  const {
  return (u[indexPart]);
}
/** get v (Y-velocity) of particle with label indexPart */
double Particles3Dcomm::getV(long long indexPart)  const {
  return (v[indexPart]);
}
/**get w (Z-velocity) of particle with label indexPart */
double Particles3Dcomm::getW(long long indexPart)  const {
  return (w[indexPart]);
}
/**get ID of particle with label indexPart */
long long Particles3Dcomm::getParticleID(long long indexPart)  const {
  return (partID[indexPart]);
}
/**get ID of particle with label indexPart */
long long Particles3Dcomm::getParticleRank(long long indexPart)  const {
  return (partRank[indexPart]);
}
/**get charge of particle with label indexPart */
double Particles3Dcomm::getQ(long long indexPart)  const {
  return (q[indexPart]);
}
/** return the number of particles */
long long Particles3Dcomm::getNOP()  const {
  return (nop);
}
/** return the Kinetic energy */
double Particles3Dcomm::getKe() {
  double localKe = 0.0;
  double totalKe = 0.0;
  if (Lagrangian == 0) {
    for (register long long i = 0; i < nop; i++) {
      localKe += .5 * (q[i] / qom) * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    }
  } else {
    for (register long long i = 0; i < nop; i++) {
      double vx = u[i] + u_bulk;
      double vy = v[i] + v_bulk;
      double vz = w[i] + w_bulk;
      localKe += .5 * (q[i] / qom) * (vx*vx + vy*vy + vz*vz);
    }

  }
  MPI_Allreduce(&localKe, &totalKe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalKe);
}
/** return the total momentum */
double Particles3Dcomm::getP() {
  double localP = 0.0;
  double totalP = 0.0;
  if (Lagrangian == 0) {
    for (register long long i = 0; i < nop; i++)
      localP += (q[i] / qom) * sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  } else {
    for (register long long i = 0; i < nop; i++) {
      double vx = u[i] + u_bulk;
      double vy = v[i] + v_bulk;
      double vz = w[i] + w_bulk;
      localP += (q[i] / qom) * sqrt(vx*vx + vy*vy + vz*vz);
    }

  }
  MPI_Allreduce(&localP, &totalP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalP);
}

/** return the Kinetic energy removed from the system*/
double Particles3Dcomm::getKremoved() {
  double localK = K_removed;
  double totalK = 0.0;
  MPI_Allreduce(&localK, &totalK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalK);
}

/** return the highest kinetic energy */
double Particles3Dcomm::getMaxVelocity() {
  double localVel = 0.0;
  double maxVel = 0.0;
  for (long long i = 0; i < nop; i++)
    localVel = max(localVel, sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
  MPI_Allreduce(&localVel, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return (maxVel);
}


/** get energy spectrum */
unsigned long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel) {
  unsigned long *f = new unsigned long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (long long i = 0; i < nop; i++) {
    Vel = sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    bin = int (floor(Vel / dv));
    if (bin >= nBins)
      f[nBins - 1] += 1;
    else
      f[bin] += 1;
  }

  MPI_Allreduce(MPI_IN_PLACE, f, nBins, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  return f;
}

/** return the total charge */
double Particles3Dcomm::getTotalCharge() {
  double localQ = 0.0;
  double totalQ = 0.0;
  for (register long long i = 0; i < nop; i++)
    localQ += q[i];
  MPI_Allreduce(&localQ, &totalQ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalQ);
}

/** return the total number of particles */
long long Particles3Dcomm::getTotalNumPart() {
  long long localN = 0.0;
  long long totalN = 0.0;
  localN = nop;
  MPI_Allreduce(&localN, &totalN, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  return (totalN);
}


/** print particles info */
void Particles3Dcomm::Print(VirtualTopology3D * ptVCT) const {
  cout << endl;
  cout << "Number of Particles: " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << ns << endl;
  for (long long i = 0; i < nop; i++)
    cout << "Particles #" << i << " x=" << x[i] << " y=" << y[i] << " z=" << z[i] << " u=" << u[i] << " v=" << v[i] << " w=" << w[i] << endl;
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp(VirtualTopology3D * ptVCT)  const {
  cout << endl;
  cout << "Number of Particles of species " << ns << ": " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << endl;
}
