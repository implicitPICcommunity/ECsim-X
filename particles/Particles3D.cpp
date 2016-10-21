/*******************************************************************************************
  Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/


#include <iostream>
#include <math.h>

#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Collective.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "MPIdata.h"

#include "Particles3D.h"


#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
// particles processed together
#define P_SAME_TIME 2

/**
 * 
 * Class for particles of the same species
 * @date Fri Jun 4 2009
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */


/**
 *  Auxiliar routine for random numbers (thread-safe)
 */
inline double Drand48(long long& s)
{
  s = (((long long) 0x5deece66d)*s + ((long long) 0xb))
        &((long long) 0xffffffffffff);
  return ((double) s)/((double) (0x1000000000000));
}


/** constructor */
Particles3D::Particles3D() {
  // see allocate(int species, Collective* col, VirtualTopology3D* vct, Grid* grid)
  //

}
/** deallocate particles */
Particles3D::~Particles3D() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
}

/** particles are uniformly distributed with zero velocity   */
void Particles3D::uniform_background(Grid * grid, Field * EMf) {
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              u[counter] = 0.0;
              v[counter] = 0.0;
              w[counter] = 0.0;
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
              counter++;
            }


  cout << "Velocity Maxwellian Distribution " << endl;
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
  <ul>
  <li> dim = 0 --> constant velocity on X direction </li>
  <li> dim = 1 --> constant velocity on Y direction </li>
  <li> dim = 2 --> constant velocity on Z direction </li>
  </ul>

*/
void Particles3D::constantVelocity(double vel, int dim, Grid * grid, Field * EMf) {
  switch (dim) {
    case 0:
      for (long long i = 0; i < nop; i++)
        u[i] = vel, v[i] = 0.0, w[i] = 0.0;
      break;
    case 1:
      for (register long long i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = vel, w[i] = 0.0;
      break;
    case 2:
      for (register long long i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = 0.0, w[i] = vel;
      break;

  }

}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + Drand48(seed)) * (dx / npcelx) + grid->getXN(i, j, k);  
              y[counter] = (jj + Drand48(seed)) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + Drand48(seed)) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());

              // u
              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
              theta = 2.0 * M_PI * Drand48(seed);
              u[counter] = u0 + uth * prob * cos(theta);
              // v
              v[counter] = v0 + vth * prob * sin(theta);
              // w
              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
              theta = 2.0 * M_PI * Drand48(seed);
              w[counter] = w0 + wth * prob * cos(theta);

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }
}


/** Maxellian random velocity and uniform spatial distribution, goes with Double GEM */
void Particles3D::MaxwellianDoubleGEM(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  const double coarsedy= grid->getDY();
  double globaly;
  double shaperz;
  
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
	      
	      globaly= grid->getYN(i,j,k)+ coarsedy;
	      shaperz= -tanh((globaly - Ly/2)/delta) ;
	      
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);  
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
              // u
              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
              theta = 2.0 * M_PI * Drand48(seed);
              u[counter] = u0 + uth * prob * cos(theta);
              // v
              v[counter] = v0 + vth * prob * sin(theta);
              // w
              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
              theta = 2.0 * M_PI * Drand48(seed);
              w[counter] = w0*shaperz + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}

void Particles3D::get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]){

  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  const double ixd = floor((xp - xstart) * inv_dx);
  const double iyd = floor((yp - ystart) * inv_dy);
  const double izd = floor((zp - zstart) * inv_dz);

  ix = 2 + int (ixd);
  iy = 2 + int (iyd);
  iz = 2 + int (izd);

  if (ix < 1)
    ix = 1;
  if (iy < 1)
    iy = 1;
  if (iz < 1)
    iz = 1;
  if (ix > nxn - 1)
    ix = nxn - 1;
  if (iy > nyn - 1)
    iy = nyn - 1;
  if (iz > nzn - 1)
    iz = nzn - 1;

  double xi  [2];
  double eta [2];
  double zeta[2];

  xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
  eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
  zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
  xi  [1] = grid->getXN(ix,iy,iz) - xp;
  eta [1] = grid->getYN(ix,iy,iz) - yp;
  zeta[1] = grid->getZN(ix,iy,iz) - zp;

  for (int ii = 0; ii < 2; ii++)
    for (int jj = 0; jj < 2; jj++)
      for (int kk = 0; kk < 2; kk++)
        weights[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
}

void Particles3D::get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez){

  Exl = 0.0;
  Eyl = 0.0;
  Ezl = 0.0;

  int l = 0;
  for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++) {
        Exl += weights[i][j][k] * Ex[ix-i][iy-j][iz-k];
        Eyl += weights[i][j][k] * Ey[ix-i][iy-j][iz-k];
        Ezl += weights[i][j][k] * Ez[ix-i][iy-j][iz-k];
        l = l + 1;
      }

}

void Particles3D::get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext){

  Bxl = 0.0;
  Byl = 0.0;
  Bzl = 0.0;

  int l = 0;
  for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++) {
        Bxl += weights[i][j][k] * (Bx[ix-i][iy-j][iz-k] + Fext*Bx_ext[ix-i][iy-j][iz-k]);
        Byl += weights[i][j][k] * (By[ix-i][iy-j][iz-k] + Fext*By_ext[ix-i][iy-j][iz-k]);
        Bzl += weights[i][j][k] * (Bz[ix-i][iy-j][iz-k] + Fext*Bz_ext[ix-i][iy-j][iz-k]);
        l = l + 1;
      }
}

// Particles mover for the ECIMM algorithm
int Particles3D::mover_ECIMM(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** ECIMM MOVER species " << ns << " ***"<< endl;
  }
  double start_mover_PC = MPI_Wtime();
  double ***Exth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEFxth());
  double ***Eyth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEFyth());
  double ***Ezth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEFzth());
  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());
  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());
  double Fext = EMf->getFext();

  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  // don't bother trying to push any particles simultaneously;
  // MIC already does vectorization automatically, and trying
  // to do it by hand only hurts performance.
#pragma omp parallel for
#pragma simd                    // this just slows things down (why?)
  for (long long rest = 0; rest < nop; rest++) {
    // Copy the velocities
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];

    // copy the positions
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double uptilde;
    double vptilde;
    double wptilde;
    const double ixd = floor((xp - xstart) * inv_dx);
    const double iyd = floor((yp - ystart) * inv_dy);
    const double izd = floor((zp - zstart) * inv_dz);
    int ix = 2 + int (ixd);
    int iy = 2 + int (iyd);
    int iz = 2 + int (izd);
    double xi[2], eta[2], zeta[2];
    xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
    eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
    zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
    xi  [1] = grid->getXN(ix,iy,iz) - xp;
    eta [1] = grid->getYN(ix,iy,iz) - yp;
    zeta[1] = grid->getZN(ix,iy,iz) - zp;

    double Exl = 0.0;
    double Eyl = 0.0;
    double Ezl = 0.0;
    double Bxl = 0.0;
    double Byl = 0.0;
    double Bzl = 0.0;

    // MIC refuses to vectorize this ...
    // 
    // double weight[2][2][2];
    // for (int ii = 0; ii < 2; ii++)
    // for (int jj = 0; jj < 2; jj++)
    // for (int kk = 0; kk < 2; kk++)
    // weight[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
    // for (int ii = 0; ii < 2; ii++)
    // for (int jj = 0; jj < 2; jj++)
    // for (int kk = 0; kk < 2; kk++) {
    // const double Exlp = weight[ii][jj][kk] * Ex.get(ix - ii, iy - jj, iz - kk);
    // const double Eylp = weight[ii][jj][kk] * Ey.get(ix - ii, iy - jj, iz - kk);
    // const double Ezlp = weight[ii][jj][kk] * Ez.get(ix - ii, iy - jj, iz - kk);
    // const double Bxlp = weight[ii][jj][kk] * Bx.get(ix - ii, iy - jj, iz - kk);
    // const double Bylp = weight[ii][jj][kk] * By.get(ix - ii, iy - jj, iz - kk);
    // const double Bzlp = weight[ii][jj][kk] * Bz.get(ix - ii, iy - jj, iz - kk);
    // Exl += Exlp;
    // Eyl += Eylp;
    // Ezl += Ezlp;
    // Bxl += Bxlp;
    // Byl += Bylp;
    // Bzl += Bzlp;
    // }

    // ... so we expand things out instead
    // 
    const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
    const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
    const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
    const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
    const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
    const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
    const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
    const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
    // 
    Bxl += weight000 * (Bx[ix][iy][iz]             + Fext*Bx_ext[ix][iy][iz]);
    Bxl += weight001 * (Bx[ix][iy][iz - 1]         + Fext*Bx_ext[ix][iy][iz-1]);
    Bxl += weight010 * (Bx[ix][iy - 1][iz]         + Fext*Bx_ext[ix][iy-1][iz]);
    Bxl += weight011 * (Bx[ix][iy - 1][iz - 1]     + Fext*Bx_ext[ix][iy-1][iz-1]);
    Bxl += weight100 * (Bx[ix - 1][iy][iz]         + Fext*Bx_ext[ix-1][iy][iz]);
    Bxl += weight101 * (Bx[ix - 1][iy][iz - 1]     + Fext*Bx_ext[ix-1][iy][iz-1]);
    Bxl += weight110 * (Bx[ix - 1][iy - 1][iz]     + Fext*Bx_ext[ix-1][iy-1][iz]);
    Bxl += weight111 * (Bx[ix - 1][iy - 1][iz - 1] + Fext*Bx_ext[ix-1][iy-1][iz-1]);
    // 
    Byl += weight000 * (By[ix][iy][iz]             + Fext*By_ext[ix][iy][iz]);
    Byl += weight001 * (By[ix][iy][iz - 1]         + Fext*By_ext[ix][iy][iz-1]);
    Byl += weight010 * (By[ix][iy - 1][iz]         + Fext*By_ext[ix][iy-1][iz]);
    Byl += weight011 * (By[ix][iy - 1][iz - 1]     + Fext*By_ext[ix][iy-1][iz-1]);
    Byl += weight100 * (By[ix - 1][iy][iz]         + Fext*By_ext[ix-1][iy][iz]);
    Byl += weight101 * (By[ix - 1][iy][iz - 1]     + Fext*By_ext[ix-1][iy][iz-1]);
    Byl += weight110 * (By[ix - 1][iy - 1][iz]     + Fext*By_ext[ix-1][iy-1][iz]);
    Byl += weight111 * (By[ix - 1][iy - 1][iz - 1] + Fext*By_ext[ix-1][iy-1][iz-1]);
    // 
    Bzl += weight000 * (Bz[ix][iy][iz]             + Fext*Bz_ext[ix][iy][iz]);
    Bzl += weight001 * (Bz[ix][iy][iz - 1]         + Fext*Bz_ext[ix][iy][iz-1]);
    Bzl += weight010 * (Bz[ix][iy - 1][iz]         + Fext*Bz_ext[ix][iy-1][iz]);
    Bzl += weight011 * (Bz[ix][iy - 1][iz - 1]     + Fext*Bz_ext[ix][iy-1][iz-1]);
    Bzl += weight100 * (Bz[ix - 1][iy][iz]         + Fext*Bz_ext[ix-1][iy][iz]);
    Bzl += weight101 * (Bz[ix - 1][iy][iz - 1]     + Fext*Bz_ext[ix-1][iy][iz-1]);
    Bzl += weight110 * (Bz[ix - 1][iy - 1][iz]     + Fext*Bz_ext[ix-1][iy-1][iz]);
    Bzl += weight111 * (Bz[ix - 1][iy - 1][iz - 1] + Fext*Bz_ext[ix-1][iy-1][iz-1]);
    // 
    Exl += weight000 * Exth[ix][iy][iz];
    Exl += weight001 * Exth[ix][iy][iz - 1];
    Exl += weight010 * Exth[ix][iy - 1][iz];
    Exl += weight011 * Exth[ix][iy - 1][iz - 1];
    Exl += weight100 * Exth[ix - 1][iy][iz];
    Exl += weight101 * Exth[ix - 1][iy][iz - 1];
    Exl += weight110 * Exth[ix - 1][iy - 1][iz];
    Exl += weight111 * Exth[ix - 1][iy - 1][iz - 1];
    // 
    Eyl += weight000 * Eyth[ix][iy][iz];
    Eyl += weight001 * Eyth[ix][iy][iz - 1];
    Eyl += weight010 * Eyth[ix][iy - 1][iz];
    Eyl += weight011 * Eyth[ix][iy - 1][iz - 1];
    Eyl += weight100 * Eyth[ix - 1][iy][iz];
    Eyl += weight101 * Eyth[ix - 1][iy][iz - 1];
    Eyl += weight110 * Eyth[ix - 1][iy - 1][iz];
    Eyl += weight111 * Eyth[ix - 1][iy - 1][iz - 1];
    // 
    Ezl += weight000 * Ezth[ix][iy][iz];
    Ezl += weight001 * Ezth[ix][iy][iz - 1];
    Ezl += weight010 * Ezth[ix][iy - 1][iz];
    Ezl += weight011 * Ezth[ix][iy - 1][iz - 1];
    Ezl += weight100 * Ezth[ix - 1][iy][iz];
    Ezl += weight101 * Ezth[ix - 1][iy][iz - 1];
    Ezl += weight110 * Ezth[ix - 1][iy - 1][iz];
    Ezl += weight111 * Ezth[ix - 1][iy - 1][iz - 1];
//fprintf(stderr,"\n================ Local Fields          \n");
//fprintf(stderr,"El: %13.6e %13.6e %13.6e \nBl:%13.6e %13.6e %13.6e \n"
              //,Exl,Eyl,Ezl, Bxl, Byl, Bzl);
    // 
    // End of interpolation
    double beta = 0.5*qom*dt;
    
    double vbEx = up + beta*Exl;
    double vbEy = vp + beta*Eyl;
    double vbEz = wp + beta*Ezl;

    // V tilde = alpha*(v + beta*E)
    double omcx = beta * Bxl;
    double omcy = beta * Byl;
    double omcz = beta * Bzl;
    double denom = 1.0 / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
    double edotb = vbEx * omcx + vbEy * omcy + vbEz * omcz;
    uptilde = (vbEx + (vbEy * omcz - vbEz * omcy + edotb * omcx)) * denom;
    vptilde = (vbEy + (vbEz * omcx - vbEx * omcz + edotb * omcy)) * denom;
    wptilde = (vbEz + (vbEx * omcy - vbEy * omcx + edotb * omcz)) * denom;


    // update the velocity
    up = 2.0 * uptilde - u[rest];
    vp = 2.0 * vptilde - v[rest];
    wp = 2.0 * wptilde - w[rest];
    u[rest] = up;
    v[rest] = vp;
    w[rest] = wp;

    // Update the positions with the old velocity
    x[rest] += up * dt;
    y[rest] += vp * dt;
    z[rest] += wp * dt;



  }                             // END OF ALL THE PARTICLES

  // ********************//
  // COMMUNICATION 
  // *******************//
  // timeTasks.start_communicate();
  const int avail = communicate(vct);
  if (avail < 0)
    return (-1);
  MPI_Barrier(MPI_COMM_WORLD);
  // communicate again if particles are not in the correct domain
  while (isMessagingDone(vct) > 0) {
    // COMMUNICATION
    const int avail = communicate(vct);
    if (avail < 0)
      return (-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // timeTasks.addto_communicate();
 
  return (0);                   // exit succcesfully (hopefully) 
}

int Particles3D::particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is){

  /* -- NOTE: Hardcoded option -- */
  enum {LINEAR,INITIAL,FFIELD};
  int rtype = INITIAL;
  /* -- END NOTE -- */

  if (vct->getCartesian_rank()==0){
    cout << "*** Repopulator species " << ns << " ***" << endl;
  }
  double  FourPI =16*atan(1.0);
  int avail;
  long long store_nop=nop;

  ////////////////////////
  // INJECTION FROM XLEFT
  ////////////////////////
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft == 2){ // use Field topology in this case
    long long particles_index=0;
    long long nplast = nop-1;

    while (particles_index < nplast+1) {
      if (x[particles_index] < 3.0*dx ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }

    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< 4;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);

                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


                particles_index++ ;
              }
    nop = particles_index;
  }

  store_nop = nop;

  ////////////////////////
  // INJECTION FROM YLEFT
  ////////////////////////
  if (vct->getYleft_neighbor() == MPI_PROC_NULL  && bcPfaceYleft == 2)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] < 3.0*dy ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< 4;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZLEFT
  ////////////////////////
  if (vct->getZleft_neighbor() == MPI_PROC_NULL  && bcPfaceZleft == 2)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] < 3.0*dz ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< 4;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM XRIGHT
  ////////////////////////
  if (vct->getXright_neighbor() == MPI_PROC_NULL  && bcPfaceXright == 2){
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (x[particles_index] > (Lx-3.0*dx) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=(grid->getNXC()-4); i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM YRIGHT
  ////////////////////////
  if (vct->getYright_neighbor() == MPI_PROC_NULL  && bcPfaceYright == 2)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] > (Ly-3.0*dy) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=(grid->getNYC()-4); j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZRIGHT
  ////////////////////////
  if (vct->getZright_neighbor() == MPI_PROC_NULL  && bcPfaceZright == 2)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] > (Lz-3.0*dz) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=(grid->getNZC()-4); k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                prob  = sqrt(-2.0*log(1.0-.999999*Drand48(seed)));
                theta = 2.0*M_PI*Drand48(seed);
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  if (vct->getCartesian_rank()==0){
    cout << "*** number of particles " << nop << " ***" << endl;
  }

  //********************//
  // COMMUNICATION
  // *******************//
  avail = communicate(vct);
  if (avail < 0) return(-1);

  MPI_Barrier(MPI_COMM_WORLD);

  // communicate again if particles are not in the correct domain
  while(isMessagingDone(vct) >0){
    // COMMUNICATION
    avail = communicate(vct);
    if (avail < 0)
      return(-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return(0); // exit succcesfully (hopefully)
}
