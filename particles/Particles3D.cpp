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
#include <time.h>       /* time */

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

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  long long seed = time(NULL) + (rank+1)*20 + ns;
  srand(seed);
  srand48(seed);

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
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOLc(i,j,k));
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

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + 0.5) * (dx / npcelx) + grid->getXN(i, j, k);  
              y[counter] = (jj + 0.5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + 0.5) * (dz / npcelz) + grid->getZN(i, j, k);
              //x[counter] = dx*(rand() / (double) RAND_MAX) + grid->getXN(i, j, k);  
              //y[counter] = dy*(rand() / (double) RAND_MAX) + grid->getYN(i, j, k);
              //z[counter] = dz*(rand() / (double) RAND_MAX) + grid->getZN(i, j, k);
              //x[counter] = Drand48(seed)*dx + grid->getXN(i,j,k);
              //y[counter] = Drand48(seed)*dy + grid->getYN(i,j,k);
              //z[counter] = Drand48(seed)*dz + grid->getZN(i,j,k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOLc(i,j,k));

              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = u0 + uth * prob * cos(theta);
              v[counter] = v0 + vth * prob * sin(theta);
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
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
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOLc(i,j,k));
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

/** Monoenergetic random velocity and uniform spatial distribution */
int Particles3D::monoenergetic_box(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

  /* initialize random generator with different seed on different processor */
  double harvest;
  double prob, theta, sign;
  long long counter=0;
  for (int i=1; i< grid->getNXC()-1;i++)
    for (int j=1; j< grid->getNYC()-1;j++)
      for (int k=1; k< grid->getNZC()-1;k++)
        for (int ii=0; ii < npcelx; ii++)
          for (int jj=0; jj < npcely; jj++)
	    for (int kk=0; kk < npcelz; kk++){
              x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
              y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
              z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
              if (fabs(x[counter] - x_center) < L_square/2 && fabs(y[counter] - y_center) < L_square/2 && fabs(z[counter] - z_center) < L_square/2){
                // q = charg
                q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
                if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
                	cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
                // u
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[counter] = multiple*u0 + uth*cos(theta);
                // v
                v[counter] = multiple*v0 + vth*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[counter] = multiple*w0 + wth*cos(theta);
                if (TrackParticleID)
                	ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
                
                counter++ ;
              }
            }
  nop = counter;
  return nop;
}

/** Maxellian random velocity and uniform spatial distribution */
int Particles3D::maxwell_box(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 

	double harvest;
	double prob, theta, sign;
	long long counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							if (fabs(x[counter] - x_center) < L_square/2 && fabs(y[counter] - y_center) < L_square/2 && fabs(z[counter] - z_center) < L_square/2){
							// q = charg
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
							if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							// u
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							u[counter] = multiple*u0 + uth*prob*cos(theta);
							// v
							v[counter] = multiple*v0 + vth*prob*sin(theta);
							// w
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							w[counter] = multiple*w0 + wth*prob*cos(theta);
							if (TrackParticleID)
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							counter++ ;
							}
						}

	nop = counter;
	return nop;
}

/** Maxellian random velocity and uniform spatial distribution */
int Particles3D::maxwell_sphere(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 
	double harvest;
	double prob, theta, sign;
	long long counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
              double dstx = x[counter] - x_center;
              double dsty = y[counter] - y_center;
              double dstz = z[counter] - z_center;
              double dst = sqrt(dstx*dstx + dsty*dsty + dstz*dstz);
              if (dst < L_square/2) {
							  // q = charg
							  q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
							  if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								  cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							  // u
							  harvest =   rand()/(double)RAND_MAX;
							  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  u[counter] = multiple*u0 + uth*prob*cos(theta);
							  // v
							  v[counter] = multiple*v0 + vth*prob*sin(theta);
							  // w
							  harvest =   rand()/(double)RAND_MAX;
							  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  w[counter] = multiple*w0 + wth*prob*cos(theta);
							  if (TrackParticleID)
								  ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							  counter++ ;
							}
						}

	nop = counter;
	return nop;
}

/** Maxellian random velocity and uniform spatial distribution */
int Particles3D::inject_maxwell_sphere(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */

  long long counter=nop;

	double harvest;
	double prob, theta, sign;
	double inject_probability;
	inject_probability = Vinj;
	inject_probability /= npcelx*npcely*npcelz ;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							double dstx = x[counter] - x_center;
							double dsty = y[counter] - y_center;
							double dstz = z[counter] - z_center;
							double dst = sqrt(dstx*dstx + dsty*dsty + dstz*dstz);
							if (dst < L_square/2) {
								harvest =   rand()/(double)RAND_MAX;
								//cout << harvest << "   " << inject_probability << endl;
								if(harvest < inject_probability){
							  // q = charg
							  q[counter] =  (qom/fabs(qom))* rhoINJECT *(1.0/grid->getInvVOLc(i,j,k));
							  if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								  cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							  // u
							  harvest =   rand()/(double)RAND_MAX;
							  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  u[counter] = multiple*u0 + uth*prob*cos(theta);
							  // v
							  v[counter] = multiple*v0 + vth*prob*sin(theta);
							  // w
							  harvest =   rand()/(double)RAND_MAX;
							  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  w[counter] = multiple*w0 + wth*prob*cos(theta);
							  if (TrackParticleID)
								  ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							  counter++ ;
							}
							}
						}

	nop = counter;
	return nop;
}
/** Maxellian random velocity and uniform spatial distribution */
int Particles3D::monoenergetic_sphere(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 
	double harvest;
	double prob, theta, sign;
	long long counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
              double dstx = x[counter] - x_center;
              double dsty = y[counter] - y_center;
              double dstz = z[counter] - z_center;
              double dst = sqrt(dstx*dstx + dsty*dsty + dstz*dstz);
              if (dst < L_square/2) {
							  // q = charg
							  q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
							  if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								  cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							  // u
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  u[counter] = multiple*u0 + uth*cos(theta);
							  // v
							  v[counter] = multiple*v0 + vth*sin(theta);
							  // w
							  harvest =   rand()/(double)RAND_MAX;
							  theta = 2.0*M_PI*harvest;
							  w[counter] = multiple*w0 + wth*cos(theta);
							  if (TrackParticleID)
								  ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							  counter++ ;
							}
						}

	nop = counter;
	return nop;
}





double Particles3D::ReturnToCenterCircle(double L_outer){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	double r;
	double udotr;
	double ukick = 3.0 *uth;

  double xC = (Cylindrical == 1 ? 0 : Lx/2);
  double yC = Ly/2;
  double zC = Lz/2;

	while (np_current < nplast+1) {
    double xp = x[np_current]-xC;
    double yp = y[np_current]-yC;
    double zp = z[np_current]-zC;
    r = 1e-10 + sqrt(xp*xp + yp*yp + zp*zp);
		if (r> L_outer) {
			u[np_current] = -ukick * (x[np_current]-xC)/r;
			v[np_current] = -ukick * (y[np_current]-yC)/r;
			w[np_current] = -ukick * (z[np_current]-zC)/r;
			x[np_current] += u[np_current]*dt ;
			y[np_current] += v[np_current]*dt ;
			z[np_current] += w[np_current]*dt ;
      

			np_current++;
		} else {
			np_current++;
		}
	}
	nop = nplast +1;
	return(0.0);
}

void Particles3D::MaxwellianFromFieldFile(Grid* grid,Field* EMf,VirtualTopology3D* vct){

  if (vct->getCartesian_rank() == 0) {
    cout << "PARTICLES INITIAZLIZED FROM FIELDS WITH MAXWELLIAN DISTRIBUTION" << endl;
  }
  /*cout << "Maxwellian From Field File" << endl;*/
  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++) {

        double M_rho = EMf->getRHOcs(i, j, k, ns);
        double M_u0  = EMf->getJxcs (i,j,k,ns)/M_rho;
        double M_v0  = EMf->getJycs (i,j,k,ns)/M_rho;
        double M_w0  = EMf->getJzcs (i,j,k,ns)/M_rho;
        // The module of rho to be used for computing the charge
        double M_rho_abs = fabs(EMf->getRHOcs(i, j, k, ns));

        for (int ii = 0; ii < npcelx; ii++){
          for (int jj = 0; jj < npcely; jj++){
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);

              // q = charge
              q[counter] = (qom / fabs(qom)) * (M_rho_abs / npcel) * (1.0 / grid->getInvVOLc(i,j,k));
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = M_u0 + uth * prob * cos(theta);
              // v
              v[counter] = M_v0 + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = M_w0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
              counter++;
            }
          }
        }
      }
}

void Particles3D::empty(Grid* grid,Field* EMf,VirtualTopology3D* vct){
        nop = 0;
printf("EMPTY\n");
}

/* Fix the charge of each particles by using the mass matrices */
void Particles3D::fixCharge(Grid* grid,Field* EMf,VirtualTopology3D* vct, Collective *col){

  computeMassMatrices(EMf, grid, vct);
  EMf->communicateGhostP2G(vct);
  EMf->computeRhotilde(vct, col, grid, ns);

  double ***rhotilde = EMf->getRhotilde(ns); 
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  for (int ip=0; ip<nop; ip++) {
    const double ixd = floor((x[ip] - xstart) * inv_dx);
    const double iyd = floor((y[ip] - ystart) * inv_dy);
    const double izd = floor((z[ip] - zstart) * inv_dz);
    int ix = 2 + int (ixd);
    int iy = 2 + int (iyd);
    int iz = 2 + int (izd);
    double xi[2], eta[2], zeta[2];
    xi  [0] = x[ip] - grid->getXN(ix-1,iy  ,iz  );
    eta [0] = y[ip] - grid->getYN(ix  ,iy-1,iz  );
    zeta[0] = z[ip] - grid->getZN(ix  ,iy  ,iz-1);
    xi  [1] = grid->getXN(ix,iy,iz) - x[ip];
    eta [1] = grid->getYN(ix,iy,iz) - y[ip];
    zeta[1] = grid->getZN(ix,iy,iz) - z[ip];

    double qp = 0.0;
    const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
    const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
    const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
    const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
    const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
    const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
    const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
    const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
    // 
    qp += weight000 *(rhotilde[ix][iy][iz]            )/grid->getInvVOLn(ix, iy, iz);
    qp += weight001 *(rhotilde[ix][iy][iz - 1]        )/grid->getInvVOLn(ix, iy, iz-1);
    qp += weight010 *(rhotilde[ix][iy - 1][iz]        )/grid->getInvVOLn(ix, iy-1, iz);
    qp += weight011 *(rhotilde[ix][iy - 1][iz - 1]    )/grid->getInvVOLn(ix, iy-1, iz-1);
    qp += weight100 *(rhotilde[ix - 1][iy][iz]        )/grid->getInvVOLn(ix-1, iy, iz);
    qp += weight101 *(rhotilde[ix - 1][iy][iz - 1]    )/grid->getInvVOLn(ix-1, iy, iz-1);
    qp += weight110 *(rhotilde[ix - 1][iy - 1][iz]    )/grid->getInvVOLn(ix-1, iy-1, iz);
    qp += weight111 *(rhotilde[ix - 1][iy - 1][iz - 1])/grid->getInvVOLn(ix-1, iy-1, iz-1);
	
    q[ip] = qp;
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
int Particles3D::updateVelocity(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** ECSIM MOVER species " << "***" << endl;
    //cout << "*** ECSIM MOVER species " << ns << "( " << nop << " part ) ***"<< endl;
  }
  double start_mover_PC = MPI_Wtime();
  double ***Exth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getExth());
  double ***Eyth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEyth());
  double ***Ezth= asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEzth());
  double ***Ex_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx_ext());
  double ***Ey_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy_ext());
  double ***Ez_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz_ext());
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
    Exl += weight000 * (Exth[ix][iy][iz]             + Fext*Ex_ext[ix][iy][iz]);
    Exl += weight001 * (Exth[ix][iy][iz - 1]         + Fext*Ex_ext[ix][iy][iz - 1]);
    Exl += weight010 * (Exth[ix][iy - 1][iz]         + Fext*Ex_ext[ix][iy - 1][iz]);
    Exl += weight011 * (Exth[ix][iy - 1][iz - 1]     + Fext*Ex_ext[ix][iy - 1][iz - 1]);
    Exl += weight100 * (Exth[ix - 1][iy][iz]         + Fext*Ex_ext[ix - 1][iy][iz]);
    Exl += weight101 * (Exth[ix - 1][iy][iz - 1]     + Fext*Ex_ext[ix - 1][iy][iz - 1]);
    Exl += weight110 * (Exth[ix - 1][iy - 1][iz]     + Fext*Ex_ext[ix - 1][iy - 1][iz]);
    Exl += weight111 * (Exth[ix - 1][iy - 1][iz - 1] + Fext*Ex_ext[ix - 1][iy - 1][iz - 1]);
    //                                                                                    
    Eyl += weight000 * (Eyth[ix][iy][iz]             + Fext*Ey_ext[ix][iy][iz]);
    Eyl += weight001 * (Eyth[ix][iy][iz - 1]         + Fext*Ey_ext[ix][iy][iz - 1]);
    Eyl += weight010 * (Eyth[ix][iy - 1][iz]         + Fext*Ey_ext[ix][iy - 1][iz]);
    Eyl += weight011 * (Eyth[ix][iy - 1][iz - 1]     + Fext*Ey_ext[ix][iy - 1][iz - 1]);
    Eyl += weight100 * (Eyth[ix - 1][iy][iz]         + Fext*Ey_ext[ix - 1][iy][iz]);
    Eyl += weight101 * (Eyth[ix - 1][iy][iz - 1]     + Fext*Ey_ext[ix - 1][iy][iz - 1]);
    Eyl += weight110 * (Eyth[ix - 1][iy - 1][iz]     + Fext*Ey_ext[ix - 1][iy - 1][iz]);
    Eyl += weight111 * (Eyth[ix - 1][iy - 1][iz - 1] + Fext*Ey_ext[ix - 1][iy - 1][iz - 1]);
    //                                                                                    
    Ezl += weight000 * (Ezth[ix][iy][iz]             + Fext*Ez_ext[ix][iy][iz]);
    Ezl += weight001 * (Ezth[ix][iy][iz - 1]         + Fext*Ez_ext[ix][iy][iz - 1]);
    Ezl += weight010 * (Ezth[ix][iy - 1][iz]         + Fext*Ez_ext[ix][iy - 1][iz]);
    Ezl += weight011 * (Ezth[ix][iy - 1][iz - 1]     + Fext*Ez_ext[ix][iy - 1][iz - 1]);
    Ezl += weight100 * (Ezth[ix - 1][iy][iz]         + Fext*Ez_ext[ix - 1][iy][iz]);
    Ezl += weight101 * (Ezth[ix - 1][iy][iz - 1]     + Fext*Ez_ext[ix - 1][iy][iz - 1]);
    Ezl += weight110 * (Ezth[ix - 1][iy - 1][iz]     + Fext*Ez_ext[ix - 1][iy - 1][iz]);
    Ezl += weight111 * (Ezth[ix - 1][iy - 1][iz - 1] + Fext*Ez_ext[ix - 1][iy - 1][iz - 1]);

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

    #ifdef __JDOTE_TEST__
    ubar[rest] = uptilde;
    vbar[rest] = vptilde;
    wbar[rest] = wptilde;
    #endif
    u[rest] = up;
    v[rest] = vp;
    w[rest] = wp;
  }

  return (0);                   // exit succcesfully (hopefully) 
}

int Particles3D::updatePosition(Grid * grid, VirtualTopology3D * vct, Field * EMf) {

#pragma omp parallel for
#pragma simd                    // this just slows things down (why?)
  for (long long rest = 0; rest < nop; rest++) {
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    // Update the positions with the new velocity
    xp = x[rest] + up * dt;
    yp = y[rest] + vp * dt;
    zp = z[rest] + wp * dt;

    if (Cylindrical == 1) {
       double gamma = atan2(zp,xp);
       xp=xp*cos(gamma)+zp*sin(gamma);
       zp=0.0;
       double uprot=up*cos(gamma)+wp*sin(gamma);
       wp=-up*sin(gamma)+wp*cos(gamma);
       up=uprot;
 
       if (xp > Lx) {
          xp = 2*Lx - xp;
          up = -up;
        }
     }

     x[rest] = xp;
     y[rest] = yp;
     z[rest] = zp;
     u[rest] = up;
     v[rest] = vp;
     w[rest] = wp;
  }                             // END OF ALL THE PARTICLES

  fixPosition();


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

int Particles3D::mover_PC_sub(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER P-C species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  double weights[2][2][2];
  double ***Ex = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx());
  double ***Ey = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy());
  double ***Ez = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz());
  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());

  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());

  double Fext = EMf->getFext();

  // const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
  // don't bother trying to push any particles simultaneously;
  // MIC already does vectorization automatically, and trying
  // to do it by hand only hurts performance.
//#pragma omp parallel for
//#pragma simd                    // this just slows things down (why?)

  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    const double xptilde = x[rest];
    const double yptilde = y[rest];
    const double zptilde = z[rest];
    double uptilde;
    double vptilde;
    double wptilde;

    double Exl = 0.0;
    double Eyl = 0.0;
    double Ezl = 0.0;
    double Bxl = 0.0;
    double Byl = 0.0;
    double Bzl = 0.0;
    int ix;
    int iy;
    int iz;

    // BEGIN OF SUBCYCLING LOOP

    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

    const double B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
    double       dt_sub     = M_PI*c/(4*abs(qom)*B_mag);
    const int    sub_cycles = int(dt/dt_sub) + 1;

    dt_sub = dt/double(sub_cycles);
    
    const double dto2 = .5 * dt_sub, qomdt2 = qom * dto2 / c;
    
    // if (sub_cycles>1) cout << " >> sub_cycles = " << sub_cycles << endl;

    for (int cyc_cnt = 0; cyc_cnt < sub_cycles; cyc_cnt++) {

      // calculate the average velocity iteratively
      int nit = NiterMover;
      if (sub_cycles > 2*NiterMover) nit = 1;

      for (int innter = 0; innter < nit; innter++) {
        // interpolation G-->P

        get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
        get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);
        get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);

        // end interpolation
        const double omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
        const double denom = 1.0 / (1.0 + omdtsq);
        // solve the position equation
        const double ut = up + qomdt2 * Exl;
        const double vt = vp + qomdt2 * Eyl;
        const double wt = wp + qomdt2 * Ezl;
        const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
        // solve the velocity equation 
        uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
        vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
        wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
        // update position
        xp = xptilde + uptilde * dto2;
        yp = yptilde + vptilde * dto2;
        zp = zptilde + wptilde * dto2;
      }                           // end of iteration
      // update the final position and velocity
      up = 2.0 * uptilde - u[rest];
      vp = 2.0 * vptilde - v[rest];
      wp = 2.0 * wptilde - w[rest];
      xp = xptilde + uptilde * dt_sub;
      yp = yptilde + vptilde * dt_sub;
      zp = zptilde + wptilde * dt_sub;
      x[rest] = xp;
      y[rest] = yp;
      z[rest] = zp;
      u[rest] = up;
      v[rest] = vp;
      w[rest] = wp;
    } // END  OF SUBCYCLING LOOP
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

  int ncr = 3;

  ////////////////////////
  // INJECTION FROM XLEFT
  ////////////////////////
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft == 2){ // use Field topology in this case
    long long particles_index=0;
    long long nplast = nop-1;

    while (particles_index < nplast+1) {
      if (x[particles_index] < ncr*dx ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }

    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i<=ncr ;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                //x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                //y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                //z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                x[particles_index] = Drand48(seed)*dx + grid->getXN(i,j,k);
                y[particles_index] = Drand48(seed)*dy + grid->getYN(i,j,k);
                z[particles_index] = Drand48(seed)*dz + grid->getZN(i,j,k);
                //x[particles_index] = (ii + 0.5)*(dx/npcelx) + grid->getXN(i,j,k);
                //y[particles_index] = (jj + 0.5)*(dy/npcely) + grid->getYN(i,j,k);
                //z[particles_index] = (kk + 0.5)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }
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
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }

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
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }
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
      if (x[particles_index] > (Lx-ncr*dx) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double prob, theta, sign;
    //particles_index;
    for (int i=(grid->getNXC()-ncr-1); i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                //x[particles_index] = (ii + Drand48(seed))*(dx/npcelx) + grid->getXN(i,j,k);
                //y[particles_index] = (jj + Drand48(seed))*(dy/npcely) + grid->getYN(i,j,k);
                //z[particles_index] = (kk + Drand48(seed))*(dz/npcelz) + grid->getZN(i,j,k);
                x[particles_index] = Drand48(seed)*dx + grid->getXN(i,j,k);
                y[particles_index] = Drand48(seed)*dy + grid->getYN(i,j,k);
                z[particles_index] = Drand48(seed)*dz + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }
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
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }
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
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOLc(i,j,k));
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

                if (Lagrangian == 1) {
                  u[particles_index] -= u_bulk;
                  v[particles_index] -= v_bulk;
                  w[particles_index] -= w_bulk;
                }
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

/** apply a linear perturbation to particle distribution */
void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  double value1 = 0.0, value2 = 0.0, max_value = 0.0, min_value = 0.0, phi, n;
  long long counter = 0, total_generated = 0;
  bool rejected;
  double harvest, prob, theta;
  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  double integral = 0.0;

  alpha = deltaBoB * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0) + EMf->getBz(1, 1, 0) * EMf->getBz(1, 1, 0)) / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;



  // find the maximum value of f=1+delta_f/f0
  for (register double vpar = -2 * uth; vpar <= 2 * uth; vpar += 0.0005)
    for (register double vperp = 1e-10; vperp <= 2 * vth; vperp += 0.0005)
      for (register double X = xstart; X <= xend; X += 2 * grid->getDX())
        for (register double Y = ystart; Y <= yend; Y += 2 * grid->getDY()) {
          value1 = 1 + delta_f(vpar, vperp, 0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(vpar, vperp);

          if (value1 > max_value)
            max_value = value1;


        }



  max_value *= 3.2;
  phi = 1.48409;
  n = 2.948687;                 // security factor...
  if (ns == 1) {
    max_value *= 3.0;
    phi = -1.65858;
    n = 2.917946;
  }                             // security factor...
  cout << "max-value=" << max_value << " min-value=" << min_value << endl;

  /* initialize random generator */

  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int ii = 0; ii < npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))); ii++)
        for (int jj = 0; jj < npcely; jj++) {
          x[counter] = (ii + .5) * (dx / (npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))))) + grid->getXN(i, j, 0);
          y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, 0);
          q[counter] = (qom / fabs(qom)) * ((0.19635) / npcel) * (1.0 / invVOL);

          // apply rejection method in velocity space
          rejected = true;
          while (rejected) {
            total_generated++;
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            // u
            u[counter] = u0 + uth * prob * cos(theta);
            // v
            v[counter] = v0 + vth * prob * sin(theta);
            // w
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            w[counter] = w0 + wth * prob * cos(theta);

            // test: if rand < (1+delta_f/f0)/max_value --> accepted
            if (rand() / (double) RAND_MAX <= (1 + delta_f(u[counter], v[counter], w[counter], x[counter], y[counter], kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(u[counter], sqrt(v[counter] * v[counter] + w[counter] * w[counter]))) / max_value)
              rejected = false;

          }
          if (TrackParticleID)
            ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
          counter++;
        }
  nop = counter + 1;
  // if (vct->getCartesian_rank()==0)
  cout << "Rejection method: " << (counter + 1) / double (total_generated) * 100 << " % of particles are accepted for species " << ns << " counter=" << counter << endl;
}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double theta, Field * EMf) {
  const complex < double >I(0.0, 1.0);
  const double vperp = sqrt(v * v + w * w);
  const double vpar = u;
  const double kpar = kx;
  double kperp;
  if (ky == 0.0)                // because this formula is not valid for exactly parallel
    kperp = 1e-9;
  else
    kperp = ky;
  const double om_c = qom / c * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0)) / 2 / M_PI;
  const double phi = atan2(w, v);
  const double lambda = kperp * vperp / om_c;
  const complex < double >omega(omega_re, omega_i);

  const int lmax = 5;           // sum from -lmax to lmax

  double bessel_Jn_array[lmax + 2];
  double bessel_Jn_prime_array[lmax + 1];
  complex < double >a1[2 * lmax + 1], a2[2 * lmax + 1], a3[2 * lmax + 1];
  complex < double >factor, deltaf;

  // rotation of x,y
  double temp;
  temp = x;
  x = x * cos(theta) - y * sin(theta);
  y = temp * sin(theta) + y * cos(theta);


  /** for compilation issues comment this part: PUT in the math stuff */
  // calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
  factor = (kpar * vperp / omega * df0_dvpar(vpar, vperp) + (1.0 - (kpar * vpar / omega)) * df0_dvperp(vpar, vperp));
  for (register int l = -lmax; l < 0; l++) {  // negative index
    a1[l + lmax] = factor / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * 0.5 * pow(-1.0, -l) * (bessel_Jn_array[-l - 1] - bessel_Jn_array[-l + 1]);
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * pow(-1.0, -l) * bessel_Jn_array[-l];
  }

  for (register int l = 0; l < lmax + 1; l++) { // positive index
    a1[l + lmax] = factor / lambda * bessel_Jn_array[l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * bessel_Jn_prime_array[l];
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * bessel_Jn_array[l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * bessel_Jn_array[l];
  }

  deltaf = (0.0, 0.0);
  for (register int l = -lmax; l < lmax + 1; l++) {
    deltaf += (a3[l + lmax] * Ex_mod * exp(I * Ex_phase) + a1[l + lmax] * Ey_mod * exp(I * Ey_phase) + a2[l + lmax] * Ez_mod * exp(I * Ez_phase)) / (kpar * vpar + l * om_c - omega) * exp(-I * phi * (double) l);
  }
  deltaf *= I * qom * exp(I * lambda * sin(phi)) * exp(I * (2 * M_PI * kx * x + 2 * M_PI * ky * y));

  return (real(deltaf));
}

double Particles3D::df0_dvpar(double vpar, double vperp) {
  double result;
  result = -2 * (vpar - u0) / uth / uth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::df0_dvperp(double vpar, double vperp) {
  double result;
  result = -2 * (vperp) / vth / vth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::f0(double vpar, double vperp) {
  double result;
  result = exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

void Particles3D::RotatePlaneXY(double theta) {
  double temp, temp2;
  for (register long long s = 0; s < nop; s++) {
    temp = u[s];
    temp2 = v[s];
    u[s] = temp * cos(theta) + v[s] * sin(theta);
    v[s] = -temp * sin(theta) + temp2 * cos(theta);
  }
}

/*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
double Particles3D::deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center){

  long long np_current = 0;
  long long nplast     = nop-1;

  bool twodim = false;
  if ( nzn == 4) twodim = true;

  while (np_current < nplast+1){

    double xd = x[np_current] - x_center;
    double yd = y[np_current] - y_center;
    double zd = z[np_current] - z_center;
    if (twodim) zd = 0.0;

    if ( (xd*xd+yd*yd+zd*zd) < R*R ){
      Q_removed += q[np_current];
      del_pack(np_current,&nplast);

    } else {
      np_current++;
    }
  }
  nop = nplast +1;
  return(Q_removed);
}


/** Delete the particles outer multx*dx, multy*dy, multz*dz celles */
double Particles3D::deleteParticlesOuterFrame(double multx, double multy, double multz){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
        if (Cylindrical == 1) {
	  while (np_current < nplast+1) {
            if (x[np_current] > (Lx-multx*dx) ||
                y[np_current] < multy*dy || y[np_current] > (Ly-multy*dx)) {
                Q_removed += q[np_current];
                del_pack(np_current,&nplast);
            } else {
              np_current++;
            }
	  }
        } else {
	  while (np_current < nplast+1) {
            if (x[np_current] < multx*dx || x[np_current] > (Lx-multx*dx) ||
                y[np_current] < multy*dy || y[np_current] > (Ly-multy*dy) ) {
                //z[np_current] < multz*dz || z[np_current] > (Lz-multz*dz)) {
                Q_removed += q[np_current];
                del_pack(np_current,&nplast);
            } else {
              np_current++;
            }
	  }
        }

	nop = nplast +1;
	return(Q_removed);
}

/** Delete the particles inside the sphere with radius R and center x_center y_center */
double Particles3D::deleteParticlesOutsideSphere(double R, double x_center, double y_center, double z_center){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	while (np_current < nplast+1){
		if (sqrt( pow(x[np_current] - x_center,2) + pow(y[np_current] - y_center,2) + pow(z[np_current] - z_center,2) ) > R){
			// delete the particle and pack the particle array, the value of nplast changes
			Q_removed += q[np_current];
			del_pack(np_current,&nplast);

		} else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}
	}
	nop = nplast +1;
	return(Q_removed);
}

//** Delete the particles outside the cube with dimension L */
double Particles3D::deleteParticlesOutsideBox(double L){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	while (np_current < nplast+1){
		if (x[np_current] > L || x[np_current] < 0 ||
			y[np_current] > L || y[np_current] < 0 ||
			z[np_current] > L || z[np_current] < 0){
			// delete the particle and pack the particle array, the value of nplast changes
			Q_removed += q[np_current];
			del_pack(np_current,&nplast);

		} else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}
	}
	nop = nplast +1;
	return(Q_removed);
}

/* Computes overalp of two intervals */
double Particles3D::interval_overlap(double xd0, double xd1, double xp0, double xp1, double& xstart_interval)
{
  double dx_inject;
  xstart_interval = max(xp0, xd0);
  double xend_interval = min(xp1,xd1);
  if(xstart_interval < xend_interval )
    dx_inject = xend_interval-xstart_interval;
  else
    dx_inject = 0.0;
  return(dx_inject);
}


/* Inject from a box, checking overal with the local processor */
int Particles3D::injector_rand_box(Grid* grid,VirtualTopology3D* vct, Field* EMf, double x_center_inject, double y_center_inject, double z_center_inject, double L_inject)
{
	double harvest;
	double prob, theta, sign;
	double r;
	int avail;
	long long store_nop=nop;
	long long counter=nop;

	//Vinj from the input file has hte number of particles to be injected
	//Ninj is the RHOinject from the input file
	double Iinj=Vinj;
	long long npinject=1;
	if (Iinj > 1)
		npinject = (int)floor(Iinj);

	double xstart_inject;
	double dx_inject = interval_overlap(x_center_inject-L_inject,x_center_inject+L_inject, xstart, xend, xstart_inject);
	double ystart_inject;
	double dy_inject = interval_overlap(y_center_inject-L_inject,y_center_inject+L_inject, ystart, yend, ystart_inject);
	double zstart_inject;
	double dz_inject = interval_overlap(z_center_inject-L_inject,z_center_inject+L_inject, zstart, zend, zstart_inject);

//	cout << xstart_inject  << "  " << ystart_inject << "   " << zstart_inject << endl;
	if (dx_inject * dy_inject * dz_inject > 0.0)
	{

//cout << dx_inject * dy_inject * dz_inject  << "  " << npinject << endl;

		for (int inject = 0; inject < npinject; inject++)
		{
			harvest =   rand()/(double)RAND_MAX ;
			x[nop] = xstart_inject+(harvest) * dx_inject;
			harvest =   rand()/(double)RAND_MAX ;
			y[nop] = ystart_inject+(harvest) * dy_inject;
			harvest =   rand()/(double)RAND_MAX ;
			z[nop] = zstart_inject+(harvest) * dz_inject;
			r = 1e-10+sqrt((x[nop]-Lx/2.0)*(x[nop]-Lx/2.0) +
					(y[nop]-Ly/2.0)*(y[nop]-Ly/2.0) +
					(z[nop]-Lz/2.0)*(z[nop]-Lz/2.0));

			q[nop] =  (qom/fabs(qom)) * rhoINJECT * dx_inject * dy_inject * dz_inject / npinject;
			// u
			harvest =   rand()/(double)RAND_MAX;
			prob  = sqrt(-2.0*log(1.0-.999999*harvest));
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;

			u[nop] = - v0 * (x[nop]-Lx/2.0)/r + uth*prob*cos(theta);
			v[nop] = - v0 * (y[nop]-Ly/2.0)/r + vth*prob*sin(theta);
			harvest =   rand()/(double)RAND_MAX;
			prob  = sqrt(-2.0*log(1.0-.999999*harvest));
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;
			w[nop] = - v0 * (z[nop]-Lz/2.0)/r + wth*prob*cos(theta);;


			if (TrackParticleID)
				ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

			//printf("*** DEBUG: Added a particle at (x,y,z) = (%f,%f,%f) with (vx,vy,vz) = (%f,%f,%f)\n", x[nop], y[nop], z[nop], u[nop], v[nop], w[nop]);

			nop++;
//			cout << nop << endl;
		}

	}
		// move particles inside the sim box
		for (int i=store_nop; i < nop; i++){
			x[i] += u[i]*dt;
			y[i] += v[i]*dt;
			z[i] += w[i]*dt;
		}


		// ******************** //
		// COMMUNICATION
		// ******************* //
	    avail = communicate(vct);
	    if (avail < 0){
		cout <<"Cannot communicate! Failing!!!!\n";
		return(-1);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
//	    cout << "Past barrier...\n";
	    // communicate again if particles are not in the correct domain
	    while(isMessagingDone(vct) >0){
			// COMMUNICATION
			avail = communicate(vct);

			if (avail < 0)
		        return(-1);
			MPI_Barrier(MPI_COMM_WORLD);
		}


//	cout << "Communicated!!\n";

	return(0); // exit succcesfully (hopefully)

}


/* Inject from a box, checking overal with the local processor */
int Particles3D::injector_rand_box_mono(Grid* grid,VirtualTopology3D* vct, Field* EMf, double x_center_inject, double y_center_inject, double z_center_inject, double L_inject)
{
	double harvest;
	double prob, theta, sign;
	double r;
	int avail;
	long long store_nop=nop;
	long long counter=nop;

	//Vinj from the input file has hte number of particles to be injected
	//Ninj is the RHOinject from the input file
	double Iinj=Vinj;
	long long npinject=1;
	if (Iinj > 1)
		npinject = (int)floor(Iinj);

	double xstart_inject;
	double dx_inject = interval_overlap(x_center_inject-L_inject,x_center_inject+L_inject, xstart, xend, xstart_inject);
	double ystart_inject;
	double dy_inject = interval_overlap(y_center_inject-L_inject,y_center_inject+L_inject, ystart, yend, ystart_inject);
	double zstart_inject;
	double dz_inject = interval_overlap(z_center_inject-L_inject,z_center_inject+L_inject, zstart, zend, zstart_inject);

//	cout << xstart_inject  << "  " << ystart_inject << "   " << zstart_inject << endl;
	if (dx_inject * dy_inject * dz_inject > 0.0)
	{

//cout << dx_inject * dy_inject * dz_inject  << "  " << npinject << endl;

		for (int inject = 0; inject < npinject; inject++)
		{
			harvest =   rand()/(double)RAND_MAX ;
			x[nop] = xstart_inject+(harvest) * dx_inject;
			harvest =   rand()/(double)RAND_MAX ;
			y[nop] = ystart_inject+(harvest) * dy_inject;
			harvest =   rand()/(double)RAND_MAX ;
			z[nop] = zstart_inject+(harvest) * dz_inject;
			r = 1e-10+sqrt((x[nop]-Lx/2.0)*(x[nop]-Lx/2.0) +
					(y[nop]-Ly/2.0)*(y[nop]-Ly/2.0) +
					(z[nop]-Lz/2.0)*(z[nop]-Lz/2.0));

			q[nop] =  (qom/fabs(qom)) * rhoINJECT * dx_inject * dy_inject * dz_inject / npinject;
			// u
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;

			u[nop] = - v0 * (x[nop]-Lx/2.0)/r + uth*cos(theta);
			v[nop] = - v0 * (y[nop]-Ly/2.0)/r + vth*sin(theta);
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;
			w[nop] = - v0 * (z[nop]-Lz/2.0)/r + wth*cos(theta);;


			if (TrackParticleID)
				ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

			//printf("*** DEBUG: Added a particle at (x,y,z) = (%f,%f,%f) with (vx,vy,vz) = (%f,%f,%f)\n", x[nop], y[nop], z[nop], u[nop], v[nop], w[nop]);

			nop++;
//			cout << nop << endl;
		}

	}
		// move particles inside the sim box
		for (int i=store_nop; i < nop; i++){
			x[i] += u[i]*dt;
			y[i] += v[i]*dt;
			z[i] += w[i]*dt;
		}


		// ******************** //
		// COMMUNICATION
		// ******************* //
	    avail = communicate(vct);
	    if (avail < 0){
		cout <<"Cannot communicate! Failing!!!!\n";
		return(-1);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
//	    cout << "Past barrier...\n";
	    // communicate again if particles are not in the correct domain
	    while(isMessagingDone(vct) >0){
			// COMMUNICATION
			avail = communicate(vct);

			if (avail < 0)
		        return(-1);
			MPI_Barrier(MPI_COMM_WORLD);
		}


//	cout << "Communicated!!\n";

	return(0); // exit succcesfully (hopefully)

}


/*
int Particles3D::injector_rand_box_mono(Grid* grid,VirtualTopology3D* vct, Field* EMf)
{
	double harvest;
	double prob, theta, sign, fi;
	int avail;
	long long store_nop=nop;
	long long counter=nop;

	//double Iinj=0.002;
	double Iinj=Vinj;
	long long npinject=1;
	if (Iinj > 1)
		npinject = (int)floor(Iinj);

	if (vct->getCartesian_rank() == 0){
			cout << "*** Injector Rand Box species " << ns << " ***" << NiterMover <<" INJECTING " << npinject << " particles per cell, each with a charge of " << (qom/fabs(qom))*(EMf->getRHOcs(1,1,1,ns)/npcel)*(1.0/grid->getInvVOLc(1,1,1)) << " C  ****" << endl;
        }


                  for (int i = 1; i < grid->getNXC() - 1; i++)
                    for (int j = 1; j < grid->getNYC() - 1; j++)
                      for (int k = 1; k < grid->getNZC() - 1; k++) {
			//if (nop > 10000) continue;
		        for (long long inject = 0; inject < npinject; inject++)
		        {
		           double xp = grid->getXN(i, j, k) + dx*Drand48(seed);	
		           double yp = grid->getYN(i, j, k) + dy*Drand48(seed);	
		           double zp = grid->getZN(i, j, k) + dz*Drand48(seed);	
			   if (xp< 10*dx) {
			   //if (fabs(xp - x_center) < L_square && fabs(yp - y_center) < L_square && fabs(zp - z_center) < L_square) {
                              x[nop] = xp;
                              y[nop] = yp;
                              z[nop] = zp;
                              // q = charge
                              q[nop] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOLc(i,j,k));
                 
                              //nop
                              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
                              theta = 2.0 * M_PI * Drand48(seed);
                              u[nop] = u0 + uth * prob * cos(theta);
                              // v
                              v[nop] = v0 + vth * prob * sin(theta);
                              // w
                              prob = sqrt(-2.0 * log(1.0 - .999999 * Drand48(seed)));
                              theta = 2.0 * M_PI * Drand48(seed);
                              w[nop] = w0 + wth * prob * cos(theta);
                 
                              if (TrackParticleID)
                                ParticleID[nop] = nop * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
                 
			     nop++;
                           }
                         }
		      }

	return(0); // exit succcesfully (hopefully)

}
*/

/*! Set the z coordintaes equal zero for cylindrical coordinates */
void Particles3D::fixZcylindrical()
{
  for (int i=0; i<nop; i++) {
    z[i] = 0;
  }
}

/*! Set to zero the corresponding components of the position of the particles for 1D and 2D cases */
void Particles3D::fixPosition()
{
  if (ndim == 1) {
    for (int i=0; i<nop; i++) {
      z[i] = 0;
      y[i] = 0;
    }
  } else if (ndim == 2) {
    for (int i=0; i<nop; i++) z[i] = 0;
  }
}


void Particles3D::setLagrangianFrame(bool mode) {
#pragma omp parallel for
  if (mode) { 
    for (long long rest = 0; rest < nop; rest++) {
      // Update the positions with the new velocity
      u[rest] -= u_bulk;
      v[rest] -= v_bulk;
      w[rest] -= w_bulk;
    }                            
  } else {
    for (long long rest = 0; rest < nop; rest++) {
      // Update the positions with the new velocity
      u[rest] += u_bulk;
      v[rest] += v_bulk;
      w[rest] += w_bulk;
    }                            
  }

}
int Particles3D::updateLagrangianFrame(VirtualTopology3D * vct) {

#pragma omp parallel for
#pragma simd                    // this just slows things down (why?)
  for (long long rest = 0; rest < nop; rest++) {
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    // Update the positions with the new velocity
    xp = x[rest] + u_bulk * dt;
    yp = y[rest] + v_bulk * dt;
    zp = z[rest] + w_bulk * dt;

    if (Cylindrical == 1) {
       double gamma = atan2(zp,xp);
       xp=xp*cos(gamma)+zp*sin(gamma);
       zp=0.0;
       double uprot=up*cos(gamma)+wp*sin(gamma);
       wp=-up*sin(gamma)+wp*cos(gamma);
       up=uprot;
 
       if (xp > Lx) {
          xp = 2*Lx - xp;
          up = -up;
        }
     }

     x[rest] = xp;
     y[rest] = yp;
     z[rest] = zp;
     u[rest] = up;
     v[rest] = vp;
     w[rest] = wp;
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

