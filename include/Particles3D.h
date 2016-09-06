/*******************************************************************************************
  Particles3D.h  -  Class for particles of the same species, in a 2D space and 3 component velocity
  -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
 ********************************************************************************************/

#ifndef Part2D_H
#define Part2D_H

#include "Particles3Dcomm.h"

/**
 * 
 * Class for particles of the same species, in a 2D space and 3 component velocity
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Particles3D:public Particles3Dcomm {
  public:
    /** constructor */
    Particles3D();
    /** destructor */
    ~Particles3D();
    /** Initial condition: uniform in space and motionless */
    void uniform_background(Grid * grid, Field * EMf);
    /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
      <ul>
      <li> dim = 0 --> constant velocity on X direction </li>
      <li> dim = 1 --> constant velocity on Y direction </li>
      </ul>
      */
    void constantVelocity(double vel, int dim, Grid * grid, Field * EMf);
    /** Initial condition: uniform in space and maxwellian in velocity */
    void maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
    /** mover for the ECIMM algorith */
    int mover_ECIMM(Grid * grid, VirtualTopology3D * vct, Field * EMf);

    void get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext);
    void get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez);
    void get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]);
    /** particle repopulator */
    int particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is);
};


#endif
