
//#define ELECTROSTATIC
#include "EMfields3D.h"
//#include "MaxwellECIMM.h"

#define NE_MASS 14
/*! constructor */
EMfields3D::EMfields3D(Collective * col, Grid * grid) {
  nxc = grid->getNXC();
  nxn = grid->getNXN();
  nyc = grid->getNYC();
  nyn = grid->getNYN();
  nzc = grid->getNZC();
  nzn = grid->getNZN();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  invVOL = grid->getInvVOL();
  xStart = grid->getXstart();
  xEnd = grid->getXend();
  yStart = grid->getYstart();
  yEnd = grid->getYend();
  zStart = grid->getZstart();
  zEnd = grid->getZend();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  ns = col->getNs();
  c = col->getC();
  dt = col->getDt();
  th = col->getTh();
  ue0 = col->getU0(0);
  ve0 = col->getV0(0);
  we0 = col->getW0(0);

  delt = c * th * dt;
  GMREStol = col->getGMREStol();
  qom = new double[ns];
  for (int i = 0; i < ns; i++)
    qom[i] = col->getQOM(i);

  bcEMfaceXright = col->getBcEMfaceXright();
  bcEMfaceXleft = col->getBcEMfaceXleft();
  bcEMfaceYright = col->getBcEMfaceYright();
  bcEMfaceYleft = col->getBcEMfaceYleft();
  bcEMfaceZright = col->getBcEMfaceZright();
  bcEMfaceZleft = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  // Earth Simulation
  B1x = col->getB1x();
  B1y = col->getB1y();
  B1z = col->getB1z();
  delta = col->getDelta();
  Smooth = col->getSmooth();
  // get the density background for the gem Challange
  rhoINIT   = new double[ns];
  rhoINJECT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i = 0; i < ns; i++) {
    rhoINIT  [i] = col->getRHOinit(i);
    rhoINJECT[i] = col->getRHOinject(i);
    if ((fabs(col->getW0(i)) != 0) || (fabs(col->getU0(i)) != 0)) // GEM and LHDI
      DriftSpecies[i] = true;
    else
      DriftSpecies[i] = false;
  }
  /*! parameters for GEM challenge */
  FourPI = 16 * atan(1.0);
  /*! Restart */
  restart1 = col->getRestart_status();
  RestartDirName = col->getRestartDirName();
  Case = col->getCase();

  // arrays allocation: nodes
  Ex = newArr3(double, nxn, nyn, nzn);
  Ey = newArr3(double, nxn, nyn, nzn);
  Ez = newArr3(double, nxn, nyn, nzn);
  Exth = newArr3(double, nxn, nyn, nzn);
  Eyth = newArr3(double, nxn, nyn, nzn);
  Ezth = newArr3(double, nxn, nyn, nzn);
  Bxn = newArr3(double, nxn, nyn, nzn);
  Byn = newArr3(double, nxn, nyn, nzn);
  Bzn = newArr3(double, nxn, nyn, nzn);
  rhon = newArr3(double, nxn, nyn, nzn);
  Jxh = newArr3(double, nxn, nyn, nzn);
  Jyh = newArr3(double, nxn, nyn, nzn);
  Jzh = newArr3(double, nxn, nyn, nzn);

  J_partx = newArr3(double, nxn, nyn, nzn);
  J_party = newArr3(double, nxn, nyn, nzn);
  J_partz = newArr3(double, nxn, nyn, nzn);

  Mxx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mxy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mxz = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myz = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzz = newArr4(double, NE_MASS, nxn, nyn, nzn);

  // External imposed fields
  Bx_ext = newArr3(double,nxn,nyn,nzn);
  By_ext = newArr3(double,nxn,nyn,nzn);
  Bz_ext = newArr3(double,nxn,nyn,nzn);
  // involving species
  rhons = newArr4(double, ns, nxn, nyn, nzn);
  rhocs = newArr4(double, ns, nxc, nyc, nzc);
  // arrays allocation: central points 
  Bxc = newArr3(double, nxc, nyc, nzc);
  Byc = newArr3(double, nxc, nyc, nzc);
  Bzc = newArr3(double, nxc, nyc, nzc);
  rhoc = newArr3(double, nxc, nyc, nzc);
  rhoh = newArr3(double, nxc, nyc, nzc);

  // temporary arrays
  tempXC = newArr3(double, nxc, nyc, nzc);
  tempYC = newArr3(double, nxc, nyc, nzc);
  tempZC = newArr3(double, nxc, nyc, nzc);

  tempXN = newArr3(double, nxn, nyn, nzn);
  tempYN = newArr3(double, nxn, nyn, nzn);
  tempZN = newArr3(double, nxn, nyn, nzn);
  tempC = newArr3(double, nxc, nyc, nzc);
  tempX = newArr3(double, nxn, nyn, nzn);
  tempY = newArr3(double, nxn, nyn, nzn);
  tempZ = newArr3(double, nxn, nyn, nzn);
  temp2X = newArr3(double, nxn, nyn, nzn);
  temp2Y = newArr3(double, nxn, nyn, nzn);
  temp2Z = newArr3(double, nxn, nyn, nzn);
  imageEX = newArr3(double, nxc, nyc, nzc);
  imageEY = newArr3(double, nxc, nyc, nzc);
  imageEZ = newArr3(double, nxc, nyc, nzc);
  imageBX = newArr3(double, nxn, nyn, nzn);
  imageBY = newArr3(double, nxn, nyn, nzn);
  imageBZ = newArr3(double, nxn, nyn, nzn);
  vectX = newArr3(double, nxn, nyn, nzn);
  vectY = newArr3(double, nxn, nyn, nzn);
  vectZ = newArr3(double, nxn, nyn, nzn);
  arr = newArr3(double,nxn,nyn,nzn);

}


/* Prepare the memory for E calculation */
void EMfields3D::startEcalc(VirtualTopology3D * vct) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;
}

void EMfields3D::endEcalc(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid)
{
  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, tempXC, tempYC, tempZC, xkrylov, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Bxc, tempXC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Byc, tempYC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Bzc, tempZC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);

  // communicate so the interpolation can have good values
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bxc,   col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Byc,   col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc,   col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

}


/*! Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
void EMfields3D::calculateEB(Grid * grid, VirtualTopology3D * vct, Collective *col) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;

  int dim = 3*(nxn - 2) * (nyn - 2) * (nzn - 2) + 3*(nxc - 2) * (nyc - 2) * (nzc - 2);

  double *xkrylov = new double[dim];  
  double *bkrylov = new double[dim];  

  // set to zero all the stuff 
  eqValue(0.0, xkrylov, dim);
  eqValue(0.0, bkrylov, dim);

  if (vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER ***" << endl;
  // prepare the source 
  MaxwellSource(bkrylov, grid, vct, col);
  phys2solver(xkrylov, Ex, Ey, Ez, Bxn, Byn, Bzn, nxn, nyn, nzn);
  // solver
  GMRES(&Field::MaxwellImage, xkrylov, dim, bkrylov, 20, 200, GMREStol, grid, vct, this);
  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, tempXC, tempYC, tempZC, xkrylov, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Bxc, tempXC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Byc, tempYC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Bzc, tempZC, nxc, nyc, nzc);
  addscale(1 / th, -(1.0 - th) / th, Ex, Exth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ey, Eyth, nxn, nyn, nzn);
  addscale(1 / th, -(1.0 - th) / th, Ez, Ezth, nxn, nyn, nzn);

  // communicate so the interpolation can have good values
  communicateNodeBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ex,   col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ey,   col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Ez,   col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;
}

/*! Calculate sorgent for Maxwell solver */
void EMfields3D::MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  // 
  #ifdef ELECTROSTATIC
  eqValue(0, Bxc, nxc, nyc, nzc);
  eqValue(0, Byc, nxc, nyc, nzc);
  eqValue(0, Bzc, nxc, nyc, nzc);
  #else
  scale(tempXC, Bxc, 1, nxc, nyc, nzc);
  scale(tempYC, Byc, 1, nxc, nyc, nzc);
  scale(tempZC, Bzc, 1, nxc, nyc, nzc);
  #endif

  // E - 4 pi dt th Jhat
  scale(tempX, Jxh, -4*M_PI*th*dt, nxn, nyn, nzn);
  scale(tempY, Jyh, -4*M_PI*th*dt, nxn, nyn, nzn);
  scale(tempZ, Jzh, -4*M_PI*th*dt, nxn, nyn, nzn);
  addscale(1, tempX, Ex, nxn, nyn, nzn);
  addscale(1, tempY, Ey, nxn, nyn, nzn);
  addscale(1, tempZ, Ez, nxn, nyn, nzn);


  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, tempXC, tempYC, tempZC, nxn, nyn, nzn);

}


/*! Mapping of Maxwell image to give to solver */
void EMfields3D::MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct) {
  int dim = 3*(nxn - 2) * (nyn - 2) * (nzn - 2) + 3*(nxc - 2) * (nyc - 2) * (nzc - 2);
  eqValue(0.0, im, dim);
  //eqValue(0.0, im, 6 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  // move from krylov space to physical space: E -> vect   B -> temp
  solver2phys(tempX, tempY, tempZ, tempXC, tempYC, tempZC, vector, nxn, nyn, nzn);


  communicateCenterBC(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);
  //communicateNodeBC(nxn, nyn, nzn, tempX, 2, 2, 2, 2, 2, 2, vct);
  //communicateNodeBC(nxn, nyn, nzn, tempY, 2, 2, 2, 2, 2, 2, vct);
  //communicateNodeBC(nxn, nyn, nzn, tempZ, 2, 2, 2, 2, 2, 2, vct);
  // WARNING: this in provisional: communicateNodeBC_P does not copy the ghost edges
  // correctly
  communicateNode_P(nxn, nyn, nzn, tempX, vct);
  communicateNode_P(nxn, nyn, nzn, tempY, vct);
  communicateNode_P(nxn, nyn, nzn, tempZ, vct);


  #ifdef ELECTROSTATIC
    for (int i=1; i<nxn-1; i++)
      for (int j=1; j<nyn-1; j++) 
        for (int k=1; k<nzn-1; k++) {
          imageEX[i][j][k] = tempXC[i][j][k];
          imageEY[i][j][k] = tempYC[i][j][k];
          imageEZ[i][j][k] = tempZC[i][j][k];
        }
  #else
  // c*th*dt*curl(Eth) + Bth
  grid->curlN2C(imageEX, imageEY, imageEZ, tempX, tempY, tempZ);
  scale(imageEX, c*th*dt, nxc, nyc, nzc);
  scale(imageEY, c*th*dt, nxc, nyc, nzc);
  scale(imageEZ, c*th*dt, nxc, nyc, nzc);
  addscale(1, imageEX, tempXC, nxc, nyc, nzc);
  addscale(1, imageEY, tempYC, nxc, nyc, nzc);
  addscale(1, imageEZ, tempZC, nxc, nyc, nzc);
  #endif

  // c*dt*th*curl(Bth) - Eth - dt*th*4*pi*Sum(M Eth)
  #ifdef ELECTROSTATIC
  eqValue(0, imageBX, nxn, nyn, nzn);
  eqValue(0, imageBY, nxn, nyn, nzn);
  eqValue(0, imageBZ, nxn, nyn, nzn);
  #else
  grid->curlC2N(imageBX, imageBY, imageBZ, tempXC, tempYC, tempZC);
  scale(imageBX, -c*th*dt, nxn, nyn, nzn);
  scale(imageBY, -c*th*dt, nxn, nyn, nzn);
  scale(imageBZ, -c*th*dt, nxn, nyn, nzn);
  #endif

  addscale(1, imageBX, tempX, nxn, nyn, nzn);
  addscale(1, imageBY, tempY, nxn, nyn, nzn);
  addscale(1, imageBZ, tempZ, nxn, nyn, nzn);

  double fac = dt*th*FourPI;
  for (int i=1; i<nxn-1; i++) 
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        double MEx, MEy, MEz;
        productMassE(&MEx, &MEy, &MEz, tempX, tempY, tempZ, i, j, k);
        imageBX[i][j][k] += fac*MEx;
        imageBY[i][j][k] += fac*MEy;
        imageBZ[i][j][k] += fac*MEz;
      }

  // Boundary condition in the known term
  // boundary condition: Xleft
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcEMfaceXleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 0);
  // boundary condition: Xright
  if (vct->getXright_neighbor() == MPI_PROC_NULL && bcEMfaceXright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 0);
  // boundary condition: Yleft
  if (vct->getYleft_neighbor() == MPI_PROC_NULL && bcEMfaceYleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 1);
  // boundary condition: Yright
  if (vct->getYright_neighbor() == MPI_PROC_NULL && bcEMfaceYright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 1);
  // boundary condition: Zleft
  if (vct->getZleft_neighbor() == MPI_PROC_NULL && bcEMfaceZleft == 0)  // perfect conductor
    perfectConductorLeftS(tempX, tempY, tempZ, 2);
  // boundary condition: Zright
  if (vct->getZright_neighbor() == MPI_PROC_NULL && bcEMfaceZright == 0)  // perfect conductor
    perfectConductorRightS(tempX, tempY, tempZ, 2);

  // move from physical space to krylov space
  phys2solver(im, imageBX, imageBY, imageBZ, imageEX, imageEY, imageEZ, nxn, nyn, nzn);
}


// Claculate the product of the mass matrix with the vector v
void EMfields3D::productMassE(double* MEx, double* MEy, double* MEz, double ***vX, double ***vY, double ***vZ, int i, int j, int k)
{
  double resX = 0;
  double resY = 0;
  double resZ = 0;

  // First the nodes that have been computed in the node i,j,k
  for (int g=0; g<NE_MASS; g++) {
    int i1 = i + NeNo.getX(g);
    int j1 = j + NeNo.getY(g);
    int k1 = k + NeNo.getZ(g);

    resX += vX[i1][j1][k1]*Mxx[g][i][j][k] 
          + vY[i1][j1][k1]*Myx[g][i][j][k] 
          + vZ[i1][j1][k1]*Mzx[g][i][j][k];
                              
    resY += vX[i1][j1][k1]*Mxy[g][i][j][k] 
          + vY[i1][j1][k1]*Myy[g][i][j][k] 
          + vZ[i1][j1][k1]*Mzy[g][i][j][k];
                              
    resZ += vX[i1][j1][k1]*Mxz[g][i][j][k] 
          + vY[i1][j1][k1]*Myz[g][i][j][k] 
          + vZ[i1][j1][k1]*Mzz[g][i][j][k];

    if (g == 0) continue;
    int i2 = i - NeNo.getX(g);
    int j2 = j - NeNo.getY(g);
    int k2 = k - NeNo.getZ(g);
    resX += vX[i2][j2][k2]*Mxx[g][i2][j2][k2] 
          + vY[i2][j2][k2]*Myx[g][i2][j2][k2] 
          + vZ[i2][j2][k2]*Mzx[g][i2][j2][k2];
                              
    resY += vX[i2][j2][k2]*Mxy[g][i2][j2][k2] 
          + vY[i2][j2][k2]*Myy[g][i2][j2][k2] 
          + vZ[i2][j2][k2]*Mzy[g][i2][j2][k2];
                              
    resZ += vX[i2][j2][k2]*Mxz[g][i2][j2][k2] 
          + vY[i2][j2][k2]*Myz[g][i2][j2][k2] 
          + vZ[i2][j2][k2]*Mzz[g][i2][j2][k2];


  }

  *MEx = resX;
  *MEy = resY;
  *MEz = resZ;
}


/*! interpolate charge density and pressure density from node to center */
void EMfields3D::interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid) {
  // do we need communication or not really?
  grid->interpN2C(rhoc, rhon);
}

/*! communicate ghost for grid -> Particles interpolation */
void EMfields3D::communicateGhostP2G(VirtualTopology3D * vct) {

  // interpolate adding common nodes among processors
  communicateInterp(nxn, nyn, nzn, rhon, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jxh, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jyh, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jzh, 0, 0, 0, 0, 0, 0, vct);
  for (int i=0; i<NE_MASS; i++) { 
    communicateInterp(nxn, nyn, nzn, i, Mxx, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Mxy, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Mxz, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Myx, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Myy, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Myz, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzx, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzy, 2, 2, 2, 2, 2, 2, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzz, 2, 2, 2, 2, 2, 2, vct);
  }

  // put the correct values on ghost cells
  communicateNode_P(nxn, nyn, nzn, rhon, vct);
  communicateNode_P(nxn, nyn, nzn, Jxh, vct);
  communicateNode_P(nxn, nyn, nzn, Jyh, vct);
  communicateNode_P(nxn, nyn, nzn, Jzh, vct);
  for (int i=0; i<NE_MASS; i++) {
    communicateNode_P(nxn, nyn, nzn, Mxx, i, vct);
    communicateNode_P(nxn, nyn, nzn, Mxy, i, vct);
    communicateNode_P(nxn, nyn, nzn, Mxz, i, vct);
    communicateNode_P(nxn, nyn, nzn, Myx, i, vct);
    communicateNode_P(nxn, nyn, nzn, Myy, i, vct);
    communicateNode_P(nxn, nyn, nzn, Myz, i, vct);
    communicateNode_P(nxn, nyn, nzn, Mzx, i, vct);
    communicateNode_P(nxn, nyn, nzn, Mzy, i, vct);
    communicateNode_P(nxn, nyn, nzn, Mzz, i, vct);
  }

}


/* Interpolate the magnetic fields from centers to nodes */
void EMfields3D::centers2nodesB(VirtualTopology3D *vct, Grid *grid, Collective *col) {
  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);
  communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

}


/*! add an amount of charge density to charge density field at node X,Y */
void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        rhons[is][X - i][Y - j][Z - k] += temp * invVOL;
      }
}

/** add an amount of current density to current density - direction X to current density field on the node*/
void EMfields3D::addJxh(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jxh[X - i][Y - j][Z - k] += temp*invVOL;
      }
}
/** add an amount of current density - direction Y to current density field on the node */
void EMfields3D::addJyh(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jyh[X - i][Y - j][Z - k] += temp*invVOL;
      }
}
/** add an amount of current density - direction Z to current density field on the node */
void EMfields3D::addJzh(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jzh[X - i][Y - j][Z - k] += temp*invVOL;
      }
}

/** add an amount of current density to mass matrix field at node X,Y */
void EMfields3D::addMxx(double val, int X, int Y, int Z, int ind) {
  Mxx[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMxy(double val, int X, int Y, int Z, int ind) {
  Mxy[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMxz(double val, int X, int Y, int Z, int ind) {
  Mxz[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMyx(double val, int X, int Y, int Z, int ind) {
  Myx[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMyy(double val, int X, int Y, int Z, int ind) {
  Myy[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMyz(double val, int X, int Y, int Z, int ind) {
  Myz[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMzx(double val, int X, int Y, int Z, int ind) {
  Mzx[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMzy(double val, int X, int Y, int Z, int ind) {
  Mzy[ind][X][Y][Z] += val*invVOL;
}
void EMfields3D::addMzz(double val, int X, int Y, int Z, int ind) {
  Mzz[ind][X][Y][Z] += val*invVOL;
}

/*! set to 0 all the densities fields */
void EMfields3D::setZeroDensities() {
  for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
        Jxh [i][j][k] = 0.0;
        Jyh [i][j][k] = 0.0;
        Jzh [i][j][k] = 0.0;
        rhon[i][j][k] = 0.0;
      }

  for (register int c = 0; c < NE_MASS; c++) 
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
         Mxx[c][i][j][k] = 0;
         Mxy[c][i][j][k] = 0;
         Mxz[c][i][j][k] = 0;
         Myx[c][i][j][k] = 0;
         Myy[c][i][j][k] = 0;
         Myz[c][i][j][k] = 0;
         Mzx[c][i][j][k] = 0;
         Mzy[c][i][j][k] = 0;
         Mzz[c][i][j][k] = 0;
        }

  for (register int i = 0; i < nxc; i++)
    for (register int j = 0; j < nyc; j++)
      for (register int k = 0; k < nzc; k++) {
        rhoc[i][j][k] = 0.0;
        rhoh[i][j][k] = 0.0;
      }
  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          rhons[kk][i][j][k] = 0.0;
        }

}
/*!SPECIES: Sum the charge density of different species on NODES */
void EMfields3D::sumOverSpecies(VirtualTopology3D * vct) {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++)
          rhon[i][j][k] += rhons[is][i][j][k];
}


/*! initialize Magnetic and Electric Field with initial configuration */
void EMfields3D::init(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  if (restart1 == 0) {

    eqValue(0.0, Exth, nxn, nyn, nzn);
    eqValue(0.0, Eyth, nxn, nyn, nzn);
    eqValue(0.0, Ezth, nxn, nyn, nzn);
    eqValue(0.0, Ex, nxn, nyn, nzn);
    eqValue(0.0, Ey, nxn, nyn, nzn);
    eqValue(0.0, Ez, nxn, nyn, nzn);
    eqValue(0.0, tempXC, nxc, nyc, nzc);
    eqValue(0.0, tempYC, nxc, nyc, nzc);
    eqValue(0.0, tempZC, nxc, nyc, nzc);
    eqValue(0.0, Bxc, nxc, nyc, nzc);
    eqValue(0.0, Byc, nxc, nyc, nzc);
    eqValue(0.0, Bzc, nxc, nyc, nzc);
    for (int i = 0; i < nxn; i++) {
      for (int j = 0; j < nyn; j++) {
        for (int k = 0; k < nzn; k++) {
          for (int is = 0; is < ns; is++) {
            rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          Bxn[i][j][k] = B0x;
          Byn[i][j][k] = B0y;
          Bzn[i][j][k] = B0z;
        }
      }
    }

    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {                        // READING FROM RESTART
    if (vct->getCartesian_rank() == 0)
      cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[3];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << name_file << endl;
      cout << "RESTART NOT POSSIBLE" << endl;
    }

    dataset_id = H5Dopen2(file_id, "/fields/Bx/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);



    // Bxn
    double *temp_storage = new double[dims_out[0] * dims_out[1] * dims_out[2]];
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    int k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bxn[i][j][jj] = temp_storage[k++];


    status = H5Dclose(dataset_id);

    // Byn
    dataset_id = H5Dopen2(file_id, "/fields/By/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Byn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Bzn
    dataset_id = H5Dopen2(file_id, "/fields/Bz/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Bzn[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ex
    dataset_id = H5Dopen2(file_id, "/fields/Ex/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ex[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);


    // Ey 
    dataset_id = H5Dopen2(file_id, "/fields/Ey/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ey[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // Ez 
    dataset_id = H5Dopen2(file_id, "/fields/Ez/cycle_0", H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
    k = 0;
    for (int i = 1; i < nxn - 1; i++)
      for (int j = 1; j < nyn - 1; j++)
        for (int jj = 1; jj < nzn - 1; jj++)
          Ez[i][j][jj] = temp_storage[k++];

    status = H5Dclose(dataset_id);

    // open the charge density for species

    stringstream *species_name = new stringstream[ns];
    for (int is = 0; is < ns; is++) {
      species_name[is] << is;
      string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
      dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
      status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storage);
      k = 0;
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int jj = 1; jj < nzn - 1; jj++)
            rhons[is][i][j][jj] = temp_storage[k++];

      communicateNode_P(nxn, nyn, nzn, rhons, is, vct);
      status = H5Dclose(dataset_id);

    }

    if (col->getCase()=="Dipole") {
      cout << "Dipole case not implemented" << endl;
    }

    // ConstantChargeOpenBC(grid, vct);

    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // communicate E
    communicateNodeBC(nxn, nyn, nzn, Ex, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ey, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Ez, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
    // close the hdf file
    status = H5Fclose(file_id);
    delete[]temp_storage;
    delete[]species_name;
  }
}

/*! initiliaze EM for GEM challange */
void EMfields3D::initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  double pertX = 0.4;
  double xpert, ypert, exp_pert;
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize GEM Challenge with Pertubation" << endl;
      cout << "------------------------------------------" << endl;
      cout << "B0x                              = " << B0x << endl;
      cout << "B0y                              = " << B0y << endl;
      cout << "B0z                              = " << B0z << endl;
      cout << "Delta (current sheet thickness) = " << delta << endl;
      for (int i = 0; i < ns; i++) {
        cout << "rho species " << i << " = " << rhoINIT[i];
        if (DriftSpecies[i])
          cout << " DRIFTING " << endl;
        else
          cout << " BACKGROUND " << endl;
      }
      cout << "-------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is])
              rhons[is][i][j][k] = ((rhoINIT[is] / (cosh((grid->getYN(i, j, k) - Ly / 2) / delta) * cosh((grid->getYN(i, j, k) - Ly / 2) / delta)))) / FourPI;
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * tanh((grid->getYN(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly );
          Byn[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXN(i, j, k) - Lx / 2;
          ypert = grid->getYN(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxn[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byn[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzn[i][j][k] = B0z;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = B0x * tanh((grid->getYC(i, j, k) - Ly / 2) / delta);
          // add the initial GEM perturbation
          // Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly );
          Byc[i][j][k] = B0y;   // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
          // add the initial X perturbation
          xpert = grid->getXC(i, j, k) - Lx / 2;
          ypert = grid->getYC(i, j, k) - Ly / 2;
          exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));
          Bxc[i][j][k] += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * ypert / delta - cos(M_PI * xpert / 10.0 / delta) * sin(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          Byc[i][j][k] += (B0x * pertX) * exp_pert * (cos(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * 2.0 * xpert / delta + sin(M_PI * xpert / 10.0 / delta) * cos(M_PI * ypert / 10.0 / delta) * M_PI / 10.0);
          // guide field
          Bzc[i][j][k] = B0z;

        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}



/*! Calculate the susceptibility on the boundary leftX */
void EMfields3D::sustensorLeftX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[1][j][k] + Fext*Bx_ext[1][j][k]);
        omcy = beta * (Byn[1][j][k] + Fext*By_ext[1][j][k]);
        omcz = beta * (Bzn[1][j][k] + Fext*Bz_ext[1][j][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][1][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary rightX */
void EMfields3D::sustensorRightX(double **susxx, double **susyx, double **suszx) {
  double beta, omcx, omcy, omcz, denom;
  for (int j = 0; j < nyn; j++)
    for (int k = 0; k < nzn; k++) {
      susxx[j][k] = 1.0;
      susyx[j][k] = 0.0;
      suszx[j][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[nxn - 2][j][k] + Fext*Bx_ext[nxn - 2][j][k]);
        omcy = beta * (Byn[nxn - 2][j][k] + Fext*By_ext[nxn - 2][j][k]);
        omcz = beta * (Bzn[nxn - 2][j][k] + Fext*Bz_ext[nxn - 2][j][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][nxn - 2][j][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxx[j][k] += (  1.0 + omcx * omcx) * denom;
        susyx[j][k] += (-omcz + omcx * omcy) * denom;
        suszx[j][k] += ( omcy + omcx * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[i][1][k] + Fext*Bx_ext[i][1][k]);
        omcy = beta * (Byn[i][1][k] + Fext*By_ext[i][1][k]);
        omcz = beta * (Bzn[i][1][k] + Fext*Bz_ext[i][1][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][1][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightY(double **susxy, double **susyy, double **suszy) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int k = 0; k < nzn; k++) {
      susxy[i][k] = 0.0;
      susyy[i][k] = 1.0;
      suszy[i][k] = 0.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int k = 0; k < nzn; k++) {
        omcx = beta * (Bxn[i][nyn - 2][k] + Fext*Bx_ext[i][nyn - 2][k]);
        omcy = beta * (Byn[i][nyn - 2][k] + Fext*By_ext[i][nyn - 2][k]);
        omcz = beta * (Bzn[i][nyn - 2][k] + Fext*Bz_ext[i][nyn - 2][k]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][nyn - 2][k] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxy[i][k] += ( omcz + omcx * omcy) * denom;
        susyy[i][k] += (  1.0 + omcy * omcy) * denom;
        suszy[i][k] += (-omcx + omcy * omcz) * denom;
      }
  }
}

/*! Calculate the susceptibility on the boundary left */
void EMfields3D::sustensorLeftZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * (Bxn[i][j][1] + Fext*Bx_ext[i][j][1]);
        omcy = beta * (Byn[i][j][1] + Fext*By_ext[i][j][1]);
        omcz = beta * (Bzn[i][j][1] + Fext*Bz_ext[i][j][1]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][1] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }

}
/*! Calculate the susceptibility on the boundary right */
void EMfields3D::sustensorRightZ(double **susxz, double **susyz, double **suszz) {
  double beta, omcx, omcy, omcz, denom;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      susxz[i][j] = 0.0;
      susyz[i][j] = 0.0;
      suszz[i][j] = 1.0;
    }
  for (int is = 0; is < ns; is++) {
    beta = .5 * qom[is] * dt / c;
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++) {
        omcx = beta * (Bxn[i][j][nzn - 2] + Fext*Bx_ext[i][j][nzn - 2]);
        omcy = beta * (Byn[i][j][nzn - 2] + Fext*By_ext[i][j][nzn - 2]);
        omcz = beta * (Bzn[i][j][nzn - 2] + Fext*Bz_ext[i][j][nzn - 2]);
        denom = FourPI / 2 * delt * dt / c * qom[is] * rhons[is][i][j][nyn - 2] / (1.0 + omcx * omcx + omcy * omcy + omcz * omcz);
        susxz[i][j] += (-omcy + omcx * omcz) * denom;
        susyz[i][j] += ( omcx + omcy * omcz) * denom;
        suszz[i][j] += (  1.0 + omcz * omcz) * denom;
      }
  }
}

/*! Perfect conductor boundary conditions: LEFT wall */
void EMfields3D::perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0:  // boundary condition on X-DIRECTION 
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorLeftX(susxx, susyx, suszx);
      for (int i=1; i <  nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[1][i][j] = vectorX[1][i][j] - (Ex[1][i][j] - susyx[i][j]*vectorY[1][i][j] - suszx[i][j]*vectorZ[1][i][j] - Jxh[1][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[1][i][j] = vectorY[1][i][j] - 0.0*vectorY[2][i][j];
          imageZ[1][i][j] = vectorZ[1][i][j] - 0.0*vectorZ[2][i][j];
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorLeftY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[i][1][j] = vectorX[i][1][j] - 0.0*vectorX[i][2][j];
          imageY[i][1][j] = vectorY[i][1][j] - (Ey[i][1][j] - susxy[i][j]*vectorX[i][1][j] - suszy[i][j]*vectorZ[i][1][j] - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][1][j] = vectorZ[i][1][j] - 0.0*vectorZ[i][2][j];
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorLeftZ(susxz, susyz, suszz);
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          imageX[i][j][1] = vectorX[i][j][1];
          imageY[i][j][1] = vectorX[i][j][1];
          imageZ[i][j][1] = vectorZ[i][j][1] - (Ez[i][j][1] - susxz[i][j]*vectorX[i][j][1] - susyz[i][j]*vectorY[i][j][1] - Jzh[i][j][1]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions: RIGHT wall */
void EMfields3D::perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid) {
  double beta, omcx, omcy, omcz, denom;
  double** susxy;
  double** susyy;
  double** suszy;
  double** susxx;
  double** susyx;
  double** suszx;
  double** susxz;
  double** susyz;
  double** suszz;
  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      susxx = newArr2(double,nyn,nzn);
      susyx = newArr2(double,nyn,nzn);
      suszx = newArr2(double,nyn,nzn);
      sustensorRightX(susxx, susyx, suszx);
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j <  nzn-1;j++){
          imageX[nxn-2][i][j] = vectorX[nxn-2][i][j] - (Ex[nxn-2][i][j] - susyx[i][j]*vectorY[nxn-2][i][j] - suszx[i][j]*vectorZ[nxn-2][i][j] - Jxh[nxn-2][i][j]*dt*th*FourPI)/susxx[i][j];
          imageY[nxn-2][i][j] = vectorY[nxn-2][i][j] - 0.0 * vectorY[nxn-3][i][j];
          imageZ[nxn-2][i][j] = vectorZ[nxn-2][i][j] - 0.0 * vectorZ[nxn-3][i][j];
        }
      delArr2(susxx,nxn);
      delArr2(susyx,nxn);       
      delArr2(suszx,nxn);
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      susxy = newArr2(double,nxn,nzn);
      susyy = newArr2(double,nxn,nzn);
      suszy = newArr2(double,nxn,nzn);
      sustensorRightY(susxy, susyy, suszy);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          imageX[i][nyn-2][j] = vectorX[i][nyn-2][j] - 0.0*vectorX[i][nyn-3][j];
          imageY[i][nyn-2][j] = vectorY[i][nyn-2][j] - (Ey[i][nyn-2][j] - susxy[i][j]*vectorX[i][nyn-2][j] - suszy[i][j]*vectorZ[i][nyn-2][j] - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
          imageZ[i][nyn-2][j] = vectorZ[i][nyn-2][j] - 0.0*vectorZ[i][nyn-3][j];
        }
      delArr2(susxy,nxn);
      delArr2(susyy,nxn);
      delArr2(suszy,nxn);
      break;
    case 2: // boundary condition on Z-DIRECTION RIGHT
      susxz = newArr2(double,nxn,nyn);
      susyz = newArr2(double,nxn,nyn);
      suszz = newArr2(double,nxn,nyn);
      sustensorRightZ(susxz, susyz, suszz);
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nyn-1;j++){
          imageX[i][j][nzn-2] = vectorX[i][j][nzn-2];
          imageY[i][j][nzn-2] = vectorY[i][j][nzn-2];
          imageZ[i][j][nzn-2] = vectorZ[i][j][nzn-2] - (Ez[i][j][nzn-2] - susxz[i][j]*vectorX[i][j][nzn-2] - susyz[i][j]*vectorY[i][j][nzn-2] - Jzh[i][j][nzn-2]*dt*th*FourPI)/suszz[i][j];
        }
      delArr2(susxz,nxn);
      delArr2(susyz,nxn);       
      delArr2(suszz,nxn);
      break;
  }
}

/*! Perfect conductor boundary conditions for source: LEFT WALL */
void EMfields3D::perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION LEFT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[1][i][j] = 0.0;
          vectorY[1][i][j] = ebc[1];
          vectorZ[1][i][j] = ebc[2];
          //+//          vectorX[1][i][j] = 0.0;
          //+//          vectorY[1][i][j] = 0.0;
          //+//          vectorZ[1][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][1][j] = ebc[0];
          vectorY[i][1][j] = 0.0;
          vectorZ[i][1][j] = ebc[2];
          //+//          vectorX[i][1][j] = 0.0;
          //+//          vectorY[i][1][j] = 0.0;
          //+//          vectorZ[i][1][j] = 0.0;
        }
      break;
    case 2: // boundary condition on Z-DIRECTION LEFT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][1] = ebc[0];
          vectorY[i][j][1] = ebc[1];
          vectorZ[i][j][1] = 0.0;
          //+//          vectorX[i][j][1] = 0.0;
          //+//          vectorY[i][j][1] = 0.0;
          //+//          vectorZ[i][j][1] = 0.0;
        }
      break;
  }
}

/*! Perfect conductor boundary conditions for source: RIGHT WALL */
void EMfields3D::perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir) {

  double ebc[3];

  // Assuming E = - ve x B
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
      for (int i=1; i < nyn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[nxn-2][i][j] = 0.0;
          vectorY[nxn-2][i][j] = ebc[1];
          vectorZ[nxn-2][i][j] = ebc[2];
          //+//          vectorX[nxn-2][i][j] = 0.0;
          //+//          vectorY[nxn-2][i][j] = 0.0;
          //+//          vectorZ[nxn-2][i][j] = 0.0;
        }
      break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
      for (int i=1; i < nxn-1;i++)
        for (int j=1; j < nzn-1;j++){
          vectorX[i][nyn-2][j] = ebc[0];
          vectorY[i][nyn-2][j] = 0.0;
          vectorZ[i][nyn-2][j] = ebc[2];
          //+//          vectorX[i][nyn-2][j] = 0.0;
          //+//          vectorY[i][nyn-2][j] = 0.0;
          //+//          vectorZ[i][nyn-2][j] = 0.0;
        }
      break;
    case 2:
      for (int i=1; i <  nxn-1;i++)
        for (int j=1; j <  nyn-1;j++){
          vectorX[i][j][nzn-2] = ebc[0];
          vectorY[i][j][nzn-2] = ebc[1];
          vectorZ[i][j][nzn-2] = 0.0;
          //+//          vectorX[i][j][nzn-2] = 0.0;
          //+//          vectorY[i][j][nzn-2] = 0.0;
          //+//          vectorZ[i][j][nzn-2] = 0.0;
        }
      break;
  }
}


void EMfields3D::UpdateRHOcs(Grid * grid){

  double r = 1.0;

  double xmin = 0.0;
  double xmax = Lx;
  double rmin = 1.0;
  double rmax = 1.0;

  for (int is = 0; is < ns; is++)
    for (int i=0; i<nxc-1; i++)
      for (int j=0; j<nyc-1; j++)
        for (int k=0; k<nzc-1; k++){
          double x = grid->getXN(i,j,k);
          r = rmin + (rmax-rmin) * (xmax - x) / (xmax - xmin);
          if (r<xmin) r = rmax;
          if (r>xmax) r = rmin;
          rhocs[is][i][j][k] = (qom[is]/fabs(qom[is]))*(r/FourPI);
        }
    //grid->interpN2C(rhocs, is, rhons);
}

/*! get Fext */
double EMfields3D::getFext(){
  return (Fext);
}

/*! get Ex(X,Y,Z) */
double &EMfields3D::getEx(int indexX, int indexY, int indexZ) const {
  return (Ex[indexX][indexY][indexZ]);
}
/*! get Electric field component X array */
double ***EMfields3D::getEx() {
  return (Ex);
}
double ***EMfields3D::getExth() {
  return (Exth);
}
/*! get Electric Field component X array cell without the ghost cells */
double ***EMfields3D::getExc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
        arr[i][j][k] = .125 * (Ex[i][j][k]     + 
                               Ex[i + 1][j][k] + 
                               Ex[i][j + 1][k] + 
                               Ex[i][j][k + 1] + 
                               Ex[i + 1][j + 1][k] + 
                               Ex[i + 1][j][k + 1] + 
                               Ex[i][j + 1][k + 1] + 
                               Ex[i + 1][j + 1][k + 1]);

  return arr;
}
/*! get Ey(X,Y,Z) */
double &EMfields3D::getEy(int indexX, int indexY, int indexZ) const {
  return (Ey[indexX][indexY][indexZ]);
}
/*! get Electric field component Y array */
double ***EMfields3D::getEy() {
  return (Ey);
}
double ***EMfields3D::getEyth() {
  return (Eyth);
}
/*! get Electric Field component Y array cell without the ghost cells */
double ***EMfields3D::getEyc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
        arr[i][j][k] = .125 * (Ey[i][j][k]     + 
                               Ey[i + 1][j][k] + 
                               Ey[i][j + 1][k] + 
                               Ey[i][j][k + 1] + 
                               Ey[i + 1][j + 1][k] + 
                               Ey[i + 1][j][k + 1] + 
                               Ey[i][j + 1][k + 1] + 
                               Ey[i + 1][j + 1][k + 1]);
  return arr;
}
/*! get Ez(X,Y,Z) */
double &EMfields3D::getEz(int indexX, int indexY, int indexZ) const {
  return (Ez[indexX][indexY][indexZ]);
}
/*! get Electric field component Z array */
double ***EMfields3D::getEz() {
  return (Ez);
}
double ***EMfields3D::getEzth() {
  return (Ezth);
}
/*! get Electric Field component Z array cell without the ghost cells */
double ***EMfields3D::getEzc() {

  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
        arr[i][j][k] = .125 * (Ez[i][j][k]     +
                               Ez[i + 1][j][k] +
                               Ez[i][j + 1][k] +
                               Ez[i][j][k + 1] +
                               Ez[i + 1][j + 1][k] +
                               Ez[i + 1][j][k + 1] +
                               Ez[i][j + 1][k + 1] +
                               Ez[i + 1][j + 1][k + 1]);

  return arr;
}
/*! get Bx(X,Y,Z) */
double &EMfields3D::getBx(int indexX, int indexY, int indexZ) const {
  return (Bxn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component X array */
double ***EMfields3D::getBx() {
  return (Bxn);
}
/*! get Magnetic Field component X array cell without the ghost cells */
double EMfields3D::getBxc(int i, int j, int k) {
  return (Bxc[i][j][k]);
}
/*! get Magnetic Field component X array cell without the ghost cells */
double ***EMfields3D::getBxc() {
  return (Bxc);
}
/*! get By(X,Y,Z) */
double &EMfields3D::getBy(int indexX, int indexY, int indexZ) const {
  return (Byn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Y array */
double ***EMfields3D::getBy() {
  return (Byn);
}
/*! get Magnetic Field component Y array cell without the ghost cells */
double EMfields3D::getByc(int i, int j, int k) {
  return (Byc[i][j][k]);
}
/*! get Magnetic Field component Y array cell without the ghost cells */
double ***EMfields3D::getByc() {
  return (Byc);
}
/*! get Bz(X,Y,Z) */
double &EMfields3D::getBz(int indexX, int indexY, int indexZ) const {
  return (Bzn[indexX][indexY][indexZ]);
}
/*! get Magnetic Field component Z array */
double ***EMfields3D::getBz() {
  return (Bzn);
}
/*! get Magnetic Field component Z array cell without the ghost cells */
double EMfields3D::getBzc(int i, int j, int k) {
  return (Bzc[i][j][k]);
}
/*! get Magnetic Field component Z array cell without the ghost cells */
double ***EMfields3D::getBzc() {
  return (Bzc);
}

/*! get rhoc(X,Y,Z) */
double &EMfields3D::getRHOc(int indexX, int indexY, int indexZ) const {
  return (rhoc[indexX][indexY][indexZ]);
} double ***EMfields3D::getRHOc() {
  return (rhoc);
}
/*! get density on node(indexX,indexY,indexZ) */
double &EMfields3D::getRHOn(int indexX, int indexY, int indexZ) const {
  return (rhon[indexX][indexY][indexZ]);
}
/*! get density array defined on nodes */
double ***EMfields3D::getRHOn() {
  return (rhon);
}
/*! get rhos(X,Y,Z) : density for species */
double &EMfields3D::getRHOns(int indexX, int indexY, int indexZ, int is) const {
  return (rhons[is][indexX][indexY][indexZ]);
}
/*! SPECIES: get density array defined on center cells */
double &EMfields3D::getRHOcs(int indexX, int indexY, int indexZ, int is) const {
  return (rhocs[is][indexX][indexY][indexZ]);
}
/*! SPECIES: get density array defined on center cells */
double ****& EMfields3D::getRHOcs() {
  return (rhocs);
}
/*! get density array defined on nodes */
double ****EMfields3D::getRHOns() {
  return (rhons);
}
/*! get density array defined on nodes for one species */
double ***& EMfields3D::getRHOns(int is) {
  return (rhons[is]);
}
/*! get species density component X array cell without the ghost cells */
double ***EMfields3D::getRHOcs(int is) {
  for (int i = 0; i < nxc; i++)
    for (int j = 0; j < nyc; j++)
      for (int k = 0; k < nzc; k++)
        arr[i][j][k] = .125 * (rhons[is][i][j][k]     +
                               rhons[is][i + 1][j][k] +
                               rhons[is][i][j + 1][k] +
                               rhons[is][i][j][k + 1] +
                               rhons[is][i + 1][j + 1][k] +
                               rhons[is][i + 1][j][k + 1] +
                               rhons[is][i][j + 1][k + 1] +
                               rhons[is][i + 1][j + 1][k + 1]);

  return arr;
}

double ***& EMfields3D::getRHOcs(int is, int dummy) {
  return (rhocs[is]);
}

/*! get Bx_ext(X,Y,Z)  */
double &EMfields3D::getBx_ext(int indexX, int indexY, int indexZ) const{
  return(Bx_ext[indexX][indexY][indexZ]);
}
/*!  get By_ext(X,Y,Z) */
double &EMfields3D::getBy_ext(int indexX, int indexY, int indexZ) const{
  return(By_ext[indexX][indexY][indexZ]);
}
/*!  get Bz_ext(X,Y,Z) */
double &EMfields3D::getBz_ext(int indexX, int indexY, int indexZ) const{
  return(Bz_ext[indexX][indexY][indexZ]);
}

/*! get Bx_ext  */
double ***EMfields3D::getBx_ext() {
  return(Bx_ext);
}
/*!  get By_ext */
double ***EMfields3D::getBy_ext() {
  return(By_ext);
}
/*!  get Bz_ext */
double ***EMfields3D::getBz_ext() {
  return(Bz_ext);
}

double ***&EMfields3D::getBxTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
        arr[i][j][k] = Bxn[i][j][k] + Fext * Bx_ext[i][j][k];

  return arr;
}
double ***&EMfields3D::getByTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
        arr[i][j][k] = Byn[i][j][k] + Fext * By_ext[i][j][k];

  return arr;
}
double ***&EMfields3D::getBzTot(){
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++)
        arr[i][j][k] = Bzn[i][j][k] + Fext * Bz_ext[i][j][k];

  return arr;
}

/*! get current -Direction X */
double &EMfields3D::getJxh(int indexX, int indexY, int indexZ) const {
  return (Jxh[indexX][indexY][indexZ]);
}
/*! get current array X component * */
double ***EMfields3D::getJxh() {
  return (Jxh);
}
/*! get current -Direction Y */
double &EMfields3D::getJyh(int indexX, int indexY, int indexZ) const {
  return (Jyh[indexX][indexY][indexZ]);
}
/*! get current array Y component * */
double ***EMfields3D::getJyh() {
  return (Jyh);
}
/*! get current -Direction Z */
double &EMfields3D::getJzh(int indexX, int indexY, int indexZ) const {
  return (Jzh[indexX][indexY][indexZ]);
}
/*! get current array Z component * */
double ***EMfields3D::getJzh() {
  return (Jzh);
}

/*! get the electric field energy */
double EMfields3D::getEenergy(void) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;

  double err; 
  for (int i = 1; i < nxn-1 ; i++) {
  }


  for (int i = 1; i < nxn-2 ; i++)
    for (int j = 1; j < nyn-2 ; j++)
      for (int k = 1; k < nzn-2 ; k++) {
        localEenergy += .5 * dx * dy * dz * (Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI);
      }

  MPI_Allreduce(&localEenergy, &totalEenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalEenergy);

}
/*! get the magnetic field energy */
double EMfields3D::getBenergy(void) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;
  for (int i = 1; i < nxc-1; i++)
    for (int j = 1; j < nyc-1; j++)
      for (int k = 1; k < nzc-1; k++){
        Bxt = Bxc[i][j][k];
        Byt = Byc[i][j][k];
        Bzt = Bzc[i][j][k];
        localBenergy += .5*dx*dy*dz*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI);
      }

  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalBenergy);
}


/*! Print info about electromagnetic field */
void EMfields3D::print(void) const {
}

/*! destructor: deallocate arrays */
EMfields3D::~EMfields3D() {
  // nodes
  delArr3(Ex, nxn, nyn);
  delArr3(Ey, nxn, nyn);
  delArr3(Ez, nxn, nyn);
  delArr3(Exth, nxn, nyn);
  delArr3(Eyth, nxn, nyn);
  delArr3(Ezth, nxn, nyn);
  delArr3(Bxn, nxn, nyn);
  delArr3(Byn, nxn, nyn);
  delArr3(Bzn, nxn, nyn);
  delArr3(rhon, nxn, nyn);
  delArr3(Jxh, nxn, nyn);
  delArr3(Jyh, nxn, nyn);
  delArr3(Jzh, nxn, nyn);
  delArr4(Mxx, NE_MASS, nxn, nyn);
  delArr4(Mxy, NE_MASS, nxn, nyn);
  delArr4(Mxz, NE_MASS, nxn, nyn);
  delArr4(Myx, NE_MASS, nxn, nyn);
  delArr4(Myy, NE_MASS, nxn, nyn);
  delArr4(Myz, NE_MASS, nxn, nyn);
  delArr4(Mzx, NE_MASS, nxn, nyn);
  delArr4(Mzy, NE_MASS, nxn, nyn);
  delArr4(Mzz, NE_MASS, nxn, nyn);
  // nodes and species
  delArr4(rhons, ns, nxn, nyn);
  // central points
  delArr3(Bxc, nxc, nyc);
  delArr3(Byc, nxc, nyc);
  delArr3(Bzc, nxc, nyc);
  delArr3(rhoc, nxc, nyc);
  delArr3(rhoh, nxc, nyc);
  // various stuff needs to be deallocated too
  delArr3(tempXC, nxc, nyc);
  delArr3(tempYC, nxc, nyc);
  delArr3(tempZC, nxc, nyc);
  delArr3(tempXN, nxn, nyn);
  delArr3(tempYN, nxn, nyn);
  delArr3(tempZN, nxn, nyn);
  delArr3(tempC, nxc, nyc);
  delArr3(tempX, nxn, nyn);
  delArr3(tempY, nxn, nyn);
  delArr3(tempZ, nxn, nyn);
  delArr3(temp2X, nxn, nyn);
  delArr3(temp2Y, nxn, nyn);
  delArr3(temp2Z, nxn, nyn);
  delArr3(imageEX, nxn, nyn);
  delArr3(imageEY, nxn, nyn);
  delArr3(imageEZ, nxn, nyn);
  delArr3(imageBX, nxn, nyn);
  delArr3(imageBY, nxn, nyn);
  delArr3(imageBZ, nxn, nyn);
  delArr3(vectX, nxn, nyn);
  delArr3(vectY, nxn, nyn);
  delArr3(vectZ, nxn, nyn);
}
