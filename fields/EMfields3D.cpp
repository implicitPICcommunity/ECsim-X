
//#define ELECTROSTATIC
#include "EMfields3D.h"
#include "MyClock.h"
//#include "MaxwellECIMM.h"
extern MyClock *clocks;


long double factorial(int k) {
  if (k == 0) return 1;
  else if (k == 1) return k;
  else return k*factorial(k-1);
}

double besselJ(int alpha, double x) {
  double val = 0;
  int iter_max = 1000;
  double TOL = 1e-16;
  int sign = -1;

  double err;
  for (int k=0; k<iter_max; k++) {
    sign *= -1;
    double v = sign*pow(x/2,2*k+alpha)/(factorial(k)*factorial(alpha+k));
    val += v;
    err = fabs(v);
    if (err < TOL) {
      return val;
    }
  }
  return val;
}



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
  x_center = col->getx_center();
  y_center = col->gety_center();
  z_center = col->getz_center();
  L_square = col->getL_square();
  coilD = col->getcoilD();
  coilSpacing = col->getcoilSpacing();
  Cylindrical = col->getCylindrical();

  invVOLn = grid->getInvVOLn();
  invVOLc = grid->getInvVOLc();

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

  Fext = 0;

  u_bulk = col->getUbulk();
  v_bulk = col->getVbulk();
  w_bulk = col->getWbulk();

  // arrays allocation: nodes
  Ex = newArr3(double, nxn, nyn, nzn);
  Ey = newArr3(double, nxn, nyn, nzn);
  Ez = newArr3(double, nxn, nyn, nzn);
  //Ex_prev= newArr3(double, nxn, nyn, nzn);
  //Ey_prev= newArr3(double, nxn, nyn, nzn);
  //Ez_prev= newArr3(double, nxn, nyn, nzn);
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
  Jxhs = newArr4(double, ns, nxn, nyn, nzn);
  Jyhs = newArr4(double, ns, nxn, nyn, nzn);
  Jzhs = newArr4(double, ns, nxn, nyn, nzn);
  Jzhs = newArr4(double, ns, nxn, nyn, nzn);

  rhotilde = newArr4(double, ns, nxn, nyn, nzn);

  divE = newArr3(double, nxn, nyn, nzn);

  Jxs = newArr4(double, ns, nxn, nyn, nzn);
  Jys = newArr4(double, ns, nxn, nyn, nzn);
  Jzs = newArr4(double, ns, nxn, nyn, nzn);


  Mxx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mxy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mxz = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Myz = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzx = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzy = newArr4(double, NE_MASS, nxn, nyn, nzn);
  Mzz = newArr4(double, NE_MASS, nxn, nyn, nzn);

  pXXsn = newArr4(double, ns, nxn, nyn, nzn);
  pXYsn = newArr4(double, ns, nxn, nyn, nzn);
  pXZsn = newArr4(double, ns, nxn, nyn, nzn);
  pYYsn = newArr4(double, ns, nxn, nyn, nzn);
  pYZsn = newArr4(double, ns, nxn, nyn, nzn);
  pZZsn = newArr4(double, ns, nxn, nyn, nzn);


  // External imposed fields
  Bx_ext = newArr3(double,nxn,nyn,nzn);
  By_ext = newArr3(double,nxn,nyn,nzn);
  Bz_ext = newArr3(double,nxn,nyn,nzn);
  Ex_ext = newArr3(double,nxn,nyn,nzn);
  Ey_ext = newArr3(double,nxn,nyn,nzn);
  Ez_ext = newArr3(double,nxn,nyn,nzn);
  // involving species
  rhons = newArr4(double, ns, nxn, nyn, nzn);
  rhocs = newArr4(double, ns, nxc, nyc, nzc);
  // arrays allocation: central points 
  Bxc = newArr3(double, nxc, nyc, nzc);
  Byc = newArr3(double, nxc, nyc, nzc);
  Bzc = newArr3(double, nxc, nyc, nzc);
  rhoc = newArr3(double, nxc, nyc, nzc);

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

  zeroCurrent = (col->getZeroCurrent() == 1 ? 1 : 0);
  Jx_ext = newArr3(double, nxn, nyn, nzn);
  Jy_ext = newArr3(double, nxn, nyn, nzn);
  Jz_ext = newArr3(double, nxn, nyn, nzn);

  if (col->getCase() == "Dipole" ||
      col->getCase() == "TwoCoils" ||
      col->getCase() == "TwoCoilsINJ" ||
      col->getCase() == "FluxRope" ||
      col->getCase() == "GEM" ) {
    Lambda = newArr3(double, nxn, nyn, nzn);
    damping = 1;
  } else damping = 0;

  // Set to zero all the memory allocated
  setAllzero();

  #ifdef _TEST_CURL_
    testCurl(grid);
    exit(0);
  #endif

  #ifdef _TEST_DIV_
    testDivergence(grid);
    exit(0);
  #endif

}

void EMfields3D::setAllzero()
{
   Fext = 0;
   eqValue(0.0, Ex, nxn, nyn, nzn);
   eqValue(0.0, Ey, nxn, nyn, nzn);
   eqValue(0.0, Ez, nxn, nyn, nzn);
   eqValue(0.0, Exth, nxn, nyn, nzn);
   eqValue(0.0, Eyth, nxn, nyn, nzn);
   eqValue(0.0, Ezth, nxn, nyn, nzn);
   eqValue(0.0, Bxn, nxn, nyn, nzn);
   eqValue(0.0, Byn, nxn, nyn, nzn);
   eqValue(0.0, Bzn, nxn, nyn, nzn);
   eqValue(0.0, rhon, nxn, nyn, nzn);
   eqValue(0.0, Jxh, nxn, nyn, nzn);
   eqValue(0.0, Jyh, nxn, nyn, nzn);
   eqValue(0.0, Jzh, nxn, nyn, nzn);
   eqValue(0.0, Jxhs, ns, nxn, nyn, nzn);
   eqValue(0.0, Jyhs, ns, nxn, nyn, nzn);
   eqValue(0.0, Jzhs, ns, nxn, nyn, nzn);

   eqValue(0.0, divE, nxn, nyn, nzn);
   eqValue(0.0, rhotilde, ns, nxn, nyn, nzn);
   eqValue(0.0, Jxs, ns, nxn, nyn, nzn);
   eqValue(0.0, Jys, ns, nxn, nyn, nzn);
   eqValue(0.0, Jzs, ns, nxn, nyn, nzn);
   eqValue(0.0, Mxx, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Mxy, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Mxz, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Myx, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Myy, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Myz, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Mzx, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Mzy, NE_MASS, nxn, nyn, nzn);
   eqValue(0.0, Mzz, NE_MASS, nxn, nyn, nzn);

   eqValue(0.0, pXXsn, ns, nxn, nyn, nzn);
   eqValue(0.0, pXYsn, ns, nxn, nyn, nzn);
   eqValue(0.0, pXZsn, ns, nxn, nyn, nzn);
   eqValue(0.0, pYYsn, ns, nxn, nyn, nzn);
   eqValue(0.0, pYZsn, ns, nxn, nyn, nzn);
   eqValue(0.0, pZZsn, ns, nxn, nyn, nzn);

   eqValue(0.0, Bx_ext, nxn, nyn, nzn);
   eqValue(0.0, By_ext, nxn, nyn, nzn);
   eqValue(0.0, Bz_ext, nxn, nyn, nzn);
   eqValue(0.0, Ex_ext, nxn, nyn, nzn);
   eqValue(0.0, Ey_ext, nxn, nyn, nzn);
   eqValue(0.0, Ez_ext, nxn, nyn, nzn);

   eqValue(0.0,rhons, ns, nxn, nyn, nzn);
   eqValue(0.0,rhocs, ns, nxc, nyc, nzc);

   eqValue(0.0, Bxc, nxc, nyc, nzc);
   eqValue(0.0, Byc, nxc, nyc, nzc);
   eqValue(0.0, Bzc, nxc, nyc, nzc);
   eqValue(0.0, rhoc, nxc, nyc, nzc);

   eqValue(0.0, tempXC, nxc, nyc, nzc);
   eqValue(0.0, tempYC, nxc, nyc, nzc);
   eqValue(0.0, tempZC, nxc, nyc, nzc);
   eqValue(0.0, tempXN, nxn, nyn, nzn);
   eqValue(0.0, tempYN, nxn, nyn, nzn);
   eqValue(0.0, tempZN, nxn, nyn, nzn);
   eqValue(0.0, tempC, nxc, nyc, nzc);
   eqValue(0.0, tempX, nxn, nyn, nzn);
   eqValue(0.0, tempY, nxn, nyn, nzn);
   eqValue(0.0, tempZ, nxn, nyn, nzn);
   eqValue(0.0, temp2X, nxn, nyn, nzn);
   eqValue(0.0, temp2Y, nxn, nyn, nzn);
   eqValue(0.0, temp2Z, nxn, nyn, nzn);
   eqValue(0.0, imageEX, nxc, nyc, nzc);
   eqValue(0.0, imageEY, nxc, nyc, nzc);
   eqValue(0.0, imageEZ, nxc, nyc, nzc);
   eqValue(0.0, imageBX, nxn, nyn, nzn);
   eqValue(0.0, imageBY, nxn, nyn, nzn);
   eqValue(0.0, imageBZ, nxn, nyn, nzn);
   eqValue(0.0, vectX, nxn, nyn, nzn);
   eqValue(0.0, vectY, nxn, nyn, nzn);
   eqValue(0.0, vectZ, nxn, nyn, nzn);
   eqValue(0.0, arr, nxn, nyn, nzn);

   eqValue(0.0, Jx_ext, nxn, nyn, nzn);
   eqValue(0.0, Jy_ext, nxn, nyn, nzn);
   eqValue(0.0, Jz_ext, nxn, nyn, nzn);

   if (damping == 1) eqValue(0.0, Lambda, nxn, nyn, nzn);
}

void EMfields3D::fixGhostB() {

  if (Cylindrical == 1) {
    if (xStart == 0) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          Bxc[0][j][k] = -Bxc[1][j][k];
          Byc[0][j][k] = Byc[1][j][k];
          Bzc[0][j][k] = -Bzc[1][j][k];
        }
    }
    if (xEnd == Lx) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          Bxc[nxc-1][j][k] = Bxc[nxc-2][j][k];
          Byc[nxc-1][j][k] = Byc[nxc-2][j][k];
          Bzc[nxc-1][j][k] = Bzc[nxc-2][j][k];
        }
    }
  }
}



/* Prepare the memory for E calculation */
void EMfields3D::startEcalc(VirtualTopology3D * vct, Collective *col) {
  if (vct->getCartesian_rank() == 0)
    cout << "*** E CALCULATION ***" << endl;

}

void EMfields3D::endEcalc(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid)
{
  #ifdef __PROFILING__
  clocks->start(14);
  #endif

  //eq(Ex_prev, Ex, nxn, nyn, nzn);
  //eq(Ey_prev, Ey, nxn, nyn, nzn);
  //eq(Ez_prev, Ez, nxn, nyn, nzn);

  // move from krylov space to physical space
  solver2phys(Exth, Eyth, Ezth, tempXC, tempYC, tempZC, xkrylov, nxn, nyn, nzn);
  smoothEth(Smooth, vct, col);
  smoothEth(Smooth, vct, col);
  smoothEth(Smooth, vct, col);
  //smoothBth(Smooth, tempXC, tempYC, tempZC, vct, col);
  //smoothBth(Smooth, tempXC, tempYC, tempZC, vct, col);
  //smoothBth(Smooth, tempXC, tempYC, tempZC, vct, col);

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
  communicateCenterBC(nxc, nyc, nzc, Bxc, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, Byc, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, Bzc, 2, 2, 2, 2, 2, 2, vct);

  #ifdef __PROFILING__
  clocks->stop(14);
  #endif
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
  phys2solver(xkrylov, Ex, Ey, Ez, Bxc, Byc, Bzc, nxn, nyn, nzn);
  // solver
  GMRES(&Field::MaxwellImage, xkrylov, dim, bkrylov, 20, 200, GMREStol, grid, vct, col, this);

  endEcalc(xkrylov, vct, col, grid);

  // deallocate temporary arrays
  delete[]xkrylov;
  delete[]bkrylov;
}

/*! Calculate sorgent for Maxwell solver */
void EMfields3D::MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  #ifdef __PROFILING__
  clocks->start(11);
  #endif

  #ifdef ELECTROSTATIC
  eqValue(0, Bxc, nxc, nyc, nzc);
  eqValue(0, Byc, nxc, nyc, nzc);
  eqValue(0, Bzc, nxc, nyc, nzc);
  #else
  scale(tempXC, Bxc, 1, nxc, nyc, nzc);
  scale(tempYC, Byc, 1, nxc, nyc, nzc);
  scale(tempZC, Bzc, 1, nxc, nyc, nzc);
  #endif

  //smoothJh(Smooth, vct, col);

  for (int i=0; i<nxn-0; i++)
    for (int j=0; j<nyn-0; j++) 
      for (int k=0; k<nzn-0; k++) {
	 double Jx_tot = Jxh[i][j][k] + zeroCurrent*Jx_ext[i][j][k];
	 double Jy_tot = Jyh[i][j][k] + zeroCurrent*Jy_ext[i][j][k];
	 double Jz_tot = Jzh[i][j][k] + zeroCurrent*Jz_ext[i][j][k];
         tempX[i][j][k] = -4*M_PI*th*dt*Jx_tot*grid->getInvVOLn(i,j,k);
         tempY[i][j][k] = -4*M_PI*th*dt*Jy_tot*grid->getInvVOLn(i,j,k);
         tempZ[i][j][k] = -4*M_PI*th*dt*Jz_tot*grid->getInvVOLn(i,j,k);
      }

  addscale(1, tempX, Ex, nxn, nyn, nzn);
  addscale(1, tempY, Ey, nxn, nyn, nzn);
  addscale(1, tempZ, Ez, nxn, nyn, nzn);

  
  if (Cylindrical == 1)  fixBC_Cyl_Source(tempX, tempY, tempZ);
  else { 
     fixBC_Source (vct, col, tempX, tempY, tempZ);
     if (col->getCase() == "Dipole"){
       fixBC_SourceC(vct, col, tempXC, tempYC, tempZC);
       fixBC_PlanetSource (vct, col, grid, tempX , tempY , tempZ );
       fixBC_PlanetSourceC(vct, col, grid, tempXC, tempYC, tempZC);
     }
  }

  // physical space -> Krylov space
  phys2solver(bkrylov, tempX, tempY, tempZ, tempXC, tempYC, tempZC, nxn, nyn, nzn);
  #ifdef __PROFILING__
  clocks->stop(11);
  #endif
}


/*! Mapping of Maxwell image to give to solver */
void EMfields3D::MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  #ifdef __PROFILING__
  clocks->start(12);
  #endif
  int dim = 3*(nxn - 2) * (nyn - 2) * (nzn - 2) + 3*(nxc - 2) * (nyc - 2) * (nzc - 2);
  eqValue(0.0, im, dim);
  //eqValue(0.0, im, 6 * (nxn - 2) * (nyn - 2) * (nzn - 2));
  // move from krylov space to physical space: E -> vect   B -> temp
  solver2phys(tempX, tempY, tempZ, tempXC, tempYC, tempZC, vector, nxn, nyn, nzn);

  #ifdef __PROFILING__
  clocks->start(13);
  #endif
  communicateCenterBC(nxc, nyc, nzc, tempXC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempYC, 2, 2, 2, 2, 2, 2, vct);
  communicateCenterBC(nxc, nyc, nzc, tempZC, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempX, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempY, 2, 2, 2, 2, 2, 2, vct);
  communicateNodeBC(nxn, nyn, nzn, tempZ, 2, 2, 2, 2, 2, 2, vct);
  #ifdef __PROFILING__
  clocks->stop(13);
  #endif

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
  double cthdt = c*th*dt;

  if (Cylindrical == 1) {
    grid->curlN2C_cyl(imageEX, imageEY, imageEZ, tempX, tempY, tempZ);
  } else {
    grid->curlN2C(imageEX, imageEY, imageEZ, tempX, tempY, tempZ);
  }

  scale(imageEX, cthdt, nxc, nyc, nzc);
  scale(imageEY, cthdt, nxc, nyc, nzc);
  scale(imageEZ, cthdt, nxc, nyc, nzc);
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
  if (Cylindrical == 1) {
    grid->curlC2N_cyl(imageBX, imageBY, imageBZ, tempXC, tempYC, tempZC);
  } else {
    grid->curlC2N(imageBX, imageBY, imageBZ, tempXC, tempYC, tempZC);
  }
  scale(imageBX, -cthdt, nxn, nyn, nzn);
  scale(imageBY, -cthdt, nxn, nyn, nzn);
  scale(imageBZ, -cthdt, nxn, nyn, nzn);
  addscale(1, imageBX, tempX, nxn, nyn, nzn);
  addscale(1, imageBY, tempY, nxn, nyn, nzn);
  addscale(1, imageBZ, tempZ, nxn, nyn, nzn);

  double fac = dt*th*FourPI;
  for (int i=1; i<nxn-1; i++) 
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        double MEx, MEy, MEz;
        double invV = grid->getInvVOLn(i, j, k);
        productMassE(&MEx, &MEy, &MEz, tempX, tempY, tempZ, i, j, k);
        //imageBX[i][j][k] += invV*fac*MEx;
        //imageBY[i][j][k] += invV*fac*MEy;
        //imageBZ[i][j][k] += invV*fac*MEz;
        temp2X[i][j][k] = invV*fac*MEx;
        temp2Y[i][j][k] = invV*fac*MEy;
        temp2Z[i][j][k] = invV*fac*MEz;
      }

  //smooth(Smooth, temp2X, temp2Y, temp2Z, vct, col);

  for (int i=1; i<nxn-1; i++) 
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        imageBX[i][j][k] += temp2X[i][j][k];
        imageBY[i][j][k] += temp2Y[i][j][k];
        imageBZ[i][j][k] += temp2Z[i][j][k];
      }
  #endif

  // Temporal damping
  if (damping == 1) {
    sumscalprod(imageBX, 1.0, tempX, Lambda, nxn, nyn, nzn);
    sumscalprod(imageBY, 1.0, tempY, Lambda, nxn, nyn, nzn);
    sumscalprod(imageBZ, 1.0, tempZ, Lambda, nxn, nyn, nzn);
  }

  if (Cylindrical == 1)  fixBC_Cyl_Image(grid, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
  else {
    fixBC_Image (vct, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
    if (col->getCase() == "Dipole"){
      fixBC_ImageC(vct, imageEX, imageEY, imageEZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
      fixBC_PlanetImage (vct, col, grid, imageBX, imageBY, imageBZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
      fixBC_PlanetImageC(vct, col, grid, imageEX, imageEY, imageEZ, tempX, tempY, tempZ, tempXC, tempYC, tempZC);
    }
  }

  // move from physical space to krylov space
  phys2solver(im, imageBX, imageBY, imageBZ, imageEX, imageEY, imageEZ, nxn, nyn, nzn);

  #ifdef __PROFILING__
  clocks->stop(12);
  #endif
}

/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smoothEth(double value, VirtualTopology3D * vct, Collective *col) {

  int nvolte = 6;
  if (value != 1.0) {
    for (int icount = 1; icount < nvolte + 1; icount++) {
      double alpha;
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Exth, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Eyth, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Ezth, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

      if (Cylindrical == 1) {
        if (xStart == 0) {
          for (int j = 0; j < nyn ; j++)
            for (int k = 0; k < nzn ; k++) {
              Exth[0][j][k] = -Exth[2][j][k];
              Eyth[0][j][k] = Eyth[2][j][k];
              Ezth[0][j][k] = -Ezth[2][j][k];
            }
          }

      }

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Exth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Exth[i][j][k] + alpha * (Exth[i - 1][j][k] + Exth[i + 1][j][k] + 
                                                             Exth[i][j - 1][k] + Exth[i][j + 1][k] + 
                                                             Exth[i][j][k - 1] + Exth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Exth[i][j][k] = temp[i][j][k];
      // Eyth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Eyth[i][j][k] + alpha * (Eyth[i - 1][j][k] + Eyth[i + 1][j][k] + 
                                                             Eyth[i][j - 1][k] + Eyth[i][j + 1][k] + 
                                                             Eyth[i][j][k - 1] + Eyth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Eyth[i][j][k] = temp[i][j][k];
      // Ezth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Ezth[i][j][k] + alpha * (Ezth[i - 1][j][k] + Ezth[i + 1][j][k] + 
                                                             Ezth[i][j - 1][k] + Ezth[i][j + 1][k] + 
                                                             Ezth[i][j][k - 1] + Ezth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Ezth[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
}


void EMfields3D::smooth(double value, double ***Vx, double ***Vy, double ***Vz, VirtualTopology3D * vct, Collective *col) {

  int nvolte = 6;
  if (value != 1.0) {
    for (int icount = 1; icount < nvolte + 1; icount++) {
      double alpha;
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Vx, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Vy, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Vz, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Vx
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Vx[i][j][k] + alpha * (Exth[i - 1][j][k] + Exth[i + 1][j][k] + 
                                                             Vx[i][j - 1][k] + Exth[i][j + 1][k] + 
                                                             Vx[i][j][k - 1] + Exth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Vx[i][j][k] = temp[i][j][k];
      // Vy
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Vy[i][j][k] + alpha * (Eyth[i - 1][j][k] + Eyth[i + 1][j][k] + 
                                                             Vy[i][j - 1][k] + Eyth[i][j + 1][k] + 
                                                             Vy[i][j][k - 1] + Eyth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Vy[i][j][k] = temp[i][j][k];
      // Vz
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Vz[i][j][k] + alpha * (Ezth[i - 1][j][k] + Ezth[i + 1][j][k] + 
                                                             Vz[i][j - 1][k] + Ezth[i][j + 1][k] + 
                                                             Vz[i][j][k - 1] + Ezth[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Vz[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
}
/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smoothBth(double value, double ***Bxth, double ***Byth, double ***Bzth, VirtualTopology3D * vct, Collective *col) {

  int nvolte = 6;
  if (value != 1.0) {
    for (int icount = 1; icount < nvolte + 1; icount++) {
      double alpha;
      communicateCenterBoxStencilBC(nxc, nyc, nzc, Bxth, 2, 2, 2 ,2 ,2 ,2 ,vct);
      communicateCenterBoxStencilBC(nxc, nyc, nzc, Byth, 2, 2, 2 ,2 ,2 ,2 ,vct);
      communicateCenterBoxStencilBC(nxc, nyc, nzc, Bzth, 2, 2, 2 ,2 ,2 ,2 ,vct);
      if (Cylindrical == 1) {
        if (xStart == 0) {
          for (int j = 0; j < nyn ; j++)
            for (int k = 0; k < nzn ; k++) {
              Bxth[0][j][k] = -Bxth[1][j][k];
              Byth[0][j][k] =  Byth[1][j][k];
              Bzth[0][j][k] = -Bzth[1][j][k];
            }
          }

      }

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Bxth
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            temp[i][j][k] = value * Bxth[i][j][k] + alpha * (Bxth[i - 1][j][k] + Bxth[i + 1][j][k] + 
                                                             Bxth[i][j - 1][k] + Bxth[i][j + 1][k] + 
                                                             Bxth[i][j][k - 1] + Bxth[i][j][k + 1]);
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            Bxth[i][j][k] = temp[i][j][k];
      // Byth
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            temp[i][j][k] = value * Byth[i][j][k] + alpha * (Byth[i - 1][j][k] + Byth[i + 1][j][k] + 
                                                             Byth[i][j - 1][k] + Byth[i][j + 1][k] + 
                                                             Byth[i][j][k - 1] + Byth[i][j][k + 1]);
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            Byth[i][j][k] = temp[i][j][k];
      // Bzth
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            temp[i][j][k] = value * Bzth[i][j][k] + alpha * (Bzth[i - 1][j][k] + Bzth[i + 1][j][k] + 
                                                             Bzth[i][j - 1][k] + Bzth[i][j + 1][k] + 
                                                             Bzth[i][j][k - 1] + Bzth[i][j][k + 1]);
      for (int i = 1; i < nxc - 1; i++)
        for (int j = 1; j < nyc - 1; j++)
          for (int k = 1; k < nzc - 1; k++)
            Bzth[i][j][k] = temp[i][j][k];


      delArr3(temp, nxc, nyc);
    }
  }
}
/* Interpolation smoothing: Smoothing (vector must already have ghost cells) TO MAKE SMOOTH value as to be different from 1.0 type = 0 --> center based vector ; type = 1 --> node based vector ; */
void EMfields3D::smoothJh(double value, VirtualTopology3D * vct, Collective *col) {

  int nvolte = 6;
  if (value != 1.0) {
    for (int icount = 1; icount < nvolte + 1; icount++) {
      double alpha;
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Jxh, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Jyh, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
      communicateNodeBoxStencilBC(nxn, nyn, nzn, Jzh, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);

      if (Cylindrical == 1) {
        if (xStart == 0) {
          for (int j = 0; j < nyn ; j++)
            for (int k = 0; k < nzn ; k++) {
              Jxh[0][j][k] = -Jxh[2][j][k];
              Jyh[0][j][k] =  Jyh[2][j][k];
              Jzh[0][j][k] = -Jzh[2][j][k];
            }
          }

      }

      double ***temp = newArr3(double, nxn, nyn, nzn);
      if (icount % 2 == 1) {
        value = 0.;
      }
      else {
        value = 0.5;
      }
      alpha = (1.0 - value) / 6;
      // Jxhth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Jxh[i][j][k] + alpha * (Jxh[i - 1][j][k] + Jxh[i + 1][j][k] + Jxh[i][j - 1][k] + Jxh[i][j + 1][k] + Jxh[i][j][k - 1] + Jxh[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Jxh[i][j][k] = temp[i][j][k];
      // Jyhth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Jyh[i][j][k] + alpha * (Jyh[i - 1][j][k] + Jyh[i + 1][j][k] + Jyh[i][j - 1][k] + Jyh[i][j + 1][k] + Jyh[i][j][k - 1] + Jyh[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Jyh[i][j][k] = temp[i][j][k];
      // Jzhth
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            temp[i][j][k] = value * Jzh[i][j][k] + alpha * (Jzh[i - 1][j][k] + Jzh[i + 1][j][k] + Jzh[i][j - 1][k] + Jzh[i][j + 1][k] + Jzh[i][j][k - 1] + Jzh[i][j][k + 1]);
      for (int i = 1; i < nxn - 1; i++)
        for (int j = 1; j < nyn - 1; j++)
          for (int k = 1; k < nzn - 1; k++)
            Jzh[i][j][k] = temp[i][j][k];


      delArr3(temp, nxn, nyn);
    }
  }
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

double EMfields3D::productMassEV(double ***vX, int i, int j, int k, Grid* grid)
{
  double res = 0;

  // First the nodes that have been computed in the node i,j,k
  for (int g=0; g<NE_MASS; g++) {
    int i1 = i + NeNo.getX(g);
    int j1 = j + NeNo.getY(g);
    int k1 = k + NeNo.getZ(g);

    double V1 = 1.0/grid->getInvVOLn(i1,j1,k1);
    res += V1*vX[i1][j1][k1]*Mxx[g][i][j][k];
                              
    if (g == 0) continue;
    int i2 = i - NeNo.getX(g);
    int j2 = j - NeNo.getY(g);
    int k2 = k - NeNo.getZ(g);
    double V2 = 1.0/grid->getInvVOLn(i2,j2,k2);
    res += V2*vX[i2][j2][k2]*Mxx[g][i2][j2][k2];

  }
  return res;
}


/*! interpolate charge density and pressure density from node to center */
void EMfields3D::interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid) {
  // do we need communication or not really?
  for (int is=0; is<ns; is++) grid->interpN2C(rhocs[is], rhons[is]);
}

/*! communicate ghost for grid -> Particles interpolation */
void EMfields3D::communicateGhostP2G(VirtualTopology3D * vct) {

  int rank = vct->getCartesian_rank();

  // interpolate adding common nodes among processors
  communicateInterp(nxn, nyn, nzn, Jxh, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jyh, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jzh, 0, 0, 0, 0, 0, 0, vct);
  for (int i=0; i<ns; i++) {
    communicateInterp(nxn, nyn, nzn, i, Jxhs, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Jyhs, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Jzhs, 0, 0, 0, 0, 0, 0, vct);
  }
  for (int i=0; i<NE_MASS; i++) { 
    communicateInterp(nxn, nyn, nzn, i, Mxx, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Mxy, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Mxz, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Myx, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Myy, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Myz, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzx, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzy, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Mzz, 0, 0, 0, 0, 0, 0, vct);
  }

  // put the correct values on ghost cells
  communicateNode_P(nxn, nyn, nzn, Jxh, vct);
  communicateNode_P(nxn, nyn, nzn, Jyh, vct);
  communicateNode_P(nxn, nyn, nzn, Jzh, vct);
  for (int i=0; i<ns; i++) {
    communicateNode_P(nxn, nyn, nzn, Jxhs, i, vct);
    communicateNode_P(nxn, nyn, nzn, Jyhs, i, vct);
    communicateNode_P(nxn, nyn, nzn, Jzhs, i, vct);
  }
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


  // If in cylindrical coordinates, set to zero the mass matrices of the ghost node
  // at the left of the axis
  if (Cylindrical == 1 && xStart == 0) {
    for (int j = 0; j < nyn; j++)
      for (int k = 0; k < nzn; k++) 
        for (int g=0; g<NE_MASS; g++) {
				  Mxx[g][0][j][k] = 0;
				  Mxy[g][0][j][k] = 0;
				  Mxz[g][0][j][k] = 0;
				  Myx[g][0][j][k] = 0;
				  Myy[g][0][j][k] = 0;
				  Myz[g][0][j][k] = 0;
				  Mzx[g][0][j][k] = 0;
				  Mzy[g][0][j][k] = 0;
				  Mzz[g][0][j][k] = 0;
        }
  }

}

/*! communicate ghost for grid -> Particles interpolation: Only curernt and density*/
void EMfields3D::communicateGhostP2G_JP(VirtualTopology3D * vct) {

  int rank = vct->getCartesian_rank();

  // interpolate adding common nodes among processors
  for (int i=0; i<ns; i++) { 
    communicateInterp(nxn, nyn, nzn, i, Jxs, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Jys, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, Jzs, 0, 0, 0, 0, 0, 0, vct);

    communicateInterp(nxn, nyn, nzn, i, pXXsn, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, pXYsn, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, pXZsn, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, pYYsn, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, pYZsn, 0, 0, 0, 0, 0, 0, vct);
    communicateInterp(nxn, nyn, nzn, i, pZZsn, 0, 0, 0, 0, 0, 0, vct);
  }

  // put the correct values on ghost cells
  for (int i=0; i<ns; i++) { 
    communicateNode_P(nxn, nyn, nzn, Jxs, i, vct);
    communicateNode_P(nxn, nyn, nzn, Jys, i, vct);
    communicateNode_P(nxn, nyn, nzn, Jzs, i, vct);

    communicateNode_P(nxn, nyn, nzn, pXXsn, i, vct);
    communicateNode_P(nxn, nyn, nzn, pXYsn, i, vct);
    communicateNode_P(nxn, nyn, nzn, pXZsn, i, vct);
    communicateNode_P(nxn, nyn, nzn, pYYsn, i, vct);
    communicateNode_P(nxn, nyn, nzn, pYZsn, i, vct);
    communicateNode_P(nxn, nyn, nzn, pZZsn, i, vct);
  }
}


/*! communicate ghost for grid -> Particles interpolation: Only curernt and density*/
void EMfields3D::communicateGhostP2G_Rho(VirtualTopology3D * vct) {

  int rank = vct->getCartesian_rank();

  // interpolate adding common nodes among processors
  for (int i=0; i<ns; i++) { 
    communicateInterp(nxn, nyn, nzn, i, rhons, 2, 2, 2, 2, 2, 2, vct);
  }

  // put the correct values on ghost cells
  for (int i=0; i<ns; i++) { 
    communicateNode_P(nxn, nyn, nzn, rhons, i, vct);
  }
}


/* Interpolate the magnetic fields from centers to nodes */
void EMfields3D::centers2nodesB(VirtualTopology3D *vct, Grid *grid, Collective *col) {
  if (Cylindrical == 1) fixGhostB();
  grid->interpC2N(Bxn,Bxc);
  grid->interpC2N(Byn,Byc);
  grid->interpC2N(Bzn,Bzc);
  communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

/*
for (register int j = 0; j < nyn; j++)
for (register int k = 0; k < nzn; k++) {
printf("%d %d %13.6e %13.6e %13.6e\n",j,k,Bxn[1][j][k], Byn[1][j][k], Bzn[1][j][k]);
}
getchar();
*/

}


/*! add an amount of charge density to charge density field at node X,Y */
void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        rhons[is][X - i][Y - j][Z - k] += temp * invVOLn[X - i][Y - j][Z - k];
      }
}

/** add an amount of current density to current density - direction X to current density field on the node*/
void EMfields3D::addJxh(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jxhs[is][X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Y to current density field on the node */
void EMfields3D::addJyh(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jyhs[is][X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Z to current density field on the node */
void EMfields3D::addJzh(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jzhs[is][X - i][Y - j][Z - k] += temp;
      }
}

/** add an amount of current density to current density - direction X to current density field on the node*/
void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int ns) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jxs[ns][X - i][Y - j][Z - k] += temp*invVOLn[X - i][Y - j][Z - k];
      }
}
/** add an amount of current density - direction Y to current density field on the node */
void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int ns) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jys[ns][X - i][Y - j][Z - k] += temp*invVOLn[X - i][Y - j][Z - k];
      }
}
/** add an amount of current density - direction Z to current density field on the node */
void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int ns) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k];
        Jzs[ns][X - i][Y - j][Z - k] += temp*invVOLn[X - i][Y - j][Z - k];
      }
}

/** add an amount of current density to mass matrix field at node X,Y */
void EMfields3D::addMass(double val[3][3], int X, int Y, int Z, int ind) {
  Mxx[ind][X][Y][Z] += val[0][0];
  Mxy[ind][X][Y][Z] += val[0][1];
  Mxz[ind][X][Y][Z] += val[0][2];
  Myx[ind][X][Y][Z] += val[1][0];
  Myy[ind][X][Y][Z] += val[1][1];
  Myz[ind][X][Y][Z] += val[1][2];
  Mzx[ind][X][Y][Z] += val[2][0];
  Mzy[ind][X][Y][Z] += val[2][1];
  Mzz[ind][X][Y][Z] += val[2][2];
}
void EMfields3D::addMassXX(double val, int X, int Y, int Z, int ind) {
  Mxx[ind][X][Y][Z] += val;
}

/*! add an amount of pressure density - direction XX to current density field on the node */
void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXXsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}
/*! add an amount of pressure density - direction XY to current density field on the node */
void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}
/*! add an amount of pressure density - direction XZ to current density field on the node */
void EMfields3D::addPxz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pXZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}
/*! add an amount of pressure density - direction YY to current density field on the node */
void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYYsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}
/*! add an amount of pressure density - direction YZ to current density field on the node */
void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pYZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}
/*! add an amount of pressure density - direction ZZ to current density field on the node */
void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        pZZsn[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOLn[X - i][Y - j][Z - k];
}

void EMfields3D::setZeroRho() {
  eqValue(0.0, rhons, ns, nxn, nyn, nzn);
  eqValue(0.0, rhocs, ns, nxc, nyc, nzc);

}

/*! set to 0 all the densities fields expect from rho*/
void EMfields3D::setZeroDensities() {
  for (register int i = 0; i < nxn; i++)
    for (register int j = 0; j < nyn; j++)
      for (register int k = 0; k < nzn; k++) {
        Jxh [i][j][k] = 0.0;
        Jyh [i][j][k] = 0.0;
        Jzh [i][j][k] = 0.0;
      }

  setZeroMass();

  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          Jxs[kk][i][j][k] = 0.0;
          Jys[kk][i][j][k] = 0.0;
          Jzs[kk][i][j][k] = 0.0;
          Jxhs[kk][i][j][k] = 0.0;
          Jyhs[kk][i][j][k] = 0.0;
          Jzhs[kk][i][j][k] = 0.0;
          pXXsn[kk][i][j][k] = 0.0;
          pXYsn[kk][i][j][k] = 0.0;
          pXZsn[kk][i][j][k] = 0.0;
          pYYsn[kk][i][j][k] = 0.0;
          pYZsn[kk][i][j][k] = 0.0;
          pZZsn[kk][i][j][k] = 0.0;
        }
}

void EMfields3D::setZeroMass() {
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
}
/*!SPECIES: Sum the charge density of different species on NODES */
void EMfields3D::sumOverSpecies(VirtualTopology3D * vct) {
  for (int is = 0; is < ns; is++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
          rhon[i][j][k] += rhons[is][i][j][k];
          Jxh[i][j][k]  += Jxhs[is][i][j][k];
          Jyh[i][j][k]  += Jyhs[is][i][j][k];
          Jzh[i][j][k]  += Jzhs[is][i][j][k];
        }
}


/*! initialize Magnetic and Electric Field with initial configuration */
void EMfields3D::init(VirtualTopology3D * vct, Grid * grid, Collective *col) {

  
  double kx = 2*M_PI/Lx;
  double ky = 2*M_PI/Ly;
  double kz = 2*M_PI/Lz;
  if (restart1 == 0) {
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
 //     	Adds damping region padded near the edge of the y axis.
        	double external_radius = Ly/20.;
        	double scale_decay = external_radius/3.0;
        	for (int i=0; i < nxn; i++)
        		for (int j=0; j < nyn; j++)
        			for (int k=0; k < nzn; k++){
        				Lambda[i][j][k]  = 0.0;
        				double r = sqrt( pow(grid->getYN(i,j,k)-Ly/2.0,2.0) );
        				if(r>Ly/2-external_radius){
        					Lambda[i][j][k]  = 1.0* tanh((r-(Ly/2-external_radius))/scale_decay);
        				}
        			}

}


void EMfields3D::initDoubleGEM(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double pertGEM = 0.0;

  double globalx;
  double globaly;
  double globalz;

  const double coarsedx= grid->getDX();
  const double coarsedy= grid->getDY();  
  const double coarsedz= grid->getDZ();

  const double deltax= Lx/2.0;
  const double deltay= Ly/2.0;
 
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize DOUBLE GEM Challenge with Pertubation" << endl;
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
	  globalx= grid->getXN(i, j, k) + coarsedx;
	  globaly= grid->getYN(i, j, k) + coarsedy;
	  globalz= grid->getZN(i, j, k) + coarsedz;
         
          double yB = globaly - .25 * Ly;
	  double yT = globaly - .75 * Ly;
          double yBd = yB / delta;
          double yTd = yT / delta;

	  double xB = globalx - .25 * Lx;
	  double xT = globalx - .75 * Lx;
	  double xBd = xB / delta;
	  double xTd = xT / delta;
	  
	  double xpert;
	  double ypert;
	  
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
              double sech_yBd = 1. / cosh(yBd);
              double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) + tanh(-yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
          
	  xpert = globalx- Lx/4;
	  ypert = globaly- Ly/4;
	  double deltax= Lx/2.0;
	  double deltay= Ly/2.0;
	  if (xpert < Lx/2 and ypert < Ly/2)
	    {
	      Bxn[i][j][k] +=(B0x*pertGEM)*(M_PI/deltay)*cos(2*M_PI*xpert/deltax)*sin(M_PI*ypert/deltay  );
	      Byn[i][j][k] = B0y -(B0x*pertGEM)*(2*M_PI/deltax)*sin(2*M_PI*xpert/deltax)*cos(M_PI*ypert/deltay);
	    }
	  // add the second initial GEM perturbation                                                                     
	  xpert = globalx- 3*Lx/4;
	  ypert = globaly- 3*Ly/4;
	  if (xpert > Lx/2 and ypert > Ly/2)
	    {
	      Bxn[i][j][k] +=(B0x*pertGEM)*(M_PI/deltay)*cos(2*M_PI*xpert/deltax)*sin(M_PI*ypert/deltay  );
	      Byn[i][j][k] = B0y -(B0x*pertGEM)*(2*M_PI/deltax)*sin(2*M_PI*xpert/deltax)*cos(M_PI*ypert/deltay);
	    }
	  double exp_pert;
	  // add the initial X perturbation                                                                              
	  xpert = globalx- Lx/4;
	  ypert = globaly- Ly/4;
	  exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

	  Bxn[i][j][k] +=(B0x*pertX)*exp_pert*(
					       -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
					       -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
					       );

	  Byn[i][j][k] +=(B0x*pertX)*exp_pert*(
					       cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
					       +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
					       );
	  // add the second initial X perturbation                                                                       
	  xpert = globalx- 3*Lx/4;
	  ypert = globaly- 3*Ly/4;
	  exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

	  Bxn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						-cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
						-cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
						);

	  Byn[i][j][k] +=(-B0x*pertX)*exp_pert*(
						cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
						+sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
						);
	  // guide field                                                                                                 
	  Bzn[i][j][k] = B0z;

          // External field
	  Bx_ext[i][j][k] = 0;
	  By_ext[i][j][k] = 0;
	  Bz_ext[i][j][k] = 0;
          
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers; same thing as on nodes but on centers

    grid->interpN2C(Bxc, Bxn);
    grid->interpN2C(Byc, Byn);
    grid->interpN2C(Bzc, Bzn);
     
    // end initialize B on centers 
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


void EMfields3D::initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  const double pertX = 0.4;
  const double deltax = 8. * delta;
  const double deltay = 4. * delta;
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
           double xM = grid->getXN(i, j, k) - .25 * Lx;
           double xMshift = grid->getXN(i, j, k) - .75 * Lx;

           double yB = grid->getYN(i, j, k) - .25 * Ly;
           double yT = grid->getYN(i, j, k) - .75 * Ly;
           double yBd = yB / delta;
           double yTd = yT / delta;
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
            if (DriftSpecies[is]) {
               double sech_yBd = 1. / cosh(yBd);
               double sech_yTd = 1. / cosh(yTd);
              rhons[is][i][j][k] = rhoINIT[is] * sech_yBd * sech_yBd / FourPI;
              rhons[is][i][j][k] += rhoINIT[is] * sech_yTd * sech_yTd / FourPI;
            }
            else
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          // Magnetic field
          Bxn[i][j][k] = B0x * (-1.0 + tanh(yBd) - tanh(yTd));
          // add the initial GEM perturbation
          Bxn[i][j][k] += 0.;
          Byn[i][j][k] = B0y;
          // add the initial X perturbation
           double xMdx = xM / deltax;
           double xMshiftdx = xMshift / deltax;
           double yBdy = yB / deltay;
           double yTdy = yT / deltay;
           double humpB = exp(-xMdx * xMdx - yBdy * yBdy);
          Bxn[i][j][k] -= (B0x * pertX) * humpB * cos(2 * M_PI * xM / Lx) * sin(M_PI * yB / Ly);
          Byn[i][j][k] += (B0x * pertX) * humpB * sin(2 * M_PI * xM / Lx) * cos(M_PI * yB / Ly);

          // add the second initial X perturbation
           double humpT = exp(-xMshiftdx * xMshiftdx - yTdy * yTdy);
          Bxn[i][j][k] += (B0x * pertX) * humpT * cos(2 * M_PI * xMshift / Lx) * sin(M_PI * yT / Ly);
          Byn[i][j][k] -= (B0x * pertX) * humpT * sin(2 * M_PI * xMshift / Lx) * cos(M_PI * yT / Ly);

          // guide field
          Bzn[i][j][k] = B0z;

	  Bx_ext[i][j][k] = 0;
	  By_ext[i][j][k] = 0;
	  Bz_ext[i][j][k] = 0;
        }
    // communicate ghost
    communicateNodeBC(nxn, nyn, nzn, Bxn, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Byn, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateNodeBC(nxn, nyn, nzn, Bzn, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    // initialize B on centers
    // initialize B on centers
	grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
    // communicate ghost
    communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
    communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

//Initializes the vacuum fields for EMC2.
void EMfields3D::initWB8(VirtualTopology3D *vct, Grid *grid, Collective *col){
  if (vct->getCartesian_rank() == 0) {
    cout << "------------------------------------------" << endl;
    cout << "Initialize two coils" << endl;
    cout << "------------------------------------------" << endl;
    cout << "B0x                              = " << B0x << endl;
    cout << "B0y                              = " << B0y << endl;
    cout << "B0z                              = " << B0z << endl;
  }
  Fext = 1;
  double distance;

  //char BlogName[256];
  //sprintf(BlogName, "/home/aws/VacuumB.%d.csv", vct->getCartesian_rank());
  //FILE *Blog = fopen(BlogName, "w");
  //fprintf(Blog, "x,y,z,Bx,By,Bz\n");

   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
       for (int k=0; k < nzn; k++){
         for (int is=0; is < ns; is++){
           rhons[is][i][j][k] = rhoINIT[is] / FourPI;
         }
         Ex[i][j][k] = 0.0;
         Ey[i][j][k] = 0.0;
         Ez[i][j][k] = 0.0;
         double blp[3];
         double a=coilD/2.0;
         double xc=x_center;
         double yc=y_center;
         double zc=z_center;
         double deltax=coilSpacing/2.0;
         double deltay=coilSpacing/2.0;
         double deltaz=coilSpacing/2.0;
         double x = grid->getXN(i,j,k);
         double y = grid->getYN(i,j,k);
         double z = grid->getZN(i,j,k);
         /*
         loopZ(blp, x, y, z, a, xc, yc, zc+deltaz, -B0z);
         Bxn[i][j][k] = blp[0];
         Byn[i][j][k] = blp[1];
         Bzn[i][j][k] = blp[2];
         loopZ(blp, x, y, z, a, xc, yc, zc-deltaz, B0z);
         Bxn[i][j][k] += blp[0];
         Byn[i][j][k] += blp[1];
         Bzn[i][j][k] += blp[2];
         loopX(blp, x, y, z, a, xc+deltax, yc, zc, -B0x);
         Bxn[i][j][k] += blp[0];
         Byn[i][j][k] += blp[1];
         Bzn[i][j][k] += blp[2];
         loopX(blp, x, y, z, a, xc-deltax, yc, zc, B0x);
         Bxn[i][j][k] += blp[0];
         Byn[i][j][k] += blp[1];
         Bzn[i][j][k] += blp[2];
         */
         loopY(blp, x, y, z, a, xc, yc+deltay, zc, -B0y);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
         loopY(blp, x, y, z, a, xc, yc-deltay, zc, B0y);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];

        //fprintf(Blog, "%f,%f,%f,%E,%E,%E\n", x, y, z, Bxn[i][j][k], Byn[i][j][k], Bzn[i][j][k]);
      }
    }
  }

  // communicate ghost
  communicateNodeBC(nxn, nyn, nzn, Bx_ext, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, By_ext, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bz_ext, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  MPI_Barrier(MPI_COMM_WORLD);
  //fclose(Blog);

  for (int is=0 ; is<ns; is++)
    grid->interpN2C(rhocs,is,rhons);

  if (restart1 ==0){
    //  do nothing
  } else { // EM initialization from RESTART
    init(vct,grid,col);  // use the fields from restart file
  }

  double dst = col->getDelta()+dx;
	double scale_decay = dst/10.0;
  for (int i=0; i < nxn; i++)
  	for (int j=0; j < nyn; j++)
  		for (int k=0; k < nzn; k++){
  			Lambda[i][j][k]  = 0.0;
        double r;
        if (Cylindrical == 1) r = sqrt(pow(grid->getXN(i,j,k),2.0) + pow(grid->getYN(i,j,k)-Ly/2.0,2.0));
        else r = sqrt(pow(grid->getXN(i,j,k)-Lx/2.0,2.0) + pow(grid->getYN(i,j,k)-Ly/2.0,2.0) + pow(grid->getZN(i,j,k)-Lz/2.0,2.0));
  			if(r>dst){
  				Lambda[i][j][k]  = 1.0* tanh((r-dst)/scale_decay);
  			}
  		}
}

//Initializes the vacuum fields for EMC2. 6 coils
void EMfields3D::initWB8_6(VirtualTopology3D *vct, Grid *grid, Collective *col){
  if (vct->getCartesian_rank() == 0) {
    cout << "------------------------------------------" << endl;
    cout << "Initialize six coils" << endl;
    cout << "------------------------------------------" << endl;
    cout << "B0x                              = " << B0x << endl;
    cout << "B0y                              = " << B0y << endl;
    cout << "B0z                              = " << B0z << endl;
  }
  Fext = 1;
  double distance;

  //char BlogName[256];
  //sprintf(BlogName, "/home/aws/VacuumB.%d.csv", vct->getCartesian_rank());
  //FILE *Blog = fopen(BlogName, "w");
  //fprintf(Blog, "x,y,z,Bx,By,Bz\n");

   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
       for (int k=0; k < nzn; k++){
         for (int is=0; is < ns; is++){
           rhons[is][i][j][k] = rhoINIT[is] / FourPI;
         }
         Ex[i][j][k] = 0.0;
         Ey[i][j][k] = 0.0;
         Ez[i][j][k] = 0.0;
         double blp[3];
         double a=coilD/2.0;
         double xc=x_center;
         double yc=y_center;
         double zc=z_center;
         double deltax=coilSpacing/2.0;
         double deltay=coilSpacing/2.0;
         double deltaz=coilSpacing/2.0;
         double x = grid->getXN(i,j,k);
         double y = grid->getYN(i,j,k);
         double z = grid->getZN(i,j,k);
         loopZ(blp, x, y, z, a, xc, yc, zc+deltaz, -B0z);
         Bx_ext[i][j][k] = blp[0];
         By_ext[i][j][k] = blp[1];
         Bz_ext[i][j][k] = blp[2];
         loopZ(blp, x, y, z, a, xc, yc, zc-deltaz, B0z);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
         loopX(blp, x, y, z, a, xc+deltax, yc, zc, -B0x);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
         loopX(blp, x, y, z, a, xc-deltax, yc, zc, B0x);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
         loopY(blp, x, y, z, a, xc, yc+deltay, zc, -B0y);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
         loopY(blp, x, y, z, a, xc, yc-deltay, zc, B0y);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];

        //fprintf(Blog, "%f,%f,%f,%E,%E,%E\n", x, y, z, Bxn[i][j][k], Byn[i][j][k], Bzn[i][j][k]);
      }
    }
  }

  if (nzn  == 4) {
    if (vct->getCartesian_rank() == 0) {
      eqValue(0.0, Bz_ext, nxn, nyn, nzn);
      cout << "-> 2D case: Bz_ext = 0" << endl;
    }
  }

  // communicate ghost
  communicateNodeBC(nxn, nyn, nzn, Bx_ext, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, By_ext, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bz_ext, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  MPI_Barrier(MPI_COMM_WORLD);
  //fclose(Blog);

  for (int is=0 ; is<ns; is++)
    grid->interpN2C(rhocs,is,rhons);

  if (restart1 ==0){
    //  do nothing
  } else { // EM initialization from RESTART
    init(vct,grid,col);  // use the fields from restart file
  }
}

//Flux Rope based on pressure equilibrium
void EMfields3D::initFluxRope(VirtualTopology3D *vct, Grid *grid, Collective *col)
{

	if (restart1 ==0){
		if (vct->getCartesian_rank() ==0){
			cout << "----------------------------------------" << endl;
			cout << "       Initialize 3D Flux Rope" << endl;
			cout << "----------------------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i] << endl;
			}
			cout << "Smoothing Factor = " << Smooth << endl;
			cout << "-------------------------" << endl;
		}

    for (int i=0; i < nxn; i++)
		  for (int j=0; j < nyn; j++)
		    for (int k=0; k < nzn; k++){

          double r, teta;
          if (Cylindrical == 1) {
			      r = sqrt(pow(grid->getXN(i,j,k),2.0) + pow(grid->getZN(i,j,k),2.0));
			      teta = 0;
          } else {
			      r = sqrt(pow(grid->getXN(i,j,k)-Lx/2.0,2.0) + pow(grid->getYN(i,j,k)-Ly/2.0,2.0));
			      teta = atan2(grid->getYN(i,j,k)-Ly/2.0,grid->getXN(i,j,k)-Lx/2.0);
          }

		       // initialize the density for species
		       for (int is=0; is < ns; is++)
			       rhons[is][i][j][k] = rhoINIT[is] /FourPI * pow(delta,4.0) / pow(r*r + delta *delta, 2.0);

			    // electric field
			    Ex[i][j][k] =  0.0;
			    Ey[i][j][k] =  0.0;
			    Ez[i][j][k] =  0.0;
			    // Magnetic field

			    double Bth = B0x * r * delta /(r*r+ delta*delta);
          if (Cylindrical == 1) {
			      Bxn[i][j][k] = -Bth * sin(teta);
			      Byn[i][j][k] = B0z;
			      Bzn[i][j][k] = Bth * cos (teta);
          } else {
			      Bxn[i][j][k] = -Bth * sin(teta);
			      Byn[i][j][k] = Bth * cos (teta);
			      Bzn[i][j][k] = B0z;
          }
		    }

    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          double r, teta;
          if (Cylindrical == 1) {
    			  r = sqrt(pow(grid->getXC(i,j,k),2.0) + pow(grid->getZC(i,j,k),2.0));
    			  teta = 0;
          } else {
    			  r = sqrt(pow(grid->getXC(i,j,k)-Lx/2.0,2.0) + pow(grid->getYC(i,j,k)-Ly/2.0,2.0));
    			  teta = atan2(grid->getYC(i,j,k)-Ly/2.0,grid->getXC(i,j,k)-Lx/2.0);
          }

    			double Bth = B0x * r * delta /(r*r+ delta*delta);
          if (Cylindrical == 1) {
    			  Bxc[i][j][k] = -Bth * sin(teta);
    			  Byc[i][j][k] = B0z;
    			  Bzc[i][j][k] = Bth * cos (teta);
          } else {
    			  Bxc[i][j][k] = -Bth * sin(teta);
    			  Byc[i][j][k] = Bth * cos (teta);
    			  Bzc[i][j][k] = B0z;
          }

        }
     	// communicate ghost
      communicateCenterBC(nxc, nyc, nzc, Bxc, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
      communicateCenterBC(nxc, nyc, nzc, Byc, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
      communicateCenterBC(nxc, nyc, nzc, Bzc, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);


    // currents are used to calculate in the Maxwell's solver
    // The ion current is equal to 0 (all current is on electrons)
    for (int i=0; i < nxn; i++)
		  for (int j=0; j < nyn; j++)
		    for (int k=0; k < nzn; k++){
			    Jxs[1][i][j][k] = 0.0; // ion species is species 1
			    Jys[1][i][j][k] = 0.0; // ion species is species 1
			    Jzs[1][i][j][k] = 0.0; // ion species is species 1
		    }

    // calculate the electron current from
    eqValue(0.0,tempXN,nxn,nyn,nzn);
    eqValue(0.0,tempYN,nxn,nyn,nzn);
    eqValue(0.0,tempZN,nxn,nyn,nzn);
    grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc); // here you calculate curl(B)
    // all current is on electrons, calculated from Ampere's law
    for (int i=0; i < nxn; i++)
      for (int j=0; j < nyn; j++)
        for (int k=0; k < nzn; k++){  // electrons are species 0
			    Jxs[0][i][j][k] = c*tempXN[i][j][k]/FourPI; // ion species is species 1
			    Jys[0][i][j][k] = c*tempYN[i][j][k]/FourPI; // ion species is species 1
			    Jzs[0][i][j][k] = c*tempZN[i][j][k]/FourPI; // ion species is species 1
		    }

		for (int is=0 ; is<ns; is++) grid->interpN2C(rhocs,is,rhons);
	} else {
	    init(vct, grid, col);;  // use the fields from restart file
	}
	double scale_decay = L_square/10.0;
	for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
			for (int k=0; k < nzn; k++){
				Lambda[i][j][k]  = 0.0;
        double r;
        if (Cylindrical == 1) r = sqrt(pow(grid->getXN(i,j,k),2.0) + pow(grid->getYN(i,j,k)-Ly/2.0,2.0));
        else r = sqrt(pow(grid->getXN(i,j,k)-Lx/2.0,2.0) + pow(grid->getYN(i,j,k)-Ly/2.0,2.0) +  pow(grid->getZN(i,j,k)-Lz/2.0,2.0));
				if(r>L_square){
					Lambda[i][j][k]  = 1.0* tanh((r-L_square)/scale_decay);
				}
			}

}


//Initializes Dipole
void EMfields3D::initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col){

  if (vct->getCartesian_rank() == 0) {
    cout << "------------------------------------------" << endl;
    cout << "Initialize dipole" << endl;
    cout << "------------------------------------------" << endl;
    cout << "B1x                              = " << B1x << endl;
    cout << "B1y                              = " << B1y << endl;
    cout << "B1z                              = " << B1z << endl;
  }
  Fext = 1;
  double distance;

  SetLambda(grid);

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);
  
  double a=0.5*L_square;
  double xc=x_center;
  double yc=y_center;
  double zc=z_center;
  
  bool twodim = false;
  if ( nzn == 4) twodim = true;

   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
       for (int k=0; k < nzn; k++){
         for (int is=0; is < ns; is++){
           rhons[is][i][j][k] = rhoINIT[is] / FourPI;
         }
         Ex[i][j][k] = ebc[0];
         Ey[i][j][k] = ebc[1];
         Ez[i][j][k] = ebc[2];
         double blp[3];
         double x = grid->getXN(i,j,k);
         double y = grid->getYN(i,j,k);
         double z = grid->getZN(i,j,k);
/*
         loopY(blp, x, y, z, a, xc, yc, zc, B1y);
         Bx_ext[i][j][k] += blp[0];
         By_ext[i][j][k] += blp[1];
         Bz_ext[i][j][k] += blp[2];
*/
         double a=0.75*L_square; // 0.75 x To avoid problems with de-centered dipoles

         double xc=x_center;
         double yc=y_center;
         double zc=z_center;
  
         double rx = x-xc;
         double ry = y-yc;
         double rz = z-zc;
         if (twodim) rz = 0.0;
  
         double r      = sqrt(rx*rx + ry*ry + rz*rz);
  
         double Mx = B1x;
         double My = B1y;
         double Mz = B1z;
  
         if (r < 1.0*a) {
           ry = sqrt(a*a - (rx*rx + rz*rz));
           double one_r3 = 1.0/(a*a*a);
           double rhx = rx/a;
           double rhy = ry/a;
           double rhz = rz/a;
  
           double Bxe = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhx - Mx );
           double Bye = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhy - My );
           double Bze = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhz - Mz );
           double Bme = sqrt(Bxe*Bxe + Bye*Bye + Bze*Bze);
  
           Bx_ext[i][j][k] = 0.0;
           By_ext[i][j][k] = (B1y == 0 ? 0.0 : (B1y/fabs(B1y)) * Bme);
           Bz_ext[i][j][k] = 0.0;
         }
         else {
  
           double one_r3 = 1.0/(r*r*r);
  
           double rhx = rx/r;
           double rhy = ry/r;
           double rhz = rz/r;
  
           Bx_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhx - Mx );
           By_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhy - My );
           Bz_ext[i][j][k] = one_r3 * ( 3 * (Mx*rhx+My*rhy+Mz*rhz) * rhz - Mz );
         }
         Bxn[i][j][k] = B0x;
         Byn[i][j][k] = B0y;
         Bzn[i][j][k] = B0z;
      }
    }
  }

  for (int is=0 ; is<ns; is++)
    grid->interpN2C(rhocs,is,rhons);

  // initialize B on centers
  for (int i=0; i < nxc; i++){
    for (int j=0; j < nyc; j++){
      for (int k=0; k < nzc; k++){
        Bxc[i][j][k] = B0x;
        Byc[i][j][k] = B0y;
        Bzc[i][j][k] = B0z;
      }
    }
  }
}

void EMfields3D::initHalfPlanet(VirtualTopology3D *vct, Grid *grid, Collective *col){
  if (vct->getCartesian_rank() == 0) {
    cout << "------------------------------------------" << endl;
    cout << "Initialize half planet" << endl;
    cout << "------------------------------------------" << endl;
    cout << "B1x                              = " << B1x << endl;
    cout << "B1y                              = " << B1y << endl;
    cout << "B1z                              = " << B1z << endl;
  }
  Fext = 1;
  double distance;

  double ebc[3];
  cross_product(ue0,ve0,we0,B0x,B0y,B0z,ebc);
  scale(ebc,-1.0,3);

  double a=0.75*L_square;
  double xc=Lx;
  double yc=y_center;
  double zc=z_center;
  
   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
       for (int k=0; k < nzn; k++){
         for (int is=0; is < ns; is++){
           rhons[is][i][j][k] = rhoINIT[is] / FourPI;
         }
         //Ex[i][j][k] = ebc[0];
         //Ey[i][j][k] = ebc[1];
         //Ez[i][j][k] = ebc[2];
         double blp[3];
         double x = grid->getXN(i,j,k);
         double y = grid->getYN(i,j,k);
         double z = grid->getZN(i,j,k);

         loopY(blp, x, y, z, a, xc, yc, zc, -B1y);
         Bx_ext[i][j][k] = blp[0];
         By_ext[i][j][k] = blp[1];
         Bz_ext[i][j][k] = blp[2];

         Bxn[i][j][k] = B0x;
         Byn[i][j][k] = B0y;
         Bzn[i][j][k] = B0z;
      }
    }
  }

  // communicate ghost
  communicateNodeBC(nxn, nyn, nzn, Bx_ext, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, By_ext, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, Bz_ext, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  // initialize B on centers
  for (int i=0; i < nxc; i++){
    for (int j=0; j < nyc; j++){
      for (int k=0; k < nzc; k++){
        Bxc[i][j][k] = B0x;
        Byc[i][j][k] = B0y;
        Bzc[i][j][k] = B0z;
      }
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  //fclose(Blog);

  for (int is=0 ; is<ns; is++)
    grid->interpN2C(rhocs,is,rhons);

  if (restart1 ==0){
    //  do nothing
  } else { // EM initialization from RESTART
    init(vct,grid,col);  // use the fields from restart file
  }
}

#define MODE_A

/*! initiliaze EM for Bessel test */
void EMfields3D::initBessel1D(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize Bessel function 1D            " << endl;
      cout << "------------------------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
#ifdef MODE_A
          Ex[i][j][k] = 0;
          Ey[i][j][k] = 0;
          Ez[i][j][k] = besselJ(1, x*ZERO_J1/Lx); //besselJ(1, x*ZERO_J1/Lx);
          // Magnetic field
          Bxn[i][j][k] = 0;
          Byn[i][j][k] = 0;
          Bzn[i][j][k] = 0;
#else

          Ex[i][j][k] = 0;
          Ey[i][j][k] = 0;
          Ez[i][j][k] = 0;
          // Magnetic field
          Bxn[i][j][k] = 0;
          Byn[i][j][k] = 0;
          Bzn[i][j][k] = besselJ(1, x*ZERO_J0/Lx); 
#endif
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          double x = grid->getXC(i,j,k);
          double y = grid->getYC(i,j,k);
#ifdef MODE_A
          Bxc[i][j][k] = 0;
          Byc[i][j][k] = 0;
          Bzc[i][j][k] = 0;
#else
          // Magnetic field
          Bxc[i][j][k] = 0;
          Byc[i][j][k] = 0;
          Bzc[i][j][k] = besselJ(1, x*ZERO_J0/Lx);
#endif
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initBessel2D(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize Bessel function 2D            " << endl;
      cout << "------------------------------------------" << endl;
    }
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
          double omega = sqrt((ZERO_J0/Lx)*(ZERO_J1/Lx) + (2*M_PI/Ly)*(2*M_PI/Ly));
          double kz = 2*M_PI/Ly;
#ifdef MODE_A
          Ex[i][j][k] = 0;
          Ey[i][j][k] = besselJ(0, x*ZERO_J0/Lx);
          Ez[i][j][k] = sin(y*2*M_PI/Ly)*besselJ(1, x*ZERO_J1/Lx); 
          // Magnetic field
          Bxn[i][j][k] = 0;
          Byn[i][j][k] = 0;
          Bzn[i][j][k] = 0;
#else
          Ex[i][j][k] = 0;
          Ey[i][j][k] = 0;
          Ez[i][j][k] = 0;
          // Magnetic field
          Bxn[i][j][k] = 0;
          Byn[i][j][k] = 0;
          Bzn[i][j][k] = sin(y*2*M_PI/Ly)*besselJ(1, x*ZERO_J1/Lx);
#endif
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          double x = grid->getXC(i,j,k);
          double y = grid->getYC(i,j,k);
          double kz = 2*M_PI/Ly;
#ifdef MODE_A
          Bxc[i][j][k] = 0;
          Byc[i][j][k] = 0;
          Bzc[i][j][k] = 0;
#else
          Bxc[i][j][k] = 0;
          Byc[i][j][k] = 0;
          Bzc[i][j][k] = sin(y*2*M_PI/Ly)*besselJ(1, x*ZERO_J1/Lx);
#endif
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initTEmode(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize TE mode                       " << endl;
      cout << "------------------------------------------" << endl;
    }

    double Kz = 2*M_PI/Ly;
    double Kr = ZERO_J1/Lx;
    double omega = sqrt(Kz*Kz + Kr*Kr);

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
          Ex[i][j][k] = 0;
          Ey[i][j][k] = 0;
          Ez[i][j][k] = -(omega/Kr)*besselJ(1, x*Kr)*sin(y*Kz);
          // Magnetic field
          Bxn[i][j][k] = (Kz/Kr)*besselJ(1, x*Kr)*sin(y*Kz);
          Byn[i][j][k] = besselJ(0, x*Kr)*cos(y*Kz);
          Bzn[i][j][k] = 0;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          double x = grid->getXC(i,j,k);
          double y = grid->getYC(i,j,k);
          Bxc[i][j][k] = (Kz/Kr)*besselJ(1, x*Kr)*sin(y*Kz);
          Byc[i][j][k] = besselJ(0, x*Kr)*cos(y*Kz);
          Bzc[i][j][k] = 0;
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}

void EMfields3D::initTMmode(VirtualTopology3D * vct, Grid * grid, Collective *col) {
  // perturbation localized in X
  if (restart1 == 0) {
    // initialize
    if (vct->getCartesian_rank() == 0) {
      cout << "------------------------------------------" << endl;
      cout << "Initialize TM mode                       " << endl;
      cout << "------------------------------------------" << endl;
    }
    double Kz = 2*M_PI/Ly;
    double Kr = ZERO_J0/Lx;
    double omega = sqrt(Kz*Kz + Kr*Kr);

    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          // initialize the density for species
          for (int is = 0; is < ns; is++) {
              rhons[is][i][j][k] = rhoINIT[is] / FourPI;
          }
          // electric field
          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
          Ex[i][j][k] = Kz*besselJ(1, x*Kr)*sin(Kz*y);
          Ey[i][j][k] = besselJ(0, x*Kr)*cos(Kz*y);
          Ez[i][j][k] = 0;
          // Magnetic field
          Bxn[i][j][k] = 0;
          Byn[i][j][k] = 0;
          Bzn[i][j][k] = omega*besselJ(1, x*Kr)*sin(Kz*y);
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          double x = grid->getXC(i,j,k);
          double y = grid->getYC(i,j,k);
          Bxc[i][j][k] = 0;
          Byc[i][j][k] = 0;
          Bzc[i][j][k] = omega*besselJ(1, x*Kr)*sin(Kz*y);
        }
    for (int is = 0; is < ns; is++)
      grid->interpN2C(rhocs, is, rhons);
  }
  else {
    init(vct, grid, col);            // use the fields from restart file
  }
}


void EMfields3D::initPulse(VirtualTopology3D * vct, Grid * grid, Collective *col) {
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          double x = grid->getXN(i,j,k);
          double y = grid->getYN(i,j,k);
          double z = grid->getZN(i,j,k);
          double xl = x-x_center;
          double yl = y-y_center;
          double zl = z-z_center;
          double r = sqrt(xl*xl + yl*yl + zl*zl);
          if (r < L_square) Ez[i][j][k] = 1.0;

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

double ***EMfields3D::getRhotilde(int is) {
  return rhotilde[is];
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

/*! get Ex_ext  */
double ***EMfields3D::getEx_ext() {
  return(Ex_ext);
}
/*!  get Ey_ext */
double ***EMfields3D::getEy_ext() {
  return(Ey_ext);
}
/*!  get Ez_ext */
double ***EMfields3D::getEz_ext() {
  return(Ez_ext);
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

double*** EMfields3D::getDivE() {
  return divE;
}

/*! SPECIES: get pressure tensor component XX defined on nodes */
double ****EMfields3D::getpXXsn() {
  return (pXXsn);
}
/*! SPECIES: get pressure tensor component XY defined on nodes */
double ****EMfields3D::getpXYsn() {
  return (pXYsn);
}
/*! SPECIES: get pressure tensor component XZ defined on nodes */
double ****EMfields3D::getpXZsn() {
  return (pXZsn);
}
/*! SPECIES: get pressure tensor component YY defined on nodes */
double ****EMfields3D::getpYYsn() {
  return (pYYsn);
}
/*! SPECIES: get pressure tensor component YZ defined on nodes */
double ****EMfields3D::getpYZsn() {
  return (pYZsn);
}
/*! SPECIES: get pressure tensor component ZZ defined on nodes */
double ****EMfields3D::getpZZsn() {
  return (pZZsn);
}

/*! get current -Direction X */
double &EMfields3D::getJxh(int indexX, int indexY, int indexZ) const {
  return (Jxh[indexX][indexY][indexZ]);
}
/*! get current array X component * */
double ***EMfields3D::getJxh() {
  return (Jxh);
}
double ***EMfields3D::getJx_ext() {
  return (Jx_ext);
}
/*! get current -Direction Y */
double &EMfields3D::getJyh(int indexX, int indexY, int indexZ) const {
  return (Jyh[indexX][indexY][indexZ]);
}
/*! get current array Y component * */
double ***EMfields3D::getJyh() {
  return (Jyh);
}
double ***EMfields3D::getJy_ext() {
  return (Jy_ext);
}
/*! get current -Direction Z */
double &EMfields3D::getJzh(int indexX, int indexY, int indexZ) const {
  return (Jzh[indexX][indexY][indexZ]);
}
/*! get current array Z component * */
double ***EMfields3D::getJzh() {
  return (Jzh);
}
double ***EMfields3D::getJz_ext() {
  return (Jz_ext);
}

/*! get current -Direction X */
double &EMfields3D::getJxs(int indexX, int indexY, int indexZ, int is) const {
  return (Jxs[is][indexX][indexY][indexZ]);
}
/*! get current array X component * */
double ****EMfields3D::getJxs() {
  return (Jxs);
}
double ***EMfields3D::getJxs(int is) {
  return (Jxs[is]);
}
/*! get current -Direction Y */
double &EMfields3D::getJys(int indexX, int indexY, int indexZ, int is) const {
  return (Jys[is][indexX][indexY][indexZ]);
}
/*! get current array Y component * */
double ****EMfields3D::getJys() {
  return (Jys);
}
double ***EMfields3D::getJys(int is) {
  return (Jys[is]);
}
/*! get current -Direction Z */
double &EMfields3D::getJzs(int indexX, int indexY, int indexZ, int is) const {
  return (Jzs[is][indexX][indexY][indexZ]);
}
/*! get current array Z component * */
double ****EMfields3D::getJzs() {
  return (Jzs);
}
double ***EMfields3D::getJzs(int is) {
  return (Jzs[is]);
}

/* get the pressure tensor */
double ***EMfields3D::getPxx(int is) { return pXXsn[is]; }
double ***EMfields3D::getPxy(int is) { return pXYsn[is]; }
double ***EMfields3D::getPxz(int is) { return pXZsn[is]; }
double ***EMfields3D::getPyy(int is) { return pYYsn[is]; }
double ***EMfields3D::getPyz(int is) { return pYZsn[is]; }
double ***EMfields3D::getPzz(int is) { return pZZsn[is]; }


double EMfields3D::getJxcs (int i,int j,int k,int is) {
  double ret = .125 * (Jxs[is][i][j][k]     +
                       Jxs[is][i + 1][j][k] +
                       Jxs[is][i][j + 1][k] +
                       Jxs[is][i][j][k + 1] +
                       Jxs[is][i + 1][j + 1][k] +
                       Jxs[is][i + 1][j][k + 1] +
                       Jxs[is][i][j + 1][k + 1] +
                       Jxs[is][i + 1][j + 1][k + 1]);

  return ret;
}

double EMfields3D::getJycs (int i,int j,int k,int is) {
  double ret = .125 * (Jys[is][i][j][k]     +
                       Jys[is][i + 1][j][k] +
                       Jys[is][i][j + 1][k] +
                       Jys[is][i][j][k + 1] +
                       Jys[is][i + 1][j + 1][k] +
                       Jys[is][i + 1][j][k + 1] +
                       Jys[is][i][j + 1][k + 1] +
                       Jys[is][i + 1][j + 1][k + 1]);

  return ret;
}

double EMfields3D::getJzcs (int i,int j,int k,int is) {
  double ret = .125 * (Jzs[is][i][j][k]     +
                       Jzs[is][i + 1][j][k] +
                       Jzs[is][i][j + 1][k] +
                       Jzs[is][i][j][k + 1] +
                       Jzs[is][i + 1][j + 1][k] +
                       Jzs[is][i + 1][j][k + 1] +
                       Jzs[is][i][j + 1][k + 1] +
                       Jzs[is][i + 1][j + 1][k + 1]);

  return ret;
}

/*! get the electric field energy */
double EMfields3D::getEenergy(void) {
  double localEenergy = 0.0;
  double totalEenergy = 0.0;

  double err; 

  for (int i = 1; i < nxn-2 ; i++)
    for (int j = 1; j < nyn-2 ; j++)
      for (int k = 1; k < nzn-2 ; k++) {
        localEenergy += .5*(Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]) / (FourPI*invVOLn[i][j][k]);
        //localEenergy += dx * dy * dz * (Ex[i][j][k]*Exth[i][j][k] + Ey[i][j][k]*Eyth[i][j][k] + Ez[i][j][k]*Ezth[i][j][k]) / (FourPI);
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
        localBenergy += .5*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI*invVOLc[i][j][k]);
      }



  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalBenergy);
}

/*! get the external magnetic field energy */
double EMfields3D::getBextenergy(void) {
  double localBenergy = 0.0;
  double totalBenergy = 0.0;
  double Bxt = 0.0;
  double Byt = 0.0;
  double Bzt = 0.0;

  for (int i = 1; i < nxn-2; i++)
    for (int j = 1; j < nyn-2; j++)
      for (int k = 1; k < nzn-2; k++){
        Bxt = Bx_ext[i][j][k];
        Byt = By_ext[i][j][k];
        Bzt = Bz_ext[i][j][k];
        localBenergy += .5*(Bxt*Bxt + Byt*Byt + Bzt*Bzt)/(FourPI*invVOLn[i][j][k]);
      }



  MPI_Allreduce(&localBenergy, &totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalBenergy);
}


void EMfields3D::computeRhotilde(VirtualTopology3D * vct, Collective *col, Grid* grid, int is)
{
  // IMPORTANT: rhons has to have the desired value of the density
    cout << "*** RHO CALCULATION ***" << endl;

  int dim = (nxn - 2) * (nyn - 2) * (nzn - 2);
  double *xkrylov = new double[dim];  
  double *bkrylov = new double[dim];  

  // set to zero all the stuff 
  eqValue(0.0, xkrylov, dim);
  eqValue(0.0, bkrylov, dim);

  // prepare the source 
  int ind = 0;
  for (int i = 1; i < nxn-1; i++)
    for (int j = 1; j < nyn-1; j++)
      for (int k = 1; k < nzn-1; k++, ind++) {
        bkrylov[ind] = rhons[is][i][j][k]/grid->getInvVOLn(i,j,k);
        xkrylov[ind] = rhons[is][i][j][k]/grid->getInvVOLn(i,j,k);
      }
  // solver
  GMRES(&Field::RhoImage, xkrylov, dim, bkrylov, 20, 200, 1e-16, grid, vct, col, this);

  ind = 0;
  for (int i = 1; i < nxn-1; i++)
    for (int j = 1; j < nyn-1; j++)
      for (int k = 1; k < nzn-1; k++, ind++) {
        rhotilde[is][i][j][j] = xkrylov[ind];
      }
  delete [] xkrylov;
  delete [] bkrylov;
}

void EMfields3D::RhoImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col) {
  int dim = (nxn - 2) * (nyn - 2) * (nzn - 2);
  eqValue(0.0, im, dim);
  // move from krylov space to physical space: E -> vect   B -> temp
  int ind = 0;
  for (int i = 1; i < nxn-1; i++)
    for (int j = 1; j < nyn-1; j++)
      for (int k = 1; k < nzn-1; k++, ind++) {
        tempX[i][j][k] = vector[ind];
      }

  communicateNodeBC(nxn, nyn, nzn, tempX, 2, 2, 2, 2, 2, 2, vct);

  for (int i=1; i<nxn-1; i++) 
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        double MEx, MEy, MEz;
        double VrhoM = productMassEV(tempX, i, j, k, grid);
        tempY[i][j][k] = VrhoM;
      }

  ind = 0;
  for (int i = 1; i < nxn-1; i++)
    for (int j = 1; j < nyn-1; j++)
      for (int k = 1; k < nzn-1; k++, ind++) {
        im[ind] = tempY[i][j][k];
      }

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
  //delArr3(Ex_prev, nxn, nyn);
  //delArr3(Ey_prev, nxn, nyn);
  //delArr3(Ez_prev, nxn, nyn);
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
  delArr4(Jxs, nxn, nyn, nzn);
  delArr4(Jys, nxn, nyn, nzn);
  delArr4(Jzs, nxn, nyn, nzn);
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
  delArr4(rhotilde, ns, nxn, nyn);
  delArr3(divE, nxn, nyn);
  // central points
  delArr3(Bxc, nxc, nyc);
  delArr3(Byc, nxc, nyc);
  delArr3(Bzc, nxc, nyc);
  delArr3(rhoc, nxc, nyc);
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
  delArr3(Bx_ext, nxn, nyn);
  delArr3(By_ext, nxn, nyn);
  delArr3(Bz_ext, nxn, nyn);
  delArr3(Ex_ext, nxn, nyn);
  delArr3(Ey_ext, nxn, nyn);
  delArr3(Ez_ext, nxn, nyn);
}




/*! Test of the curl in cylindrical coordinates */
//#define TEST_NODES
#ifdef TEST_NODES
void EMfields3D::testCurl(Grid* grid)
{
  printf("Test Curl Nodes -> Centers\n");
  FILE* fd = fopen("testCurl.txt","w");

  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) {
        double r = grid->getXN(i,j,k);
        double z = grid->getYN(i,j,k);
        double p = grid->getZN(i,j,k);
        //Ex[i][j][k] = sin(z);
        //Ey[i][j][k] = sin(z);
        //Ez[i][j][k] = sin(z);
        //Ex[i][j][k] = cos(z);
        //Ey[i][j][k] = cos(z);
        //Ez[i][j][k] = cos(z);
        //Ex[i][j][k] = cos(r);
        //Ey[i][j][k] = cos(r);
        //Ez[i][j][k] = cos(r);
        Ex[i][j][k] = cos(r*z);
        Ey[i][j][k] = cos(r*z);
        Ez[i][j][k] = cos(r*z);
        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e\n", r,z,p,Ex[i][j][k],Ey[i][j][k],Ez[i][j][k]);
      } 

  fprintf(fd,"\n\n\n\n");

  grid->curlN2C_cyl(Bxc, Byc, Bzc, Ex, Ey, Ez);

  int ntc = (nxc-2)*(nyc-2)*(nzc-2);
  double errR = 0;
  double errZ = 0;
  double errP = 0;
  for (int i=1; i<nxc-1; i++)
    for (int j=1; j<nyc-1; j++) 
      for (int k=1; k<nzc-1; k++) {
        double r = grid->getXC(i,j,k);
        double z = grid->getYC(i,j,k);
        double p = grid->getZC(i,j,k);

        // Sin(z)
        //double cR = -cos(z);
        //double cZ = sin(z)/r;
        //double cP = -cos(z);
        // Cos(z)
        //double cR = sin(z);
        //double cZ = cos(z)/r;
        //double cP = -sin(z);
        // Cos(r)
        //double cR = 0.0;
        //double cZ = (cos(r)-r*sin(r))/r;
        //double cP = sin(r);
        // Cos(r*z)
        double cR = r*sin(r*z);
        double cZ = (cos(r*z)/r-z*sin(r*z));
        double cP = z*sin(r*z) - r*sin(r*z);

        double nR = Bxc[i][j][k];
        double nZ = Byc[i][j][k];
        double nP = Bzc[i][j][k];

        double eR = fabs(cR - nR);
        double eZ = fabs(cZ - nZ);
        double eP = fabs(cP - nP);
        errR += eR;
        errZ += eZ;
        errP += eP;
        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e   %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", r,z,p,nR,nZ,nP, cR, cZ, cP, eR, eZ, eP);
      } 
  fprintf(stderr,"%3d %3d  %13.6e %13.6e %13.6e  %13.6e \n",nxc-2, nyc-2, sqrt(errR)/ntc, sqrt(errZ)/ntc, errP/ntc, sqrt(errR+errZ+errP)/ntc);
  fclose(fd);
}
#else
void EMfields3D::testCurl(Grid* grid)
{
  printf("Test Curl Centers -> Nodes\n");
  FILE* fd = fopen("testCurl.txt","w");

  for (int i=0; i<nxc; i++)
    for (int j=0; j<nyc; j++) 
      for (int k=0; k<nzc; k++) {
        double r = grid->getXC(i,j,k);
        double z = grid->getYC(i,j,k);
        double p = grid->getZC(i,j,k);
        //Bxc[i][j][k] = cos(z);
        //Byc[i][j][k] = cos(z);
        //Bzc[i][j][k] = cos(z);
        //Bxc[i][j][k] = cos(r);
        //Byc[i][j][k] = cos(r);
        //Bzc[i][j][k] = cos(r);
        Bxc[i][j][k] = cos(r*z);
        Byc[i][j][k] = cos(r*z);
        Bzc[i][j][k] = cos(r*z);
        //Bxc[i][j][k] = 0;
        //Byc[i][j][k] = 0;
        //Bzc[i][j][k] = r*cos(r*z);
        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e\n", r,z,p,Bxc[i][j][k],Byc[i][j][k],Bzc[i][j][k]);
      } 

  fprintf(fd,"\n\n\n\n");

  grid->curlC2N_cyl(Ex, Ey, Ez, Bxc, Byc, Bzc);

  int ntc = (nxn-4)*(nyn-2)*(nzn-2);
  double errR = 0;
  double errZ = 0;
  double errP = 0;
  for (int i=2; i<nxn-2; i++)
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        double r = grid->getXN(i,j,k);
        double z = grid->getYN(i,j,k);
        double p = grid->getZN(i,j,k);

        // Cos(z)
        //double cR = sin(z);
        //double cZ = cos(z)/r;
        //double cP = -sin(z);
        // Cos(r)
        //double cR = 0.0;
        //double cZ = (cos(r)-r*sin(r))/r;
        //double cP = sin(r);
        // Cos(r*z)
        double cR = r*sin(r*z);
        double cZ = (cos(r*z)/r-z*sin(r*z));
        double cP = z*sin(r*z) - r*sin(r*z);
        // Bx = 0 By = 0 Bz=r*cos(r*z);
        //double cR = r*r*sin(r*z);
        //double cZ = (2*cos(r*z)-r*z*sin(r*z));
        //double cP = 0;
      
        //cR /= grid->getInvVOLn(i,j,k);
        //cZ /= grid->getInvVOLn(i,j,k);
        //cP /= grid->getInvVOLn(i,j,k);

        double eR = fabs(cR - Ex[i][j][k]);
        double eZ = fabs(cZ - Ey[i][j][k]);
        double eP = fabs(cP - Ez[i][j][k]);
        errR += eR;
        errZ += eZ;
        errP += eP;
        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e   %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n", r,z,p,Ex[i][j][k],Ey[i][j][k],Ez[i][j][k], cR, cZ, cP, eR, eZ, eP);
      } 
  fprintf(stderr,"%3d %3d  %13.6e %13.6e %13.6e  %13.6e \n",nxc-2, nyc-2, sqrt(errR)/ntc, sqrt(errZ)/ntc, errP/ntc, sqrt(errR+errZ+errP)/ntc);
  fclose(fd);
}
#endif

/*! Test of the divergence in cylindrical coordinates */
void EMfields3D::testDivergence(Grid* grid)
{
  printf("Test Divergence Nodes -> Centers\n");
  FILE* fd = fopen("testDivergence.txt","w");

  for (int i=0; i<nxn; i++)
    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) {
        double r = grid->getXN(i,j,k);
        double z = grid->getYN(i,j,k);
        double p = grid->getZN(i,j,k);
        Ex[i][j][k] = r*r*z;
        Ey[i][j][k] = 5*cos(2*z);
        Ez[i][j][k] = 0;
        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e\n", r,z,p,Ex[i][j][k],Ey[i][j][k],Ez[i][j][k]);
      } 

  fprintf(fd,"\n\n\n\n");

  grid->divN2C_cyl(tempX, Ex, Ey, Ez);

  int ntc = (nxc-2)*(nyc-2)*(nzc-2);
  double err = 0;
  for (int i=1; i<nxc-1; i++)
    for (int j=1; j<nyc-1; j++) 
      for (int k=1; k<nzc-1; k++) {
        double r = grid->getXC(i,j,k);
        double z = grid->getYC(i,j,k);
        double p = grid->getZC(i,j,k);

        double ref = 3*r*z - sin(2*z)*10;
        double val = tempX[i][j][k];

        fprintf(fd,"%8.4f %8.4f %8.4f  %13.6e  %13.6e  %13.6e \n", r,z,p, ref, val, ref-val);
				err += (ref-val)*(ref-val);
      } 
  fprintf(stderr,"%3d %3d  %13.6e %13.6e \n",nxc-2, nyc-2, sqrt(err)/ntc, err/ntc);
  fclose(fd);
}

/*! Check if JdotE is equal from the particles and from the fields.
 * It also checks if the variation of the magnetic energy equals the
 * variation of the kinetic energy (and equals to JdotE) */
void EMfields3D::testJdotE(Grid* grid, double deltaK) 
{
  double Een = getEenergy();
  double Ben = getBenergy();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char name[128];
  sprintf(name, "testJdotE_%d.txt",rank);
  FILE* fd = fopen(name,"w");

  double delE = 0;
  double maxErr = 0;
  for (int i=1; i<nxn-1; i++)
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        double Jxbar = 0;
        double Jybar = 0;
        double Jzbar = 0;
        for (int s=0; s<ns; s++) { 
          Jxbar += Jxs[s][i][j][k];
          Jybar += Jys[s][i][j][k];
          Jzbar += Jzs[s][i][j][k];
        }
        Jxbar /= grid->getInvVOLn(i, j, k);
        Jybar /= grid->getInvVOLn(i, j, k);
        Jzbar /= grid->getInvVOLn(i, j, k);
        double JE_part   = Jxbar*Exth[i][j][k] + Jybar*Eyth[i][j][k] + Jzbar*Ezth[i][j][k];
        double mex, mey, mez;
        productMassE(&mex, &mey, &mez, Exth, Eyth, Ezth, i, j, k);
        double Jxfield = Jxh[i][j][k] + mex;
        double Jyfield = Jyh[i][j][k] + mey;
        double Jzfield = Jzh[i][j][k] + mez;
        //double Jxfield = mex;
        //double Jyfield = mey;
        //double Jzfield = mez;
        double JE_fields =  Jxfield*Exth[i][j][k] + Jyfield*Eyth[i][j][k] + Jzfield*Ezth[i][j][k];
        double err = fabs(JE_part - JE_fields);
        if (err > maxErr) maxErr = err;
        fprintf(fd, "%2d %2d %2d %13.6e %13.6e %13.6e\n", i,j,k,JE_part, JE_fields, err);
        if (i < nxn-2 && j<nyn-2 && k < nzn-2) delE += JE_fields;
      }
  fclose(fd);

  delE *= dt;
  fprintf(stderr,"TEST JdotE :rank [%d] Max. error %13.6e\n", rank, maxErr);
  fprintf(stderr,"     Delta Energy: %13.6e(JdotE)  %13.6e(EM) %13.6e (Kinetic)   (%13.6e) \n", delE, (Een+Ben - Eenergy - Benergy), deltaK, (Een+Ben - Eenergy - Benergy)+deltaK);

  getchar();
  Eenergy = Een;
  Benergy = Ben;
}

void EMfields3D::PoyntingFluxTest()
{

}

void EMfields3D::updateLagrangianFrame(VirtualTopology3D * vct, Collective *col)
{
  int i0, i1, j0, j1, k0, k1;

  double xn = u_bulk*dt;
  double yn = v_bulk*dt;
  double zn = w_bulk*dt;
 
  if (u_bulk >= 0) {
    i0 = 0;
    i1 = 1;
  } else {
    i0 = -1;
    i1 = 0;
    xn += dx;
  }

  if (v_bulk >= 0) {
    j0 = 0;
    j1 = 1;
  } else {
    j0 = -1;
    j1 = 0;
    yn += dy;
  }

  if (w_bulk >= 0) {
    k0 = 0;
    k1 = 1;
  } else {
    k0 = -1;
    k1 = 0;
    zn += dz;
  }

  xn /= dx;
  yn /= dy;
  zn /= dz;

  double weight[2][2][2];
  double xi[2], eta[2], zeta[2];
  xi[0]   = xn;
  eta[0]  = yn;
  zeta[0] = zn;
  xi[1]   = 1 - xn;
  eta[1]  = 1 - yn;
  zeta[1] = 1 - zn;

/*
  for (int k=0; k<nzn; k++) 
    for (int j=0; j<nyn; j++) 
      for (int i=0; i<nxn; i++) {
        double x = i*dx;
        double y = j*dy;
        double z = k*dz;
        Ex[i][j][k] = cos(0.5*x) + 1.0/(y +1);
        Ey[i][j][k] = 1.0/(y +1);
        Ez[i][j][k] = 1.0/(z +1);
  }
  printf("weights: %g %g    -    %g %g    -     %g %g \n",xi[0], xi[1], eta[0], eta[1], zeta[0], zeta[1]);
*/

  for (int i=1; i<nxn-1; i++) 
    for (int j=1; j<nyn-1; j++) 
      for (int k=1; k<nzn-1; k++) {
        tempX[i][j][k] = xi[0]*eta[0]*zeta[0]*Ex[i+i1][j+j1][k+k1] 
                       + xi[1]*eta[0]*zeta[0]*Ex[i+i0][j+j1][k+k1]
                       + xi[0]*eta[1]*zeta[0]*Ex[i+i1][j+j0][k+k1] 
                       + xi[1]*eta[1]*zeta[0]*Ex[i+i0][j+j0][k+k1]
                       + xi[0]*eta[0]*zeta[1]*Ex[i+i1][j+j1][k+k0] 
                       + xi[1]*eta[0]*zeta[1]*Ex[i+i0][j+j1][k+k0]
                       + xi[0]*eta[1]*zeta[1]*Ex[i+i1][j+j0][k+k0] 
                       + xi[1]*eta[1]*zeta[1]*Ex[i+i0][j+j0][k+k0];

        tempY[i][j][k] = xi[0]*eta[0]*zeta[0]*Ey[i+i1][j+j1][k+k1] 
                       + xi[1]*eta[0]*zeta[0]*Ey[i+i0][j+j1][k+k1]
                       + xi[0]*eta[1]*zeta[0]*Ey[i+i1][j+j0][k+k1] 
                       + xi[1]*eta[1]*zeta[0]*Ey[i+i0][j+j0][k+k1]
                       + xi[0]*eta[0]*zeta[1]*Ey[i+i1][j+j1][k+k0] 
                       + xi[1]*eta[0]*zeta[1]*Ey[i+i0][j+j1][k+k0]
                       + xi[0]*eta[1]*zeta[1]*Ey[i+i1][j+j0][k+k0] 
                       + xi[1]*eta[1]*zeta[1]*Ey[i+i0][j+j0][k+k0];

        tempZ[i][j][k] = xi[0]*eta[0]*zeta[0]*Ez[i+i1][j+j1][k+k1] 
                       + xi[1]*eta[0]*zeta[0]*Ez[i+i0][j+j1][k+k1]
                       + xi[0]*eta[1]*zeta[0]*Ez[i+i1][j+j0][k+k1] 
                       + xi[1]*eta[1]*zeta[0]*Ez[i+i0][j+j0][k+k1]
                       + xi[0]*eta[0]*zeta[1]*Ez[i+i1][j+j1][k+k0] 
                       + xi[1]*eta[0]*zeta[1]*Ez[i+i0][j+j1][k+k0]
                       + xi[0]*eta[1]*zeta[1]*Ez[i+i1][j+j0][k+k0] 
                       + xi[1]*eta[1]*zeta[1]*Ez[i+i0][j+j0][k+k0];

        temp2X[i][j][k] = xi[0]*eta[0]*zeta[0]*Bxn[i+i1][j+j1][k+k1] 
                        + xi[1]*eta[0]*zeta[0]*Bxn[i+i0][j+j1][k+k1]
                        + xi[0]*eta[1]*zeta[0]*Bxn[i+i1][j+j0][k+k1] 
                        + xi[1]*eta[1]*zeta[0]*Bxn[i+i0][j+j0][k+k1]
                        + xi[0]*eta[0]*zeta[1]*Bxn[i+i1][j+j1][k+k0] 
                        + xi[1]*eta[0]*zeta[1]*Bxn[i+i0][j+j1][k+k0]
                        + xi[0]*eta[1]*zeta[1]*Bxn[i+i1][j+j0][k+k0] 
                        + xi[1]*eta[1]*zeta[1]*Bxn[i+i0][j+j0][k+k0];

        temp2Y[i][j][k] = xi[0]*eta[0]*zeta[0]*Byn[i+i1][j+j1][k+k1] 
                        + xi[1]*eta[0]*zeta[0]*Byn[i+i0][j+j1][k+k1]
                        + xi[0]*eta[1]*zeta[0]*Byn[i+i1][j+j0][k+k1] 
                        + xi[1]*eta[1]*zeta[0]*Byn[i+i0][j+j0][k+k1]
                        + xi[0]*eta[0]*zeta[1]*Byn[i+i1][j+j1][k+k0] 
                        + xi[1]*eta[0]*zeta[1]*Byn[i+i0][j+j1][k+k0]
                        + xi[0]*eta[1]*zeta[1]*Byn[i+i1][j+j0][k+k0] 
                        + xi[1]*eta[1]*zeta[1]*Byn[i+i0][j+j0][k+k0];

        temp2Z[i][j][k] = xi[0]*eta[0]*zeta[0]*Bzn[i+i1][j+j1][k+k1] 
                        + xi[1]*eta[0]*zeta[0]*Bzn[i+i0][j+j1][k+k1]
                        + xi[0]*eta[1]*zeta[0]*Bzn[i+i1][j+j0][k+k1] 
                        + xi[1]*eta[1]*zeta[0]*Bzn[i+i0][j+j0][k+k1]
                        + xi[0]*eta[0]*zeta[1]*Bzn[i+i1][j+j1][k+k0] 
                        + xi[1]*eta[0]*zeta[1]*Bzn[i+i0][j+j1][k+k0]
                        + xi[0]*eta[1]*zeta[1]*Bzn[i+i1][j+j0][k+k0] 
                        + xi[1]*eta[1]*zeta[1]*Bzn[i+i0][j+j0][k+k0];

  
      }

  communicateNodeBC(nxn, nyn, nzn, tempX, col->bcEx[0],col->bcEx[1],col->bcEx[2],col->bcEx[3],col->bcEx[4],col->bcEx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, tempY, col->bcEy[0],col->bcEy[1],col->bcEy[2],col->bcEy[3],col->bcEy[4],col->bcEy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, tempZ, col->bcEz[0],col->bcEz[1],col->bcEz[2],col->bcEz[3],col->bcEz[4],col->bcEz[5], vct);
  communicateNodeBC(nxn, nyn, nzn, temp2X, col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(nxn, nyn, nzn, temp2Y, col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(nxn, nyn, nzn, temp2Z, col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  for (int k=0; k<nzn; k++) 
    for (int j=0; j<nyn; j++) 
      for (int i=0; i<nxn; i++) {
        Ex[i][j][k] = tempX[i][j][k];
        Ey[i][j][k] = tempY[i][j][k];
        Ez[i][j][k] = tempZ[i][j][k];
        Bxn[i][j][k] = temp2X[i][j][k];
        Byn[i][j][k] = temp2Y[i][j][k];
        Bzn[i][j][k] = temp2Z[i][j][k];
      }

#ifdef DEBUG_LAGRANGIAN
  FILE* fd1 = fopen("original.txt","w");
  for (int k=1; k<nzn-1; k++) {
    for (int j=1; j<nyn-1; j++) {
      for (int i=1; i<nxn-1; i++) {
        double x = i*dx;
        double y = j*dy;
        double z = k*dz;
        fprintf(fd1,"%.3f %.3f %.3f %13.6e %13.6e %13.6e\n",x,y,z,Ex[i][j][k], Ey[i][j][k], Ez[i][j][k]);
      }
    }
    fprintf(fd1,"\n\n\n");
  }
  fclose(fd1);


  FILE* fd2 = fopen("shifted.txt","w");
  for (int k=1; k<nzn-1; k++) {
    for (int j=1; j<nyn-1; j++) {
      for (int i=1; i<nxn-1; i++) {
        double x = i*dx + u_bulk*dt;
        double y = j*dy + v_bulk*dt;
        double z = k*dz + w_bulk*dt;
        fprintf(fd2,"%.3f %.3f %.3f %13.6e %13.6e %13.6e\n",x,y,z,tempX[i][j][k], tempY[i][j][k], tempZ[i][j][k]);
      }
    }
    fprintf(fd2,"\n\n\n");
  }
  fclose(fd2);
  getchar();
#endif
}

/*! Compute the external electric and magnetic field in the lagrangian frame */
void EMfields3D::setLagrangianFrame(VirtualTopology3D * vct, Collective *col, bool mode)
{
  double tmpEx, tmpEy, tmpEz, tmpBx, tmpBy, tmpBz;
  double vx, vy, vz; 
  double c = col->getC();
  double c2 = c*c;

  if (mode) {
    vx = u_bulk;
    vy = v_bulk;
    vz = w_bulk;
  } else {
    vx = -u_bulk;
    vy = -v_bulk;
    vz = -w_bulk;
  }

  for (int k=0; k<nzn; k++) 
    for (int j=0; j<nyn; j++) 
      for (int i=0; i<nxn; i++) {
        // External fields -> v_bulk
        tmpEx = Ex_ext[i][j][k];
        tmpEy = Ey_ext[i][j][k];
        tmpEz = Ez_ext[i][j][k];
        tmpBx = Bx_ext[i][j][k];
        tmpBy = By_ext[i][j][k];
        tmpBz = Bz_ext[i][j][k];

        Ex_ext[i][j][k] = tmpEx + (vy*tmpBz - vz*tmpBy);
        Ey_ext[i][j][k] = tmpEy + (vz*tmpBx - vx*tmpBz);
        Ez_ext[i][j][k] = tmpEz + (vx*tmpBy - vy*tmpBx);

        Bx_ext[i][j][k] = tmpBx - (vy*tmpEz - vz*tmpEy)/c2;
        By_ext[i][j][k] = tmpBy - (vz*tmpEx - vx*tmpEz)/c2;
        Bz_ext[i][j][k] = tmpBz - (vx*tmpEy - vy*tmpEx)/c2;

        // Local fields ->  v_bulk
        tmpEx = Ex[i][j][k];
        tmpEy = Ey[i][j][k];
        tmpEz = Ez[i][j][k];
        tmpBx = Bxn[i][j][k];
        tmpBy = Byn[i][j][k];
        tmpBz = Bzn[i][j][k];

        Ex[i][j][k] = tmpEx + (vy*tmpBz - vz*tmpBy);
        Ey[i][j][k] = tmpEy + (vz*tmpBx - vx*tmpBz);
        Ez[i][j][k] = tmpEz + (vx*tmpBy - vy*tmpBx);

        Bxn[i][j][k] = tmpBx - (vy*tmpEz - vz*tmpEy)/c2;
        Byn[i][j][k] = tmpBy - (vz*tmpEx - vx*tmpEz)/c2;
        Bzn[i][j][k] = tmpBz - (vx*tmpEy - vy*tmpEx)/c2;

      }
}

void EMfields3D::setZeroCurrent(VirtualTopology3D * vct)
{
  for (int k=1; k<nzn-1; k++) 
    for (int j=1; j<nyn-1; j++) 
      for (int i=1; i<nxn-1; i++) {
        Jx_ext[i][j][k] = 0;
        Jy_ext[i][j][k] = 0;
        Jz_ext[i][j][k] = 0;
        for (int s=0; s<ns; s++) {
          Jx_ext[i][j][k] -= Jxs[s][i][j][k];
          Jy_ext[i][j][k] -= Jys[s][i][j][k];
          Jz_ext[i][j][k] -= Jzs[s][i][j][k];
        }
      }
  communicateInterp(nxn, nyn, nzn, Jx_ext, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jy_ext, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, Jz_ext, 0, 0, 0, 0, 0, 0, vct);
  communicateNode_P(nxn, nyn, nzn, Jx_ext, vct);
  communicateNode_P(nxn, nyn, nzn, Jy_ext, vct);
  communicateNode_P(nxn, nyn, nzn, Jz_ext, vct);

  /*
  for (int k=1; k<nzn-1; k++) 
    for (int j=1; j<nyn-1; j++) 
      for (int i=1; i<nxn-1; i++) 
        printf("%g %g %g \n", Jx_ext[i][j][k], Jy_ext[i][j][k], Jz_ext[i][j][k]);

  getchar();
  */
}

void EMfields3D::SetLambda(Grid *grid){

  for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){

        double x = grid->getXN(i,j,k);
        double y = grid->getYN(i,j,k);
        double z = grid->getZN(i,j,k);

        double xmin_r = Lx - 75.0 * dx;
        double xmax_r = Lx - 25.0  * dx;

        Lambda[i][j][k] = 0.0;

        if (x > xmin_r) {
          if (x < xmax_r) Lambda[i][j][k] = ((x - xmin_r) /  (xmax_r - xmin_r)) * 4.0 * M_PI / dx;
          else            Lambda[i][j][k] = 4.0 * M_PI / dx;
        }

      }
    }
  }

}

double*** EMfields3D::GetLambda(){
  return(Lambda);
}

// Compute the charge density from the 
void EMfields3D::chargeFromCurrent(Grid* grid, VirtualTopology3D * vct){
  for (int is=0; is<ns; is++) {
    for (int i=1; i<nxn-1; i++)
      for (int j=1; j<nyn-1; j++) 
        for (int k=1; k<nzn-1; k++) {
          double mex, mey, mez;
          productMassE(&mex, &mey, &mez, Exth, Eyth, Ezth, i, j, k);
          tempX[i][j][k]  = (Jxhs[is][i][j][k] + mex)*grid->getInvVOLn(i,j,k);
          tempY[i][j][k]  = (Jyhs[is][i][j][k] + mey)*grid->getInvVOLn(i,j,k);
          tempZ[i][j][k]  = (Jzhs[is][i][j][k] + mez)*grid->getInvVOLn(i,j,k);
        }
    if (Cylindrical == 1) grid->divN2C_cyl(tempC, tempX, tempY, tempZ);
    else grid->divN2C(tempC, tempX, tempY, tempZ);


    for (int i=1; i<nxc-1; i++)
      for (int j=1; j<nyc-1; j++) 
        for (int k=1; k<nzc-1; k++) {
          rhocs[is][i][j][k] -= dt* tempC[i][j][k];
        }

    communicateCenterBC(nxc, nyc, nzc, rhocs[is], 2, 2, 2, 2, 2, 2, vct);
    grid->interpC2N(rhons[is], rhocs[is]);
  }
}
 
// Compute the divergence of E
void EMfields3D::divergenceOfE(Grid* grid, VirtualTopology3D * vct) {
  if (Cylindrical == 1) grid->divN2C_cyl(tempC, Ex, Ey, Ez);
  else grid->divN2C(tempC, Ex, Ey, Ez);

  communicateCenterBC(nxc, nyc, nzc, tempC, 2, 2, 2, 2, 2, 2, vct);
  grid->interpC2N(divE, tempC); 
}
