
#include "EMfields3D.h"

#define NUM_NODE_BC 6
#define NUM_CELL_BC 5

void EMfields3D::fixBC_Image(VirtualTopology3D *vct, double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz) {

  // X left
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[1][j][k] = nEx[1][j][k];
          ImY[1][j][k] = nEy[1][j][k];
          ImZ[1][j][k] = nEz[1][j][k];
        }
    } else if (bcEMfaceXleft == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[1][j][k] = nEx[1][j][k];
          ImY[1][j][k] = nEy[1][j][k] - nEy[2][j][k];
          ImZ[1][j][k] = nEz[1][j][k];
        }
    } else if (bcEMfaceXleft == 2) {
      for (int i=0; i<NUM_NODE_BC; i++) 
        for (int j=0; j<nyn; j++) 
          for (int k=0; k<nzn; k++) {
            ImX[i][j][k] = nEx[i][j][k];
            ImY[i][j][k] = nEy[i][j][k];
            ImZ[i][j][k] = nEz[i][j][k];
          }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // X right
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXright == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[nxn-2][j][k] = nEx[nxn-2][j][k];
          ImY[nxn-2][j][k] = nEy[nxn-2][j][k];
          ImZ[nxn-2][j][k] = nEz[nxn-2][j][k];
        }
    } else if (bcEMfaceXright == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[nxn-2][j][k] = nEx[nxn-2][j][k];
          ImY[nxn-2][j][k] = nEy[nxn-2][j][k] - nEy[nxn-3][j][k];
          ImZ[nxn-2][j][k] = nEz[nxn-2][j][k];
        }
    } else if (bcEMfaceXright == 2) {
      for (int i=0; i<NUM_NODE_BC; i++) 
        for (int j=0; j<nyn; j++) 
          for (int k=0; k<nzn; k++) {
            ImX[nxn-1-i][j][k] = nEx[nxn-1-i][j][k];
            ImY[nxn-1-i][j][k] = nEy[nxn-1-i][j][k];
            ImZ[nxn-1-i][j][k] = nEz[nxn-1-i][j][k];
          }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y left
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYleft == 0) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          ImX[i][1][k] = nEx[i][1][k];
          ImY[i][1][k] = nEy[i][1][k];
          ImZ[i][1][k] = nEz[i][1][k];
        }
    } else if (bcEMfaceYleft == 1) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          ImX[i][1][k] = nEx[i][1][k];
          ImY[i][1][k] = nEy[i][1][k] - nEy[i][2][k];
          ImZ[i][1][k] = nEz[i][1][k];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y right
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYright == 0) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          ImX[i][nyn-2][k] = nEx[i][nyn-2][k];
          ImY[i][nyn-2][k] = nEy[i][nyn-2][k];
          ImZ[i][nyn-2][k] = nEz[i][nyn-2][k];
        }
    } else if (bcEMfaceYright == 1) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          ImX[i][nyn-2][k] = nEx[i][nyn-2][k];
          ImY[i][nyn-2][k] = nEy[i][nyn-2][k] - nEy[i][nyn-3][k];
          ImZ[i][nyn-2][k] = nEz[i][nyn-2][k];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z left
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZleft == 0) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          ImX[i][j][1] = nEx[i][j][1];
          ImY[i][j][1] = nEy[i][j][1];
          ImZ[i][j][1] = nEz[i][j][1];
        }
    } else if (bcEMfaceZleft == 1) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          ImX[i][j][1] = nEx[i][j][1];
          ImY[i][j][1] = nEy[i][j][1] - nEy[i][j][2];
          ImZ[i][j][1] = nEz[i][j][1];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z right
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZright == 0) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          ImX[i][j][nzn-2] = nEx[i][j][nzn-2];
          ImY[i][j][nzn-2] = nEy[i][j][nzn-2];
          ImZ[i][j][nzn-2] = nEz[i][j][nzn-2];
        }
    } else if (bcEMfaceZright == 1) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          ImX[i][j][nzn-2] = nEx[i][j][nzn-2];
          ImY[i][j][nzn-2] = nEy[i][j][nzn-2] - nEy[i][j][nzn-3];
          ImZ[i][j][nzn-2] = nEz[i][j][nzn-2];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

void EMfields3D::fixBC_Source(VirtualTopology3D *vct, Collective *col, double ***SoX, double ***SoY, double ***SoZ) {

  double ut0 = col->getU0(0);
  double vt0 = col->getV0(0);
  double wt0 = col->getW0(0);

  double tmpBx = B0x;
  double tmpBy = B0y;
  double tmpBz = B0z;
  double tmpEx = wt0*B0y-vt0*B0z;
  double tmpEy = ut0*B0z-wt0*B0x;
  double tmpEz = vt0*B0x-ut0*B0y;

  double valx, valy, valz;
  if (col->getLagrangian() == 1) {
    double vx = u_bulk;
    double vy = v_bulk;
    double vz = w_bulk;
    valx = tmpEx + (vy*tmpBz - vz*tmpBy);
    valy = tmpEy + (vz*tmpBx - vx*tmpBz);
    valz = tmpEz + (vx*tmpBy - vy*tmpBx);
  } else {
    valx = tmpEx;
    valy = tmpEy;
    valz = tmpEz;
  }

  // X left
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[1][j][k] = valx;
          SoY[1][j][k] = valy;
          SoZ[1][j][k] = valz;
        }
    } else if (bcEMfaceXleft == 2) {
      for (int i=0; i<NUM_NODE_BC; i++) 
        for (int j=0; j<nyn; j++) 
          for (int k=0; k<nzn; k++) {
            SoX[i][j][k] = valx;
            SoY[i][j][k] = valy;
            SoZ[i][j][k] = valz;
          }

    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // X right
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXright == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[nxn-2][j][k] = valx;
          SoY[nxn-2][j][k] = valy;
          SoZ[nxn-2][j][k] = valz;
        }
    } else if (bcEMfaceXright == 2) {
      for (int i=0; i<NUM_NODE_BC; i++) 
        for (int j=0; j<nyn; j++) 
          for (int k=0; k<nzn; k++) {
            SoX[nxn-1-i][j][k] = valx;
            SoY[nxn-1-i][j][k] = valy;
            SoZ[nxn-1-i][j][k] = valz;
          }

    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y left
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYleft == 0) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          SoX[i][1][k] = valx;
          SoY[i][1][k] = valy;
          SoZ[i][1][k] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y right
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYright == 0) {
      for (int i=0; i<nxn; i++) 
        for (int k=0; k<nzn; k++) {
          SoX[i][nyn-2][k] = valx;
          SoY[i][nyn-2][k] = valy;
          SoZ[i][nyn-2][k] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z left
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZleft == 0) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          SoX[i][j][1] = valx;
          SoY[i][j][1] = valy;
          SoZ[i][j][1] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z right
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZright == 0 || bcEMfaceZright == 1) {
      for (int i=0; i<nxn; i++) 
        for (int j=0; j<nyn; j++) {
          SoX[i][j][nzn-2] = valx;
          SoY[i][j][nzn-2] = valy;
          SoZ[i][j][nzn-2] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

}

/********************************************************************************************************
 * CYLINDRICAL COORDINATES
 * Fix the boundary conditions for r=0. BC are only needed for the electric field.

 */
void EMfields3D::fixBC_Cyl_Image(double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz) {

  if (xStart == 0) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[1][j][k] = nEx[1][j][k];
          ImY[1][j][k] = nEy[1][j][k];
          ImZ[1][j][k] = nEz[1][j][k];
        }
    } else if (bcEMfaceXleft == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[1][j][k] = nEx[1][j][k] - nEx[2][j][k];
          ImY[1][j][k] = nEy[1][j][k] - nEy[2][j][k];
          ImZ[1][j][k] = nEz[1][j][k] - nEz[2][j][k];
        }
    } else if (bcEMfaceXleft == 2) {
      for (int j=1; j<nyn; j++) 
        for (int k=1; k<nzn; k++) {
          ImX[1][j][k] = nEx[1][j][k] - nEx[2][j][k];
          ImY[1][j][k] = nEy[1][j][k] - (4*dt*th/dx)*(cBz[1][j-1][k] + cBz[1][j][k])*0.5;
          ImZ[1][j][k] = nEz[1][j][k] - 0.5*(2*dt*th/dx)*(cBy[1][j-1][k] + cBy[1][j][k]) + 0.5*(dt*th/dy)*(cBx[1][j][k] - cBx[1][j-1][k]);
        }
    }
  }

  if (xEnd == Lx) {
    if (bcEMfaceXright == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[nxn-2][j][k] = nEx[nxn-2][j][k];
          ImY[nxn-2][j][k] = nEy[nxn-2][j][k];
          ImZ[nxn-2][j][k] = nEz[nxn-2][j][k];
        }
    } else if (bcEMfaceXright == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          ImX[nxn-2][j][k] = nEx[nxn-2][j][k];
          ImY[nxn-2][j][k] = nEy[nxn-2][j][k]-nEy[nxn-3][j][k];
          ImZ[nxn-2][j][k] = nEz[nxn-2][j][k]-nEz[nxn-3][j][k];
        }
    }
  }
}

void EMfields3D::fixBC_Cyl_Source(double ***SoX, double ***SoY, double ***SoZ) 
{
  if (xStart == 0) {
    if (bcEMfaceXleft == 0 || bcEMfaceXleft == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[1][j][k] = 0.0;
          SoY[1][j][k] = 0.0;
          SoZ[1][j][k] = 0.0;
        }
    } else if (bcEMfaceXleft == 2) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[1][j][k] = 0;
          SoY[1][j][k] = Ey[1][j][k];
          SoZ[1][j][k] = Ez[1][j][k];
        }
    }
  }

  if (xEnd == Lx) {
    if (bcEMfaceXright == 0 || bcEMfaceXright == 1) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[nxn-2][j][k] = 0.0;
          SoY[nxn-2][j][k] = 0.0;
          SoZ[nxn-2][j][k] = 0.0;
        }
    }
  }
}

void EMfields3D::fixBC_ImageC(VirtualTopology3D *vct, double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz) {

  // X left
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          ImX[1][j][k] = cBx[1][j][k];
          ImY[1][j][k] = cBy[1][j][k];
          ImZ[1][j][k] = cBz[1][j][k];
        }
    } else if (bcEMfaceXleft == 1) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          ImX[1][j][k] = cBx[1][j][k];
          ImY[1][j][k] = cBy[1][j][k] - cBy[2][j][k];
          ImZ[1][j][k] = cBz[1][j][k];
        }
    } else if (bcEMfaceXleft == 2) {
      for (int i=0; i<NUM_CELL_BC; i++) 
        for (int j=0; j<nyc; j++) 
          for (int k=0; k<nzc; k++) {
            ImX[i][j][k] = cBx[i][j][k];
            ImY[i][j][k] = cBy[i][j][k];
            ImZ[i][j][k] = cBz[i][j][k];
          }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // X right
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXright == 0) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          ImX[nxc-2][j][k] = cBx[nxc-2][j][k];
          ImY[nxc-2][j][k] = cBy[nxc-2][j][k];
          ImZ[nxc-2][j][k] = cBz[nxc-2][j][k];
        }
    } else if (bcEMfaceXright == 1) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          ImX[nxc-2][j][k] = cBx[nxc-2][j][k];
          ImY[nxc-2][j][k] = cBy[nxc-2][j][k] - cBy[nxc-3][j][k];
          ImZ[nxc-2][j][k] = cBz[nxc-2][j][k];
        }
    } else if (bcEMfaceXright == 2) {
      for (int i=0; i<NUM_CELL_BC; i++) 
        for (int j=0; j<nyc; j++) 
          for (int k=0; k<nzc; k++) {
            ImX[nxc-1-i][j][k] = cBx[nxc-1-i][j][k];
            ImY[nxc-1-i][j][k] = cBy[nxc-1-i][j][k];
            ImZ[nxc-1-i][j][k] = cBz[nxc-1-i][j][k];
          }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y left
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYleft == 0) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          ImX[i][1][k] = cBx[i][1][k];
          ImY[i][1][k] = cBy[i][1][k];
          ImZ[i][1][k] = cBz[i][1][k];
        }
    } else if (bcEMfaceYleft == 1) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          ImX[i][1][k] = cBx[i][1][k];
          ImY[i][1][k] = cBy[i][1][k] - cBy[i][2][k];
          ImZ[i][1][k] = cBz[i][1][k];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y right
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYright == 0) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          ImX[i][nyc-2][k] = cBx[i][nyc-2][k];
          ImY[i][nyc-2][k] = cBy[i][nyc-2][k];
          ImZ[i][nyc-2][k] = cBz[i][nyc-2][k];
        }
    } else if (bcEMfaceYright == 1) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          ImX[i][nyc-2][k] = cBx[i][nyc-2][k];
          ImY[i][nyc-2][k] = cBy[i][nyc-2][k] - cBy[i][nyc-3][k];
          ImZ[i][nyc-2][k] = cBz[i][nyc-2][k];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z left
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZleft == 0) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          ImX[i][j][1] = cBx[i][j][1];
          ImY[i][j][1] = cBy[i][j][1];
          ImZ[i][j][1] = cBz[i][j][1];
        }
    } else if (bcEMfaceZleft == 1) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          ImX[i][j][1] = cBx[i][j][1];
          ImY[i][j][1] = cBy[i][j][1] - cBy[i][j][2];
          ImZ[i][j][1] = cBz[i][j][1];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z right
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZright == 0) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          ImX[i][j][nzc-2] = cBx[i][j][nzc-2];
          ImY[i][j][nzc-2] = cBy[i][j][nzc-2];
          ImZ[i][j][nzc-2] = cBz[i][j][nzc-2];
        }
    } else if (bcEMfaceZright == 1) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          ImX[i][j][nzc-2] = cBx[i][j][nzc-2];
          ImY[i][j][nzc-2] = cBy[i][j][nzc-2] - cBy[i][j][nzc-3];
          ImZ[i][j][nzc-2] = cBz[i][j][nzc-2];
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

void EMfields3D::fixBC_SourceC(VirtualTopology3D *vct, Collective *col, double ***SoX, double ***SoY, double ***SoZ) {

  double ut0 = col->getU0(0);
  double vt0 = col->getV0(0);
  double wt0 = col->getW0(0);

  double tmpBx = B0x;
  double tmpBy = B0y;
  double tmpBz = B0z;
  double tmpEx = wt0*B0y-vt0*B0z;
  double tmpEy = ut0*B0z-wt0*B0x;
  double tmpEz = vt0*B0x-ut0*B0y;

  double valx, valy, valz;
  if (col->getLagrangian() == 1) {
    double vx = u_bulk;
    double vy = v_bulk;
    double vz = w_bulk;
    valx = tmpBx - (vy*tmpEz - vz*tmpEy);
    valy = tmpBy - (vz*tmpEx - vx*tmpEz);
    valz = tmpBz - (vx*tmpEy - vy*tmpEx);
  } else {
    valx = tmpBx;
    valy = tmpBy;
    valz = tmpBz;
  }

  // X left
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          SoX[1][j][k] = valx;
          SoY[1][j][k] = valy;
          SoZ[1][j][k] = valz;
        }
    } else if (bcEMfaceXleft == 2) {
      for (int i=0; i<NUM_CELL_BC; i++) 
        for (int j=0; j<nyc; j++) 
          for (int k=0; k<nzc; k++) {
            SoX[i][j][k] = valx;
            SoY[i][j][k] = valy;
            SoZ[i][j][k] = valz;
          }

    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // X right
  if (vct->getXright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXright == 0) {
      for (int j=0; j<nyc; j++) 
        for (int k=0; k<nzc; k++) {
          SoX[nxc-2][j][k] = valx;
          SoY[nxc-2][j][k] = valy;
          SoZ[nxc-2][j][k] = valz;
        }
    } else if (bcEMfaceXright == 2) {
      for (int i=0; i<NUM_CELL_BC; i++) 
        for (int j=0; j<nyc; j++) 
          for (int k=0; k<nzc; k++) {
            SoX[nxc-1-i][j][k] = valx;
            SoY[nxc-1-i][j][k] = valy;
            SoZ[nxc-1-i][j][k] = valz;
          }

      } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceXright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y left
  if (vct->getYleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYleft == 0) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          SoX[i][1][k] = valx;
          SoY[i][1][k] = valy;
          SoZ[i][1][k] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Y right
  if (vct->getYright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceYright == 0) {
      for (int i=0; i<nxc; i++) 
        for (int k=0; k<nzc; k++) {
          SoX[i][nyc-2][k] = valx;
          SoY[i][nyc-2][k] = valy;
          SoZ[i][nyc-2][k] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceYright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z left
  if (vct->getZleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZleft == 0) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          SoX[i][j][1] = valx;
          SoY[i][j][1] = valy;
          SoZ[i][j][1] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZleft);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  // Z right
  if (vct->getZright_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceZright == 0 || bcEMfaceZright == 1) {
      for (int i=0; i<nxc; i++) 
        for (int j=0; j<nyc; j++) {
          SoX[i][j][nzc-2] = valx;
          SoY[i][j][nzc-2] = valy;
          SoZ[i][j][nzc-2] = valz;
        }
    } else {
      fprintf(stderr,"ERROR: Boundary contitions (%d) not implemented for the fields\n", bcEMfaceZright);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

/* ====================================== */
/* FIX SOURCE AND IMAGE INSIDE THE PLANET */
/* ====================================== */

void EMfields3D::fixBC_PlanetImage(VirtualTopology3D *vct, Collective *col, Grid *grid, double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz) {

  double xc = col->getx_center();
  double yc = col->gety_center();
  double zc = col->getz_center();
  double L2 = col->getL_square() * col->getL_square();

  bool twodim = false;
  if ( nzn == 4) twodim = true;

  for (int i=0; i<nxn; i++) 
    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) {
        double x = grid->getXN(i,j,k);
        double y = grid->getYN(i,j,k);
        double z = grid->getZN(i,j,k);
        double rx = x-xc;
        double ry = y-yc;
        double rz = z-zc;
        if (twodim) rz = 0.0;
  
        double r = rx*rx + ry*ry + rz*rz;

        if (r<L2) {
          ImX[i][j][k] = nEx[i][j][k];
          ImY[i][j][k] = nEy[i][j][k];
          ImZ[i][j][k] = nEz[i][j][k];
        }
      }

}

void EMfields3D::fixBC_PlanetSource(VirtualTopology3D *vct, Collective *col, Grid *grid, double ***SoX, double ***SoY, double ***SoZ) {

  double ut0 = col->getU0(0);
  double vt0 = col->getV0(0);
  double wt0 = col->getW0(0);

  double tmpBx = B0x;
  double tmpBy = B0y;
  double tmpBz = B0z;
  double tmpEx = wt0*B0y-vt0*B0z;
  double tmpEy = ut0*B0z-wt0*B0x;
  double tmpEz = vt0*B0x-ut0*B0y;

  double valx, valy, valz;
  if (col->getLagrangian() == 1) {
    double vx = u_bulk;
    double vy = v_bulk;
    double vz = w_bulk;
    valx = tmpEx + (vy*tmpBz - vz*tmpBy);
    valy = tmpEy + (vz*tmpBx - vx*tmpBz);
    valz = tmpEz + (vx*tmpBy - vy*tmpBx);
  } else {
    valx = tmpEx;
    valy = tmpEy;
    valz = tmpEz;
  }

  double xc = col->getx_center();
  double yc = col->gety_center();
  double zc = col->getz_center();
  double L2 = col->getL_square() * col->getL_square();

  bool twodim = false;
  if ( nzn == 4) twodim = true;

  for (int i=0; i<nxn; i++) 
    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) {
        double x  = grid->getXN(i,j,k);
        double y  = grid->getYN(i,j,k);
        double z  = grid->getZN(i,j,k);
        double rx = x-xc;
        double ry = y-yc;
        double rz = z-zc;
        if (twodim) rz = 0.0;
  
        double r = rx*rx + ry*ry + rz*rz;

        if (r<L2) {
          SoX[i][j][k] = valx;
          SoY[i][j][k] = valy;
          SoZ[i][j][k] = valz;
        }
      }
}

void EMfields3D::fixBC_PlanetImageC(VirtualTopology3D *vct, Collective *col, Grid *grid, double ***ImX, double ***ImY, double ***ImZ, double ***nEx, double ***nEy, double ***nEz, double ***cBx, double ***cBy, double ***cBz) {

  double xc = col->getx_center();
  double yc = col->gety_center();
  double zc = col->getz_center();
  double L2 = col->getL_square() * col->getL_square();

  bool twodim = false;
  if ( nzn == 4) twodim = true;

  for (int i=0; i<nxc; i++) 
    for (int j=0; j<nyc; j++) 
      for (int k=0; k<nzc; k++) {
        double x = grid->getXC(i,j,k);
        double y = grid->getYC(i,j,k);
        double z = grid->getZC(i,j,k);
        double rx = x-xc;
        double ry = y-yc;
        double rz = z-zc;
        if (twodim) rz = 0.0;
  
        double r = rx*rx + ry*ry + rz*rz;

        if (r<L2) {
          ImX[i][j][k] = cBx[i][j][k];
          ImY[i][j][k] = cBy[i][j][k];
          ImZ[i][j][k] = cBz[i][j][k];
        }
      }
}

void EMfields3D::fixBC_PlanetSourceC(VirtualTopology3D *vct, Collective *col, Grid *grid, double ***SoX, double ***SoY, double ***SoZ) {

  double ut0 = col->getU0(0);
  double vt0 = col->getV0(0);
  double wt0 = col->getW0(0);

  double tmpBx = B0x;
  double tmpBy = B0y;
  double tmpBz = B0z;
  double tmpEx = wt0*B0y-vt0*B0z;
  double tmpEy = ut0*B0z-wt0*B0x;
  double tmpEz = vt0*B0x-ut0*B0y;

  double valx, valy, valz;
  if (col->getLagrangian() == 1) {
    double vx = u_bulk;
    double vy = v_bulk;
    double vz = w_bulk;
    valx = tmpBx - (vy*tmpEz - vz*tmpEy);
    valy = tmpBy - (vz*tmpEx - vx*tmpEz);
    valz = tmpBz - (vx*tmpEy - vy*tmpEx);
  } else {
    valx = tmpBx;
    valy = tmpBy;
    valz = tmpBz;
  }

  double xc = col->getx_center();
  double yc = col->gety_center();
  double zc = col->getz_center();
  double L2 = col->getL_square() * col->getL_square();

  bool twodim = false;
  if ( nzn == 4) twodim = true;

  for (int i=0; i<nxc; i++) 
    for (int j=0; j<nyc; j++) 
      for (int k=0; k<nzc; k++) {
        double x  = grid->getXC(i,j,k);
        double y  = grid->getYC(i,j,k);
        double z  = grid->getZC(i,j,k);
        double rx = x-xc;
        double ry = y-yc;
        double rz = z-zc;
        if (twodim) rz = 0.0;
  
        double r = rx*rx + ry*ry + rz*rz;

        if (r<L2) {
          SoX[i][j][k] = valx;
          SoY[i][j][k] = valy;
          SoZ[i][j][k] = valz;
        }
      }
}

