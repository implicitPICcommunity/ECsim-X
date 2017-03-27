
#include "EMfields3D.h"

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

  double ut0 = col->getU0(0) - zeroCurrent*u_bulk;
  double vt0 = col->getV0(0) - zeroCurrent*v_bulk;
  double wt0 = col->getW0(0) - zeroCurrent*w_bulk;
 
  double valx = wt0*B0y-vt0*B0z;
  double valy = ut0*B0z-wt0*B0x;
  double valz = vt0*B0x-ut0*B0y;

  // X left
  if (vct->getXleft_neighbor() == MPI_PROC_NULL) {
    if (bcEMfaceXleft == 0) {
      for (int j=0; j<nyn; j++) 
        for (int k=0; k<nzn; k++) {
          SoX[1][j][k] = valx;
          SoY[1][j][k] = valy;
          SoZ[1][j][k] = valz;
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

/*
    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) 
        for (int c=0; c<NE_MASS; c++) {
          Mxx[c][1][j][k] = 0;
          Mxy[c][1][j][k] = 0;
          Mxz[c][1][j][k] = 0;
          Myx[c][1][j][k] = 0;
          Myy[c][1][j][k] = 0;
          Myz[c][1][j][k] = 0;
          Mzx[c][1][j][k] = 0;
          Mzy[c][1][j][k] = 0;
          Mzz[c][1][j][k] = 0;

    for (int j=0; j<nyn; j++) 
      for (int k=0; k<nzn; k++) 
        for (int c=0; c<NE_MASS; c++) {
          Mxx[c][nxn-2][j][k] = 0;
          Mxy[c][nxn-2][j][k] = 0;
          Mxz[c][nxn-2][j][k] = 0;
          Myx[c][nxn-2][j][k] = 0;
          Myy[c][nxn-2][j][k] = 0;
          Myz[c][nxn-2][j][k] = 0;
          Mzx[c][nxn-2][j][k] = 0;
          Mzy[c][nxn-2][j][k] = 0;
          Mzz[c][nxn-2][j][k] = 0;
        }
  }
*/
}

