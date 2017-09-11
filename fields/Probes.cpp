#include "Probes.h"

Probes::Probes(Collective* col, Grid* grid)
{
  // 9 values: Ex, Ey, Ez, Bx, By, Bz, Btx, Bty, Btz
  nval = 9;

  probesfile = col->getProbesFile().c_str();  // File with the position of the probes
  if (strcmp(probesfile,"no") == 0) {
    active = false;
    return ;
  }

  FILE* fd = fopen(probesfile, "r");
  if (fd != 0) {
    active = true;
    nptmax = 1000;   // Number of time step to store in memory
    npmax = 20;      // Initial allocation number of probes
    npt = 0;
    np = 0;
    cycle = 1;
    dt = col->getDt();
    id = (int*)    malloc(npmax*sizeof(int));
    x  = (double*) malloc(npmax*sizeof(double));
    y  = (double*) malloc(npmax*sizeof(double));
    z  = (double*) malloc(npmax*sizeof(double));
  
    sprintf(savedir,"%s", col->getSaveDirName().c_str());

    double xmin = grid->getXstart();
    double xmax = grid->getXend();
    double ymin = grid->getYstart();
    double ymax = grid->getYend();
    double zmin = grid->getZstart();
    double zmax = grid->getZend();

    size_t lenmax = 256;
    char* txt;
    char* saveptr;
    char* buffer = (char*) malloc(lenmax*sizeof(char));

    while (getline(&buffer, &lenmax, fd) != -1) {
      if (buffer[0] == '#' || buffer[0] == '\n') continue;
      txt = strtok_r(buffer, " ", &saveptr);
      int idg  = atoi(txt);
      txt = strtok_r(NULL, " ", &saveptr);
      double xg = atof(txt);
      txt = strtok_r(NULL, " ", &saveptr);
      double yg = atof(txt);
      txt = strtok_r(NULL, " ", &saveptr);
      double zg = atof(txt);
      if (xg >= xmin && xg < xmax && yg >= ymin && yg < ymax && zg >= zmin && zg < zmax) {
        if (np == npmax) {
          npmax += 20;
          id = (int*)    realloc(id, npmax*sizeof(int));
          x  = (double*) realloc(x, npmax*sizeof(double));
          y  = (double*) realloc(y, npmax*sizeof(double));
          z  = (double*) realloc(z, npmax*sizeof(double));
        }
        id[np] = idg;
        x[np] = xg;
        y[np] = yg;
        z[np] = zg;
        np++;
      }
    }
    fclose(fd);
    free(buffer);

    data = newArr3(double, np, nval, nptmax);
		W = newArr4(double, np, 2, 2, 2);
    indx  = (int*) malloc(np*sizeof(int));
    indy  = (int*) malloc(np*sizeof(int));
    indz  = (int*) malloc(np*sizeof(int));
    init(grid);
    writeHeader();
  } else {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("ERROR while creating probes: file %s cannot be read\n", probesfile);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
}

Probes::~Probes()
{
  if (active) {
    free(x);
    free(y);
    free(z);
    delArr3(data, np, nval);
    delArr4(W, np, 2, 2);
    free(data);
  }
}

void Probes::setValues(EMfields3D* EMf)
{
  if (!active) return;

  if (npt == nptmax) {
    writeFiles();
  }

  double*** Ex = EMf->getEx();
  double*** Ey = EMf->getEy();
  double*** Ez = EMf->getEz();
  double*** Bx = EMf->getBx();
  double*** By = EMf->getBy();
  double*** Bz = EMf->getBz();
  double*** Bx_ext = EMf->getBx_ext();
  double*** By_ext = EMf->getBy_ext();
  double*** Bz_ext = EMf->getBz_ext();
  for(int p=0; p<np; p++) {
      data[p][0][npt] = interpolate(p, Ex);
      data[p][1][npt] = interpolate(p, Ey);
      data[p][2][npt] = interpolate(p, Ez);
      data[p][3][npt] = interpolate(p, Bx);
      data[p][4][npt] = interpolate(p, By);
      data[p][5][npt] = interpolate(p, Bz);
      data[p][6][npt] = interpolate(p, Bx_ext);
      data[p][7][npt] = interpolate(p, By_ext);
      data[p][8][npt] = interpolate(p, Bz_ext);
  }
  npt++;
}

void Probes::writeHeader()
{
  char fileres[300]; 
  for (int p=0; p<np; p++) {
    sprintf(fileres,"%s/Probe%03d.txt",savedir,id[p]);
    FILE* fd = fopen(fileres,"a");
    fprintf(fd,"# Probe %d\n"
               "# Location: \n"
               "#  x = %.4f\n" 
               "#  y = %.4f\n" 
               "#  z = %.4f\n" 
               "#   t           Ex            Ey            Ez            Bx            By            Bz            Bx_ext        By_ext        Bz_ext   \n"
               "#----------------------------------------------------------------------------------------------------------------------------------------\n"
               ,p, x[p], y[p], z[p]);
    fclose(fd);
  }
}

void Probes::writeFiles()
{
  if (!active) return;
  char fileres[300]; 
  for (int p=0; p<np; p++) {
    sprintf(fileres,"%s/Probe%03d.txt",savedir,id[p]);
    FILE* fd = fopen(fileres,"a");
  	for (int t=0; t<npt; t++) {
      fprintf(fd,"%10.3e", dt*(cycle+t));
  	  for (int i=0; i<nval; i++) fprintf(fd," %13.6e", data[p][i][t]);
  		fprintf(fd,"\n");
    }
    fclose(fd);
  }
  cycle += npt;
  npt = 0;
}

void Probes::init(Grid* grid)
{
  for (int i=0; i<np; i++) {
    for (int j=0; j<nval; j++)
      for (int k=0; k<nptmax; k++)
      data[i][j][k] = 0;
  }

  // Compute the weights of each probe for interpolation
  double xi[2], eta[2], zeta[2];
  int ix,iy,iz;

  double dx = grid->getDX();
  double dy = grid->getDY();
  double dz = grid->getDZ();
  double invVOL = 1.0/(dx*dy*dz);
  double xstart = grid->getXstart();
  double ystart = grid->getYstart();
  double zstart = grid->getZstart();

  for (int p=0; p<np; p++) {
    int ix = 2 + int (floor((x[p] - xstart) / dx));
    int iy = 2 + int (floor((y[p] - ystart) / dy));
    int iz = 2 + int (floor((z[p] - zstart) / dz));
    xi[0]   = x[p] - grid->getXN(ix - 1, iy, iz);
    eta[0]  = y[p] - grid->getYN(ix, iy - 1, iz);
    zeta[0] = z[p] - grid->getZN(ix, iy, iz - 1);
    xi[1]   = grid->getXN(ix, iy, iz) - x[p];
    eta[1]  = grid->getYN(ix, iy, iz) - y[p];
    zeta[1] = grid->getZN(ix, iy, iz) - z[p];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++)
          W[p][i][j][k] = xi[i] * eta[j] * zeta[k] * invVOL;
    indx[p] = ix;
    indy[p] = iy;
    indz[p] = iz;
  }

}

double Probes::interpolate(int p, double*** field)
{
  double ret = 0.0;
  int ix = indx[p];
  int iy = indy[p];
  int iz = indz[p];
	for (int i=0; i<=1; i++)
		for (int j=0; j<=1; j++)
			for (int k=0; k<=1; k++) 
				ret += W[p][i][j][k] * field[ix-i][iy-j][iz-k];
  return ret;
}


void Probes::print()
{
  printf("np = %d\n", np);
  for (int i=0; i<np; i++) {
    printf("%.3f %.3f %.3f ", x[i], y[i], z[i]);
    for (int j=0; j<nval; j++) printf(" %.3f", data[i][j][0]);
    printf("\n");
  }
}

