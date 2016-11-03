
#include "ECsim.h"
#include "MyClock.h"

using namespace ECsim;

extern MyClock *clocks;

int c_Solver::Init(int argc, char **argv) {

  // initialize MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;

  clocks = new MyClock(15);
  clocks->start(10);
  
  #ifdef __PROFILING__
  clocks->start(0);
  #endif

  col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  verbose = col->getVerbose();
  restart_cycle = col->getRestartOutputCycle();
  SaveDirName = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart = col->getRestart_status();
  ns = col->getNs();            // get the number of particle species involved in simulation
  first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  // initialize the virtual cartesian topology 
  vct = new VCtopology3D(col);
  // Check if we can map the processes into a matrix ordering defined in Collective.cpp
  if (nprocs != vct->getNprocs()) {
    if (myrank == 0) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << "x" << vct->getZLEN() << " matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()" << endl;
      mpi->finalize_mpi();
      return (1);
    }
  }
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology(MPI_COMM_WORLD);
  // initialize the central cell index
 
  nx0 = col->getNxc() / vct->getXLEN(); // get the number of cells in x for each processor
  ny0 = col->getNyc() / vct->getYLEN(); // get the number of cells in y for each processor
  nz0 = col->getNzc() / vct->getZLEN(); // get the number of cells in z for each processor
  // Print the initial settings to stdout and a file
  if (myrank == 0) {
    mpi->Print();
    vct->Print();
    col->Print();
    col->save();
  }
  // Create the local grid
  MPI_Barrier(MPI_COMM_WORLD);
  grid = new Grid3DCU(col, vct);  // Create the local grid
  EMf = new EMfields3D(col, grid);  // Create Electromagnetic Fields Object
  
  if (col->getSolInit()) {
    /* -------------------------------------------- */
    /* If using parallel H5hut IO read initial file */
    /* -------------------------------------------- */
    ReadFieldsH5hut(ns, false, EMf, col, vct, grid);

  }
  else {
    /* --------------------------------------------------------- */
    /* If using 'default' IO initialize fields depending on case */
    /* --------------------------------------------------------- */
    if      (col->getCase()=="GEM")         EMf->initGEM(vct, grid,col);
    else if (col->getCase()=="DoubleGEM")   EMf->initDoubleGEM(vct, grid,col);
    else if (col->getCase()=="DoubleHarris")EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
    else if (col->getCase()=="Whistler")    EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);

    else {
      if (myrank==0) {
        cout << " =========================================================== " << endl;
        cout << " WARNING: The case '" << col->getCase() << "' was not recognized. " << endl;
        cout << "          Runing simulation with the default initialization. " << endl;
        cout << " =========================================================== " << endl;
      }
      EMf->init(vct,grid,col);
    }
  }

  // PETSc solver:
#ifdef __WITH_PETSC__
  petscSolver = new PetscSolver(EMf, grid, vct, col);
#endif

  // Allocation of particles
  part = new Particles3D[ns];
  for (int i = 0; i < ns; i++) {
    part[i].allocate(i, 0, col, vct, grid);
    if (col->getCase()=="DoubleGEM") part[i].MaxwellianDoubleGEM(grid, EMf, vct);
    else                             part[i].maxwellian(grid, EMf, vct);
  }


  if (col->getWriteMethod() == "default") {
    // Initialize the output (simulation results and restart file)
    // PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    // myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
    hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
    for (int i = 0; i < ns; ++i)
      hdf5_agent.set_simulation_pointers_part(&part[i]);
    output_mgr.push_back(&hdf5_agent);  // Add the HDF5 output agent to the Output Manager's list
    if (myrank == 0 & restart < 2) {
      hdf5_agent.open(SaveDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
      hdf5_agent.open(RestartDirName + "/settings.hdf");
      output_mgr.output("collective + total_topology + proc_topology", 0);
      hdf5_agent.close();
    }
    // Restart
    num_proc << myrank;
    if (restart == 0) {           // new simulation from input file
      hdf5_agent.open(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
    else {                        // restart append the results to the previous simulation 
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("proc_topology ", 0);
      hdf5_agent.close();
    }
  }


  Eenergy, Benergy, TOTenergy = 0.0, TOTmomentum = 0.0;
  Ke = new double[ns];
  momentum = new double[ns];
  cq = SaveDirName + "/ConservedQuantities.txt";
  
  #ifdef __PROFILING__
  clocks->stop(0);
  #endif
  return 0;
}

void c_Solver::GatherMoments(){

  #ifdef __PROFILING__
  clocks->start(1);
  #endif

  #ifdef __PROFILING__
  clocks->start(2);
  #endif
  EMf->centers2nodesB(vct, grid, col);
  EMf->setZeroDensities();                  // set to zero the densities
  #ifdef __PROFILING__
  clocks->stop(2);
  #endif

  for (int i = 0; i < ns; i++)
    part[i].computeMoments(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)

  EMf->sumOverSpecies (vct);

  // Communicate everything:
  #ifdef __PROFILING__
  clocks->start(11);
  #endif
  EMf->communicateGhostP2G(vct);
  #ifdef __PROFILING__
  clocks->stop(11);
  #endif

  #ifdef __PROFILING__
  clocks->stop(1);
  #endif
}

void c_Solver::CalculateField() {

  #ifdef __PROFILING__
  clocks->start(7);
  #endif
  EMf->interpDensitiesN2C(vct, grid);       // calculate densities on centers from nodes
  MPI_Barrier(MPI_COMM_WORLD);

  // MAXWELL'S SOLVER
  #ifdef __WITH_PETSC__
    petscSolver->solveE();
  #else
    EMf->calculateEB(grid, vct, col);               // calculate the E field
  #endif

  #ifdef __PROFILING__
  clocks->stop(7);
  #endif

}

bool c_Solver::ParticlesMover() {

  #ifdef __PROFILING__
  clocks->start(8);
  #endif
  /*  -------------- */
  /*  Particle mover */
  /*  -------------- */

  for (int i = 0; i < ns; i++)  // move each species
  {
    mem_avail = part[i].mover_ECIMM(grid, vct, EMf); // use the Predictor Corrector scheme 
  }

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */
  InjectBoundaryParticles();

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  #ifdef __PROFILING__
  clocks->stop(8);
  #endif
  return (false);
}

void c_Solver::InjectBoundaryParticles(){

  for (int i=0; i < ns; i++) {
    if (col->getRHOinject(i)>0.0){
      mem_avail = part[i].particle_repopulator(grid,vct,EMf,i);
    }
  }
}

void c_Solver::WriteRestart(int cycle) {

  if (cycle % restart_cycle == 0 && cycle != first_cycle) {
    if (myrank == 0) 
      cout << " WARNING: Restart writing not implemented." << endl ;
  }
}

void c_Solver::WriteConserved(int cycle) {

  if (cycle == 0) {
    totEn0 = EMf->getEenergy() + EMf->getBenergy();
    for (int is = 0; is < ns; is++) {
      totEn0 +=  part[is].getKe();
    }
  }

  // write the conserved quantities
  if (cycle % col->getDiagnosticsOutputCycle() == 0) {
    Eenergy = EMf->getEenergy();
    Benergy = EMf->getBenergy();
    TOTenergy = 0.0;
    TOTmomentum = 0.0;
    for (int is = 0; is < ns; is++) {
      Ke[is] = part[is].getKe();
      TOTenergy += Ke[is];
      momentum[is] = part[is].getP();
      TOTmomentum += momentum[is];
    }
    if (myrank == 0) {
      FILE* fd = fopen(cq.c_str(),"a");
      fprintf(fd,"%4d %13.6e  %13.6e    %13.6e %13.6e %13.6e  %13.6e\n"
                ,cycle
                ,(Eenergy + Benergy + TOTenergy), fabs((Eenergy + Benergy + TOTenergy) - totEn0)
                ,Eenergy , Benergy, TOTenergy
                ,TOTmomentum);

      fclose(fd);
    }
  }
}


void c_Solver::WriteOutput(int cycle) {


  if (col->getWriteMethod() == "h5hut") {

    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */

    if (cycle%(col->getFieldOutputCycle())==0)        WriteFieldsH5hut(ns, grid, EMf,  col, vct, cycle);
    if (cycle%(col->getParticlesOutputCycle())==0 &&
        cycle!=col->getLast_cycle() && cycle!=0)      WritePartclH5hut(ns, grid, part, col, vct, cycle);

  }
  else
  {

    // OUTPUT to large file, called proc**
    if (cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("Eall + Ball + rhos + Jsall + pressure", cycle);
      // Pressure tensor is available
      hdf5_agent.close();
    }
    if (cycle % (col->getParticlesOutputCycle()) == 0 && col->getParticlesOutputCycle() != 1) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("position + velocity + q ", cycle, 1);
      hdf5_agent.close();
    }
  }
}


void c_Solver::WriteOutputASCII(int cycle) {
  if (cycle%(col->getFieldOutputCycle())!=0) return;


  char name[128];
  sprintf(name, "data/part_%d.txt",cycle);
  FILE* fdpart = fopen(name, "w");
  
  for (int i=0; i < ns; i++) {
    int n = part[i].getNOP();
    for (int j=0; j<n; j++) {
      double u = part[i].getU(j);
      double v = part[i].getV(j);
      double w = part[i].getW(j);
      double m = sqrt(u*u+v*v+w*w);
      fprintf(fdpart,"%13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                ,part[i].getX(j), part[i].getY(j), part[i].getZ(j), u,v,w);
    }
    fprintf(fdpart,"\n\n\n");
  }
  fclose(fdpart);

  sprintf(name, "data/rho_%d.txt",cycle);
  FILE* fdrho = fopen(name, "w");
  for (int i=1; i<grid->getNXC()-1; i++)
    for (int j=1; j<grid->getNYC()-1; j++)
      for (int k=1; k<grid->getNZC()-1; k++)
        fprintf(fdrho,"%13.6e %13.6e %13.6e  %13.6e\n"
                     ,grid->getDX()*(i+0.5), grid->getDY()*(j+0.5), grid->getDZ()*(k+0.5)
                     ,EMf->getRHOc(i,j,k));
  fclose(fdrho);

  sprintf(name,"data/fieldE_%d.txt",cycle);
  FILE* fdE = fopen(name,"w");
  sprintf(name,"data/fieldB_%d.txt",cycle);
  FILE* fdB = fopen(name,"w");

  int nxn = grid->getNXN();
  int nyn = grid->getNYN();
  int nzn = grid->getNZN();
  double dx = grid->getDX();
  double dy = grid->getDY();
  double dz = grid->getDZ();
  for (int k=1; k<nzn-1; k++) {
    for (int j=1; j<nyn-1; j++) {
      for (int i=1; i<nxn-1; i++) {
        double x = grid->getXstart() + i*dx; 
        double y = grid->getYstart() + j*dy; 
        double z = grid->getZstart() + k*dz; 
        double Ex = EMf->getEx(i,j,k);
        double Ey = EMf->getEy(i,j,k);
        double Ez = EMf->getEz(i,j,k);
        double Bx = EMf->getBx(i,j,k);
        double By = EMf->getBy(i,j,k);
        double Bz = EMf->getBz(i,j,k);
        fprintf(fdE,"%13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                    ,x,y,z,Ex,Ey,Ez);
        fprintf(fdB,"%13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                    ,x,y,z,Bx,By,Bz);
      }
    }
    fprintf(fdE,"\n\n\n");
    fprintf(fdB,"\n\n\n");
  }

  fclose(fdE);
  fclose(fdB);

}

void c_Solver::Finalize() {

  // stop profiling
  clocks->stop(10);

  // deallocate
  delete[]Ke;
  delete[]momentum;
#ifdef __WITH_PETSC__
  delete petscSolver;
#endif
  // close MPI
  mpi->finalize_mpi();
}
