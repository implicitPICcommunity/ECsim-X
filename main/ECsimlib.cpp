
#include "ECsim.h"
#include "MyClock.h"
#include "Memory.h"
#include "Probes.h"
#include "FilePartID.h"

using namespace ECsim;

extern MyClock *clocks;

int c_Solver::Init(int argc, char **argv) {

  // initialize MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;

  clocks = new MyClock(16);
  clocks->start(0);
  clocks->start(15);

  col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  verbose = col->getVerbose();
  restart_cycle = col->getRestartOutputCycle();
  SaveDirName = col->getSaveDirName();
  RestartDirName = col->getRestartDirName();
  restart = col->getRestart_status();
  ns = col->getNs();            // get the number of particle species involved in simulation
  TPns = col->getTP_Ns();          // get the number of test particle species involved in simulation
  first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart

  qom  = new double[ns];
  for (int i=0; i < ns; i++)
    qom[i] = col->getQOM(i);

  x_center = col->getx_center();
  y_center = col->gety_center();
  z_center = col->getz_center();
  L_square = col->getL_square();

  freqInj = (int) col->getVinj();



  // initialize the virtual cartesian topology 
  vct = new VCtopology3D(col);
  // Check if we can map the processes into a matrix ordering defined in Collective.cpp
  if (nprocs != vct->getNprocs()) {
    if (myrank == 0) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << "x" << vct->getZLEN() << " matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()" << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
      return (-1);
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

  probes = new Probes(col, grid);
  
  if (col->getSolInit()) {
    /* -------------------------------------------- */
    /* If using parallel H5hut IO read initial file */
    /* -------------------------------------------- */
    if (myrank==0) 
      cout << "INIT FIELDS FROM FILE" << endl;
    ReadFieldsH5hut(ns, true, true, EMf, col, vct, grid);

    if (col->getCase()=="TwoCoils" || col->getCase()=="SixCoils" || col->getCase()=="TwoCoilsINJ" || col->getCase()=="SixCoilsINJ" ||
        col->getCase()=="Mirror"){
      EMf->setRhoInit(grid, col);
      EMf->SetLambda(grid, col);
      EMf->setFext(1);
    }

  }
  else {
    /* --------------------------------------------------------- */
    /* If using 'default' IO initialize fields depending on case */
    /* --------------------------------------------------------- */
    if      (col->getCase()=="GEM")         EMf->initGEM(vct, grid,col);
    else if (col->getCase()=="DoubleGEM")   EMf->initDoubleGEM(vct, grid,col);
    else if (col->getCase()=="DoubleHarris")EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
    else if (col->getCase()=="FF")          EMf->initForceFree(vct,grid,col);
    else if (col->getCase()=="FF_kappa")    EMf->initForceFree(vct,grid,col);
    else if (col->getCase()=="TD")          EMf->tangential_discontinuity(vct,grid,col); 
    else if (col->getCase()=="mKHI")        EMf->mag_KHI(vct,grid,col);
    else if (col->getCase()=="SI")          EMf->slabs_instability(vct,grid,col);
    else if (col->getCase()=="IC")          EMf->ion_cyclotron(vct,grid,col);
    else if (col->getCase()=="Whistler")    EMf->initDoublePeriodicHarrisWithGaussianHumpPerturbation(vct,grid,col);
    else if (col->getCase()=="TwoCoils" || col->getCase()=="TwoCoilsINJ")    EMf->initWB8(vct,grid,col);
    else if (col->getCase()=="SixCoils" || col->getCase()=="SixCoilsINJ")    EMf->initWB8_6(vct,grid,col);
    else if (col->getCase()=="Mirror")    EMf->initMirror(vct,grid,col);
    else if (col->getCase()=="Dipole")      EMf->initDipole(vct,grid,col);
    else if (col->getCase()=="Bessel1D")    EMf->initBessel1D(vct,grid,col);
    else if (col->getCase()=="Bessel2D")    EMf->initBessel2D(vct,grid,col);
    else if (col->getCase()=="TEmode")      EMf->initTEmode(vct,grid,col);
    else if (col->getCase()=="TMmode")      EMf->initTMmode(vct,grid,col);
    else if (col->getCase()=="pulse")       EMf->initPulse(vct,grid,col);
    else if (col->getCase()=="HalfPlanet")  EMf->initHalfPlanet(vct,grid,col);
    else if (col->getCase()=="FluxRope")    EMf->initFluxRope(vct,grid,col);
    else if (col->getCase()=="Par_LH_circPol_mono_Alfven") EMf->initPar_LH_circPol_mono_Alfven(vct, grid, col);
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

  if (col->getSolInit() && col->getPartInit()=="File") {
    ReadPartclH5hut(ns, part, col, vct, grid);
    for (int i = 0; i < ns; i++) {
      if (col->getTrackParticleID(i)) part[i].tagParticles();
    }
  } else {
    for (int i = 0; i < ns; i++) {
		  part[i].allocate(false, i, 0, col, vct, grid);
  
      if      (col->getSolInit() && col->getPartInit() == "FromFieldFile")  part[i].MaxwellianFromFieldFile(grid, EMf, vct);
      else if (col->getCase()=="DoubleGEM")      part[i].MaxwellianDoubleGEM(grid, EMf, vct);
      else if (col->getCase()=="FF")             part[i].Maxwellianspacedist(col, grid, EMf, vct);
      else if (col->getCase()=="TD")             part[i].MaxwellianTangentialDiscontinuity(grid, EMf, vct);
      else if (col->getCase()=="mKHI")           part[i].Maxwellian_magKHI(grid, EMf, vct);
      else if (col->getCase()=="SI")             part[i].MaxwellianSlabInstability(grid, EMf, vct);
      else if (col->getCase()=="IC")             part[i].MaxwellianIonCyclotron(grid, EMf, vct);
      else if (col->getCase()=="Dipole")         part[i].maxwellian(grid, EMf, vct);
      else if (col->getCase()=="FluxRope")       part[i].maxwellian(grid, EMf, vct);
      else if (col->getCase()=="Par_LH_circPol_mono_Alfven") part[i].maxwellian_Par_LH_circPol_mono_Alfven(grid, EMf, vct, col);
      else if (col->getCase()=="FF_kappa") part[i].Kappaspacedist(col,grid,EMf,vct);
      else if (col->getPartInit()=="Kappa") part[i].kappa(grid,EMf,vct);
      else if (col->getCase()=="TwoCoils" || col->getCase()=="SixCoils" || col->getCase()=="TwoCoilsINJ" || col->getCase()=="SixCoilsINJ"
               || col->getCase()=="Mirror"){
        if (col->getRHOinit(i) > 0.0) 
          if (col->getPartInit() == "mono") {
            part[i].monoenergetic_sphere(grid,EMf,vct,L_square,x_center,y_center,z_center, 1.0);  // generated maxwellian in a box
          } else {
            part[i].maxwell_sphere(grid,EMf,vct,L_square,x_center,y_center,z_center, 1.0);  // generated maxwellian in a box
          }
        else 
          part[i].empty(grid, EMf, vct);
  
      } else part[i].maxwellian(grid, EMf, vct);
      
      // Compute the charge of the particles using the mass matrices
      //part[i].fixCharge(grid, EMf, vct, col);
  
      if (col->getCylindrical() == 1) part[i].fixZcylindrical();
  
      if (col->getLagrangian() == 1) part[i].setLagrangianFrame(true);
    }
  }

  // Fix the kinetic energy of species (it is done locally)
  if ((col->getRHOinit(0) > 0.0)) {
    if (col->getSameKineticEnergy() == "yes"){
      if (myrank == 0) printf("ALL SPECIES WITH THE SAME KINETIC ENERGY: Fixing velocities\n");
      long long np1 = part[0].getNOP();
      double Kref = part[0].getLocalKe(0, np1);
      for (int i=1; i<ns; i++) {
        if (col->getRHOinit(i) > 0.0) {
          long long np1 = part[i].getNOP();
          double K = part[i].getLocalKe(0, np1);
          part[i].fixKineticEnergy(Kref, K, 0, np1, 100);
        }
      }
    }
  }

  // Test particles
  testpart = new Particles3D[TPns];
  if (col->getTestIDfile() == "no") {
    for (int i=0; i<TPns; i++) {
  		  testpart[i].allocate(true, i, 0, col, vct, grid);
        if (col->getCase()=="DoubleGEM")           testpart[i].MaxwellianDoubleGEM(grid, EMf, vct);
        else if (col->getCase()=="FF")             testpart[i].Maxwellianspacedist(col, grid, EMf, vct);
        else if (col->getCase()=="TD")             testpart[i].MaxwellianTangentialDiscontinuity(grid, EMf, vct);
        else if (col->getCase()=="SI")             testpart[i].MaxwellianSlabInstability(grid, EMf, vct);
        else if (col->getCase()=="IC")             testpart[i].MaxwellianIonCyclotron(grid, EMf, vct);
        else if (col->getCase()=="Dipole")         testpart[i].maxwellian(grid, EMf, vct);
        else if (col->getCase()=="FluxRope")       testpart[i].maxwellian(grid, EMf, vct);
	else if (col->getCase()=="FF_kappa")       testpart[i].Kappaspacedist(col,grid,EMf,vct);
	else if (col->getPartInit()=="Kappa")      testpart[i].kappa(grid,EMf,vct);
        else if (col->getCase()=="TwoCoils" || col->getCase()=="SixCoils" || col->getCase()=="TwoCoilsINJ" || col->getCase()=="SixCoilsINJ"
                 || col->getCase()=="Mirror"){
            if (col->getPartInit() == "mono") {
              testpart[i].monoenergetic_sphere(grid,EMf,vct,col->getTP_L_square(i),col->getTP_x_center(i),col->getTP_y_center(i),col->getTP_z_center(i), 1.0);  // generated maxwellian in a box
            } else {
              testpart[i].maxwell_sphere(grid,EMf,vct, col->getTP_L_square(i),col->getTP_x_center(i),col->getTP_y_center(i),col->getTP_z_center(i) , 1.0);  // generated maxwellian in a box
            }
        } else testpart[i].maxwellian(grid, EMf, vct);
        
        if (col->getCylindrical() == 1) testpart[i].fixZcylindrical();
    
        if (col->getLagrangian() == 1) testpart[i].setLagrangianFrame(true);
    }
  } else {
    FilePartID* fileID = new FilePartID(col);
    for (int i=0; i<TPns; i++) {
      testpart[i].initFromIDfile(i, fileID, part, col, vct, grid);
    }

    delete fileID;
  }



  if (col->getLagrangian() == 1) EMf->setLagrangianFrame(vct, col, true);



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

  Qremoved = new double[ns];

  Eenergy, Benergy, TOTenergy = 0.0, TOTmomentum = 0.0;
  Ke = new double[ns];
  momentum = new double[ns];
  cq = SaveDirName + "/ConservedQuantities.txt";
  

  /* Compute the external current to have the total current equal to zero */
  if (col->getZeroCurrent() == 1) {
    EMf->setZeroDensities();                  // set to zero the densities
    for (int i = 0; i < ns; i++)
      part[i].computeCurrent(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)

    EMf->setZeroCurrent(vct);
  }

  
  EMf->setZeroRho();                  // set to zero the charge density
  #ifdef _CHARGE_FROM_CURRENT_
  for (int i = 0; i < ns; i++) {
    part[i].computeCharge(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)
  }
  EMf->communicateGhostP2G_Rho(vct);
  EMf->interpDensitiesN2C(vct, grid);       // calculate densities on centers from nodes
  #endif


  gridUpdated = false;

  clocks->stop(0);
  return 0;
}

void c_Solver::GatherMoments(){

  if (col->getStaticField()) return;


  #ifdef __PROFILING__
  clocks->start(5);
  #endif
  EMf->setZeroDensities();                  // set to zero the densities
  #ifdef __PROFILING__
  clocks->stop(5);
  #endif

  for (int i = 0; i < ns; i++)
    part[i].computeMoments(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)


  EMf->sumOverSpecies (vct);

  // Communicate everything:
  #ifdef __PROFILING__
  clocks->start(10);
  #endif
  EMf->communicateGhostP2G(vct);
  #ifdef __PROFILING__
  clocks->stop(10);
  #endif


  // Compute the charge density as the divergence of the current
  #ifdef _CHARGE_FROM_CURRENT_
  EMf->chargeFromCurrent(grid, vct);
  #endif
}

/** Compute the charge, current and pressure densities. They are not needed for the 
 * calculation, only for the output */
void c_Solver::GatherMoments_extra() {

 if (!gridUpdated) {
    // Compute the charge density and the current of each specie 
    EMf->setZeroDensities();                  // set to zero the densities
    for (int i = 0; i < ns; i++)
      part[i].computeCurrent(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)
  
    #ifndef _CHARGE_FROM_CURRENT_
    EMf->setZeroRho();                  // set to zero the densities
    for (int i = 0; i < ns; i++)
      part[i].computeCharge(EMf, grid, vct);      // interpolate Particles to Grid(Nodes)
  
    EMf->communicateGhostP2G_Rho(vct);
    #endif
    
    EMf->communicateGhostP2G_JP(vct);

    gridUpdated = true;
  }
}

void c_Solver::CalculateField() {

  if (col->getStaticField()) return;

  MPI_Barrier(MPI_COMM_WORLD);

  // MAXWELL'S SOLVER
  #ifdef __WITH_PETSC__
    petscSolver->solveE();
  #else
    EMf->calculateEB(grid, vct, col);               // calculate the E field
  #endif

}

bool c_Solver::ParticlesMover(int cycle) {

  /*  -------------- */
  /*  Particle mover */
  /*  -------------- */

  #ifdef  __JDOTE_TEST__
		double totK0 = 0;
    for (int is = 0; is < ns; is++) {
      totK0 +=  part[is].getKe();
    }
  #endif 

  // move each species
#ifdef PC_MOVER
  for (int i = 0; i < ns; i++)  
  mem_avail = part[i].mover_PC_sub(grid, vct, EMf);
  for (int i = 0; i < TPns; i++)  
  mem_avail = testpart[i].mover_PC_sub(grid, vct, EMf);
#else
  for (int i = 0; i < ns; i++) 
    mem_avail = part[i].updateVelocity(grid, vct, EMf);
  for (int i = 0; i < TPns; i++) 
    mem_avail = testpart[i].updateVelocity(grid, vct, EMf);

  #ifdef  __JDOTE_TEST__
  double totK1 = 0;
  for (int i = 0; i < ns; i++) {
    part[i].computeCurrentBar(EMf, grid, vct); 
    totK1 +=  part[i].getKe();
  }
  double deltaK = totK1 - totK0;
  EMf->communicateGhostP2G_JP(vct);
  EMf->testJdotE(grid, deltaK);
  #endif

  for (int i = 0; i < ns; i++)  {
    mem_avail = part[i].updatePosition(grid, vct, EMf);
    part[i].fixPosition();
  }
  for (int i = 0; i < TPns; i++)  {
    mem_avail = testpart[i].updatePosition(grid, vct, EMf);
    testpart[i].fixPosition();
  }
#endif

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  // Update the values of the magnetic field in the nodes: n+1
  EMf->centers2nodesB(vct, grid, col);

  /*-----------------------------------*/
  /*  Update the lagrangian frame      */
  /*-----------------------------------*/
  if (col->getLagrangian() == 1) {
    EMf->updateLagrangianFrame(vct, col);
    for (int i=0; i<ns; i++) part[i].updateLagrangianFrame(vct);
  }

  /* -------------------------------------- */
  /* Repopulate the buffer zone at the edge */
  /* -------------------------------------- */
  InjectBoundaryParticles(cycle);

  if (mem_avail < 0) {          // not enough memory space allocated for particles: stop the simulation
    if (myrank == 0) {
      cout << "*************************************************************" << endl;
      cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
      cout << "*************************************************************" << endl;
    }
    return (true);              // exit from the time loop
  }

  gridUpdated = false;

  return (false);
}

void c_Solver::InjectBoundaryParticles(int cycle){

  int np0[ns];

  for (int i=0; i < ns; i++) {
    np0[i] = part[i].getNOP();
    if (col->getRHOinject(i)>0.0){
      if (col->getCase()=="HalfPlanet") 
         part[i].injector_rand_box_mono(grid,vct,EMf, col->getx_center(),col->gety_center(),col->getz_center(), col->getL_square());
      else if (col->getCase()=="Dipole" || col->getCase()=="GEM") {
         mem_avail = part[i].particle_repopulator(grid,vct,EMf,i);
      }
    }

    if (col->getCase()=="TwoCoils") {
       Qremoved[i] = part[i].ReturnToCenterCircle(col->getDelta());
       double dst_2 = col->getcoilSpacing()/2.0;
       Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()+dst_2,col->getz_center());
       Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()-dst_2,col->getz_center());

    } else if (col->getCase()=="SixCoils") {
       Qremoved[i] = part[i].ReturnToCenterCircle(col->getDelta());

    } else if (col->getCase()=="TwoCoilsINJ") {
      Qremoved[i]  = part[i].deleteParticlesOutsideSphere(EMf, col->getDelta(),col->getx_center(),col->gety_center(),col->getz_center());
      double dst_2 = col->getcoilSpacing()/2.0;
      Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()+dst_2,col->getz_center());
      Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()-dst_2,col->getz_center());
      if ((freqInj > 0) && (cycle > 0) && ((cycle-1) % freqInj == 0) && (col->getRHOinject(i)>0.0)) {
        if (col->getPartInj() == "mono") 
          part[i].inject_mono_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);
        else
          part[i].inject_maxwell_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);  // generated maxwellian in a box
      }

    } else if (col->getCase()=="SixCoilsINJ") {
      Qremoved[i]  = part[i].deleteParticlesOutsideSphere(EMf, col->getDelta(),col->getx_center(),col->gety_center(),col->getz_center());
      if (col->getRHOinject(i)>0.0)
        if (col->getPartInj() == "mono") 
          part[i].inject_mono_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);
        else
          part[i].inject_maxwell_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);  // generated maxwellian in a box

    } else if (col->getCase()=="Mirror") {
      Qremoved[i]  = part[i].deleteParticlesOutsideSphere(EMf, col->getDelta(),col->getx_center(),col->gety_center(),col->getz_center());
      double dst_2 = col->getcoilSpacing()/2.0;
      Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()+dst_2,col->getz_center());
      Qremoved[i] += part[i].deleteParticlesInsideTorus_XZ(EMf, col->getcoilD(), col->getcoilThickness(), col->getx_center(),col->gety_center()-dst_2,col->getz_center());
      if ((freqInj > 0) && (cycle > 0) && ((cycle-1) % freqInj == 0) && (col->getRHOinject(i)>0.0)) {
        if (col->getPartInj() == "mono") 
          part[i].inject_mono_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);
        else
          part[i].inject_maxwell_sphere(grid,EMf,vct,col->getL_square_inj(i),col->getx_center_inj(i),col->gety_center_inj(i),col->getz_center_inj(i), 1.0);  // generated maxwellian in a box
      }
          
    } else if (col->getCase()=="Dipole") {
       Qremoved[i] = part[i].deleteParticlesInsideSphere(EMf, col->getL_square(),col->getx_center(),col->gety_center(),col->getz_center());
       //Qremoved[i] = part[i].deleteParticlesOuterFrame(2.0, 2.0, 0.0);
    } else if (col->getCase()=="HalfPlanet") {
       Qremoved[i] = part[i].deleteParticlesInsideSphere(EMf, col->getL_square(),col->getLx(),col->gety_center(),col->getz_center());
       Qremoved[i] = part[i].deleteParticlesOuterFrame(EMf, 6.0, 6.0, 6.0);
    }
  }

  // Fix the kinetic energy of species (it is done locally)
  if ((freqInj > 0) && (cycle > 0) && ((cycle-1) % freqInj == 0) && (col->getRHOinject(0)>0.0)) {
    if (col->getSameKineticEnergy() == "yes"){
      if (myrank== 0) printf("ALL SPECIES WITH THE SAME KINETIC ENERGY: Fixing velocities\n");
      long long np1 = part[0].getNOP();
      double Kref = part[0].getLocalKe(np0[0], np1);
      for (int i=1; i<ns; i++) {
        if (col->getRHOinject(i)>0.0) {
          long long np1 = part[i].getNOP();
          double K = part[i].getLocalKe(np0[i], np1);
          part[i].fixKineticEnergy(Kref, K, np0[i], np1, 100.0);
        }
      }
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
  if (cycle == first_cycle) {
  return;
  }

  if (probes->getNP() > 0) {
    GatherMoments_extra();
    probes->setValues(EMf);
    if (cycle == LastCycle()) probes->writeFiles();
  }


  if (cycle == first_cycle+1) {

    totEn0 = EMf->getEenergy() + EMf->getBintenergy() + EMf->getEenergyRemoved(false);
    for (int is = 0; is < ns; is++) {
      totEn0 +=  part[is].getKe();
      totEn0 += part[is].getKremoved();
    }
    if (myrank == 0) {
      FILE* fd = fopen(cq.c_str(),"a");
      fprintf(fd,"# 1 -> Cycle\n"
                 "# 2 -> Total internal energy: sum(B_plasma*B_plasma) + sum(E_plasma*E_plasma) + sum(K_part)\n"
                 "# 3 -> Variation of total internal energy\n"
                 "# 4 -> Electric energy: sum(E_plasma*E_plasma)\n"
                 "# 5 -> Local magnetic energy: sum(B_total*B_total) - sum(B_extern*B_extern)\n"
                 "# 6 -> Kinetic energy (currently in the domain plus the one removed in the past)\n"
                 "# 7 -> Momentum\n"
                 "# 8 -> Total magnetic energy: sum(B_total*B_total)\n"
                 "# 9 -> Internal magnetic energy: sum(B_local*B_local)\n"
                 "# 10-> Kinetic energy removed\n"
                 "# 11-> Electric energy removed\n"
                 "# 12-> Number of particles of species 0\n"
                 "# 13-> Total charge of species 0\n"
                 "# 14-> Kinetic energy of species 0\n"
                 "#(15-> Number of particles of species 1)\n"
                 "#(16-> Total charge of species 1)\n"
                 "#(17-> Kinetic energy of species 1)\n"
                 "#(...)\n");
      fprintf(fd,"# 1        2           3           4           5           6           7           8           9           10        11        12    \n"
                 "#----------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
      fclose(fd);
    }

  }

  // write the conserved quantities
  if (cycle % col->getDiagnosticsOutputCycle() == 0) {
    // Check the memory usage
    check_mem_usage();

    Eenergy = EMf->getEenergy();
    double ExtBenergy = EMf->getBextenergy();
    double IntBenergy = EMf->getBintenergy();
    Benergy = EMf->getBtotenergy();
    double EenRem  = EMf->getEenergyRemoved(true);

    TOTenergy = 0.0;
    TOTmomentum = 0.0;
    double Kr= 0;
    long long nop[ns];
    double tq[ns];
    for (int is = 0; is < ns; is++) {
      Ke[is] = part[is].getKe();
      double lKr = part[is].getKremoved();
      Kr     += lKr;
      TOTenergy += Ke[is] + lKr;
      momentum[is] = part[is].getP();
      TOTmomentum += momentum[is];
      nop[is] = part[is].getTotalNumPart();
      tq[is]  = part[is].getTotalCharge();
    }
    if (myrank == 0) {
      FILE* fd = fopen(cq.c_str(),"a");
      fprintf(fd,"%5d %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e"
                ,cycle
                ,(EenRem + Eenergy + IntBenergy + TOTenergy), fabs((EenRem + Eenergy + IntBenergy + TOTenergy) - totEn0)
                ,Eenergy , Benergy-ExtBenergy, TOTenergy
                ,TOTmomentum
                ,Benergy
                ,IntBenergy
                ,Kr
                ,EenRem
                );
      for (int is = 0; is < ns; is++) fprintf(fd," %7ld %11.4e %11.4e", nop[is], tq[is], Ke[is]);
      fprintf(fd,"\n");

      fclose(fd);
    }
  }
}


void c_Solver::WriteOutput(int cycle) {

  // Go to the Eulerian frame 
  if (col->getLagrangian() == 1) {
    EMf->setLagrangianFrame(vct, col, false);
    for (int i = 0; i < ns; i++) part[i].setLagrangianFrame(false);
  }

  EMf->divergenceOfE(grid, vct);
  GatherMoments_extra();

  // test
  if (0){
    EMf->printDensity(0);
    cout << endl;
    EMf->printDensity(1);
  }

  
  if (col->getWriteMethod() == "ascii") {
    if (cycle%(col->getFieldOutputCycle())!=0) return;
    char place[128];
    const char* dirtxt = SaveDirName.c_str();
    sprintf(place,"%s/filedsProc%02d_%06d.txt",dirtxt, myrank, cycle);
    FILE *fd = fopen(place,"w");
    double ***Ex = EMf->getEx();
    double ***Ey = EMf->getEy();
    double ***Ez = EMf->getEz();
    double ***Bx = EMf->getBx();
    double ***By = EMf->getBy();
    double ***Bz = EMf->getBz();
    double ***Jxh= EMf->getJxh();
    double ***Jyh= EMf->getJyh();
    double ***Jzh= EMf->getJzh();
    double ***Jx_ext= EMf->getJx_ext();
    double ***Jy_ext= EMf->getJy_ext();
    double ***Jz_ext= EMf->getJz_ext();
    double ***rhon = EMf->getRHOn();
    int nxn = grid->getNXN();
    int nyn = grid->getNYN();
    int nzn = grid->getNZN();


    if (col->getCase() == "Bessel1D") {
      double Kx = ZERO_J1/col->getLx();
      double dt = col->getDt();
      double omega = Kx;
      double t = dt*cycle;
      double totErr = 0;
      static double sumTotErr = 0;
      fprintf(fd,"# 1,2,3   : x, y, z \n"
                 "# 4,5,6   : Ex, Ey, Ez\n"
                 "# 7,8,9   : Bx, By, Bz\n"
                 "# 10,11,12: Ex, Ey, Ez (Analytical solution)\n"
                 "# 13,14,15: Bx, By, Bz (Analytical solution)\n"
                 "# 16,17,18: error Ex, Ey, Ez \n"
                 "# 19,20,21: error Bx, By, Bz \n"
             );
      for (int i=1; i<nxn-1; i++) 
        for (int j=1; j<nyn-1; j++) 
          for (int k=1; k<nzn-1; k++) {
            double x = grid->getXN(i,j,k);
            double y = grid->getYN(i,j,k);
            double z = grid->getZN(i,j,k);
            double feX = 0;
            double feY = 0;
            double feZ = besselJ(1, x*Kx)*cos(t*omega); 
            double fbX = 0;
            double fbY = -besselJ(0, x*Kx)*sin(t*omega)*Kx/omega;
            double fbZ = 0;
            double errEx = fabs(Ex[i][j][k] - feX);
            double errEy = fabs(Ey[i][j][k] - feY);
            double errEz = fabs(Ez[i][j][k] - feZ);
            double errBx = fabs(Bx[i][j][k] - fbX);
            double errBy = fabs(By[i][j][k] - fbY);
            double errBz = fabs(Bz[i][j][k] - fbZ);
           
            totErr += errEx + errEy + errEz + errBx + errBy + errBz;
            fprintf(fd,"%8.4f %8.4f %8.4f   %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                      ,x,y,z
                      ,Ex[i][j][k],  Ey[i][j][k],  Ez[i][j][k]
                      ,Bx[i][j][k],  By[i][j][k],  Bz[i][j][k]
                      ,feX, feY, feZ, fbX, fbY, fbZ
                      ,errEx, errEy, errEz, errBx, errBy, errBz);
          }
      totErr /= (nxn-2)*(nyn-2)*(nzn-2);
      sumTotErr += totErr;
      if (cycle == LastCycle()) fprintf(stderr,"%3d  %13.6e  %13.6e \n", cycle, totErr, sumTotErr/cycle);

    } else if (col->getCase() == "Bessel2D") {
      double Kx = ZERO_J1/col->getLx();
      double Ky = 2*M_PI/col->getLy();
      double dt = col->getDt();
      double omega = sqrt(Ky*Ky + Kx*Kx);
      double t = dt*cycle;
      double totErr = 0;
      static double sumTotErr = 0;
      fprintf(fd,"# 1,2,3   : x, y, z \n"
                 "# 4,5,6   : Ex, Ey, Ez\n"
                 "# 7,8,9   : Bx, By, Bz\n"
                 "# 10,11,12: Ex, Ey, Ez (Analytical solution)\n"
                 "# 13,14,15: Bx, By, Bz (Analytical solution)\n"
                 "# 16,17,18: error Ex, Ey, Ez \n"
                 "# 19,20,21: error Bx, By, Bz \n"
             );
      for (int i=1; i<nxn-1; i++) 
        for (int j=1; j<nyn-1; j++) 
          for (int k=1; k<nzn-1; k++) {
            double x = grid->getXN(i,j,k);
            double y = grid->getYN(i,j,k);
            double z = grid->getZN(i,j,k);
            double feX = 0;
            double feY = 0;
            double feZ = cos(y*Ky)*besselJ(1, x*Kx)*cos(t*omega); 
            double fbX = -sin(y*Ky)*besselJ(1, x*Kx)*sin(t*omega)*Ky/omega;
            double fbY = -cos(y*Ky)*besselJ(0, x*Kx)*sin(t*omega)*Kx/omega;
            double fbZ = 0;
            double errEx = fabs(Ex[i][j][k] - feX);
            double errEy = fabs(Ey[i][j][k] - feY);
            double errEz = fabs(Ez[i][j][k] - feZ);
            double errBx = fabs(Bx[i][j][k] - fbX);
            double errBy = fabs(By[i][j][k] - fbY);
            double errBz = fabs(Bz[i][j][k] - fbZ);
           
            totErr += errEx + errEy + errEz + errBx + errBy + errBz;
            fprintf(fd,"%8.4f %8.4f %8.4f   %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                      ,x,y,z
                      ,Ex[i][j][k],  Ey[i][j][k],  Ez[i][j][k]
                      ,Bx[i][j][k],  By[i][j][k],  Bz[i][j][k]
                      ,feX, feY, feZ, fbX, fbY, fbZ
                      ,errEx, errEy, errEz, errBx, errBy, errBz);
          }
      totErr /= (nxn-2)*(nyn-2)*(nzn-2);
      sumTotErr += totErr;
      if (cycle == LastCycle()) fprintf(stderr,"%3d  %13.6e  %13.6e \n", cycle, totErr, sumTotErr/cycle);

    } else if (col->getCase() == "TEmode") {
      double Kr = ZERO_J1/col->getLx();
      double Kz = 2*M_PI/col->getLy();
      double dt = col->getDt();
      double omega = sqrt(Kz*Kz + Kr*Kr);
      double t = dt*cycle;
      double totErr = 0;
      static double sumTotErr = 0;
      fprintf(fd,"# 1,2,3   : x, y, z \n"
                 "# 4,5,6   : Ex, Ey, Ez\n"
                 "# 7,8,9   : Bx, By, Bz\n"
                 "# 10,11,12: Ex, Ey, Ez (Analytical solution)\n"
                 "# 13,14,15: Bx, By, Bz (Analytical solution)\n"
                 "# 16,17,18: error Ex, Ey, Ez \n"
                 "# 19,20,21: error Bx, By, Bz \n"
             );
      for (int i=1; i<nxn-1; i++) 
        for (int j=1; j<nyn-1; j++) 
          for (int k=1; k<nzn-1; k++) {
            double x = grid->getXN(i,j,k);
            double y = grid->getYN(i,j,k);
            double z = grid->getZN(i,j,k);
            double feX = 0;
            double feY = 0;
            double feZ = -(omega/Kr)*besselJ(1, x*Kr)*sin(y*Kz-omega*t);
            double fbX = (Kz/Kr)*besselJ(1, x*Kr)*sin(y*Kz-omega*t);
            double fbY = besselJ(0, x*Kr)*cos(y*Kz-omega*t);
            double fbZ = 0;
            double errEx = fabs(Ex[i][j][k] - feX);
            double errEy = fabs(Ey[i][j][k] - feY);
            double errEz = fabs(Ez[i][j][k] - feZ);
            double errBx = fabs(Bx[i][j][k] - fbX);
            double errBy = fabs(By[i][j][k] - fbY);
            double errBz = fabs(Bz[i][j][k] - fbZ);
           
            totErr += errEx + errEy + errEz + errBx + errBy + errBz;
            fprintf(fd,"%8.4f %8.4f %8.4f   %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                      ,x,y,z
                      ,Ex[i][j][k],  Ey[i][j][k],  Ez[i][j][k]
                      ,Bx[i][j][k],  By[i][j][k],  Bz[i][j][k]
                      ,feX, feY, feZ, fbX, fbY, fbZ
                      ,errEx, errEy, errEz, errBx, errBy, errBz);
          }
      totErr /= (nxn-2)*(nyn-2)*(nzn-2);
      sumTotErr += totErr;
      if (cycle == LastCycle()) fprintf(stderr,"%3d  %13.6e  %13.6e \n", cycle, totErr, sumTotErr/cycle);

    } else if (col->getCase() == "TMmode") {
      double Kr = ZERO_J0/col->getLx();
      double Kz = 2*M_PI/col->getLy();
      double dt = col->getDt();
      double omega = sqrt(Kz*Kz + Kr*Kr);
      double t = dt*cycle;
      double totErr = 0;
      static double sumTotErr = 0;
      fprintf(fd,"# 1,2,3   : x, y, z \n"
                 "# 4,5,6   : Ex, Ey, Ez\n"
                 "# 7,8,9   : Bx, By, Bz\n"
                 "# 10,11,12: Ex, Ey, Ez (Analytical solution)\n"
                 "# 13,14,15: Bx, By, Bz (Analytical solution)\n"
                 "# 16,17,18: error Ex, Ey, Ez \n"
                 "# 19,20,21: error Bx, By, Bz \n"
             );
      for (int i=1; i<nxn-1; i++) 
        for (int j=1; j<nyn-1; j++) 
          for (int k=1; k<nzn-1; k++) {
            double x = grid->getXN(i,j,k);
            double y = grid->getYN(i,j,k);
            double z = grid->getZN(i,j,k);
            double feX = Kz*besselJ(1, x*Kr)*sin(Kz*y - omega*t);
            double feY = besselJ(0, x*Kr)*cos(Kz*y - omega*t);
            double feZ = 0;
            double fbX = 0;
            double fbY = 0;
            //double fbZ = -besselJ(1, x*Kr)*sin(Kz*y + omega*t);
            double fbZ = omega*besselJ(1, x*Kr)*sin(Kz*y - omega*t);
            double errEx = fabs(Ex[i][j][k] - feX);
            double errEy = fabs(Ey[i][j][k] - feY);
            double errEz = fabs(Ez[i][j][k] - feZ);
            double errBx = fabs(Bx[i][j][k] - fbX);
            double errBy = fabs(By[i][j][k] - fbY);
            double errBz = fabs(Bz[i][j][k] - fbZ);
           
            totErr += errEx + errEy + errEz + errBx + errBy + errBz;
            fprintf(fd,"%8.4f %8.4f %8.4f   %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e\n"
                      ,x,y,z
                      ,Ex[i][j][k],  Ey[i][j][k],  Ez[i][j][k]
                      ,Bx[i][j][k],  By[i][j][k],  Bz[i][j][k]
                      ,feX, feY, feZ, fbX, fbY, fbZ
                      ,errEx, errEy, errEz, errBx, errBy, errBz);
          }
      totErr /= (nxn-2)*(nyn-2)*(nzn-2);
      sumTotErr += totErr;
      if (cycle == LastCycle()) fprintf(stderr,"%3d  %13.6e  %13.6e \n", cycle, totErr, sumTotErr/cycle);


    } else {
      for (int i=1; i<nxn-1; i++) 
        for (int j=1; j<nyn-1; j++) 
          for (int k=1; k<nzn-1; k++) {
            double x = grid->getXN(i,j,k);
            double y = grid->getYN(i,j,k);
            double z = grid->getZN(i,j,k);
            fprintf(fd,"%8.4f %8.4f %8.4f   %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e  %13.6e %13.6e %13.6e   %13.6e \n"
                      ,x,y,z
                      ,Ex[i][j][k],  Ey[i][j][k],  Ez[i][j][k]
                      ,Bx[i][j][k],  By[i][j][k],  Bz[i][j][k]
                      ,Jxh[i][j][k], Jyh[i][j][k], Jzh[i][j][k]
                      ,Jxh[i][j][k]+Jx_ext[i][j][k], Jyh[i][j][k]+Jy_ext[i][j][k], Jzh[i][j][k]+Jz_ext[i][j][k]
                      ,rhon[i][j][k]
                    );
          }
    }
    fclose(fd);
  } else if (col->getWriteMethod() == "h5hut") {

    /* -------------------------------------------- */
    /* Parallel HDF5 output using the H5hut library */
    /* -------------------------------------------- */
    if (col->getFieldOutputCycle() != 0) 
      //if (cycle == first_cycle || cycle%(col->getFieldOutputCycle())==0)        
      if (cycle == first_cycle || (cycle%(col->getFieldOutputCycle()))< col->getSuccessiveFieldsSaved() )
	WriteFieldsH5hut(ns, grid, EMf,  col, vct, cycle);

    if (col->getParticlesOutputCycle() != 0) 
      if (cycle == first_cycle || cycle%(col->getParticlesOutputCycle())==0)    WritePartclH5hut(ns, grid, part, col, vct, cycle, false);
    
    if (TPns > 0) 
      if (cycle % col->getTestParticlesOutputCycle() == 0) WritePartclH5hut(TPns, grid, testpart, col, vct, cycle, true);

  } else {

    // OUTPUT to large file, called proc**
    if (col->getFieldOutputCycle() != 0) {
      if ((cycle % (col->getFieldOutputCycle())) < col->getSuccessiveFieldsSaved()  ) { 
      //if (cycle % (col->getFieldOutputCycle()) == 0 ) {
        hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
        output_mgr.output("Eall + Ball + rhos + Jsall + pressure", cycle);
        // Pressure tensor is available
        hdf5_agent.close();
      }
    }
    if (col->getParticlesOutputCycle() != 0) {
      if (cycle % (col->getParticlesOutputCycle()) == 0) {
        hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
        output_mgr.output("position + velocity + q ", cycle, 1);
        hdf5_agent.close();
      }
    }
  }

  // Go back to the Lagrangian frame 
  if (col->getLagrangian() == 1) {
    EMf->setLagrangianFrame(vct, col, true);
    for (int i = 0; i < ns; i++) part[i].setLagrangianFrame(true);
  }
}


void c_Solver::Finalize() {

  // stop profiling
  clocks->stop(15);

  // deallocate
  delete[]Ke;
  delete[]momentum;
  delete[] Qremoved;
  // close MPI
  mpi->finalize_mpi();
}

