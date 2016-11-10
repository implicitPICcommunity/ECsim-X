#include <iomanip>
#include "ECsim.h"
#include "MyClock.h"

MyClock *clocks;

using namespace ECsim;

int main(int argc, char **argv) {

  ECsim::c_Solver KCode;
  bool b_err = false;

  

  /* ------------------------------ */
  /* 0- Initialize the solver class */
  /* ------------------------------ */

  KCode.Init(argc, argv);
  KCode.InjectBoundaryParticles();

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */
  for (int i = KCode.FirstCycle(); i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */
    KCode.GatherMoments();
    KCode.CalculateField();
    b_err = KCode.ParticlesMover();
    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    //KCode.WriteOutputASCII(i);
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);
  }


  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  KCode.Finalize();


  // Profiling:
  #ifdef __PROFILING__
  if (myrank == 0) {
    std::cout << "################################################" << std::endl
            << "Initialization                : " << clocks->get(0) << " s" << std::endl
            << "Moments Gathering             : " << clocks->get(1) << " s" << std::endl
            << " +- B interp. centers->nodes  : " << clocks->get(5) << " s" << std::endl
            << " +- B interp. nodes->particles: " << clocks->get(6) << " s" << std::endl
            << " +- Charge density            : " << clocks->get(7) << " s" << std::endl
            << " +- Implicit Current          : " << clocks->get(8) << " s" << std::endl
            << " +- Mass matrix               : " << clocks->get(9) << " s" << std::endl
            << " +- Communications            : " << clocks->get(10)<< " s" << std::endl
            << "Field Calculation             : " << clocks->get(2) << " s" << std::endl
            << "Particle Mover                : " << clocks->get(3) << " s" << std::endl
            << "Writing                       : " << clocks->get(4) << " s" << std::endl
            << "----------------------------------------------------------" << std::endl
            << "Total                         : " << clocks->get(11)<< " s" << std::endl;
  }
  #else
  if (myrank == 0) {
    std::cout << "################################################" << std::endl
            << "Initialization                : " << clocks->get(0) << " s" << std::endl
            << "Moments Gathering             : " << clocks->get(1) << " s" << std::endl
            << "Field Calculation             : " << clocks->get(2) << " s" << std::endl
            << "Particle Mover                : " << clocks->get(3) << " s" << std::endl
            << "Writing                       : " << clocks->get(4) << " s" << std::endl
            << "----------------------------------------------------------" << std::endl
            << "Total                         : " << clocks->get(11)<< " s" << std::endl;
  }
  #endif


  return 0;
}
