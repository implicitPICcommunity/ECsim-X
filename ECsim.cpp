#include <iomanip>
#include "ECsim.h"
#include "MyClock.h"
#include "nvtx.h"


MyClock *clocks;

using namespace ECsim;

int main(int argc, char **argv) {

  ECsim::c_Solver KCode;
  bool b_err = false;

  

  /* ------------------------------ */
  /* 0- Initialize the solver class */
  /* ------------------------------ */

  KCode.Init(argc, argv);
  KCode.InjectBoundaryParticles(0);
  KCode.WriteOutput(KCode.FirstCycle());

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */
  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */
    mynvtxstart_("momentGathering");
    clocks->start(1);
    KCode.GatherMoments();
    clocks->stop(1);
    mynvtxstop_();

    mynvtxstart_("fields");
    clocks->start(2);
    KCode.CalculateField();
    clocks->stop(2);
    mynvtxstop_();

    mynvtxstart_("particleMover");
    clocks->start(3);
    b_err = KCode.ParticlesMover(i);
    clocks->stop(3);
    mynvtxstop_();

    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */
    mynvtxstart_("IO");
    clocks->start(4);
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);
    clocks->stop(4);
    mynvtxstop_();
  }


  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  KCode.Finalize();


  // Profiling:
  #ifdef __PROFILING__
  if (myrank == 0) {
    std::cout << "################################################" << std::endl
            << "Initialization                : " << clocks->get(0) <<  std::endl
            << "Moments Gathering             : " << clocks->get(1) <<  std::endl
            << " +- B interp. centers->nodes  : " << clocks->get(5) <<  std::endl
            << " +- B interp. nodes->particles: " << clocks->get(6) <<  std::endl
            << " +- Charge density            : " << clocks->get(7) <<  std::endl
            << " +- Implicit Current          : " << clocks->get(8) <<  std::endl
            << " +- Mass matrix               : " << clocks->get(9) <<  std::endl
            << " +- Communications            : " << clocks->get(10)<<  std::endl
            << "Field Calculation             : " << clocks->get(2) <<  std::endl
            << " +- Maxwell sources           : " << clocks->get(11)<<  std::endl
            << " +- Maxwell imagen            : " << clocks->get(12)<<  std::endl
            << "     +- Communications        : " << clocks->get(13)<<  std::endl
            << " +- End E calc (comm)         : " << clocks->get(14)<<  std::endl
            << "Particle Mover                : " << clocks->get(3) <<  std::endl
            << "Writing                       : " << clocks->get(4) <<  std::endl
            << "----------------------------------------------------------" << std::endl
            << "Total                         : " << clocks->get(15)<<  std::endl;
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
            << "Total                         : " << clocks->get(15)<< " s" << std::endl;
  }
  #endif


  return 0;
}
