ECsim-X: a Particle-In-Cell code implementing the Energy Conserving Semi Implicit Method
========================================================================================

Currently there are two codebases
*CPU for running on CPU
*OpenACC for using NVIDIA GPUs


For reference on the method, please check Lapenta G. (2017), “Exactly energy conserving semi-implicit particle in cell formulation, Journal of Computational Physics, vol. 334, p. 349.
For reference on the code implementation, please check D. Gonzalez-Herrero, E. Boella, and G. Lapenta (2018), Performance analysis and implementation details of the Energy Conserving Semi-Implicit Method code (ECsim), Computer Physics Communications, vol. 229, p. 162.
For reference on the OpenACC implementation, please check E. Boella, N. Shukla, F. Spiga, M. K. Chimeh, M. Bettencourt, and M. E. Innocenti (2026), "Accelerating the Particle-In-Cell Code ECsim with OpenACC," in 2026 34th Euromicro International Conference on Parallel, Distributed, and Network-Based Processing (PDP), Cluj-Napoca, Romania, p. 224.

Dependencies
------------

The following dependencies are required to build and run the project:
- C/C++ compiler: A compiler supporting C and C++ (e.g., GCC or LLVM/Clang).
- MPI: An MPI distribution supporting MPI >= 4.0. Tested with: Open MPI and MPICH
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/): Version >= 1.10.6.
- [H5hut](https://github.com/eth-cscs/h5hut): Version >= 2.0.0.
- [PETSc](https://petsc.org/release/): Version >= 3.20.1.
- [CMake](https://cmake.org/): Version >= 2.8.5.

Building ECsim
--------------
```bash
mkdir build
cd build
cmake -DIPIC_PETSC_SOLVER=ON -DIPIC_H5HUT_OUTPUT=ON -DH5HUT_HOME=${H5HUT_HOME} ..
make -j
```
