cmake_minimum_required(VERSION 2.8.5)

#
# Set project name
#

project(ECsim)

#
# Find packages (already managed by CMake)
#

find_package(HDF5 COMPONENTS HL C REQUIRED)
find_package(MPI REQUIRED)

#find_package(CUDA_NVTX_LIBRARY REQUIRED)
# Find PETSc
#
list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")

## NSIGHTS
##set (CMAKE_CXX_FLAGS "-cuda")
set (CMAKE_CXX_FLAGS "-cuda -acc -gpu=cc70,managed,lineinfo -Minfo=accel")

#-Minfo=accel: info on compilation of parallel region
#lineinfo: compare GPU assembly and code-- profiling, kernel levle
#
# Options
#
option(PROFILING "Measure the times in each stage" OFF)
if(PROFILING)
  add_definitions( -D__PROFILING__ )
endif()

option(TEST_JDOTE "Check if the current from the particles and fields are the
same" OFF)
if(TEST_JDOTE)
  add_definitions( -D__JDOTE_TEST__)
endif()

option(CHARGE_FROM_CURRENT "Compute the charge density from the charge conservation continuity equation" OFF)
if(CHARGE_FROM_CURRENT)
  add_definitions( -D_CHARGE_FROM_CURRENT_)
endif()

option(SAME_CHARGE "All particles have the same charge. The number of particles per cell might no be the same in all cells" OFF)
if(SAME_CHARGE)
  add_definitions( -D_WITH_SAME_CHARGE_)
endif()

option(IPIC_TESTS "Set up the code tests" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Chose the type of build, options
are: None, Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
endif()

#
# Include sub-directories
#

add_subdirectory(main)

include_directories(
    include
    ${MPI_INCLUDE_PATH}
    ${HDF5_INCLUDE_DIRS}
    ${PETSC_INCLUDES}
    /home/sw-hackathons/opt/spack/linux-ubuntu18.04-broadwell/gcc-7.5.0/cuda-11.2.2-jfcwf7pdscjm6ep3t4i4fhrrqrinzu6u/include
)

add_definitions(
  ${PETSC_DEFINITIONS}
)


#
# Executable compilation
#

add_executable(
    ECsim
    ECsim.cpp
)

#
# Libraries linked
#

target_link_libraries(
    ECsim
    ECsimlib
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    ${PETSC_LIBRARIES}
    /home/sw-hackathons/opt/spack/linux-ubuntu18.04-broadwell/gcc-7.5.0/cuda-11.2.2-jfcwf7pdscjm6ep3t4i4fhrrqrinzu6u/lib64/libnvToolsExt.so
)

#
# Code testing
#

if(IPIC_TESTS)
  enable_testing()
  set(IPIC_TESTS_DIR "${CMAKE_BINARY_DIR}/tests" CACHE STRING "Location of the
source files for ECsim")

  add_test(NAME GEM-test
           COMMAND ${CMAKE_COMMAND}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -DIPIC_SOURCE_DIR=${CMAKE_SOURCE_DIR}
           -DIPICEXEC=$<TARGET_FILE:ECsim>
           -DMPIEXEC=${MPIEXEC}
           -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
           -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-GEM.txt)

  add_test(NAME uname-test
           COMMAND ${CMAKE_COMMAND}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-uname.txt)

endif(IPIC_TESTS)
