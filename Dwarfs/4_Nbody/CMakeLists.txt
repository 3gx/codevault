cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(CodeVault NONE)

set(CodeVault_VERSION_MAJOR 0)
set(CodeVault_VERSION_MINOR 1)
set(CodeVault_VERSION_PATCH 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake" "${CMAKE_SOURCE_DIR}/cmake/Modules")

enable_language (C)
enable_language (CXX)

#include(detect_fortran_compiler)
#if (${Fortran_COMPILER_FOUND})
#  enable_language (Fortran)
#else(${Fortran_COMPILER_FOUND})
#  message(WARNING "*** Fortran is not found. Will not build Fortran examples")
#endif(${Fortran_COMPILER_FOUND})


# RPATH settings
set(CMAKE_SKIP_BUILD_RPATH false) # Use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_BUILD_WITH_INSTALL_RPATH false) # When building, don't use the install RPATH already
set(CMAKE_INSTALL_RPATH "") # The RPATH to be used when installing
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH false) # Don't add the automatically determined parts

# CXX Compiler-version check (requires at least CMake 2.8.10)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    message(FATAL_ERROR "GCC version must be at least 4.8")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.3)
    message(FATAL_ERROR "Clang version must be at least 3.3")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    message(FATAL_ERROR "Clang (OS X) version must be at least 5.0")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14.0)
    message(FATAL_ERROR "ICC version must be at least 14.0")
  endif()
endif()

#if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
#  if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.8)
#    message(FATAL_ERROR "gfortran version must be at least 4.8")
#  endif()

#message("CMAKE_Fortran_COMPILER = ${CMAKE_Fortran_COMPILER}")
#message("CMAKE_Fortran_COMPILER_ID= ${CMAKE_Fortran_COMPILER_ID}")
#message("CMAKE_Fortran_COMPILER_NAMES= ${CMAKE_Fortran_COMPILER_NAMES}")
#if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  # gfortran
#  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
#elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#  # ifort (untested)
#  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
#elseif (Fortran_COMPILER_NAME MATCHES "g77")
#  # g77
#  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
#else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#  message ("No optimized Fortran compiler flags are known, we just try -O2...")
#  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
#  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
#endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")



# C++ compiler settings
set(FLAGS "-O3 -std=c++11")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(FLAGS "${FLAGS} -Wall -Wno-comment")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")

# ==================================================================================================

find_package(ISPC QUIET)
find_package(CUDA QUIET)
find_package(MPI QUIET)
find_package(OpenMP QUIET)
find_package(GPUDirectRDMA QUIET)
find_package(OpenCL QUIET)
find_package(Threads QUIET)

if (${ISPC_FOUND})
else (${ISPC_FOUND})
  message(WARNING "*** ISPC is not found. Will not build ISPC examples")
endif (${ISPC_FOUND})

if (${CUDA_FOUND})
else (${CUDA_FOUND})
  message(WARNING "*** CUDA is not found. Will not build CUDA examples")
endif (${CUDA_FOUND})

if (${MPI_FOUND})
else (${MPI_FOUND})
  message(WARNING "*** MPI is not found. Will not build MPI examples")
endif (${MPI_FOUND})

if (${OPENMP_FOUND})
else (${OPENMP_FOUND})
  message(WARNING "*** OpenMP is not found. Will not build OpenMP examples")
endif (${OPENMP_FOUND})

if (${OpenCL_FOUND})
else (${OpenCL_FOUND})
  message(WARNING "*** OpenCL is not found. Will not build OpenCL examples")
endif (${OpenCL_FOUND})

if (${GPUDirectRDMA_FOUND})
else (${GPUDirectRDMA_FOUND})
  message(WARNING "*** GPUDirect RDMA is not found. Will not build GPUDirect RDMA examples")
endif (${GPUDirectRDMA_FOUND})

if(${Threads_FOUND})
else(${Threads_FOUND})
  message(WARNING "*** [p]Threads is not found. Will not build [p]Threads examples")
endif(${Threads_FOUND})

configure_file (
  "${PROJECT_SOURCE_DIR}/include/CodeVaultConfig.h.in"
  "${PROJECT_BINARY_DIR}/CodeVaultConfig.h"
  )


add_subdirectory(hermite4)


