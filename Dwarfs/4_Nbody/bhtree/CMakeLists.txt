cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(BHTree NONE)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

enable_language(CXX)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
set(CXX11 ${COMPILER_SUPPORTS_CXX11})
set(CXX11_FLAGS -std=c++11)
 
set(example "\"BHTree\"")
find_package(OpenMP)

macro(build_openmp)
  set_property(SOURCE src/bhtree_omp.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
  add_executable(bhtree_omp src/bhtree_omp.cpp)
  set_target_properties(bhtree_omp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endmacro(build_openmp)

if (NOT CXX11)
  message(WARNING "*** C++11 is not supported by compiler. Ingoring ${example} example")
else (NOT CXX11)
  add_definitions(${CXX11_FLAGS} -Wall -Werror)
  add_executable(bhtree src/bhtree.cpp)
  add_executable(bhtree_mod src/bhtree_mod.cpp)

  find_package(OpenMP)

  if (NOT ${OPENMP_FOUND})
    message(WARNING "*** OpenMP is not support by compiler. Ignoring ${example} OpenMP examples")
  else (NOT ${OPENMP_FOUND})
    build_openmp()
  endif()

endif (NOT CXX11)


