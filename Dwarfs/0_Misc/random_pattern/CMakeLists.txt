cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(Random CXX)

find_package(OpenMP)

if (${OPENMP_FOUND})
  add_definitions("-std=c++11 -Wall -Werror ${OpenMP_CXX_FLAGS}")
  include_directories(".")
  add_executable(random random.cpp)
  set_target_properties(random PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif(${OPENMP_FOUND})


