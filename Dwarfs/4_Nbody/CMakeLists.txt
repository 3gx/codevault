cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(4_NBody NONE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

add_subdirectory(hermite4)
#add_subdirectory(bhtree)
#add_subdirectory(fmm)
#add_subdirectory(sph)
#add_subdirectory(md)


