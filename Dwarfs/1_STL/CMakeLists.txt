cmake_minimum_required(VERSION 3.0)
project(1_STL)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  
endif()
include( ${CGAL_USE_FILE} )

# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")
  return()  
endif()
include( CGAL_CreateSingleSourceCGALProgram )

add_definitions(-std=c++14 -Wall -Werror)

add_executable( weld_vtx  weld_vtx.cpp )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS weld_vtx )
# Link the executable to CGAL and third-party libraries
target_link_libraries(weld_vtx   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

find_package(OpenMP)
if (NOT ${OPENMP_FOUND})
  message("## Skipping OpenMP. OpenMP not found.")
  return()
endif()

add_executable( weld_vtx_omp  weld_vtx.cpp )
add_to_cached_list( CGAL_EXECUTABLE_TARGETS weld_vtx_omp )
target_link_libraries(weld_vtx_omp   ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )
set_target_properties(weld_vtx_omp PROPERTIES COMPILE_FLAGS "-D_GLIBCXX_PARALLEL ${OpenMP_CXX_FLAGS}")
set_target_properties(weld_vtx_omp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
