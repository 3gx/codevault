if (CODEVAULT_ENABLE_SERIAL)
  
endif()

if (CODEVAULT_ENABLE_OPENMP)
  
endif()

if (CODEVAULT_ENABLE_CUDA)
  
endif()

#add_executable(serial serial.cpp)

#if (True)
#  add_custom_command(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/serial.o
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/serial.cpp
#    COMMAND g++ -std=c++11 -o ${CMAKE_CURRENT_BINARY_DIR}/serial.o    -I${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/serial.cpp -c
#    COMMENT "Build CXX object ${CMAKE_CURRENT_BINARY_DIR}/serial.o"
#    )
#  add_custom_target(
#    my_serial ALL
#    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/serial.o
#    COMMAND g++ -o my_serial ${CMAKE_CURRENT_BINARY_DIR}/serial.o
#    COMMENT "Linking CXX object my_serial")
#endif()


#if (OPENMP_FOUND)
#  #set_source_files_properties(openmp.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
#  set_property(SOURCE openmp.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
#  add_executable(openmp openmp.cpp)
#  set_target_properties(openmp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
#endif()

# Locate IntelCompilers if available
find_package(IntelCompilers QUIET)

#if found Intel C++ compiler ,build C++ example
#if (ICPC_FOUND)
#  clone_source(openmp.cpp openmp_icc.cpp)
#  icpc_sources(openmp_icc.cpp COMPILE_FLAGS " -openmp -Wall -Werror")
#  icpc_add_executable(openmp_icc openmp_icc.cpp LINK_FLAGS "-openmp")
#endif()

#if found Intel Fortran compiler, build fortran example
#if (IFORT_FOUND)
#  ifort_sources(openmp.f90 COMPILE_FLAGS "-openmp")
#  ifort_add_executable(openmp_f90 openmp.f90 LINK_FLAGS "-openmp")
#endif()


# Attempt localte ISPC compiler
find_package(ISPC)

# if found, build ISPC Example
if (ISPC_FOUND)
  set(ispc_obj "")
  ispc_compile(ispc_obj ispc/hermite4.ispc COMPILE_FLAGS "-O2")
  include_directories("util")

  if (True)
    icpc_sources(
      ispc/omp_tasksys.cpp
      ispc/hermite4.cpp
      COMPILE_FLAGS "-O3 -g -openmp")
    #  add_definitions("-fopenmp")
    icpc_add_executable(hermite4_ispc 
      ${ispc_obj}
      ispc/hermite4.ispc
      ispc/hermite4.cpp
      ispc/omp_tasksys.cpp
      LINK_FLAGS "-openmp")
  else()
    add_definitions("-fopenmp")
    add_executable(hermite4_ispc 
      ${ispc_obj}
      ispc/hermite4.ispc
      ispc/hermite4.cpp
      ispc/omp_tasksys.cpp)
    set_target_properties(hermite4_ispc PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    #    set(path
    #  ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/hermite4_ispc.dir)
  #  set(ispc_path
  #  ${CMAKE_CURRENT_BINARY_DIR}/ispc)
  #  message("path= ${path}")
  #  add_custom_command(
  #  TARGET hermite4_ispc
  #  PRE_LINK
  #  COMMAND g++ -o hermite4_ispc.x    ${ispc_path}/hermite4.dev.o    ${path}/ispc/hermite4.cpp.o    ${path}/ispc/omp_tasksys.cpp.o    -fopenmp
  #  DEPENDS  ${ispc_path}/hermite4.dev.o   ${path}/ispc/hermite4.cpp.o  ${path}/ispc/omp_tasksys.cpp.o  
  #  COMMENT "Pre_link command")
  endif()
endif()

if (ISPC_FOUND AND XEONPHI_FOUND)
endif()

