if (CODEVAULT_ENABLE_SERIAL)
  
endif()

if (CODEVAULT_ENABLE_OPENMP)
  
endif()

if (CODEVAULT_ENABLE_CUDA)
  
endif()

add_executable(serial serial.cpp)

find_package(IntelCompilers QUIET)

if (OPENMP_FOUND)
  #set_source_files_properties(openmp.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  set_property(SOURCE openmp.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
  add_executable(openmp openmp.cpp)
  set_target_properties(openmp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

if (ICPC_FOUND)
  clone_source(openmp.cpp openmp_icc.cpp)
  #  set_property(SOURCE openmp_icc.cpp APPEND_STRING PROPERTY COMPILE_FLAGS "  -openmp -Wall -Werror")
  icpc_sources(openmp_icc.cpp COMPILE_FLAGS " -openmp -Wall -Werror")
  icpc_add_executable(openmp_icc openmp_icc.cpp LINK_FLAGS "-openmp")
  #  add_custom_command(
  #  TARGET openmp_icc
  #  POST_BUILD
  #  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/openmp_icc.cpp
  #  COMMENT "Remove openmp_icc.cpp"
  #  )
endif()

if (IFORT_FOUND)
  ifort_sources(openmp.f90 COMPILE_FLAGS "-openmp")
  ifort_add_executable(openmp_f90 openmp.f90 LINK_FLAGS "-openmp")
endif()

