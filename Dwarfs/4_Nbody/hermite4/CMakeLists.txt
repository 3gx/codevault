if (CODEVAULT_ENABLE_SERIAL)
  
endif()

if (CODEVAULT_ENABLE_OPENMP)
  
endif()

if (CODEVAULT_ENABLE_CUDA)
  
endif()

add_executable(serial serial.cpp)

if (${OPENMP_FOUND})
  set_source_files_properties(openmp.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  add_executable(openmp openmp.cpp)
  set_target_properties(openmp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()
