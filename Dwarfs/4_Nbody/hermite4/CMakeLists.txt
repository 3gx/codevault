if (CODEVAULT_ENABLE_SERIAL)
  
endif()

if (CODEVAULT_ENABLE_OPENMP)
  
endif()

if (CODEVAULT_ENABLE_CUDA)
  
endif()

add_executable(serial serial.cpp)
