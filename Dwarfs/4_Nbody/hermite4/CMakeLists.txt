if (CODEVAULT_ENABLE_SERIAL)
  
endif()

if (CODEVAULT_ENABLE_OPENMP)
  
endif()

if (CODEVAULT_ENABLE_CUDA)
  
endif()

add_executable(serial serial.cpp)

find_package(IntelCompilers QUIET)

if (${OPENMP_FOUND})
  set_source_files_properties(openmp.cpp PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")
  add_executable(openmp openmp.cpp)
  set_target_properties(openmp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif()

if (${ICPC_FOUND})
  set_source_files_properties(openmp_icc.cpp PROPERTIES COMPILE_FLAGS "${FLAGS} -openmp")
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/openmp_icc.cpp
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/openmp.cpp ${CMAKE_CURRENT_BINARY_DIR}/openmp_icc.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/openmp.cpp
    COMMENT "Copy openmp.cpp to openmp_icc.cpp"
    )
  set_source_files_properties(openmp_icc.cpp PROPERTIES COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  icpc_add_executable(openmp_icc openmp_icc.cpp)
  set_target_properties(openmp_icc PROPERTIES LINK_FLAGS "-openmp")
  #  add_custom_command(
  #  TARGET openmp_icc
  #  POST_BUILD
  #  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_CURRENT_BINARY_DIR}/openmp_icc.cpp
  #  COMMENT "Remove openmp_icc.cpp"
  #  )
endif()

