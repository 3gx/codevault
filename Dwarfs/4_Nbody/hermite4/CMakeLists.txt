cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
project(Hermite4 NONE)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

enable_language(CXX)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
set(CXX11 ${COMPILER_SUPPORTS_CXX11})
set(CXX11_FLAGS -std=c++11)

if (NOT CXX11)
  message(WARNING "*** C++11 is not supported by compiler. Ingoring \"hermite4\"  example")
else (NOT CXX11)
  add_definitions(${CXX11_FLAGS} -Wall -Werror)
  add_executable(hermite4_serial hermite4_serial.cpp)
  
  find_package(OpenMP)

  if (NOT ${OPENMP_FOUND})
    message(WARNING "*** OpenMP is not support by compiler. Ignoring \"hermite4\" OpenMP examples")
  else (NOT ${OPENMP_FOUND})

    # Build OpenMP example
    set_property(SOURCE hermite4_omp.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
    add_executable(hermite4_omp hermite4_omp.cpp)
    set_target_properties(hermite4_omp PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
  
    find_package(ISPC)
    if (${ISPC_FOUND})

      # Build ISPC example
      set (ispc_obj "")
      if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        ispc_compile(ispc_obj ispc/hermite4.ispc COMPILE_FLAGS "-O2 --opt=fast-math --math-lib=svml")
      else()
        ispc_compile(ispc_obj ispc/hermite4.ispc COMPILE_FLAGS "-O2 --opt=fast-math")
      endif()
      set_property(SOURCE ispc/hermite4.cpp    APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")
      set_property(SOURCE ispc/omp_tasksys.cpp APPEND_STRING PROPERTY COMPILE_FLAGS " ${OpenMP_CXX_FLAGS}")

      add_executable(hermite4_ispc
        ${ispc_obj} 
        ispc/typeReal.h
        ispc/hermite4.ispc
        ispc/hermite4.cpp 
        ispc/omp_tasksys.cpp)

      set_target_properties(hermite4_ispc PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    else()
      message(WARNING "*** ISPC is not found. Ignoring  \"hermite4\" ISPC example")
    endif()

  endif (NOT ${OPENMP_FOUND})

endif (NOT CXX11)


