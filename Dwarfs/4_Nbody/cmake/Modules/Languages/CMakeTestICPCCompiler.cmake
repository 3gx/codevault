# This file is used by EnableLanguage in cmGlobalGenerator to determine that
# the ICPC builder GNAT_EXECUTABLE_BUILDER = gnatmake can actually compile
# and link the most basic of programs.  If not, a fatal error is set and
# cmake stops processing commands and will not generate any makefiles or
# projects.

function(PrintTestCompilerStatus LANG MSG)
  if(CMAKE_GENERATOR MATCHES Make)
    message(STATUS "Check for working ${LANG} compiler: ${CMAKE_${LANG}_COMPILER}${MSG}")
  else()
    message(STATUS "Check for working ${LANG} compiler using: ${CMAKE_GENERATOR}${MSG}")
  endif()
endfunction()

unset(CMAKE_CXX_COMPILER_WORKS CACHE)
if (NOT CMAKE_ICPC_COMPILER_WORKS)
  PrintTestCompilerStatus("ICPC" "")
  file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testICPCCompiler.cxx
    "#ifndef __cplusplus\n"
    "# error \"The CMAKE_ICPC_COMPILER is set to a C compiler\"\n"
    "#endif\n"
    "int main(){return 0;}\n")
  try_compile(CMAKE_ICPC_COMPILER_WORKS ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/testICPCCompiler.cxx
    OUTPUT_VARIABLE __CMAKE_ICPC_COMPILER_OUTPUT)
  # Move result from cache to normal variable.
  set(CMAKE_ICPC_COMPILER_WORKS ${CMAKE_ICPC_COMPILER_WORKS})
  unset(CMAKE_ICPC_COMPILER_WORKS CACHE)
  set(ICPC_TEST_WAS_RUN 1)
endif (NOT CMAKE_ICPC_COMPILER_WORKS)

if(NOT CMAKE_ICPC_COMPILER_WORKS)
  PrintTestCompilerStatus("ICPC" " -- broken")
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining if the ICPC compiler works failed with "
    "the following output:\n${__CMAKE_ICPC_COMPILER_OUTPUT}\n\n")
  message(FATAL_ERROR "The C++ compiler \"${CMAKE_ICPC_COMPILER}\" "
    "is not able to compile a simple test program.\nIt fails "
    "with the following output:\n ${__CMAKE_ICPC_COMPILER_OUTPUT}\n\n"
    "CMake will not be able to correctly generate this project.")
else()
  if(ICPC_TEST_WAS_RUN)
    PrintTestCompilerStatus("ICPC" " -- works")
    file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
      "Determining if the ICPC compiler works passed with "
      "the following output:\n${__CMAKE_ICPC_COMPILER_OUTPUT}\n\n")
  endif()

  # Try to identify the ABI and configure it into CMakeICPCCompiler.cmake
  include(${CMAKE_ROOT}/Modules/CMakeDetermineCompilerABI.cmake)
  CMAKE_DETERMINE_COMPILER_ABI(ICPC ${CMAKE_ROOT}/Modules/CMakeCXXCompilerABI.cpp)
  # Try to identify the compiler features
  include(${CMAKE_ROOT}/Modules/CMakeDetermineCompileFeatures.cmake)
  CMAKE_DETERMINE_COMPILE_FEATURES(ICPC)

  # Re-configure to save learned information.
  #  configure_file(
  #  ${CMAKE_ROOT}/Modules/CMakeICPCCompiler.cmake.in
  #  ${CMAKE_PLATFORM_INFO_DIR}/CMakeICPCCompiler.cmake
  #  @ONLY
  #  )
  #include(${CMAKE_PLATFORM_INFO_DIR}/CMakeICPCCompiler.cmake)

  if(CMAKE_ICPC_SIZEOF_DATA_PTR)
    foreach(f ${CMAKE_ICPC_ABI_FILES})
      include(${f})
    endforeach()
    unset(CMAKE_ICPC_ABI_FILES)
  endif()
endif()

unset(__CMAKE_ICPC_COMPILER_OUTPUT)
