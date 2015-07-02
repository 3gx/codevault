
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/Languages")
enable_language(ICPC OPTIONAL)
#enable_language(IFORT OPTIONAL)

set(ICPC_FOUND False)
if (${CMAKE_ICPC_COMPILER} STREQUAL "CMAKE_ICPC_COMPILER-NOTFOUND")
else()
  function(icpc_add_executable icpc_target)
    add_executable(${icpc_target} ${ARGN})
    foreach(file ${ARGN})
      set_source_files_properties(${file} PROPERTIES LANGUAGE ICPC)
    endforeach()
    set_target_properties(${icpc_target} PROPERTIES LINKER_LANGUAGE ICPC)
  endfunction()
  set(ICPC_FOUND True)

  #  function(icpc_copy
endif()
