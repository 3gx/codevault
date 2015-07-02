
function(clone_source src dst)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${dst}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${src} ${CMAKE_CURRENT_BINARY_DIR}/${dst}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${src}
    COMMENT "Clone ${src} to ${dst}"
    )
endfunction()

macro(get_options _sources _option_list _option_name)
  set( ${_sources} )
  set( ${_option_list} )
  set( _found_options False)
  foreach(arg ${ARGN})
    if ("x${arg}" STREQUAL "x${_option_name}")
      set (_found_options True)
    else()
      if (_found_options)
        list(APPEND ${_option_list} ${arg})
      else()
        list(APPEND ${_sources} ${arg})
      endif()
    endif()
  endforeach()
endmacro()
