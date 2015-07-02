
function(clone_source src dst)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${dst}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${src} ${CMAKE_CURRENT_BINARY_DIR}/${dst}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${src}
    COMMENT "Clone ${src} to ${dst}"
    )
endfunction()
