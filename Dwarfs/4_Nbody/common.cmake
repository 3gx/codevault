if (NOT COMMON_CMAKE_SET)
  set(COMMON_CMAKE_SET True)

  # set path to basic modules
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake" "${CMAKE_CURRENT_LIST_DIR}/cmake/Modules")

  # users can add their paths here
  # ....
endif()
