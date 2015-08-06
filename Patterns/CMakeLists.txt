cmake_minimum_required(VERSION 2.8.10 FATAL_ERROR)
project(Patterns CXX)

set(FLAGS "-O3 -std=c++1y -Wall -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAGS}")
include_directories("../include")

add_executable(copy-and-swap copy-and-swap.cpp)
add_executable(singleton singleton.cpp)
add_executable(multiple_factories multiple_factories.cpp)
add_executable(type_erasure type_erasure.cpp)

