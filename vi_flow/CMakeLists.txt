cmake_minimum_required(VERSION 3.2)
project(imn_libs)
# Version number
set(vi_flow_VERSION_MAJOR 1)
set(vi_flow_VERSION_MINOR 0)

add_library(Vi_Flow vi_flow.cxx)
target_include_directories(Vi_Flow PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Vi_Flow DESTINATION bin)
install(FILES vi_flow.hpp DESTINATION include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")
