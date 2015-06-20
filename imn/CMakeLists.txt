cmake_minimum_required(VERSION 3.2)
project(imn_libs)
# Version
set(Imn_VERSION_MAJOR 2)
set(Imn_VERSION_MINOR 1)

add_library(Imn diff_schemes.cxx vi_flow.cxx)
target_include_directories(Imn PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Imn DESTINATION bin)
install(FILES diff_schemes.hpp vi_flow.hpp DESTINATION include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")

