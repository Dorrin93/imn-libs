cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -lm -O3")

project(imn_libs CXX)
# Version
set(Imn_VERSION_MAJOR 3)
set(Imn_VERSION_MINOR 0)

add_library(Imn diff_schemes.cxx grid.cxx iter_schemes.cxx ogrid.cxx p_flow.cxx utils.cxx vi_flow.cxx wave_equation.cxx)
target_include_directories(Imn PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Imn DESTINATION bin)
install(FILES diff_schemes.hpp grid.hpp iter_schemes.hpp utils.hpp ogrid.hpp p_flow.hpp vi_flow.hpp wave_equation.hpp DESTINATION include)


