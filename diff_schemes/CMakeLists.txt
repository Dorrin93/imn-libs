cmake_minimum_required(VERSION 3.2)
project(imn_libs)
# Version number
set(diff_schemes_VERSION_MAJOR 1)
set(diff_schemes_VERSION_MINOR 1)

add_library(Diff_Schemes diff_schemes.cxx)
target_include_directories(Diff_Schemes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS Diff_Schemes DESTINATION bin)
install(FILES diff_schemes.hpp DESTINATION include)

