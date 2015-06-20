cmake_minimum_required(VERSION 3.2)
project(imn_libs)
# Version number
set(imn_libs_VERSION_MAJOR 0)
set(imn_libs_VERSION_MINOR 2)

# Differential Schemes lib
include_directories("${PROJECT_SOURCE_DIR}/diff_schemes")
add_subdirectory(diff_schemes)

# Viscous Flux lib
include_directories("${PROJECT_SOURCE_DIR}/vi_flow")
add_subdirectory(vi_flow)


include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")

set(SOURCE_FILES main.cpp)
add_executable(imn ${SOURCE_FILES})

# Differential Schemes lib
target_link_libraries(imn LINK_PUBLIC Diff_Schemes)
# Viscous Flux link
target_link_libraries(imn LINK_PUBLIC Vi_Flow)

install(TARGETS imn DESTINATION bin)

