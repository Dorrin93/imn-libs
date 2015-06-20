cmake_minimum_required(VERSION 3.2)
project(imn_libs)
# Version number
set(imn_libs_VERSION_MAJOR 0)
set(imn_libs_VERSION_MINOR 2)

include_directories("${PROJECT_SOURCE_DIR}/imn")
add_subdirectory(imn)

include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lm")

set(SOURCE_FILES main.cpp)
add_executable(imn_exe ${SOURCE_FILES})

target_link_libraries(imn_exe LINK_PUBLIC Imn)

install(TARGETS imn_exe DESTINATION bin)
