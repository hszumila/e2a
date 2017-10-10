# CMakeLists for the e2a analysis repository

# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

# Project name
project (e2a)

# ADD SUBDIRECTORIES
add_subdirectory(calibration_data)

include_directories(${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(skim skim.cxx)
target_link_libraries(skim ${ROOT_LIBRARIES})

# Add executable for acceptance corrections
add_executable(skim_tree skim_tree.cpp)
target_link_libraries(skim_tree ${ROOT_LIBRARIES})