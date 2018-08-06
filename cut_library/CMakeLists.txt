# CMakeLists for the e2a analysis repository

SET(LIBNAME "e2aCuts")

include_directories( ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR} )

# Library for all of the cuts
add_library(${LIBNAME} Fiducial.cpp Run_dependent.cpp)
target_link_libraries(${LIBNAME} ${ROOT_LIBRARIES} )
