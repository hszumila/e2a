# CMakeLists for e2a maps

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(e2a_maps e2a_maps.cpp ${CMAKE_SOURCE_DIR}/Fiducial.cpp ${CMAKE_SOURCE_DIR}/Run_dependent.cpp) 
target_link_libraries(e2a_maps ${ROOT_LIBRARIES})
