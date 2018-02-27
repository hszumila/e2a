# CMakeLists for e2a maps

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(sim_pid_electron sim_pid_electron.cpp ${CMAKE_SOURCE_DIR}/Fiducial.cpp ${CMAKE_SOURCE_DIR}/Run_dependent.cpp) 
target_link_libraries(sim_pid_electron ${ROOT_LIBRARIES})
