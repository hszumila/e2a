# CMakeLists for Adin's fiducial cuts

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add executable for acceptance corrections
add_executable(count count.cpp ) 
target_link_libraries(count ${ROOT_LIBRARIES} e2aCuts)

# Add executable for acceptance corrections
add_executable(fiducial_histmaker fiducial_histmaker.cpp)
target_link_libraries(fiducial_histmaker ${ROOT_LIBRARIES})

# Add executable for acceptance corrections
add_executable(fiducial_analysis fiducial_analysis.cpp)
target_link_libraries(fiducial_analysis ${ROOT_LIBRARIES})

