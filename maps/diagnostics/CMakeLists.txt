# CMakeLists for diagnostic tools for testing maps

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

add_executable(basic_test basic_test.cpp)
target_link_libraries(basic_test ${ROOT_LIBRARIES})
