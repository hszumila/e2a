# CMakeLists for creating maps and other programs
#   that require CLASTOOL


set (CLASTOOL /u/home/baraks/packages/ClasTool)
set (ANALYZER /u/home/baraks/packages/Analyser/analysis_lib )
set (OSNAME Linux64RHEL7 )

include_directories( ${CLASTOOL}/include ${ANALYZER}/include ${ROOT_INCLUDE_DIR} )
link_directories( ${ANALYZER}/slib ${CLASTOOL}/slib/${OSNAME} )

add_executable(convert_uniform convert_uniform.cxx)
target_link_libraries(convert_uniform ${ROOT_LIBRARIES})

add_executable(mapmaker mapmaker.cpp)
target_link_libraries(mapmaker ${ROOT_LIBRARIES} e2aCuts)

add_executable(generator generator.cxx)
target_link_libraries(generator ${ROOT_LIBRARIES})
