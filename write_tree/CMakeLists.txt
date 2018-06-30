# CMakeLists for write_tree_e2a and other programs
#   that require CLASTOOL

set (CLASTOOL $ENV{CLASTOOL} )
set (ANALYZER /u/home/baraks/packages/Analyser/analysis_lib )
set (OSNAME $ENV{OS_NAME} )

include_directories( ${CLASTOOL}/include ${ANALYZER}/include ${ROOT_INCLUDE_DIR} )
link_directories( ${ANALYZER}/slib ${CLASTOOL}/slib/${OSNAME} )

add_executable(write_tree_e2a write_tree_e2a.cxx)
target_link_libraries(write_tree_e2a ${ROOT_LIBRARIES} ClasTool ClasBanks TIdentificator)
