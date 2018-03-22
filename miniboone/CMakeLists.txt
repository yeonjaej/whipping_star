add_executable (minirepo miniboone.cxx)
add_executable (minirepoL minibooneL.cxx)

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (minirepo SBNfitlib)
target_link_libraries (minirepo tinyxmllib)
target_link_libraries (minirepoL SBNfitlib)
target_link_libraries (minirepoL tinyxmllib)

