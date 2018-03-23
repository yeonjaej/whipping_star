add_executable (exampledemo examples.cxx)

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (exampledemo SBNfitlib)
target_link_libraries (exampledemo tinyxmllib)


