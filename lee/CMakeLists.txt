add_executable (DLSens DLSens.cc)
add_executable (DLCovar DLCovariance.cc)
include_directories ( ${CMAKE_SOURCE_DIR}/inc)

target_link_libraries (DLSens SBNfitlib)
target_link_libraries (DLSens tinyxmllib)

target_link_libraries (DLCovar SBNfitlib)
target_link_libraries (DLCovar tinyxmllib)
