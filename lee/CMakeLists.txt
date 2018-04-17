
add_executable (LEEfit LEEfit.cxx)
add_executable (LEEmudis LEEmudis.cxx)
add_executable (LEEcovariance LEEcovariance.cxx)

include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (LEEfit SBNfitlib)
target_link_libraries (LEEfit tinyxmllib)

target_link_libraries (LEEcovariance SBNfitlib)
target_link_libraries (LEEcovariance tinyxmllib)

target_link_libraries (LEEmudis SBNfitlib)
target_link_libraries (LEEmudis tinyxmllib)

