# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.

add_executable (LEEfit LEEfit.cxx)
add_executable (buildcovar build_microboone_covariance.cxx)


# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
include_directories ( ${CMAKE_SOURCE_DIR}/inc)
target_link_libraries (LEEfit SBNfitlib)
target_link_libraries (LEEfit tinyxmllib)

target_link_libraries (buildcovar SBNfitlib)
target_link_libraries (buildcovar tinyxmllib)

