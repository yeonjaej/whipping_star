


add_library (SBNfitlib SBNgenerate.cxx SBNmultisim.cxx SBNconfig.cxx SBNspec.cxx SBNfit.cxx SBNchi.cxx prob.cxx SBNfit3pN.cxx SBNosc.cxx SBNprob.cxx )

include_directories ( ${CMAKE_SOURCE_DIR}/inc/)
include_directories ( ${CMAKE_SOURCE_DIR}/src/)
include_directories ( ${CMAKE_SOURCE_DIR}/tinyxml/)

#ROOT_GENERATE_DICTIONARY(LINKDEF mdict.h )
target_link_libraries(SBNfitlib ${ROOT_LIBRARIES})
target_link_libraries(SBNfitlib ${GSL_LIBRARIES})
target_link_libraries(SBNfitlib tinyxmllib)
