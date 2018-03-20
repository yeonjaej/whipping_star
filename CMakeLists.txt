cmake_minimum_required (VERSION 2.6)
project (SBNfit)
# The version number.
set (SBNfit_VERSION_MAJOR 2)
set (SBNfit_VERSION_MINOR 0)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
          set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
    else ()
	    set (CMAKE_CXX_STANDARD 11)
endif ()


FIND_PACKAGE(ROOT)# COMPONENTS Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread)

if(ROOT_FOUND) 
	message("Found ROOT using CMAKE FIND_PACKAGE")
	include(${ROOT_USE_FILE})
ELSE()
	message("ROOT not compiled using CMAKE, use root-config")
	#Get our ROOT libraries and such
	FIND_PROGRAM(ROOT_CONFIG root-config)
	EXEC_PROGRAM(${ROOT_CONFIG} ARGS --cflags --glibs OUTPUT_VARIABLE ROOT_CXX_FLAGS)
	message("ROOT cxx flags set to: ${ROOT_CXX_FLAGS}")
ENDIF()



#Get our GSL libraries
FIND_PROGRAM(GSL_CONFIG gsl-config)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --cflags OUTPUT_VARIABLE GSL_CXX_FLAGS)
EXEC_PROGRAM(${GSL_CONFIG} ARGS --libs OUTPUT_VARIABLE GSL_LIBS)
set (CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} ${GSL_LIBS} ${GSL_CXX_FLAGS}")


# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories(${PROJECT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc/)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/tinyxml/)
#ROOT_GENERATE_DICTIONARY(G__Event MCEventWeight.h LINKDEF MCEventWeightLinkDef.h)

#add_subdirectory(inc)
add_subdirectory(tinyxml)
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(miniboone)
add_subdirectory(lee)
