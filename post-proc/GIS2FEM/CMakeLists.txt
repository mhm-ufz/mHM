##-------------------------------------------------------
##          
##    Tool to convert GIS shape data into FEM data
##                By Wenqing Wang
##
##                     03.2010        
##-------------------------------------------------------

cmake_minimum_required(VERSION 2.4.5)

#IF (WIN32)
	SET ( CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE TYPE INTERNAL FORCE )
#ENDIF (WIN32)

IF (UNIX)
	string( LENGTH "$ENV{DESTDIR}" DESTDIRLEN )
	IF ( ${DESTDIRLEN} GREATER 0 )
		message( STATUS "DESTDIR specified ($ENV{DESTDIR}). Overriding 
		CMAKE_INSTALL_PREFIX (${CMAKE_INSTALL_PREFIX})" )
		set( CMAKE_INSTALL_PREFIX "" CACHE PATH "Use DESTDIR instead of CMAKE_INSTALL_PREFIX" FORCE )
	ENDIF()
ENDIF (UNIX)

# Project name
SET(PROJECT_NAME GIS2FEM)
PROJECT( ${PROJECT_NAME}  )
INCLUDE( Cmake/base.cmake )

# Set build directories
# Binaries are created in /bin and libraries in /lib
SET( EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin )
#SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )
IF (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/release/g2f)
ELSE (MSVC)
	SET(OGS_EXECUTABLE ${EXECUTABLE_OUTPUT_PATH}/g2f)
ENDIF (MSVC)

# Relative paths are used in the created projects
# May not work correcty due to the CMake docs
#SET( CMAKE_USE_RELATIVE_PATHS TRUE )

IF (MSVC)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
	ADD_DEFINITIONS(-D_CRT_NONSTDC_NO_WARNINGS)
ENDIF (MSVC)
#IF(UNIX)
##	ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
##	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-deprecated -Wall -Wextra")
#ENDIF (UNIX)
IF(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -Wall -Wextra")
        ADD_DEFINITIONS( -std=c++11 )
ENDIF()

# Create Configure.h
CONFIGURE_FILE (Cmake/Configure.h.in ${PROJECT_BINARY_DIR}/Configure.h)



SET(OGS_PATH ../..)
SET(OGS_FEM_PATH ${OGS_PATH}/FEM)
SET(OGS_MSH_PATH ${OGS_PATH}/MSH)
SET(OGS_GEO_PATH ${OGS_PATH}/GEO)
SET(OGS_MATH_PATH ${OGS_PATH}/MathLib)
SET(OGS_BASE_PATH ${OGS_PATH}/Base)

INCLUDE_DIRECTORIES(
.
${OGS_FEM_PATH}
${OGS_MSH_PATH}
${OGS_GEO_PATH}
${OGS_MATH_PATH}
${OGS_BASE_PATH}
# This is for Configure.h 
${PROJECT_BINARY_DIR}
)

ADD_DEFINITIONS(-DNONFLUID_MOMENTUM 
		-DNON_GEO
		-DNON_PROCESS
		-DNEW_EQS)

SET(GIS2FEM_SOURCES main.cpp misc.cpp geo.cpp
   ${OGS_FEM_PATH}/fem_ele.cpp ${OGS_FEM_PATH}/femlib.cpp
   ${OGS_FEM_PATH}/mathlib.cpp ${OGS_FEM_PATH}/matrix_class.cpp ${OGS_MSH_PATH}/msh_core.cpp  ${OGS_MSH_PATH}/msh_edge.cpp
   ${OGS_MSH_PATH}/msh_node.cpp ${OGS_MSH_PATH}/msh_elem.cpp ${OGS_MSH_PATH}/msh_mesh.cpp ${OGS_MSH_PATH}/MSHEnums.cpp
   misc.h geo.h
   ${OGS_FEM_PATH}/matrix_class.h  ${OGS_FEM_PATH}/fem_ele.h  ${OGS_FEM_PATH}/mathlib.h ${OGS_MSH_PATH}/msh_core.h  ${OGS_MSH_PATH}/msh_edge.h
   ${OGS_MSH_PATH}/msh_node.h ${OGS_MSH_PATH}/msh_elem.h ${OGS_MSH_PATH}/msh_mesh.h ${OGS_MSH_PATH}/MSHEnums.h
  )

  
add_executable(${PROJECT_NAME}  ${GIS2FEM_SOURCES})
