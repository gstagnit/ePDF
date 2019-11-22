cmake_minimum_required (VERSION 2.6)
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# the project name
project(ePDF)

# activating some global properties for the project
set(ePDF_VERSION 1.3.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/inc/ePDF/version.h.in"
  "${PROJECT_BINARY_DIR}/inc/ePDF/version.h"
)

find_package(PkgConfig REQUIRED)
pkg_search_module(YAML REQUIRED yaml-cpp)

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file(
  "${PROJECT_SOURCE_DIR}/bin/ePDF-config.in"
  "${PROJECT_SOURCE_DIR}/bin/ePDF-config"
)

include_directories(${PROJECT_SOURCE_DIR}/inc)

# GSL
find_program(GSL_CONFIG gsl-config REQUIRED)
if (GSL_CONFIG)
  exec_program(${GSL_CONFIG}
    ARGS --cflags
    OUTPUT_VARIABLE GSL_CXX_FLAGS
  )
  set(GSL_CXX_FLAGS ${GSL_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${GSL_CONFIG}
    ARGS --libs
    OUTPUT_VARIABLE GSL_LIBRARIES
  )
  set(GSL_LIBRARIES ${GSL_LIBRARIES} CACHE STRING INTERNAL)
endif(GSL_CONFIG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS} ${YAML_CFLAGS} -std=c++11 -Wall -fPIC -Wunused")

# build target library
add_subdirectory(src)
add_subdirectory(run)

install(FILES ${PROJECT_SOURCE_DIR}/bin/ePDF-config
	   DESTINATION bin
	   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
