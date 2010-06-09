# CMake script for CoMap
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.4)
PROJECT(CoMap CXX)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF()

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF()

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF()

#static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE("-- Static linkage requested.")
  SET(CMAKE_CXX_FLAGS "-static -static-libgcc ${CMAKE_CXX_FLAGS}")
ENDIF()

#build info?
IF(NOT DEFINED INFO)
  SET(INFO UNIX)
ENDIF(NOT DEFINED INFO)

#find executables for documentation
FIND_PROGRAM(MAKEINFO_EXE NAMES makeinfo)
FIND_PROGRAM(PDFTEX_EXE NAMES pdftex)

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  UNSET(${lib_name}_LIBRARY CACHE) #Reset the cache.

  IF(BUILD_STATIC)
    SET(${lib_name}_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  ELSE(BUILD_STATIC)
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  ENDIF(BUILD_STATIC)

  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes : ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   libraries: ${${lib_name}_LIBRARY}")
  ELSE(${lib_name}_LIBRARY)
    MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
  ENDIF(${lib_name}_LIBRARY)
  
  #add the dependency:
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

#Find the libraries. The order is very important for staitc linkage, it won't
#link if you change it!
IMPROVED_FIND_LIBRARY(LIBS bpp-phyl Phyl/Tree.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-seq Seq/Alphabet.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-numcalc NumCalc/NumTools.h)
IMPROVED_FIND_LIBRARY(LIBS bpp-utils Utils/Clonable.h)

IF(MAKEINFO_EXE)
  MESSAGE("-- Found makeinfo here: ${MAKEINFO_EXE}")
  MESSAGE("   Adding targets: info, html")

  IF(INFO)
    SET(ADD_INFO_TO "ALL")
    MESSAGE("   Adding target info to target all")
  ENDIF(INFO)

  ADD_CUSTOM_TARGET(info
    ${ADD_INFO_TO}
    COMMAND ${MAKEINFO_EXE} comap.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )

  ADD_CUSTOM_TARGET(html
    COMMAND ${MAKEINFO_EXE} --html --css-ref=http://www.w3.org/StyleSheets/Core/Steely bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(MAKEINFO_EXE)

IF(PDFTEX_EXE)
  MESSAGE("-- Found pdftex here: ${PDFTEX_EXE}")
  MESSAGE("   Adding target: pdf")

  ADD_CUSTOM_TARGET(pdf
    COMMAND ${PDFTEX_EXE} bppsuite.texi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
    )
ENDIF(PDFTEX_EXE)

# Subdirectories
ADD_SUBDIRECTORY(CoMap)
ADD_SUBDIRECTORY(doc)

# Packager
SET(CPACK_PACKAGE_NAME "CoMap")
SET(CPACK_PACKAGE_VENDOR "CNRS - University Montpellier 2")
SET(CPACK_PACKAGE_VERSION "1.4.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "4")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "CoMap: a program for mapping substitutions and detecting molecular co-evolution.")
SET(CPACK_RESOURCE_FILE_LICENSE "COPYING")
SET(CPACK_RESOURCE_FILE_AUTHORS "AUTHORS")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.svn"
 ".*\\\\.gz"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 "\\\\..*\\\\.swp"
 "CoMap/\\\\..*"
 "CoMap/comap"
 "doc/comap/"
 "doc/comap.info"
 "debian/tmp"
 "debian/comap"
 ${CPACK_SOURCE_IGNORE_FILES}
)
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
INCLUDE(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(deb fakeroot debian/rules binary)
ADD_DEPENDENCIES(deb dist info)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist info)

#ADD_DEPENDENCIES(info install)
ENDIF(UNIX)
