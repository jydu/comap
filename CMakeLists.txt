# CMake script for CoMap
# Author: Julien Dutheil
# Created: 22/08/2009

# Global parameters
cmake_minimum_required(VERSION 3.10)
project(comap CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Compile options
add_compile_options(
    -Wall
    -Weffc++
    -Wshadow
    -Wconversion
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif()

if(NOT NO_DEP_CHECK)
    set(NO_DEP_CHECK
        FALSE
        CACHE BOOL
        "Disable dependencies check for building distribution only."
        FORCE
    )
endif(NOT NO_DEP_CHECK)

if(NO_DEP_CHECK)
    message(
        STATUS
        "Dependencies checking disabled. Only distribution can be built."
    )
else(NO_DEP_CHECK)
    #static linkage?
    if(NOT BUILD_STATIC)
        set(BUILD_STATIC FALSE CACHE BOOL "Enable static linkage." FORCE)
    endif()
    if(BUILD_STATIC)
        message(STATUS "Static linkage requested.")
    endif()

    if(CMAKE_INSTALL_PREFIX)
        set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
    endif(CMAKE_INSTALL_PREFIX)

    # Check compression program
    # COMPRESS_PROGRAM controls the choice of program
    # COMPRESS_EXT can be used to override the file extension
    if(NOT COMPRESS_PROGRAM)
        set(COMPRESS_PROGRAM
            gzip
            CACHE STRING
            "Set program for compressing documentation"
            FORCE
        )
    endif()
    find_program(
        COMPRESS_BIN
        NAMES ${COMPRESS_PROGRAM}
        DOC "${COMPRESS_PROGRAM} compression program"
    )
    if(NOT COMPRESS_BIN)
        message(
            STATUS
            "${COMPRESS_PROGRAM} program not found, text doc will not be compressed"
        )
    else()
        # Deduce COMPRESS_EXT for known compression programs if not set
        if(NOT COMPRESS_EXT)
            if(${COMPRESS_PROGRAM} STREQUAL "gzip")
                set(COMPRESS_EXT "gz")
            elseif(${COMPRESS_PROGRAM} STREQUAL "bzip2")
                set(COMPRESS_EXT "bz2")
            else()
                set(COMPRESS_EXT "${COMPRESS_PROGRAM}") # Default: program name (works for xz/lzma)
            endif()
        endif()
        # Generate command line args (always add -c to output compressed file to stdout)
        if(${COMPRESS_PROGRAM} STREQUAL "gzip")
            # -n for no timestamp in files (reproducible builds)
            set(COMPRESS_ARGS -c -n)
        else()
            set(COMPRESS_ARGS -c)
        endif()
        message(
            STATUS
            "Found ${COMPRESS_BIN} compression program, using file extension .${COMPRESS_EXT}"
        )
    endif()

    # Find dependencies (add install directory to search)
    if(CMAKE_INSTALL_PREFIX)
        set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
    endif(CMAKE_INSTALL_PREFIX)

    include(GNUInstallDirs)
    find_package(bpp-phyl 3.0.0 REQUIRED)

    # Subdirectories
    add_subdirectory(CoMap)
    add_subdirectory(doc)
    add_subdirectory(man)
endif(NO_DEP_CHECK)

# Packager
set(CPACK_PACKAGE_NAME "comap")
set(CPACK_PACKAGE_VENDOR "Julien Y. Dutheil")
set(CPACK_PACKAGE_VERSION "1.6.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "6")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "CoMap: a program for mapping substitutions and detecting molecular co-evolution."
)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
set(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
set(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
set(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/\\\\.git/"
    "/\\\\.gitignore"
    ${CPACK_SOURCE_IGNORE_FILES}
)

set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)
set(CPACK_DEBSOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig"
)
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

if(NOT NO_DEP_CHECK)
    if(UNIX)
        #This creates rpm packages:
        #ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
        #ADD_DEPENDENCIES(rpm dist info man)

        #ADD_DEPENDENCIES(info install)
    endif(UNIX)
endif(NOT NO_DEP_CHECK)
