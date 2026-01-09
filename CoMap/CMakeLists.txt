# CMake script for Bio++ CoMap
# Author: Julien Dutheil
# Created: 17/09/2009

set(comap_CPP
    CoETools.cpp
    CoMap.cpp
    Cluster.cpp
    ClusterTools.cpp
    AnalysisTools.cpp
    Domain.cpp
    IntervalData.cpp
)

set(mica_CPP Mica.cpp Domain.cpp)

add_executable(comap ${comap_CPP})
add_executable(mica ${mica_CPP})

set(comap-targets comap mica)

foreach(target ${comap-targets})
    # Link (static or shared)
    if(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_STATIC})
        set_target_properties(${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
	target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
    else(BUILD_STATIC)
        target_link_libraries(${target} ${BPP_LIBS_SHARED})
        set_target_properties(
            ${target}
            PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
        )
    endif(BUILD_STATIC)
endforeach(target)

install(TARGETS ${comap-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
