# CMake script for Bio++ CoMap
# Author: Julien Dutheil
# Created: 17/09/2009

SET(comap_CPP
 CoETools.cpp
 CoMap.cpp
 Cluster.cpp
 ClusterTools.cpp
 AnalysisTools.cpp
 Domain.cpp
 IntervalData.cpp
)

SET(mica_CPP
 Mica.cpp
 Domain.cpp
)

ADD_EXECUTABLE (comap ${comap_CPP})
ADD_EXECUTABLE (mica  ${mica_CPP})

SET (comap-targets
	comap
	mica
	)

foreach (target ${comap-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    set_target_properties (${target} LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${comap-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

