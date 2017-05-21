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

FOREACH (target ${bppsuite-targets})
  # Link (static or shared)
  IF (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_STATIC})
    SET_TARGET_PROPERTIES (${target} LINK_SEARCH_END_STATIC TRUE)
  ELSE (BUILD_STATIC)
    TARGET_LINK_LIBRARIES (${target} ${BPP_LIBS_SHARED})
  ENDIF (BUILD_STATIC)
ENDFOREACH (target)

INSTALL (TARGETS ${comap-targets} DESTINATION bin)
