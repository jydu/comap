# CMake script for Bio++ CoMap
# Author: Julien Dutheil
# Created: 17/09/2009

SET(comap_CPP
 CoETools.cpp
 CoMap.cpp
 Cluster.cpp
 ClusterTools.cpp
 AnalysisTools.cpp
)

SET(mica_CPP
 Mica.cpp
)

ADD_EXECUTABLE(comap ${comap_CPP})
TARGET_LINK_LIBRARIES(comap ${LIBS})
SET_TARGET_PROPERTIES(comap PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})

ADD_EXECUTABLE(mica ${mica_CPP})
TARGET_LINK_LIBRARIES(mica ${LIBS})
SET_TARGET_PROPERTIES(mica PROPERTIES LINK_SEARCH_END_STATIC ${BUILD_STATIC})

# Install libs
INSTALL(TARGETS comap mica DESTINATION bin)
