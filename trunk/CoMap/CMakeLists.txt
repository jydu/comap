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

ADD_EXECUTABLE(comap ${comap_CPP})
TARGET_LINK_LIBRARIES(comap ${LIBS})

# Install libs
INSTALL(TARGETS comap DESTINATION bin)
