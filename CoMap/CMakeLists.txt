# CMake script for Bio++ Program Suite
# Author: Julien Dutheil
# Created: 17/09/2009

set(comap_CPP
 CoETools.cpp
 CoMap.cpp
 Cluster.cpp
 ClusterTools.cpp
 AnalysisTools.cpp
)

add_executable(comap ${comap_CPP})
target_link_libraries(comap ${LIBS})

# Install libs
install(TARGETS comap DESTINATION bin)
