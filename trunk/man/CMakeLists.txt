# CMake script for the CoMap package.
# Author: Julien Dutheil
# Created: 22/08/2009

IF(MAN)
  INSTALL(FILES comap.1.gz DESTINATION share/man/man1)
  INSTALL(FILES mica.1.gz DESTINATION share/man/man1)
ENDIF(MAN)
