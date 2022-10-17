#! /bin/sh
arch=`uname -m`
version=1.5.5-1

strip CoMap/comap
strip CoMap/mica
tar cvzf comap-${arch}-bin-static-${version}.tar.gz CoMap/comap CoMap/mica

