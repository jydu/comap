#! /bin/bash

VERSION=1.4.0
FILTER="--exclude */Simple* --exclude */Grantham* --exclude */Polarity* --exclude */Volume* --exclude */Charge* --exclude *.svn* --exclude *.out"

cd examples/Proteins
cd GroupsCompensation
zip -r -FS MAP-$VERSION.zip MAP $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*
zip -r -FS Myoglobin-$VERSION.zip Myoglobin $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*
zip -r -FS SRK-$VERSION.zip SRK $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*

cd ../GroupsCorrelation
zip -r -FS MAP-$VERSION.zip MAP $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*
zip -r -FS Myoglobin-$VERSION.zip Myoglobin $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*
zip -r -FS SRK-$VERSION.zip SRK $FILTER --exclude *CoMap/MAP* --exclude *Results/MAP*



