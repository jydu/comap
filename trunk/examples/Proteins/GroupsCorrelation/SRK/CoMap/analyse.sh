#! /bin/sh

DATA=Myo
	
echo "Analysing dataset $DATA."

echo "Running CoMap with unweighted mapping:"
mkdir Simple
comap param=comap.bpp \
      nijt=simple \
      output.vectors.file=Simple/${DATA}.vec \
      clustering.output.tree.file=Simple/${DATA}_clust.dnd \
      clustering.output.groups.file=Simple/${DATA}_groups.csv \
      clustering.null.output.file=Simple/${DATA}_simulations.csv 

echo "Running CoMap with Grantham mapping:"
mkdir Grantham
comap param=comap.bpp \
      nijt=aadist \
      nijt_aadist.type=grantham \
      nijt_aadist.sym=yes \
      output.vectors.file=Grantham/${DATA}.vec \
      clustering.output.tree.file=Grantham/${DATA}_clust.dnd \
      clustering.output.groups.file=Grantham/${DATA}_groups.csv \
      clustering.null.output.file=Grantham/${DATA}_simulations.csv 

echo "Running CoMap with Volume mapping:"
mkdir Volume
comap param=comap.bpp \
      nijt=aadist \
      nijt_aadist.type=grantham.volume \
      nijt_aadist.sym=yes \
      output.vectors.file=Volume/${DATA}.vec \
      clustering.output.tree.file=Volume/${DATA}_clust.dnd \
      clustering.output.groups.file=Volume/${DATA}_groups.csv \
      clustering.null.output.file=Volume/${DATA}_simulations.csv

echo "Running CoMap with Polarity mapping:"
mkdir Polarity
comap param=comap.bpp \
      nijt=aadist \
      nijt_aadist.type=grantham.polarity \
      nijt_aadist.sym=yes \
      output.vectors.file=Polarity/${DATA}.vec \
      clustering.output.tree.file=Polarity/${DATA}_clust.dnd \
      clustering.output.groups.file=Polarity/${DATA}_groups.csv \
      clustering.null.output.file=Polarity/${DATA}_simulations.csv

echo "Running CoMap with Charge mapping:"
mkdir Charge
comap param=comap.bpp \
      nijt=aadist \
      nijt_aadist.type=charge \
      nijt_aadist.sym=yes \
      output.vectors.file=Charge/${DATA}.vec \
      clustering.output.tree.file=Charge/${DATA}_clust.dnd \
      clustering.output.groups.file=Charge/${DATA}_groups.csv \
      clustering.null.output.file=Charge/${DATA}_simulations.csv

echo "Done."

