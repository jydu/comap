#! /bin/bash

DATA=Myo
	
echo "Analysing dataset $DATA."

echo "Running CoMap with unweighted mapping: Naive"
/usr/bin/time -v comap param=comap.bpp \
      nijt=naive \
      output.vectors.file=${DATA}_naive.vec >& ${DATA}_naive.out&

echo "Running CoMap with unweighted mapping: Laplace"
/usr/bin/time -v comap param=comap.bpp \
      nijt=laplace \
      output.vectors.file=${DATA}_laplace.vec >& ${DATA}_laplace.out&

echo "Running CoMap with unweighted mapping: Uniformization"
/usr/bin/time -v comap param=comap.bpp \
      nijt=uniformization \
      output.vectors.file=${DATA}_unif.vec >& ${DATA}_unif.out&

echo "Running CoMap with unweighted mapping: Decomposition"
/usr/bin/time -v comap param=comap.bpp \
      nijt=decomposition \
      output.vectors.file=${DATA}_decomp.vec >& ${DATA}_decomp.out&




echo "Running CoMap with weighted mapping: Naive"
/usr/bin/time -v comap param=comap.bpp \
      "nijt=Naive(weight=AAdist(type=grantham, sym=yes))" \
      output.vectors.file=${DATA}_naive_grantham.vec >& ${DATA}_naive_grantham.out&

echo "Running CoMap with weighted mapping: Uniformization"
/usr/bin/time -v comap param=comap.bpp \
      "nijt=Uniformization(weight=AAdist(type=grantham, sym=yes))" \
      output.vectors.file=${DATA}_unif_grantham.vec >& ${DATA}_unif_grantham.out&

echo "Running CoMap with weighted mapping: Decomposition"
/usr/bin/time -v comap param=comap.bpp \
      "nijt=Decomposition(weight=AAdist(type=grantham, sym=yes))" \
      output.vectors.file=${DATA}_decomp_grantham.vec >& ${DATA}_decomp_grantham.out&

echo "Done."

