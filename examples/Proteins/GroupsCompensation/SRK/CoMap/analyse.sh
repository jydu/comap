#! /bin/sh

DATA=SRK
	
echo "Analysing dataset $DATA."

echo "Running CoMap with Grantham mapping:"
mkdir Grantham
comap param=comap.bpp \
      "nijt=Uniformization(weight=Grantham(symmetrical=no))" \
      DATA=$DATA OUTPUT=Grantham

echo "Running CoMap with Volume mapping:"
mkdir Volume
comap param=comap.bpp \
      "nijt=Uniformization(weight=Diff(index1=GranthamVolume, symmetrical=no))" \
      DATA=$DATA OUTPUT=Volume

echo "Running CoMap with Polarity mapping:"
mkdir Polarity
comap param=comap.bpp \
      "nijt=Uniformization(weight=Diff(index1=GranthamPolarity, symmetrical=no))" \
      DATA=$DATA OUTPUT=Polarity

echo "Running CoMap with Charge mapping:"
mkdir Charge
comap param=comap.bpp \
      "nijt=Uniformization(weight=Diff(index1=KleinCharge, symmetrical=no))" \
      DATA=$DATA OUTPUT=Charge

echo "Done."

