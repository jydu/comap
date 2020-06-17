<!-- Created 16/06/20 by jdutheil -->

Run MI analysis
===============


Non-parametric bootstrap
------------------------

MI is conditioned on the minimum site-entropy. No model and no phylogenetic tree is used.
```bash
mica param=mica.bpp null.method=nonparametric-bootstrap
```

Results:
```r
d<-read.table("Myglobin.MI.csv", header=T)
d[d$Bs.p.value < 0.005,]
```

Most significant pairs:
```
         Group       MI       APC        RCW   Hjoint     Hmin Bs.p.value Bs.nb
1181 [173;292] 0.598094 0.5673500 0.05265180 1.010390 0.651822 0.00432126 20132
1220 [173;338] 0.610926 0.5110380 0.04742580 0.743630 0.651822 0.00273183 20132
1224 [173;342] 0.584071 0.5371970 0.04985340 1.265080 0.651822 0.00496697 20132
4286 [241;292] 0.461762 0.4118340 0.03821940 0.994854 0.499953 0.00416408 11286
5347 [254;279] 0.493918 0.3925950 0.03643400 0.809008 0.493918 0.00194915 11286
5352 [254;292] 0.493918 0.4323050 0.04011920 0.956663 0.493918 0.00310091 11286
6662 [289;292] 0.282347 0.2042270 0.01895290 0.998252 0.323936 0.00282486 18053
6709 [289;347] 0.282347 0.1740890 0.01615600 1.121050 0.323936 0.00282486 18053
6712 [289;350] 0.268556 0.0466599 0.00433017 0.323936 0.268556 0.00398804 18053
```

Parametric bootstrap
--------------------

A phylogeny and a model are used. MI is conditionned on the minimum norm of the substitution vector.
```bash
mica param=mica.bpp null.method=parametric-bootstrap use_model=yes
```

