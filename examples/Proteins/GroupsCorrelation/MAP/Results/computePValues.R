# Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
# 06/06/2006
# 15/11/2007
# 04/08/2007
# 04/12/2009
# 17/02/2011
# 18/04/2015

# ----------------------------------------------------------------------------------------------------------------
# EDITABLE SECTION:
# ----------------------------------------------------------------------------------------------------------------

# File paths.
# Warning: Method names (within [[ ]]) must match between 1) and 2).

# 1) Clustering output

all.data.files<-list()
all.data.files[["Raw"]]       <-"../CoMap/Raw/MAP_groups.csv"
all.data.files[["Grantham"]]  <-"../CoMap/Grantham/MAP_groups.csv"
all.data.files[["Volume"]]    <-"../CoMap/Volume/MAP_groups.csv"
all.data.files[["Polarity"]]  <-"../CoMap/Polarity/MAP_groups.csv"
all.data.files[["Charge"]]    <-"../CoMap/Charge/MAP_groups.csv"


# 2) Simulations

all.sim.files<-list()
all.sim.files[["Raw"]]       <-"../CoMap/Raw/MAP_null.csv"
all.sim.files[["Grantham"]]  <-"../CoMap/Grantham/MAP_null.csv"
all.sim.files[["Volume"]]    <-"../CoMap/Volume/MAP_null.csv"
all.sim.files[["Polarity"]]  <-"../CoMap/Polarity/MAP_null.csv"
all.sim.files[["Charge"]]    <-"../CoMap/Charge/MAP_null.csv"

# 3) Output files:
# Type I test:
output.file<-"MAP_groups_predictions_pvalues.csv"

#Sliding windows sizes:
window.Nmin<-0.2

# General options:

#Maximum p-value level or groups in the output files.
#Use NA if you don't want to filter the output.
level<-0.05
#Minimum number of simulated points required for computing p-value:
min.nobs<-1000
#Log file ("" = terminal)
logFile<-"Cliques.txt"
#FDR (NA if you don't want to perform a FDR analysis):
fdr<-0.05
#Number of replicates to use for computing FDR:
nfdr<-100
#Grid preprocessing (very recommended for pairwise analysis)
gridpp<-TRUE

# Clustering options (will have no effect for pairwise analyses):

#Maximum group size to test:
maxgs<-10
#Correction for nested groups:
cng<-TRUE

# ----------------------------------------------------------------------------------------------------------------
# END OF EDITABLE SECTION
# ----------------------------------------------------------------------------------------------------------------

source("CoMapFunctions.R")

methods<-names(all.data.files)
if(any(methods != names(all.sim.files))) {
 stop("Method names do not match in data and simulations.")
} else {

all.data<-list()
for(m in methods)
{
  cat("Reading groups for method '", m, "'.\n", sep="")
  all.data[[m]]<-read.table(all.data.files[[m]], header=T, sep="\t")
}

all.sim<-list()
for(m in methods)
{
  cat("Reading simulations for method '", m, "'.\n", sep="")
  all.sim[[m]]<-read.table(all.sim.files[[m]], header=T, sep="\t")
}

all.pred<-list()
for(m in methods)
{
  cat("Testing method '", m, "'.\n", sep="")
  all.pred[[m]]<-format.pred(all.data[[m]], all.sim[[m]], 2:maxgs, window.Nmin, min.nobs, gridpp, method=m, level=level, cng=cng, logFile=paste(m,logFile,sep="_"), fdr=fdr, nfdr=nfdr)
}

# Merging all results:
cat("Merging results for all methods...\n")
detected<-character(0)
for(m in methods)
{
  if(nrow(all.pred[[m]]) > 0) detected<-c(detected, m);
}
if(length(detected) > 0)
{
  pred<-all.pred[[detected[1]]]
  if(length(detected) > 1)
  {
    for(i in 2:length(detected))
    {
      pred<-merge(pred, all.pred[[detected[i]]], all=T)
    }
  }
  write.table(pred, file=output.file, sep="\t", row.names=F, quot=F)
} else {
  cat("No group detected!\n");
}

} #END

