# Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
# 06/06/2006
# 15/11/2007

# ----------------------------------------------------------------------------------------------------------------
# Parameters to edit:
# ----------------------------------------------------------------------------------------------------------------

# File paths.
# Warning: Method names (within [[ ]]) must match between 1) and 2).

# 1) Clustering output

all.data.files<-list()
all.data.files[["Simple"]]    <-"../Coe/Simple/Myo_groups.csv"
all.data.files[["Grantham"]]  <-"../Coe/Grantham/Myo_groups.csv"
all.data.files[["Volume"]]    <-"../Coe/Volume/Myo_groups.csv"
all.data.files[["Polarity"]]  <-"../Coe/Polarity/Myo_groups.csv"
all.data.files[["Charge"]]    <-"../Coe/Charge/Myo_groups.csv"


# 2) Simulations

all.sim.files<-list()
all.sim.files[["Simple"]]    <-"../Coe/Simple/Myo_simulations.csv"
all.sim.files[["Grantham"]]  <-"../Coe/Grantham/Myo_simulations.csv"
all.sim.files[["Volume"]]    <-"../Coe/Volume/Myo_simulations.csv"
all.sim.files[["Polarity"]]  <-"../Coe/Polarity/Myo_simulations.csv"
all.sim.files[["Charge"]]    <-"../Coe/Charge/Myo_simulations.csv"

# 3) Output files:
# Type I test:
output.file1<-"Myo_predictions_pvalues.csv"
# type II test:
output.file2<-"Myo_predictions2_pvalues.csv"

#Sliding windows sizes:
window.Nmin=0.2
window.Delta=0.2

# General options:

#Maximum group size to test:
maxgs<-20
#Maximum p-value level or groups in the output files:
level<-0.01
#Correction for nested groups:
cng<-TRUE
#Log file ("" = terminal)
logFile<-"Cliques.txt"
logFile2<-"Cliques2.txt"

# ----------------------------------------------------------------------------------------------------------------
# END OF PARAMETER TO EDIT
# ----------------------------------------------------------------------------------------------------------------

source("CoMapFunctions.R")

methods=names(all.data.files)
if(any(methods != names(all.sim.files))) {
 stop("Method names do not match in data and simulations.")
} else {
  
all.data<-list()
for(m in methods)
{
  cat("Reading groups for method", m, ".\n")
  all.data[[m]]<-read.table(all.data.files[[m]], header=T, sep="\t")
}

all.sim<-list()
for(m in methods)
{
  cat("Reading simulations for method", m, ".\n")
  all.sim[[m]]<-read.table(all.sim.files[[m]], header=T, sep="\t")
}

all.pred<-list()
for(m in methods)
{
  cat("Testing method", m, "with type I test.\n")
  all.pred[[m]]<-get.pred(all.data[[m]], all.sim[[m]], 2:maxgs, window.Nmin, m, level, cng, logFile)
}

# Merging all results:
cat("Merging results for all methods (Type I tests).\n")
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
  write.table(pred, file=output.file1, sep="\t", row.names=F, quot=F)
}
else
{
  cat("No group detected with test I!\n");
}

all.pred2<-list()
for(m in methods)
{
  cat("Testing method", m, "with type II test.\n")
  all.pred2[[m]]<-get.pred2(all.data[[m]], all.sim[[m]], 2:maxgs, window.Nmin, window.Delta, m, level, cng, logFile2)
}

# Merging all results:
cat("Merging results for all methods (Type II tests).\n")
detected2<-character(0)
for(m in methods)
{
  if(nrow(all.pred2[[m]]) > 0) detected2<-c(detected2, m);
}
if(length(detected2) > 0)
{
  pred2<-all.pred2[[detected2[1]]]
  if(length(detected2) > 1)
  {
    for(i in 2:length(detected2))
    {
      pred2<-merge(pred2, all.pred2[[detected2[i]]], all=T)
    }
  }
  write.table(pred2, file=output.file2, sep="\t", row.names=F, quot=F)
}
else
{
  cat("No group detected with test II!\n");
}

} #END

