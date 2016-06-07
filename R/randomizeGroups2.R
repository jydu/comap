# Created on 24/05/16 by jdutheil
# Last update on 24/05/16 by jdutheil
# Generates a list of groups with characteristics similar to a given set of groups.
# The output groups have the same size and similar site norm (or any other variable), but sites are taken randomly.
# In this version, site norms are not discretized but a similarity threshold is used.

# Input:
sitesPath <- "Myo_sites.csv"
groupsPath <- "Myo_stats_pvalues.csv"

# Variable to codition over
cond.var <- "N" #Norm

# Number of randomizations to perform for each group:
nrep <- 100

# Similarity threshold used for norm:
sim.t <- 0.1 #10% difference: (Nsim - Nobs) / Nobs < 0.1
min.obs <- 5 #Minimum number of matching site. If not enough sites are found for a given position, a warning will be issued. You may try to increase the threshold if too many warnings are produced. 

# Output path:
outputPath <- "Myo_random.csv"

###########################################

# Read input data:
sites <- read.table(sitesPath, header = TRUE, stringsAsFactors = FALSE)
sites$Group <- substr(sites$Group, 2, nchar(sites$Group) - 1)
row.names(sites) <- sites$Group
groups <- read.table(groupsPath, header = TRUE, stringsAsFactors = FALSE)
groupsLst <- strsplit(substr(groups$Group, 2, nchar(groups$Group) - 1), ";")

# Set of sites available for each replicate:
sitesSet <- list()
for (i in 1:nrep) sitesSet[[i]] <- sites

# Now replicate each group:
x.rep <- numeric(nrow(groups) * nrep)
x.grp <- character(nrow(groups) * nrep)
x.ave <- numeric(nrow(groups) * nrep) #Average of the sampled group
x.siz <- numeric(nrow(groups) * nrep)
x.oav <- numeric(nrow(groups) * nrep) #Average the original group
i <- 1
for (grp in 1:nrow(groups)) {
  size <- groups[grp, "Size"]
  nmin <- groups[grp, "Nmin"]
  
  # Get all sites with adequate value for each position:
  gp <- groupsLst[[grp]]
  if (length(gp) != size) stop("!!! Error in input file, group size does not match number or sites!")
  gp.vals <- numeric(size)
  for (j in 1:size) {
    gp.vals[j] <- sites[gp[j], cond.var]
  }

  x.rep[i:(i + nrep - 1)] <- 1:nrep
  x.siz[i:(i + nrep - 1)] <- size
  x.grp[i:(i + nrep - 1)] <- "["
  x.ave[i:(i + nrep - 1)] <- 0
  x.oav[i:(i + nrep - 1)] <- mean(gp.vals)
  
  # Loop over each site in the group:
  for (sit in 1:size) {
    # Loop over each simulation replicate:
    for (sim in 1:nrep) {
      x <- sitesSet[[sim]][,cond.var]
      t <- abs(x - gp.vals[sit]) / gp.vals[sit]
    
      # Get all sites with sufficient variability:
      tmp <- subset(sitesSet[[sim]], t <= sim.t)
      if (nrow(tmp) == 0)
        warning(paste("No more site available for site", sit, "in group", grp, "replicate", sim))
      if (nrow(tmp) < min.obs)
        warning(paste("Minimum site frequency not matched for site", sit, "in group", grp, "replicate", sim))

      # Now sample sites:
      x <- sample(1:nrow(tmp), 1)
      x.grp[i + sim - 1] <- paste(x.grp[i + sim - 1], paste(tmp[x, "Group"], sep = ";"), sep = ifelse(x.grp[i + sim - 1] == "[", "", ";"))
      x.ave[i + sim - 1] <- x.ave[i + sim - 1] + sum(tmp[x, cond.var])

      # As we sample without replacement:
      sitesSet[[sim]] <- subset(sitesSet[[sim]], Group != tmp[x, "Group"])
    }
  }
  x.grp[i:(i + nrep - 1)] <- paste(x.grp[i:(i + nrep - 1)], "]", sep = "")
  i <- i + nrep
}

x.ave <- x.ave / x.siz
results <- data.frame(Replicate = x.rep, Group = x.grp, Size = x.siz, RandMean = x.ave, OrigMean = x.oav)

write.table(results, file = outputPath, sep = "\t")


#### DONE
