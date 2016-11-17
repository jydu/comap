# Created on 24/05/16 by jdutheil
# Last update on 24/05/16 by jdutheil
# Generates a list of groups with characteristics similar to a given set of groups.
# The output groups have the same size and similar site norm (or any other variable), but sites are taken randomly.
# In this version, site norms are not discretized but a similarity threshold is used.
# To remove bias due to skewed distribution of norms, a correction is further introduced.

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

# Now replicate each group:
x.rep <- numeric(nrow(groups) * nrep)
l.grp <- list()
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
    # Get all sites with similar rate (or any other conditional variable):
    x <- sites[,cond.var]
    t <- abs(x - gp.vals[sit]) / gp.vals[sit]
    condSites <- subset(sites, t <= sim.t)
    # Loop over each simulation replicate:
    for (sim in 1:nrep) {    
      # Remove sites already present in the group
      if (i + sim - 1 <= length(l.grp)) { # test if element already exists in list
        tmp <- subset(condSites, ! (Group %in% l.grp[[i + sim - 1]]))
      } else {
        tmp <- condSites
      }

      # Sampling correction:
      tmpl <- tmp[tmp[,cond.var] < gp.vals[sit],]
      tmpg <- tmp[tmp[,cond.var] > gp.vals[sit],]
      tmpe <- tmp[tmp[,cond.var] == gp.vals[sit],]
      n <- min(nrow(tmpl), nrow(tmpg))
      n <- max(n, min.obs) # we add min.obs there to avoid getting no replicate site when we have extreme values for the candidate site.
      tmp2 <- rbind(tmpl[sample(1:nrow(tmpl), min(n, nrow(tmpl))),], tmpe, tmpg[sample(1:nrow(tmpg), min(n, nrow(tmpg))),])
      if (nrow(tmp2) == 0) {
        warning(paste("No more similar site available for candidate site", sit, "in group", grp, "replicate", sim))
        x.grp[i + sim - 1] <- paste(x.grp[i + sim - 1], "NA", sep = ifelse(x.grp[i + sim - 1] == "[", "", ";"))
        x.ave[i + sim - 1] <- NA
      } else {
        if (nrow(tmp2) < min.obs)
          warning(paste("Minimum site frequency not matched for candidate site", sit, "in group", grp, "replicate", sim))

        # Now sample sites:
        x <- sample(1:nrow(tmp2), 1)
        if (i + sim - 1 <= length(l.grp)) {
          l.grp[[i + sim - 1]] <- append(l.grp[[i + sim - 1]], tmp2[x, "Group"])
        } else {
          l.grp[[i + sim - 1]] <- tmp2[x, "Group"]
        }
        x.grp[i + sim - 1] <- paste(x.grp[i + sim - 1], tmp2[x, "Group"], sep = ifelse(x.grp[i + sim - 1] == "[", "", ";"))
        x.ave[i + sim - 1] <- x.ave[i + sim - 1] + tmp2[x, cond.var]
      }
    }
  }
  x.grp[i:(i + nrep - 1)] <- paste(x.grp[i:(i + nrep - 1)], "]", sep = "")
  i <- i + nrep
}

x.ave <- x.ave / x.siz
results <- data.frame(Replicate = x.rep, Group = x.grp, Size = x.siz, RandMean = x.ave, OrigMean = x.oav)

write.table(results, file = outputPath, sep = "\t", row.names = FALSE)


#### DONE
