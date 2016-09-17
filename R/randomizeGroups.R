# Created on 23/03/15 by jdutheil
# Last update on 19/10/15 by jdutheil
# Generates a list of groups with characteristics similar to a given set of groups.
# The output groups have the same size and similar site norm, but sites are taken randomly.

# Input:
sitesPath <- "Myo_sites.csv"
groupsPath <- "Myo_stats_pvalues.csv"

# Number of randomizations to perform for each group:
nrep <- 10

# Number of classes for norm discretization:
nclass <- 5

# Output path:
outputPath <- "Myo_random.csv"

###########################################

# Read input data:
sites <- read.table(sitesPath, header = TRUE, stringsAsFactors = FALSE)
sites$Group <- substr(sites$Group, 2, nchar(sites$Group) - 1)
row.names(sites) <- sites$Group
groups <- read.table(groupsPath, header = TRUE, stringsAsFactors = FALSE)
groupsLst <- strsplit(substr(groups$Group, 2, nchar(groups$Group) - 1), ";")

# Get the distribution of norms in all sites:
nBounds <- quantile(sites$N, prob=seq(0, 1, len = nclass + 1))
nBounds[nclass + 1] <- Inf

# Assign all sites a rate class:
for (i in 1:nrow(sites)) {
  sites[i, "NRC"] <- max(which(sites[i, "N"] >= nBounds))
}

# Now replicate each group:
x.rep <- numeric(nrow(groups) * nrep)
x.grp <- character(nrow(groups) * nrep)
x.min <- numeric(nrow(groups) * nrep)
x.siz <- numeric(nrow(groups) * nrep)
x.omi <- numeric(nrow(groups) * nrep)
k <- 1
for (i in 1:nrow(groups)) {
  size <- groups[i, "Size"]
  nmin <- groups[i, "Nmin"]
  
  # Get all sites with adequate norm for each position:
  gp <- groupsLst[[i]]
  gpRc <- numeric(length(gp))
  for (j in 1:length(gp)) {
    gpRc[j] <- sites[gp[j], "NRC"]
  }

  x.rep[k:(k + nrep - 1)] <- 1:nrep
  x.siz[k:(k + nrep - 1)] <- size
  x.grp[k:(k + nrep - 1)] <- "["
  x.min[k:(k + nrep - 1)] <- Inf
  x.omi[k:(k + nrep - 1)] <- nmin
  #print(x.grp)
  
  tbl <- table(gpRc)
  for (j in 1:length(tbl)) {
    rc <- as.numeric(names(tbl[j]))
    rf <- tbl[j]
    # Get all sites with sufficient variability:
    tmp <- subset(sites, NRC == rc)
    if (nrow(tmp) < rf)
      stop("Error, input files are inconsistant :(\n")

    for (l in k:(k + nrep - 1)) {
      x <- sample(1:nrow(tmp), rf)
      x.grp[l] <- paste(x.grp[l], paste(tmp[x, "Group"], collapse = ";"), sep = ifelse(x.grp[l] == "[", "", ";"))
      x.min[l] <- min(x.min[l], tmp[x, "N"])
    }
  }
  x.grp[k:(k + nrep - 1)] <- paste(x.grp[k:(k + nrep - 1)], "]", sep = "")
  k <- k + nrep
}

results <- data.frame(Replicate = x.rep, Group = x.grp, Size = x.siz, Nmin = x.min)

write.table(results, file = outputPath, sep = "\t")

#rates <- data.frame(NminOrig = x.omi, Size = x.siz, NminRand = x.min)

#plot(NminRand~NminOrig, rates)
#abline(0, 1)

#### DONE
