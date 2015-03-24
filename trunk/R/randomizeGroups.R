# Created on 23/03/15 by jdutheil
# Last update on 23/03/15 by jdutheil
# Generates a list of groups with characteristics similar to a given set of groups.
# The output groups have the same size and minimum site norm, but sites are taken randomly.

# Input:
sitesPath <- "Myo_sites.csv"
groupsPath <- "Myo_stats_pvalues.csv"

# Number of randomizations to perform for each group:
nrep <- 10

# Output path:
outputPath <- "Myo_random.csv"

###########################################

# Read input data:
sites <- read.table(sitesPath, header = TRUE, stringsAsFactors = FALSE)
sites$Group <- substr(sites$Group, 2, nchar(sites$Group) - 1)
groups <- read.table(groupsPath, header = TRUE, stringsAsFactors = FALSE)

# Now replicate each group:
x.rep <- numeric(nrow(groups) * nrep)
x.gpr <- character(nrow(groups) * nrep)
x.min <- numeric(nrow(groups) * nrep)
x.siz <- numeric(nrow(groups) * nrep)

k <- 0
for (i in 1:nrow(groups)) {
  nmin <- groups[i, "Nmin"]
  size <- groups[i, "Size"]
  # Get all sites with sufficient variability:
  tmp <- subset(sites, N >= nmin)
  if (nrow(tmp) < size)
    stop("Error, input files are inconsistant :(\n")
  for (j in 1:nrep) {
    k <- k + 1
    x <- sample(1:nrow(tmp), size)
    x.rep[k] <- j
    x.gpr[k] <- paste("[", paste(tmp[x, "Group"], collapse = ";"), "]", sep = "")
    x.siz[k] <- size
    x.min[k] <- min(tmp[x, "N"])
  }
}
results <- data.frame(Replicate = x.rep, Group = x.gpr, Size = x.siz, Nmin = x.min)

write.table(results, file = outputPath, sep = "\t")

#### DONE
