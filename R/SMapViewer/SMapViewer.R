# Created on 10/07/14 by jdutheil
require(ape)     # for trees
require(seqinr)  # for sequences
require(ggplot2) # for plotting heatmaps
require(reshape2)# for data manipulation
require(plyr)    # for summary statistics
require(treeio)  # for parsing trees and attach data
require(ggtree)  # for plotting trees with ggplot
require(phangorn)# for midpoint rooting

read.mapping<-function(file) {
  mapping<-read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  class(mapping)<-c("mapping", class(mapping))
  mapping$Branches<-as.character(mapping$Branches)
  return(mapping)
}

attach.tree<-function(mapping, treefile, tlnfile) {
  tree<-read.tree(treefile)
  #We check if the tree is compatible:
  if (any(!mapping$Branches %in% c(tree$node.label, tree$tip.label))) {
    stop("Error, tree is not compatible with this mapping.")
  }
  attr(mapping, "treedata") <- tree
  tln<-read.table(tlnfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  attr(mapping, "tiptln") <- tln
  index<-c(tree$tip.label, tree$node.label)
  attr(mapping, "index") <- index
  return(mapping)
}

getIds<-function(labels, index) {
  ids<-numeric(length(labels))
  for (i in 1:length(labels)) {
    ids[i] <- which(index == labels[i])
  }
  return(ids)
}

attach.sequences<-function(mapping, seqfile, ...) {
  seqs<-read.alignment(seqfile, ...)
  # Ancestral sequence have node ids as names,
  # but extant species have their real names.
  # We therefore first need to convert extant
  # species names to node id. Another issue is
  # that potential additional extant species
  # might be included, if they were present in
  # the original alignment but not in the tree!
  tln<-attr(mapping, "tiptln") # => gets the node ids of the analyzed leaves.
  # We first get read of the putative additional sequences:
  tree<-attr(mapping, "tree")
  keepers <- seqs$nam %in% tln$Name | seqs$nam %in% tree$node.label
  cat(sum(!keepers), "sequences not in the tree were discarded.\n")
  seqs$nb  <- sum(keepers)
  seqs$seq <- seqs$seq[keepers]
  seqs$nam <- seqs$nam[keepers]
  seqnames <- seqs$nam
  # Now translates extant sequences to node ids:
  for (i in 1:length(seqs$nam)) {
    j<-which(tln$Name == seqs$nam[i])
    if (length(j) == 1) seqnames[i] <- tln[j, "Id"]
    else if (length(j) > 1) stop(paste("Duplicated sequence", seqnames[i], "in translation table.\n"))
  }
  #ids<-getIds(seqnames, attr(mapping, "index"))
  #seqs$nam<-ids
  seqs$nam <- seqnames
  attr(mapping, "sequences")<-as.matrix(seqs)
  return(mapping)
}

plot.mapping<-function(mapping, sites, drawing = "tree", top = -1,
                       show.states = TRUE, show.branch.labels = FALSE, root = FALSE,
                       nsim = 100, threshold = -1, ...) {
  counts<-mapping[,paste("Site", sites, sep = "")]
  s<-paste("Site", sites, sep="")
  vec<-subset(counts, select = s)
  vec$Branches <- mapping$Branches
  vec$BrLen <- mapping$Mean

  # We scale colors according to the selected positions only
  mx<-max(abs(range(counts)))
  # Compute compensation coefficient for each branch:
  sumvec <- data.frame(Branches = vec$Branches, Sum = apply(vec[, s], 1, sum), SumAbs = apply(abs(vec[, s]), 1, sum), stringsAsFactors = FALSE)
  cc <- with(sumvec, abs(Sum) / SumAbs)
  # Sort branches according to compensation signal
  vec <- vec[order(cc),]
  sumvec <- sumvec[order(cc),]
  vec$Branches <- ordered(vec$Branches, levels = vec$Branches)
  sumvec$Branches <- ordered(sumvec$Branches, levels = sumvec$Branches)
  # Select branches for display:
  t <- apply(abs(vec[, s]), 1, max) >= threshold;
  if (top > 0) {
    b <- vec$Branches[t][1:top]
  } else {
    b <- vec$Branches[t]
  }

  if (drawing == "tree") {
    tree <- attr(mapping, "tree")
    if (root) tree <- midpoint(tree)
    tree <- as.treedata(tree)
    # Get tree index:
    index<-attr(mapping, "index")
    # Get the corresponding vector:
    pp <- list()
    for (site in sites) {
      si <- paste("Site", site, sep="")
      cat("Mapping", si, "\n")
      vec<-subset(counts, select = si)
      dd <- data.frame(node = mapping$Branches, count = vec)
      names(dd) <- c("node", "count")
      p <- ggplot(tree) %<+% dd
      p <- p + geom_tree(aes(col = count))
      p <- p + scale_color_gradient2(low = "blue", mid = "white", high = "red", limits = c(-mx, mx))
      p <- p + ggtitle(label = si)
      pp[[si]] <- p
      ####
      # Plot ancestral states:
      #if (show.states) {
      #  k <- which(names(mapping) == s) - 2 #2 first columns in counts are "Branches" and "Mean"
      #  ids <- row.names(attr(mapping, "sequences"))
      #  states <- attr(mapping, "sequences")[, k]
      #  nodelabels(text = toupper(states), node = as.numeric(ids), frame = "none")
      #}
    }
    return(pp)
  } else if (drawing == "heatmap") {
    if (show.states) {
      # Add ancestral states reconstructions:
      tree<-attr(mapping, "tree")
      k <- which(names(mapping) %in% s) - 2 #2 first columns in counts are "Branches" and "Mean"
      ids <- row.names(attr(mapping, "sequences"))
      # The "To" states are directly given by the node ids (top nodes)
      statesTo <- as.data.frame(toupper(attr(mapping, "sequences")[as.character(vec$Branches), k]), stringsAsFactors = FALSE)
      names(statesTo) <- s
      statesTo <- melt(statesTo, measure.vars = 1:length(sites), variable.name = "Sites", value.name = "ToState")
      # The "From" states are given by the id of the parent node (bottom nodes),
      # which we first need to retrieve:
      ids <- getIds(as.character(vec$Branches), attr(mapping, "index"))
      parent.ids <- numeric(length = length(ids))
      for (i in 1:length(parent.ids)) {
        parent.ids[i] <- tree$edge[which(tree$edge[,2] == ids[i]), 1]
      }
      statesFrom <- as.data.frame(toupper(attr(mapping, "sequences")[attr(mapping, "index")[parent.ids], k]), stringsAsFactors = FALSE)
      names(statesFrom) <- s
      statesFrom <- melt(statesFrom, measure.vars = 1:length(sites), variable.name = "Sites", value.name = "FromState")
    }
    vec2 <- melt(vec, id.vars = c("Branches", "BrLen"), variable.name = "Sites", value.name = "Count")
    if (show.states) {
      vec2 <- cbind(vec2, statesTo, statesFrom)
    }
    vec3 <- subset(vec2, as.character(Branches) %in% b)

    # Build plot:
    p <- ggplot(vec3, aes(x = Sites, y = Branches))
    p <- p + geom_tile(height = vec3$BrLen, aes(fill = Count))
    p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-mx, mx))
    p <- p + scale_x_discrete(position = "top")
    p <- p + scale_y_discrete(limits = rev(b))
    p <- p + panel_border(size = 1, col = "black")
    p <- p + theme(legend.position = "left",
                   axis.title.x = element_blank())

    if (!show.branch.labels) {
      p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    }

    if (show.states) {
      p <- p + geom_label(aes(label = paste(FromState, sprintf("\u27a1"), ToState, sep = "")))
    }
    return(p)
  } else if (drawing == "compensogram") {
    sumvec2 <- subset(sumvec, as.character(Branches) %in% b)
    sumvec2$Type <- "Observed"
    # Randomization:
    pb <- txtProgressBar(1, 100, style = 3)
    for (i in 1:100) {
      setTxtProgressBar(pb, i)
      sim <- vec
      for (j in s) sim[, j] <- sample(sim[, j], replace = FALSE)
      sumsim <- data.frame(Branches = sim$Branches,
                           Sum = apply(sim[, s], 1, sum),
                           SumAbs = apply(abs(sim[, s]), 1, sum),
                           stringsAsFactors = TRUE)
      sumsim <- subset(sumsim, as.character(Branches) %in% b)
      sumsim$Type <- "Simulation"
      sumvec2 <- rbind(sumvec2, sumsim)
    }
    p <- ggplot() + xlab("Abs(Sum) / Sum(Abs)")

    simsum <- ddply(subset(sumvec2, Type == "Simulation"),
                    .variables = "Branches", .fun = summarise,
                    Mean = mean(abs(Sum)/SumAbs),
                    Q95 = quantile(abs(Sum)/SumAbs, probs = 0.05))
    p <- p + geom_errorbarh(data = simsum, aes(x = Mean, xmin = Q95, y = Branches), xmax = 1, col = "blue")
    p <- p + geom_point(data = simsum, aes(x = Mean, y = Branches), col = "blue")
    p <- p + geom_point(data = subset(sumvec2, Type == "Observed"), aes(x = abs(Sum) / SumAbs, y = Branches), col = "red")
    p <- p + scale_x_continuous(position = "top")
    p <- p + scale_y_discrete(limits = rev(b))
    p <- p + panel_border(size = 1, col = "black")
    if (!show.branch.labels) {
      p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
    }
    return(p)
  }
}

