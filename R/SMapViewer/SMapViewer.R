# Created on 10/07/14 by jdutheil
require(ape)    # for trees
require(seqinr) # for sequences

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
  attr(mapping, "tree")<-tree
  tln<-read.table(tlnfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  attr(mapping, "tiptln")<-tln
  index<-c(tree$tip.label, tree$node.label)
  attr(mapping, "index")<-index
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
  seqnames<-seqs$nam
  tln<-attr(mapping, "tiptln")
  for (i in 1:length(seqnames)) {
    j<-which(tln$Name == seqnames[i])
    if (length(j) > 0) seqnames[i]<-tln[j, "Id"]
    else if (length(j) > 1) stop(paste("Duplicated sequence", seqnames[i], "in translation table.\n"))
  }
  ids<-getIds(seqnames, attr(mapping, "index"))
  seqs$nam<-ids
  attr(mapping, "sequences")<-as.matrix(seqs)
  return(mapping)
  
}

plot.mapping<-function(mapping, sites, nrow=1, ...) {
  counts<-round(mapping[,-(1:2)])
  counts[counts > 1]<-1
  lim<-range(counts)
  f<-seq(lim[1], lim[2], length.out = 10)
  g<-rev(gray.colors(10))
  getCols<-function(values) {
    cols<-character(length(values))
    for (i in 1:length(values)) {
      cols[i]<-g[min(which(values[i] <= f))]
    }
    return(cols)
  }
  ####
  layout(matrix(1:length(sites), nrow=nrow))
  tree<-attr(mapping, "tree")
  # Get tree index:
  index<-attr(mapping, "index")
  # Get the corresponding vector:
  for (site in sites) {
    s<-paste("Site", site, sep="")
    cat("Mapping", s, "\n")
    vec<-subset(counts, select=s)
    rownames(vec)<-mapping$Branches
    vec<-vec[index[tree$edge[,2]],]
    plot(tree, edge.color=getCols(vec), ...)
  
    ####
    # Plot ancestral states:
    k <- which(names(counts) == s)
    ids <- row.names(attr(mapping, "sequences"))
    states <- attr(mapping, "sequences")[, k]
    nodelabels(text = toupper(states), node = as.numeric(ids), frame="none")
  }
}

