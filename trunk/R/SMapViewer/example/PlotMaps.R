# Created on 21/07/14 by jdutheil

source("../SMapViewer.R")

counts<-read.mapping(file="Myo_counts.txt")
counts<-attach.tree(counts, "Myo_tags.dnd", "Myo_tags_tln.txt")
counts<-attach.sequences(counts, "Myo_ancestors.fasta", format="fasta")

plot(counts, sites=c(324, 332), edge.width = 10, show.tip.label=FALSE, show.node.label=FALSE, type="cladogram")

plot(counts, sites=c(307, 338, 279), edge.width = 10, show.tip.label=FALSE, show.node.label=FALSE, type="cladogram")
plot(counts, sites=c(254, 279, 292, 243), edge.width = 10, show.tip.label=FALSE, show.node.label=FALSE, type="cladogram")
plot(counts, sites=c(199, 267), edge.width = 10, show.tip.label=FALSE, show.node.label=FALSE, type="cladogram")
plot(counts, sites=c(326, 328), edge.width = 10, show.tip.label=FALSE, show.node.label=FALSE, type="cladogram")



