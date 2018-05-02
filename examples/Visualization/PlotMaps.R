# Created on 21/07/14 by jdutheil
# Updated 21/01/16 by jdutheil

source("../../R/SMapViewer/SMapViewer.R")
require(cowplot)
require(svglite)

counts <- read.mapping(file="Myo_volume.txt")
counts <- attach.tree(counts, "Myo_tags.dnd", "Myo_tags_translation.txt")
counts <- attach.sequences(counts, "Myo_ancestors.fasta", format="fasta")

plot.candidate <- function(counts, sites, threshold, top, nrow = 1) {

  p.treelst <- plot(counts, sites, drawing = "tree", root = TRUE, show.states = FALSE, show.tip.label = FALSE)
  p.heatmap <- plot(counts, sites, drawing = "heatmap", top = top, threshold = threshold, show.branch.labels = FALSE)
  p.dotplot <- plot(counts, sites, drawing = "compensogram", top = top, threshold = threshold, show.branch.labels = FALSE)

  for (i in 1:length(p.treelst)) {
    p.treelst[[i]] <- p.treelst[[1]] + theme_tree(bgcolor = "grey30", plot.title = element_text(hjust = 0.5))
  }
  p.trees <- plot_grid(plotlist = p.treelst, nrow = 1)

  p.graphs <- plot_grid(
    p.heatmap,
    p.dotplot + theme(axis.title.y = element_blank()),
    labels = c("B", "C"), rel_widths = c(1,1), nrow = 1, align = "h")

  p <- plot_grid(p.trees, p.graphs, labels = c("A", ""), rel_widths = c(1.5,2), nrow = nrow)

  return(p)
}

# Visualizing significant groups (cf ../Simple/ProteinGroupCompensation/Myo_stats_pvalues.csv)
p <- plot.candidate(counts, sites=c(193, 260), threshold = 10, top = 20)
p

p <- plot.candidate(counts, sites=c(295, 313), threshold = 10, top = 20)
p

p <- plot.candidate(counts, sites=c(304, 341, 345), threshold = 10, top = 5, nrow = 2)
p
