#!/usr/bin/env Rscript

# Author: Robert J. Prill <rjprill@us.ibm.com>
# Copyright IBM Corp 2016

rm(list=ls())

d = read.table("output/species_precision_recall.txt", header=F, sep="\t")
names(d) = c("sample", "tool", "Precision", "Recall", "F-score", "threshold")

x = subset(d, threshold == "best")$"F-score"
z = subset(d, threshold == "last")$"F-score"
tools = as.character(subset(d, threshold == "best")$tool)

# cbind(tools, x, z)  # sanity check

o = order(z-x, x, z)  # distance between best and last points
x = x[o]
z = z[o]
tools = tools[o]

y = seq_along(x)

show <- function() {
  par(mar=c(5, 4, 4, 15) + 1)
  plot(c(x, z), c(seq_along(x), seq_along(z)), type="n", main="ds.soil", ylab="", xlab="F-score (Presence / Absence)", xlim=c(0,1), yaxt="n")
  axis(4, at=seq_along(y), las=1, labels=tools)
  abline(v=1, lty=2)
  for(i in seq_along(x)) {
    a = c(x[i], i)
    a = rbind(a, c(z[i], i))
    lines(a)
  }
  points(x, y, pch=16, cex=1.5, col=2)
  points(z, y, pch=1,  cex=1.5, lwd=2, col=1)
  legend("topleft", cex=1, pch=c(16, 1), lty=0, lwd=c(1, 2), col=c(2, 1), legend=c("Optimal threshold", "No threshold"))
}

# show()

pdf("output/optimal_threshold_figure.pdf")
show()
dev.off()

