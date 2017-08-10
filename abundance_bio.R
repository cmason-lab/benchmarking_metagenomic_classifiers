#!/usr/bin/env Rscript
library(plyr)
library(Hmisc)
library(RColorBrewer)

options(warn=1)

abd <- "../output/species/"
arb_result2matrix <- function(filename, sample) {
  
  filename <- paste0(abd, sample, filename)
  
  A <- try(read.delim(filename, header=FALSE))
  if (inherits(A, 'try-error') || nrow(A) < 1){
    A <- data.frame() 
    
  } else {
    A$V3 <- A$V3/sum(A$V3)
    A$V3 <- 100*A$V3
  }
  
  return(A)
}


samples <- c("BioPool_BioPool")
#samples <- c("Huttenhower_LC1", "Huttenhower_LC2", "Huttenhower_LC3", "Huttenhower_LC4", "Huttenhower_LC5", "Huttenhower_LC6", "Huttenhower_LC7", "Huttenhower_LC8", "Huttenhower_HC1", "Huttenhower_HC2")
Mat <- list()
readEstimators <- c("B","Bl","C","Cs","D","K","Kf","L","N")

for (s in samples) {
  T <- read.delim(paste0("../../Truth_set_IDs/",s,"_abund.txt"), header = TRUE) # abundance estimate in automated build is incorrect
  T$RelativeAbundance <- T$Abundance/(T$GenomeLength*T$Ploidy)
  T$RelativeAbundance <- 100*T$RelativeAbundance/sum(T$RelativeAbundance)
  
  toolfilenames <- c("_BlastMeganFiltered.txt", "_BlastMeganFilteredLiberal.txt", 
                     "_ClarkM1Default.txt", "_ClarkM4Spaced.txt","_DiamondMegan.txt", 
                     "_KrakenFiltered.txt","_Kraken.txt", "_LMAT.txt", "_Metaphlan.txt", 
                     "_NBC.txt", "_Gottcha.txt", "_Phylosift90pct.txt", "_Phylosift.txt", "_MetaFlow.txt","_bmfmetagotlmat-priority.txt","_diakrametagotlmat-priority.txt"
  )
  results <- lapply(toolfilenames, function (t) arb_result2matrix(t,s))  
  names(results) <- c("B", "Bl", "C", "Cs",  "D", "Kf", "K", "L", "M", "N", "G", "Pf", "P","MF","BE","DE")

  nn <- c("B", "Bl", "C", "Cs",  "D", "Kf", "K", "L", "M", "N", "G", "Pf", "P","MF","BE","DE")
  results <- results[nn]
  
  results <- lapply(nn, function (y){x <- results[[y]];   x$V3[is.na(x$V3)] <-0.000001; x$V6 <- rep(y, nrow(x)); x}  )
  
  Mat[[s]] <- do.call(rbind, results)
  Mat[[s]]$V3[is.na(Mat[[s]]$V3)] <- 0
  
  true_abundance <- mapply(function(species,tool) if (is.element(tool,readEstimators)){T$Abundance[match(species, T$taxid_species)]} else {T$RelativeAbundance[match(species, T$taxid_species)]} , Mat[[s]]$V1, Mat[[s]]$V6) 
  #Mat[[s]]$Deviation <- 100*(Mat[[s]]$V3 - true_abundance)/true_abundance
  Mat[[s]]$Deviation <- Mat[[s]]$V3 - true_abundance
  
  Mat[[s]]$TrueAbundance <- true_abundance #T$Abundance[match(Mat[[s]]$V1, T$taxid_species)]
  Mat[[s]]$sample <- rep(s, nrow(Mat[[s]]))
  
  
}

Mat <- do.call(rbind, Mat)
Mat <- rename(Mat, c("V6" = "tool"))
Mat$tool <- factor(Mat$tool, levels = unique(Mat$tool))



outfile <-  paste0("BioPool", ".boxplot.abundance.pdf")
maintitle <- paste0("BioOmics", " sample")


pal = brewer.pal(12,"Paired")
pal <- c(pal, "gray", "maroon", "darkviolet", "magenta", "chartreuse", "chartreuse3")
pal <- pal[c(1:4,13,7:8,14,15:16,10,17:18,12,5:6)]

alln <- c( "BlastMegan_filtered", "BlastMegan_liberal", "CLARK", "CLARK-S", "DiamondMegan_filtered","Kraken_filtered", "Kraken", "LMAT",  "MetaPhlAn", "NBC", "GOTTCHA",   "PhyloSift_filtered", "PhyloSift","MetaFlow","BlastEnsemble","DiamondEnsemble")


## code to sort boxplot by median
names(pal) <- unique(Mat$tool)
names(alln) <- unique(Mat$tool)
bp <- boxplot(Deviation~tool, data=Mat, plot=FALSE)

bpmedian <- abs(bp$stats[3,])
names(bpmedian) <- bp$names
bpmedian <- sort(bpmedian)
Mat$tool <- factor(Mat$tool, levels = names(bpmedian))
pal <- pal[names(bpmedian)]
alln <- alln[names(bpmedian)]
##

Mat$DevTransformed <- sapply(Mat$Deviation, function(y) if(!is.na(y)){ if (y<0){ -1*log10(1+abs(y))} else { log10(1+abs(y))}} else {NA})

pdf(outfile,width=6,height=6)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(4,0, 0, 0))
boxplot(DevTransformed~tool, data=Mat,  medcol="white", medlwd=1.5, outpch=16, frame.plot=FALSE, axes=FALSE, ylim=c(-2,2))
box(bty="l")
yticks <- c(-5:5)
axis(side = 2, at=yticks, las=2, cex.axis=0.75)

abline(h = yticks, col="gray")

boxplot(DevTransformed~tool, data=Mat,  medcol="white", medlwd=1.5, outpch=16,  frame.plot=FALSE, add=TRUE, main=maintitle, outline=TRUE, yaxt='n', xaxt='n', col=pal, border=pal,ylab='log-modulus deviation')
axis(1, labels=FALSE, at=seq_along(alln))
text(x=seq_along(alln),  y = par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),  srt = 45, labels=alln, xpd=TRUE, adj=1,cex=0.75)
dev.off()

write.table(Mat, "abundancetablealltools.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(T, "truth_abundances.txt", sep="\t", row.names=FALSE, quote=FALSE)







