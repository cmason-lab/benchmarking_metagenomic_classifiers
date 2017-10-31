#!/usr/bin/env Rscript
library(plyr)
library(Hmisc)
library(RColorBrewer)

options(warn=1)

#rmall()


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


#samples <- c("BioPool_BioPool", "Mason_NARG1")
samples <- c("Huttenhower_LC1", "Huttenhower_LC2", "Huttenhower_LC3", "Huttenhower_LC4", "Huttenhower_LC5", "Huttenhower_LC6", "Huttenhower_LC7", "Huttenhower_LC8", "Huttenhower_HC1", "Huttenhower_HC2")
Mat <- list()
readEstimators <- c("B","Bl","C","Cs","D","K","Kf","L","N")

hc_lc_readcounts <- read.table("../../Truth_set_IDs/HC-LC-abundances.csv",header=TRUE,sep=",") #separate file with read counts for each species for each dataset - available from IMMSA website in species subfolder of truth_sets
hc_lc_readcounts$count <- as.numeric(hc_lc_readcounts$count)
L1_distances <- data.frame(L1=character(),sample=character(),tool=character(),stringsAsFactors=FALSE)

for (s in samples) {
	#s <- samples[[1]]
  T <- read.delim(paste0("../../Truth_set_IDs/",s,"_abund.txt"), header = TRUE) #abundance estimate in automated build is incorrect - see truth_set_relative_abundances.zip 
  T$Abundance <- 100*T$Abundance/sum(T$Abundance)
  
  toolfilenames <- c("_BlastMeganFiltered.txt", "_BlastMeganFilteredLiberal.txt", 
                     "_ClarkM1Default.txt", "_ClarkM4Spaced.txt","_DiamondMegan.txt", 
                     "_KrakenFiltered.txt","_Kraken.txt", "_LMAT.txt", "_Metaphlan.txt", 
                     "_NBC.txt", "_Gottcha.txt", "_Phylosift90pct.txt", 
                     "_Phylosift.txt", "_MetaFlow.txt","_bmfmetagotlmat-priority.txt","_diakrametagotlmat-priority.txt"
  )
  results <- lapply(toolfilenames, function (t) arb_result2matrix(t,s))  
  names(results) <- c("B", "Bl", "C", "Cs", "D", "Kf", "K", "L", "M", "N", "G", "Pf", "P","MF","BE","DE")
#  results$Cbigf <- Cbigf
#  results$Cbigu <- Cbigu
  nn <- c("B", "Bl", "C", "Cs", "D", "Kf", "K", "L", "M", "N", "G", "Pf", "P","MF","BE","DE")
  results <- results[nn]
  
  results <- lapply(nn, function (y){x <- results[[y]];   x$V3[is.na(x$V3)] <-0.000001; x$V6 <- rep(y, nrow(x)); x}  )
  
  Mat[[s]] <- do.call(rbind, results)
  Mat[[s]]$V3[is.na(Mat[[s]]$V3)] <- 0
  T2 <- hc_lc_readcounts[which(hc_lc_readcounts$dataset == s),]
  T2$counts <- 100*T2$counts/sum(T2$counts)
  #Abundance = relative abundance for the _abund.txt truth sets
  true_abundance <- mapply(function(species,tool) if (is.element(tool,readEstimators)){T2$counts[match(species, T2$taxid)]} else {T$Abundance[match(species, T$taxid_species)]} , Mat[[s]]$V1, Mat[[s]]$V6) 
  Mat[[s]]$Deviation <- Mat[[s]]$V3 - true_abundance
  
  Mat[[s]]$TrueAbundance <- true_abundance
  Mat[[s]]$sample <- rep(s, nrow(Mat[[s]]))
  
  for (tool in unique(Mat[[s]]$V6)){
  	new_row <- c(as.character(sum(unlist(sapply(Mat[[s]]$Deviation[which(Mat[[s]]$V6 == tool)], function(x) if (! (is.na(x))){abs(x)})))),as.character(s),as.character(tool))
  	#names(new_row) <- c("L1","sample","tool")
  	L1_distances[sapply(L1_distances, is.factor)] <- lapply(L1_distances[sapply(L1_distances, is.factor)], as.character) 
  	L1_distances <- rbind(L1_distances,sapply(new_row,as.character))
  	colnames(L1_distances) <- c("L1","sample","tool")
  }
}
L1_distances$L1 <- sapply(L1_distances$L1,as.numeric)

Mat <- do.call(rbind, Mat)
Mat <- rename(Mat, c("V6" = "tool"))
Mat$tool <- factor(Mat$tool, levels = unique(Mat$tool))


pal = brewer.pal(12,"Paired")
pal <- c(pal, "gray", "maroon", "darkviolet", "magenta", "chartreuse", "chartreuse3")
pal <- pal[c(1:4,13,7:8,14,15:16,10,17:18,12,5:6)] #with filtered CosmosID: 1:6
alln <- c( "BlastMegan_filtered", "BlastMegan_liberal", "CLARK", "CLARK-S", "DiamondMegan_filtered","Kraken_filtered", "Kraken", "LMAT",  "MetaPhlAn", "NBC", "GOTTCHA",   "PhyloSift_filtered", "PhyloSift","MetaFlow","BlastEnsemble","DiamondEnsemble")


#scatter plot of expected vs. true abundances for GOTTCHA
outfile <-  paste0("GOTTCHA.LC-HC.scatter.abundance.pdf")
maintitle <- "GOTTCHA"
pdf(outfile,width=4,height=4)
plot(Mat[which(Mat$tool == 'G'),]$TrueAbundance,Mat[which(Mat$tool == 'G'),]$Deviation, col=pal[[14]], axes=FALSE, ylab='Deviation',xlab='Expected Abundance',pch=10,main=maintitle,xlim=c(0,70)) 
axis(2, mgp=c(3, .5, 0),labels=TRUE)
xticks <- 10*c(0:7)
axis(1, at=xticks, mgp=c(3, .3, 0),labels=TRUE)
dev.off()

#scatter plot of expected vs. true abundances for 
outfile <-  paste0("Diamond.LC-HC.scatter.abundance.pdf")
maintitle <- "Diamond-MEGAN"
pdf(outfile,width=4,height=4)
plot(Mat[which(Mat$tool == 'D'),]$TrueAbundance,Mat[which(Mat$tool == 'D'),]$Deviation, col=pal[[7]], axes=FALSE, ylab='Deviation',xlab='Expected Abundance',pch=5,main=maintitle,xlim=c(0,70)) 
axis(2, mgp=c(3, .5, 0),labels=TRUE)
axis(1, at=xticks, mgp=c(3, .3, 0),labels=TRUE)
dev.off()


#main boxplots 
outfile <-  paste0("LC_HC", ".boxplot.abundance.mean.pdf")
maintitle <- paste0("Difference of Estimated to True Abundance\n", "LC & HC", " samples")

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

pdf(outfile,width=6,height=8)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(4,0, 0, 0))
boxplot(DevTransformed~tool, data=Mat,  medcol="white", medlwd=1.5, outpch=16, frame.plot=FALSE, axes=FALSE, ylim=c(-2,4))
box(bty="l")
yticks <- c(-2:2)
axis(side = 2, at=yticks, las=2, cex.axis=0.75)

abline(h = yticks, col="gray")

boxplot(DevTransformed~tool, data=Mat,  medcol="white", medlwd=1.5, outpch=16, frame.plot=FALSE, add=TRUE, main=maintitle, outline=TRUE, yaxt='n',xaxt='n', col=pal, border=pal, ylab='log-modulus deviation')
axis(1, labels=FALSE, at=seq_along(alln))
text(x=seq_along(alln),  y = par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),  srt = 45, labels=alln, xpd=TRUE, adj=1,cex=0.75)

LogModX <- c(-100:100)
LogModY <- sapply(LogModX, function(y) if (y<0){ -1*log10(1+abs(y))} else { log10(1+abs(y))})

# calculate position of inset
plotdim <- par("plt")
xleft    = plotdim[2] - (plotdim[2] - plotdim[1]) * 0.25
xright   = plotdim[2]  #
ybottom  = plotdim[4] - (plotdim[4] - plotdim[3]) * 0.25  #
ytop     = plotdim[4]  #

# set position for inset
par(
  fig = c(xleft, xright, ybottom, ytop)
  , mar=c(0,0,0,0)
  , new=TRUE
  )

# add inset
plot(LogModX,LogModY, col=2, axes=FALSE, ylab='log-modulus (x)',xlab='x',pch='.') # inset bottomright
axis(2, mgp=c(3, .5, 0),cex.lab=0.7,cex.axis=0.6,labels=TRUE)
axis(1, mgp=c(3, .3, 0),cex.lab=0.7,cex.axis=0.6,labels=TRUE)
dev.off()

write.table(Mat, "abundancetablealltools_sim.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(L1_distances, "l1_abundance_results.txt", sep="\t", row.names=FALSE, quote=FALSE)

#from Steve's violin plot:
library(ggplot2)
library(tidyr)
library(dplyr)
#library(readr)
library(RColorBrewer) 

bp <- boxplot(L1~tool, data=L1_distances, plot=FALSE)
bpmedian <- abs(bp$stats[3,])
names(bpmedian) <- bp$names
bpmedian <- sort(bpmedian)
L1_distances$tool <- factor(L1_distances$tool, levels = names(bpmedian))
pal <- pal[names(bpmedian)]
alln <- alln[names(bpmedian)]

df<-rbind(L1_distances)%>%
  group_by(tool) %>%
  #mutate(N=n()) %>%
  ungroup() %>%
  #filter(N>1) %>%
  ggplot(aes(x=reorder(tool,`L1`,median),y=`L1`)) + xlab("")+scale_x_discrete(labels=alln)+

  geom_jitter(colour='black',alpha=.25,size=2) + geom_violin(aes(fill=tool),colour='black') + geom_boxplot(width=0.1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none',panel.background = element_rect(fill = "white"))+
  scale_fill_manual(values=pal)+
  ggsave('l1_species.pdf',width = 7, height = 4, dpi=500)
