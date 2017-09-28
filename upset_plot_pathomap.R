library('UpSetR')
library('RColorBrewer')
require(grid)

setwd('species')
ds <- 'Downsample_Deep.50.Mreads' #'Huttenhower_HC1' #'UnambiguouslyMapped_ds.soil'
list_tools <- c("_BlastMeganFiltered.txt", "_BlastMeganFilteredLiberal.txt", 
                     "_ClarkM1Default.txt", "_ClarkM4Spaced.txt","_DiamondMegan.txt", 
                     "_KrakenFiltered.txt","_Kraken.txt", "_LMAT.txt", "_Metaphlan.txt", 
                     "_NBC.txt", "_Gottcha.txt", "_Phylosift90pct.txt", "_Phylosift.txt", "_MetaFlow.txt","_blastensemble-priority.txt","_diamondensemble-priority.txt")
#list_tools <- c("_ClarkM1Default.txt","_Kraken.txt","_LMAT.txt","_NBC.txt","_Phylosift.txt","_TRUTH.txt")
pal = brewer.pal(12,"Paired")
pal <- c(pal, "gray", "maroon", "darkviolet", "magenta", "chartreuse", "chartreuse3")
pal <- pal[c(1:4,13,7:8,14,15:16,10,17:18,12,5:6)]
#pal <- c(pal[c(1:4,13,7:8,14,15:16,9:10,17:18,12)],"black")
#pal <- c(pal[c(3,8,14,16,18)],"black")

alln <- c( "BlastMegan_filtered", "BlastMegan_liberal", "CLARK", "CLARK_S", "DiamondMegan_filtered","Kraken_filtered", "Kraken", "LMAT",  "MetaPhlAn", "NBC", "GOTTCHA",   "PhyloSift_filtered", "PhyloSift","MetaFlow","BlastEnsemble","DiamondEnsemble")
#alln <- c("CLARK","Kraken", "LMAT","NBC", "PhyloSift","Truth")

col.classes <- rep(list(c('numeric')),length(alln))

plot_upset <- function(name,setsize){

  species.df <- read.table(text = "",nrows=1,
                 colClasses = col.classes,
                 col.names = alln)
  rows <- list()
  set_lengths <- c()
  num.tools <- length(alln)
  for (i in 1:num.tools) {
  	tool <- alln[i]
  	finame = paste(name,list_tools[i],sep='')
  	results = read.table(finame,header=FALSE,sep='\t')
  	if (setsize > 0){
  		if ( ! tool %in% c("Truth")){ nset = min(setsize,length(results$V1)) } else{ nset = length(results$V1) }
  	} else{
  		nset = length(results$V1)
  	}
  	species = results$V1[1:nset]
  	print(paste(tool,':',length(species),'species'))
  	set_lengths <- c(set_lengths,length(species))
  	for (spe in species){
  		spe <- as.character(spe)
    	if (! spe %in% rows ){
    		#print(spe)
      		rows[[length(rows)+1]] <- spe
      		#print(paste(spe,'species'))
      		#print(integer(num.tools))
      		#print(nrow(species.df)+1)
      		species.df[nrow(species.df)+1,] <- integer(num.tools) 
    	}
    	row.names(species.df) <- rows
    	#print(row.names(species.df))
    	species.df[spe,tool] = as.numeric(1)
    	#print(species.df[nrow(species.df),])
    	#print(species.df[spe,tool])
  	}
  }	
  Species <- row.names(species.df)
  species.df <- cbind(Species,species.df)
  lengthind <- order(-set_lengths)
  print('plotting...')
    
  #pdf(paste0("../upset_",name,'_',setsize,".pdf"),width=18,height=5) #, width = 1800, height = 850, units='px', pointsize=10, res = 200)
  #upset(species.df, sets = alln, nsets=50, sets.bar.color = pal[lengthind],order.by=c("freq","degree"),point.size=1.5,decreasing=c(TRUE,TRUE), group.by = "degree",nintersects=200,mb.ratio = c(0.4, 0.6), scale.intersections='log10')

  #dev.off() 
  print('finished plotting')
  species.df$sum <- rowSums(species.df[,2:ncol(species.df)])
  for (overlap in 1:16){
  	print(paste(overlap,length(row.names(species.df[which(species.df$sum >= overlap),]))))
  }
  write.csv(row.names(species.df[which(species.df$sum >= 10),]),file="thirtyfour_species.txt",sep='\n',quote=F)
}

plot_upset(ds,0)
#plot_upset(ds,100)
#plot_upset(ds,25)
setwd('..')