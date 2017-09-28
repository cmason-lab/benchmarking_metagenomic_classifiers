library('UpSetR')
library('RColorBrewer')

setwd('species')
ds <- 'UnambiguouslyMapped_ds.soil' #'Huttenhower_HC1' #'UnambiguouslyMapped_ds.soil'
#list_tools <- c("_BlastMeganFiltered.txt", "_BlastMeganFilteredLiberal.txt", 
#                     "_ClarkM1Default.txt", "_ClarkM4Spaced.txt","_DiamondMegan.txt", 
#                     "_KrakenFiltered.txt","_Kraken.txt", "_LMAT.txt", "_Metaphlan.txt", 
#                     "_NBC.txt", "_OneCodexAbundanceFiltered.txt", "_Gottcha.txt", "_Phylosift90pct.txt", "_Phylosift.txt", "_MetaFlow.txt","_TRUTH.txt"
#  )
list_tools <- c("_ClarkM1Default.txt","_Kraken.txt","_LMAT.txt","_NBC.txt","_Phylosift.txt","_TRUTH.txt")
pal = brewer.pal(12,"Paired")
pal <- c(pal, "gray", "maroon", "darkviolet", "magenta", "chartreuse", "chartreuse3")
#pal <- c(pal[c(1:4,13,7:8,14,15:16,9:10,17:18,12)],"black")
pal <- c(pal[c(3,8,14,16,18)],"black")

#alln <- c( "BlastMegan_filtered", "BlastMegan_liberal", "CLARK", "CLARK_S", "DiamondMegan_filtered","Kraken_filtered", "Kraken", "LMAT",  "MetaPhlAn", "NBC", "OneCodex_filtered",  "GOTTCHA",   "PhyloSift_filtered", "PhyloSift","MetaFlow","Truth")
alln <- c("CLARK","Kraken", "LMAT","NBC", "PhyloSift","Truth")

col.classes <- rep(list(c('numeric')),length(alln))

plot_upset <- function(name,setsize){

  species.df <- read.table(text = "",nrows=1,
                 colClasses = col.classes,
                 col.names = alln)
  rows <- list()
  set_lengths <- c()
  for (i in 1:length(alln)) {
  	tool <- alln[i]
  	finame = paste(name,list_tools[i],sep='')
  	results = read.table(finame,header=FALSE,sep='\t')
  	if ( ! tool %in% c("Truth")){ nset = min(setsize,length(results$V5)) } else{ nset = length(results$V5) }
  	species = results$V5[1:nset]
  	print(paste(tool,':',length(species),'species'))
  	set_lengths <- c(set_lengths,length(species))
  	for (spe in species){
  		#print(spe)
    	if (! spe %in% rows ){
      		rows[[length(rows)+1]] <- spe 
      		#print(paste(spe,'species'))
      		species.df[nrow(species.df)+1,] <- integer(length(alln))  
      		#print(nrow(species.df))
    	}
    	row.names(species.df) <- rows
    	#print(row.names(species.df))
    	species.df[spe,tool] = as.numeric(1)
    	#print(species.df[spe,tool])
  	}
  }	
  Species <- row.names(species.df)
  species.df <- cbind(Species,species.df)
  lengthind <- order(-set_lengths)
  print('plotting...')
  
  pdf(paste0("../upset_",name,'_',setsize,".pdf"),width=6,height=3.5,res=500) #, width = 1800, height = 850, units='px', pointsize=10, res = 200)
  upset(species.df, sets = alln, nsets=5, sets.bar.color = pal[lengthind],order.by=c("freq","degree"),point.size=1.5,decreasing=c(TRUE,TRUE), group.by = "degree",nintersects=200,mb.ratio = c(0.7, 0.3)) #, scale.intersections='log10')

  dev.off() 
  print('finished plotting')
  species.df$sum <- rowSums(species.df[,2:ncol(species.df)])
  
}

plot_upset(ds,50)
plot_upset(ds,100)
plot_upset(ds,25)
setwd('..')