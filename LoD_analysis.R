#!/usr/local/bin/Rscript
library(ggplot2)
library(reshape2)


# Read in the set of input read counts
counts <- read.table('HC-LC-input.csv', sep=',', header=TRUE, stringsAsFactors=FALSE)
# Calculate genome coverage
read_len <- 100  # Mean read length
counts$coverage <- counts$counts * read_len / counts$size

# Read in the set of results
result_folder <- '../output/species'
output_folder <- '.'
all_dat <- lapply(list.files(result_folder), function(f){
    fp <- paste(result_folder, f, sep='/')
    dataset <- strsplit(f, '_')[[1]][2]
    if(!dataset %in% counts$dataset) return(data.frame())
    method <- strsplit(f, '_')[[1]][3]
    method <- strsplit(method, '.txt')[[1]][1]
    # Check if the file is empty
    info = file.info(fp)
	cat(paste(method,dataset,info$size[1],'\n'))
    if(info$size[1] <= 1) return(data.frame())

    # Record which organisms were detected
    dat <- counts[counts$dataset == dataset, ]
    results <- read.table(fp, sep='\t', header=FALSE, stringsAsFactors=FALSE)
    dat$detected <- sapply(dat$species, function(species){species %in% results$V5})
    dat$method <- method
    dat$dataset <- dataset
    return(dat)
})
all_dat <- do.call(rbind, all_dat)

# Present a reasonable number of methods
# Remove the secondary entries for each method
for(m in c('ENSEMBLE1','MetaPallette-SmallDB-qual4-default',
           'MetaPallette-SmallDB-qual4-specific', 'OneCodexCounts','CosmosIDRelativeAbundance',
           'GeniusFiltered', 'CosmosIDFiltered', 'KrakenFiltered', 'TRUTH','diakrametagotlmat-max','diakrametagotlmat-average','diakrametagotlmat-median',
           'bmfmetagotlmat-average','bmfmetagotlmat-max','bmfmetagotlmat-median','f1ensemble-average','f1ensemble-max',
           'OneCodexCountsFiltered', 'OneCodexAbundance','OneCodexAbundanceFiltered','CosmosIDRelativeAbundanceFiltered',
           'PhyloSift90pct', 'BlastMeganFiltered','CosmosIDRelativeAbundanceFiltered','CosmosID')) all_dat <- all_dat[all_dat$method != m, ]

# Clean up names
all_dat$method[all_dat$method == 'ClarkM1Default'] <- 'CLARK'
all_dat$method[all_dat$method == 'ClarkM4Spaced'] <- 'CLARK-S'
#all_dat$method[all_dat$method == 'MetaPallette-SmallDB-qual4-default'] <- 'MetaPalette'
all_dat$method[all_dat$method == 'OneCodexAbundanceFiltered'] <- 'OneCodex_filtered'
all_dat$method[all_dat$method == 'BlastMeganFilteredLiberal'] <- 'BlastMegan_liberal'
all_dat$method[all_dat$method == 'DiamondMegan'] <- 'DiamondMegan_filt'
all_dat$method[all_dat$method == 'Metaphlan'] <- 'MetaPhlAn'
all_dat$method[all_dat$method == 'Gottcha'] <- 'GOTTCHA'
all_dat$method[all_dat$method == 'bmfmetagotlmat-priority'] <- 'BlastEnsemble'
all_dat$method[all_dat$method == 'diakrametagotlmat-priority'] <- 'DiamondEnsemble'
all_dat$method[all_dat$method == 'BlastMeganFilteredLMAT'] <- 'BlastMegan+LMAT'
all_dat$method[all_dat$method == 'DiamondMeganKrakenFiltered'] <- 'DiamondM.+Kraken'
all_dat$method[all_dat$method == 'ClarkM1DefaultGottcha'] <- 'CLARK+GOTTCHA'
#all_dat$method[all_dat$method == 'CosmosIDRelativeAbundance'] <- 'CosmosID'

# Bin the data by genome coverage
step <- 0.025
n_steps <- 6
bins <- c(1:n_steps) * step

sens <- do.call(rbind, lapply(unique(all_dat$method), function(method){
  do.call(rbind, lapply(bins, function(max_cov){
    sensitivity <- mean(all_dat[all_dat$method == method & all_dat$coverage > max_cov - step & all_dat$coverage <= max_cov, 'detected'])
    return(data.frame(method=method, bin=max_cov, sensitivity=sensitivity, stringsAsFactors=FALSE))
  }))
}))

print(do.call(rbind, lapply(bins, function(max_cov){
  method <- unique(all_dat$method)[1]
  n <- length(all_dat[all_dat$method == method & all_dat$coverage > max_cov - step & all_dat$coverage <= max_cov, 'detected'])
  return(data.frame(bin=max_cov, n_species=n, stringsAsFactors=FALSE))
})))

# Make a plot
#fp_out <- paste(output_folder, 'LoD_analysis.overlay.eps', sep='/')
#postscript(file=fp_out, width=6, height=5)
#p <- qplot(data=sens, x=bin, y=sensitivity, geom='line', color=method)+theme_bw()
#p <- p + ggtitle(paste('Limit of Detection'))
#p <- p + xlab('Genome Coverage (X)')
#p <- p + ylab('Recall')
#print(p)
#dev.off()

# Make a plot
fp_out <- paste(output_folder, 'LoD_analysis.pdf', sep='/')
pdf(file=fp_out, width=9, height=6)
p <- qplot(data=sens, x=bin, y=sensitivity, geom='line')
p <- p + facet_wrap(~method, ncol=6, scales='free_x')
p <- p + theme_bw()
p <- p + ggtitle(paste('Limit of Detection'))
p <- p + xlab('Genome Coverage (X)')
p <- p + ylab('Recall')
print(p)
dev.off()

# Make a table
wide_sens <- dcast(sens, bin ~ method)
names(wide_sens)[1] <- 'Genome Coverage (X)'
write.table(wide_sens, file=paste(output_folder, 'LoD.sensitivity.csv', sep='/'), sep=',', row.names=FALSE)
