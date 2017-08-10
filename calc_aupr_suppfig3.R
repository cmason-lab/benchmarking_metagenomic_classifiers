#!/usr/local/bin/Rscript
library(ggplot2)

# Read in the txt file, accounting for an extra tab in the second column
read.txt <- function(fp){
    dat <- read.table(fp, sep='\t', stringsAsFactors=FALSE, row.names=NULL, header=TRUE)
    if(all(is.na(dat['algorithm']))){
        dat['algorithm'] <- dat['sample_tag']
        dat['sample_tag'] <- dat['row.names']
        dat['row.names'] <- NULL
    }
    return(dat)
}

calc_aupr <- function(dat){
    # Calculate AURP using the lower trapezoid method
    # http://pages.cs.wisc.edu/~boyd/aucpr_final.pdf
    recall_vals <- rev(sort(unique(dat[['recall']])))
    aupr <- 0
    for(ix in c(1:(length(recall_vals) - 1))){
        r1 <- recall_vals[ix]
        r2 <- recall_vals[ix + 1]
        p1 <- min(dat[dat['recall'] == r1, 'precision'])
        p2 <- max(dat[dat['recall'] == r2, 'precision'])
        aupr <- aupr + ((p1 + p2)/2 * (r1 - r2))
    }
    return(aupr)
}


for(level in c('genus', 'species', 'subspecies')){
    out <- data.frame(algorithm=c(), sample_tag=c(), aupr=c(), stringsAsFactors=FALSE)
    fp_in <- paste('output/', level, '_precision_recall.txt', sep='')
    # Read in CSV
    dat <- read.txt(fp_in)
    # Calculate AUPR for each combination of sample_tag and algorithm
    for(sample in unique(dat[['sample_tag']])){
        for(alg in unique(dat[['algorithm']])){
            # Get the subset of data for a given sample and algorithm
            t_dat <- dat[dat['sample_tag'] == sample & dat['algorithm'] == alg, ]
            if(nrow(t_dat) <= 1) next
            if(length(unique(t_dat[['recall']])) == 1) next
            aupr <- calc_aupr(t_dat)
            out <- rbind(out, data.frame(algorithm=c(alg),
                                         sample_tag=c(sample),
                                         aupr=c(aupr),
                                         stringsAsFactors=FALSE))
        }
    }

    # Save results
    fp_out <- paste('output/', level, '_AUPR.txt', sep='')
    write.table(out, file=fp_out, sep='\t', quote=FALSE, row.names=FALSE)

    # Calculate summary
    total_n_samples <- length(unique(out$sample_tag))
    summ <- do.call(rbind, lapply(unique(out$algorithm), function(alg){
        alg_out <- out[out$algorithm == alg, ]
        n_missing_samples <- total_n_samples - nrow(alg_out)
        if(n_missing_samples > 0){
            note<-paste("Warning: missing data for", n_missing_samples, "samples.")
        }else{
            note<-""
        }
        return(data.frame(algorithm=c(alg),
                          num_samples=nrow(alg_out),
                          median_aupr=c(median(alg_out$aupr)),
                          mean_aupr=c(mean(alg_out$aupr)),
                          note=note, stringsAsFactors=FALSE))
    }))
    fp_out <- paste('output/', level, '_AUPR.summary.txt', sep='')
    write.table(summ, file=fp_out, sep='\t', quote=FALSE, row.names=FALSE)

    # Make summary plot
    fp_out <- paste('output/', level, '_AUPR.pdf', sep='')
    pdf(file=fp_out)
    p <- ggplot(out, aes(x=factor(algorithm), y=aupr)) + geom_boxplot()
    p <- p + ggtitle(paste('Area Under Precision Recall Curve ', level))
    p <- p + xlab('Method')
    p <- p + ylab('AUPR') + ylim(0, 1)
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    print(p)

    # Make plots for each dataset
    for(sample in unique(out$sample_tag)){
        p <- ggplot(out[out$sample_tag == sample, ], aes(x=factor(algorithm), y=aupr)) + geom_bar(stat='identity')
        p <- p + ggtitle(paste(level, '-', sample))
        p <- p + xlab('Method')
        p <- p + ylab('AUPR') + ylim(0, 1)

        # Rotate axis labels
        p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        print(p)
    }
    dev.off()
    
}


###SUPP FIGURE 3

# Read in the AUPR results
all_precision_recall <- do.call(rbind, lapply(c('genus', 'species', 'subspecies'), function(level){
    fp <- paste(level, '_precision_recall.txt', sep='')
    fp <- paste('output', fp, sep='/')
    if(!file.exists(fp)){
        fp <- paste(level, '_precision_recall.txt', sep='')
        fp <- paste(results_folder, 'scores', fp, sep='/')
        if(!file.exists(fp)) return(data.frame())
    }
    dat <- read.table(file=fp, header=TRUE)
    if(names(dat)[1] != 'sample_tag'){
        dat <- read.table(file=fp, header=FALSE)
        names(dat) <- c('sample_tag', 'algorithm', 'precision', 'recall')
    }
    # Add the implicit termini
    x <- do.call(rbind, lapply(unique(dat$sample_tag), function(sample_tag){
        do.call(rbind, lapply(unique(dat$algorithm), function(algorithm){
            data.frame(sample_tag=sample_tag,
                       algorithm=algorithm,
                       precision=c(1, 0),
                       recall=c(0, 1), stringsAsFactors=FALSE)
        }))
    }))
    dat <- rbind(dat, x)
    dat$level <- level
    dat <- dat[order(dat$algorithm, dat$sample_tag, dat$precision, -dat$recall), ]
    return(dat)
}))

n_bins <- 10
bins <- c(0:n_bins) / n_bins
all_pr <- subset(all_precision_recall, level == 'species')
all_pr$algorithm <- sapply(all_pr$algorithm, as.character)
all_pr <-do.call(rbind, lapply(unique(all_pr$level), function(level){
    do.call(rbind, lapply(unique(all_pr$algorithm), function(alg){
        do.call(rbind, lapply(unique(all_pr$sample_tag), function(d){
            dat <- subset(all_pr, algorithm == alg & sample_tag == d)
            if(nrow(dat) == 0) return(data.frame())
            recall <- sapply(bins, function(precision){
                if(precision %in% dat$precision) return(max(dat$recall[dat$precision == precision]))
                lower_p <- max(dat$precision[dat$precision < precision])
                upper_p <- min(dat$precision[dat$precision > precision])
                lower_r <- max(dat$recall[dat$precision == lower_p])
                upper_r <- max(dat$recall[dat$precision == upper_p])
                y_span <- upper_r - lower_r
                x_span <- upper_p - lower_p
                delta_x <- precision - lower_p
                return(lower_r + (y_span * delta_x / x_span ))
            })
            return(data.frame(algorithm=alg, sample_tag=d, precision=bins,
                              recall=recall, level=level, stringsAsFactors=FALSE))
        }))
    }))
}))


# Put each method in a single subplot
fp <- paste('tool_off_analysis.combined.pdf')
pdf(file=fp, width=5, height=8)

# Put each method in a single subplot
p <- ggplot(subset(all_pr, level='species'), aes(x=precision, y=recall))
p <- p + geom_smooth(se=FALSE, n=5)
p <- p + facet_wrap(~algorithm)
p <- p + xlab('Precision')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + ylab('Recall')
print(p)

dev.off()

# Plot the individual ROC curves
fp <- paste('tool_off_analysis.individual.ROC.pdf', sep='/')
pdf(file=fp, height=10, width=13)

p <- qplot(data=subset(all_pr, level='species'),
               x=precision, y=recall, color=algorithm, group=algorithm,
               geom='line')
p <- p + facet_wrap(~sample_tag)
#p <- p + ggtitle(dataset)
p <- p + xlab('Precision')
p <- p + ylab('Recall')
print(p)

dev.off()

