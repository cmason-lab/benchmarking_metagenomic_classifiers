library(MASS)
tools <- c('BlastMegan_filtered','DiamondMegan_filtered','Kraken_filtered','Kraken','PhyloSift',
           'BlastMegan_filtered_liberal','CLARK','MetaPhlAn','LMAT','NBC','PhyloSift_filtered','CLARK-S',
           'GOTTCHA','MetaFlow')[14:14]
for (tool in tools) {
  print(tool)
  medpar = read.csv(paste0(tool,'_nb.csv'))
  f = fps~n_reads+n_taxa+read_len+factor(sim)
  mod = glm.nb(f, medpar,control = glm.control(maxit = 100))
  outtab = coef(summary(mod))
  write.csv(outtab, file = paste0(tool,"_nbcoefs.csv"), quote = FALSE)
}
