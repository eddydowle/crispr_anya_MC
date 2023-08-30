library(tidyverse)
library(vegan)
## RAREFACTION CURVES: determine if sufficient sequencing is achieved.

rarefaction_df_crispr <- read.csv("crispr read counts.txt",sep='\t',header=TRUE, row.names = 1)
rarefaction_df_meta <- read.csv("meta read counts.txt",sep='\t',header=TRUE, row.names = 1)

# identify lowest number of reads for samples

raremax_rarefaction_crispr <- min(rowSums(rarefaction_df_crispr))
raremax_rarefaction_meta <- min(rowSums(rarefaction_df_meta))
# generate graph
rarecurve(rarefaction_df_crispr, step = 100000, sample = raremax_rarefaction_crispr, col = 'blue', cex = 0.6)
rarecurve(rarefaction_df_meta, step = 100, sample = raremax_rarefaction_meta, col = 'blue', cex = 0.6)


#from vegan manual
S_cr <- specnumber(rarefaction_df_crispr) # observed number of species
Srare_cr <- rarefy(rarefaction_df_crispr, raremax_rarefaction_crispr)
plot(S_cr, Srare_cr, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(rarefaction_df_crispr, step = 20, sample = raremax_rarefaction_crispr, col = "blue", cex = 0.6)


S_meta <- specnumber(rarefaction_df_meta) # observed number of species
Srare_meta <- rarefy(rarefaction_df_meta, raremax_rarefaction_meta)
plot(S_meta, Srare_meta, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(rarefaction_df_meta, step = 20000, sample = raremax_rarefaction_meta, col = "blue", cex = 0.6)

