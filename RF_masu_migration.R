library(data.table)
fread(input = "genotypes.beagle") -> genotypes
source("R/format_genotype_likelihoods.R")
Phenos <- read.table("data/sample_meta.txt", header = TRUE)
get_expected_counts(genotypes = genotypes) -> dat 

library(ranger)
fdat<-filter_expected_counts_maf(dat = dat, .1)
fdat$Phenos<- Phenos$phen
fdat$Pop<- Phenos$pop
rm(dat, genotypes)
RF <- ranger(data = dat, dependent.variable.name = "Phenos", 
             num.trees = 100000, mtry = floor(ncol(fdat)/2), 
             num.threads = 24)
saveRDS(RF, "data/RFoutput.RDS")
