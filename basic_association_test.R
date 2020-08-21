library(snpR); library(qqman);
# prep data
sample.meta <- read.table("data/sample_meta.txt", header = T, stringsAsFactors = F)
genos <- read.table("genotypes.geno", stringsAsFactors = F)
snp.meta <- genos[,1:2]
colnames(snp.meta) <- c("scaffold", "position")
dat <- import.snpR.data(genos[,-c(1:2)], snp.meta, sample.meta)
dat <- filter_snps(dat, maf = 0.05, HWE = .0000005)

# run assocition test
dat <- calc_association(dat, response = "phen", maxiter = 1000)

# plot
p <- plot_manhattan(dat, "gmmat_pval_phen", chr = "scaffold", significant = 0.00001, suggestive = 0.0001, log.p = T)