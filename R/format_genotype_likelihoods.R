#' Format beagle genotype likelihood files into wide format, one individual per row.
#' 
#' Format beagle genotype likelihood files, like those produced by ANGSD, into
#' a format where each individual is on one row and columns are marker-genotype combinations.
#' 
#' @param genotypes data.table or an object coercable into a data.table. Beagle format genotype likelihoods.
#' 
#' @return A data.table containing genotype likelihoods where each individual sample is on one row
#'   and each column is a marker-genotype pair.
format_genotype_likelihoods <- function(genotypes){
  if(!"data.table" %in% class(genotypes)){
    genotypes <- data.table::as.data.table(genotypes)
  }
  
  # melt
  colnames(genotypes)[-c(1:3)] <- paste0(c("hoM", "het", "hom"), "_", colnames(genotypes[,-c(1:3)]))
  genotypes <- data.table::melt(genotypes, id.vars = c("marker", "allele1", "allele2"))
  
  # add genotype letters
  #a.tab <- data.frame(a = c("A", "C", "G", "T"), n = 0:3)
  #genotypes$allele1 <- a.tab$a[match(genotypes$allele1, a.tab$n)]
  #genotypes$allele2 <- a.tab$a[match(genotypes$allele2, a.tab$n)]
  genotypes$genotype <- substr(genotypes$variable, 1, 3)
  #genotypes$genotype[genotypes$genotype == "hoM"] <- paste0(genotypes$allele1[genotypes$genotype == "hoM"], genotypes$allele1[genotypes$genotype == "hoM"])
  #genotypes$genotype[genotypes$genotype == "hom"] <- paste0(genotypes$allele2[genotypes$genotype == "hom"], genotypes$allele2[genotypes$genotype == "hom"])
  #genotypes$genotype[genotypes$genotype == "het"] <- paste0(genotypes$allele1[genotypes$genotype == "het"], genotypes$allele2[genotypes$genotype == "het"])
  
  # pull out individual ids
  genotypes$ind <- gsub(".+_", "", genotypes$variable)
  
  # cast and return
  return(data.table::dcast(genotypes[,-c("allele1", "allele2", "variable")], 
                           formula = ind ~ marker + genotype, value.var = "value"))
}


get_expected_counts <- function(genotypes, method = "af"){
  genotypes[genotypes == 0.333333] <- NA
  counts <- t(t(genotypes[, -c(1:3)]) * c(0, 1, 2))
  counts <- as.data.table(counts)
  data.table::set(genotypes, j = 4:ncol(genotypes), value = counts)
  
  genotypes <- data.table::melt(genotypes, id.vars = c("marker", "allele1", "allele2"))
  genotypes <- data.table::dcast(genotypes[, -c("allele1", "allele2")], variable ~ marker, fun.aggregate = sum)
  
  inds <- genotypes[,1]
  genotypes <- genotypes[,-1]
  af <- colMeans(genotypes, na.rm = T) # this is the allele frequency of the "1" allele
  
  # identify all of the NAs and the columns that they belong to
  NAs <- which(is.na(genotypes)) # cells with NA
  NA.cols <- floor(NAs/nrow(genotypes)) + 1 # figure out the column
  adj <- which(NAs%%nrow(genotypes) == 0) # adjust for anything that sits in the last row, since it'll get assigned the wrong column
  NA.cols[adj] <- NA.cols[adj] - 1
  
  # do interpolation for each missing data point
  if(method == "bernoulli"){
    ndat <- rbinom(length(NAs), 2, af[NA.cols])
  }
  else if(method == "af"){
    ndat <- af[NA.cols]
  }
  
  # replace
  genotypes <- as.matrix(genotypes)
  genotypes[NAs] <- ndat
  
  return(data.table::as.data.table(genotypes))
}

filter_expected_counts_maf <- function(dat){
  af <- colSums(dat)/(2*nrow(dat))
  low.af <- which(af < .1)
  dat <- dat[,-..low.af]
  return(dat)
}