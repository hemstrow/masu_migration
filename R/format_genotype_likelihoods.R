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
  a.tab <- data.frame(a = c("A", "C", "G", "T"), n = 0:3)
  genotypes$allele1 <- a.tab$a[match(genotypes$allele1, a.tab$n)]
  genotypes$allele2 <- a.tab$a[match(genotypes$allele2, a.tab$n)]
  genotypes$genotype <- substr(genotypes$variable, 1, 3)
  genotypes$genotype[genotypes$genotype == "hoM"] <- paste0(genotypes$allele1[genotypes$genotype == "hoM"], genotypes$allele1[genotypes$genotype == "hoM"])
  genotypes$genotype[genotypes$genotype == "hom"] <- paste0(genotypes$allele2[genotypes$genotype == "hom"], genotypes$allele2[genotypes$genotype == "hom"])
  genotypes$genotype[genotypes$genotype == "het"] <- paste0(genotypes$allele1[genotypes$genotype == "het"], genotypes$allele2[genotypes$genotype == "het"])
  
  # pull out individual ids
  genotypes$ind <- gsub(".+_", "", genotypes$variable)
  
  # cast and return
  return(data.table::dcast(genotypes[,-c("allele1", "allele2", "variable")], 
                           formula = ind ~ marker + genotype, value.var = "value"))
}