#' @export
make_vcf <- function(x, meta, missing_genotypes = -9){
  #==================convert to 0/0, 0/1, 1/1, or ./.=============
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  tab <- data.frame(gt = c(0,1,2,missing_genotypes + missing_genotypes), recode = c("0/0", "0/1", "1/1", "./."))
  ind.genos <- tab[match(ind.genos, tab$gt),2]
  ind.genos <- matrix(ind.genos, ncol = ncol(x)/2)
  
  #==================add metadata to genos============
  # NOTE:
  # (vcf expects bp, may edit later to allow a meta containing that info to be imported)
  # for now, 0s will be A, 1s will be T, -1s will be .
  # for now, going to do this via conversion to 0 1 2, then converting. In the future will have a skip for phased data.
  
  vcf <- data.table::data.table(CHROM = as.numeric(as.factor(meta[,1])),
                                POS = meta[,2],
                                ID = paste0("snp", 1:nrow(meta)),
                                REF = "A",
                                ALT = "T",
                                QUAL = ".",
                                FILTER = "PASS",
                                INFO = ".",
                                FORMAT = "GT"
  )
  colnames(ind.genos) <- paste0("SAMP_", 1:ncol(ind.genos))
  vcf <- cbind(vcf, as.data.table(ind.genos))
  colnames(vcf)[1] <- '#CHROM'
  
  writeLines("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=Integer,Description='Genotype'>\n##FORMAT=<ID=GP,Number=G,Type=Float,Description='Genotype Probabilities'>\n##FORMAT=<ID=PL,Number=G,Type=Float,Description='Phred-scaled Genotype Likelihoods'>", "data.vcf")
  data.table::fwrite(vcf, "data.vcf", sep = "\t", append = T, col.names = T, row.names = F, scipen = 999)
  return(vcf)
}

impute_and_phase_beagle <- function(x = NULL, meta = NULL,
                                    beagle_path = "/usr/bin/beagle.jar",
                                    num_threads = 1,
                                    ne = 1000000,
                                    additional_args = NULL){
  #===============sanity checks===========
  # make a vcf if it doesn't exist
  if(!file.exists("data.vcf")){
    make_vcf(x, meta)
  }
  
  #===============construct call==========
  old.scipen <- getOption("scipen")
  options(scipen = 999)
  call <- paste0("java -jar ", beagle_path, " gt=data.vcf out=data.gt nthreads=",
                 num_threads, " ne=", ne)
  if(!is.null(additional_args)){
    call <- paste0(call, " ", additional_args)
  }
  
  #===============call beagle============
  system(call)
  options(scipen = old.scipen)
  
  
  #==============parse results===========
  # read in vcf and convert back to input format
  res <- readLines("data.gt.vcf.gz")
  res <- res[-which(grepl("^#", res))]
  res <- gsub("\\|", "\t", res)
  writeLines(res, "data.gt.two_col.txt")
  res <- data.table::fread("data.gt.two_col.txt", drop = 1:9)
  
  return(res)
}