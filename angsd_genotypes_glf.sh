#!/bin/bash -l
#SBATCH -t 4-12:00:00
#SBATCH --mem 20G
#SBATCH -n 12

list=bamlist
minInd=48

echo "$list"

angsd -bam ${list} -GL 1 -out genotypes -doMaf 2 -doMajorMinor 1 -SNP_pval 0.00000001 -doGlf 2 -doPost 2 -minQ 20 -minMapQ 20 -minInd $minInd -minMaf 0.05 -P 12 
