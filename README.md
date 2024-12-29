This is a R function to perform a Zscore meta genome wide association analysis. It was designed to work on the output of summary statistics from the GCTA software. It can also work on any summary statistics with the following column 
names (in any order): Chr, SNP, bp, b and p. The input of this function is the summary statistics from the various GWAS and their sample size (in the same order).

Example
#names of the various summary statistics with their path
meta_dat = c("C:/metafile/gwas1.mlma","C:/metafile/gwas2.mlma","C:/metafile/gwas3.mlma","C:/metafile/gwas4.mlma")
#sample size of each GWAS
n = c(2006,640,2911,2949)

source("zscore_meta.R")
#To run function
zscore_sumstat <- zscore_meta(dat = meta_dat,size = n)
