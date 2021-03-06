args <- commandArgs(trailingOnly=TRUE)
require(biomaRt)

CluStrat_file <- args[1]
num_annots <- as.numeric(args[2])

ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

# getBM(attributes=c(
#     "refsnp_id", "chr_name", "chrom_start", "chrom_end",
#     "allele", "mapweight", "validated", "allele_1", "minor_allele",
#     "minor_allele_freq", "minor_allele_count", "clinical_significance",
#     "synonym_name", "ensembl_gene_stable_id"),
#     filters="snp_filter", values="rs6025",
#     mart=ensembl, uniqueRows=TRUE)

# data <- scan("CluStrat_sigSNPs.txt", what="character", sep = "\n")

data<- read.table(CluStrat_file,header=FALSE,sep=" ")

getBM(attributes=c("refsnp_id", "ensembl_gene_stable_id", "associated_gene"), filters="snp_filter", values=data[1:num_annots,2], mart=ensembl, uniqueRows=TRUE)
