
## seems to want to install these all seperately before I can install the vcf package??

#### Reload the bioconductor package ###
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")


#### Load gwasvcf
devtools::install_github("mrcieu/gwasvcftools")

library(gwasvcf)


### install gwasglue 

devtools::install_github("mrcieu/gwasglue")
library(gwasglue)


## read in the vcf file for protein 1
library(vcfR)
prot_a_1 <- read.vcfR ("/Users/hw15842/Documents/PhD_work/protein_interactions/prot-a-1.vcf")

head(prot_a_1_VCF)
head(prot_a_1@gt)
head(prot_a_1@meta)
head(prot_a_1@fix)

vcffile <- system.file("/Users/hw15842/Documents/PhD_work/protein_interactions/prot-a-1.vcf", package="gwasvcf")


library(VariantAnnotation)
prot_a_1 <- readVcf ("/Users/hw15842/Documents/PhD_work/protein_interactions/prot-a-1.vcf")



prot_a_1

header(prot_a_1)

samples(header(prot_a_1))

a <- vcf_to_granges(prot_a_1)
b <- vcf_to_granges(prot_a_1) %>% dplyr::as_tibble()

## convert to outcome data ##
prot_a_1_VCF <- gwasvcf_to_TwoSampleMR(prot_a_1, type="outcome")





### need to run MR for each of the pQTLs in zheng et al as instruments on all the VCF files from BC4 
## 2113 pQTLs
## 3282 proteins "a" on bc4 - think these are the ones from Sun
## 83 proteins "b" on bc4 - think these are the ones from Folkersen 
## https://gwas.mrcieu.ac.uk/datasets/



zheng_pQTLs <- read.csv("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/pQTL_instruments_Zheng.csv")



## read in zheng data to use as the exposures
## read in a prot-a-"number" to use as the outcome - just read it in one dataset per script, will run as multiple jobs 
## perform MR using all the pQTLs from zheng on the prot data read in from BC4 - should be 2113 results 
## need to make sure to run the pQTLs individually and as multiple instruments for a protein where there are multiple instruments
## store the results - new folder for each outcome protein 

library(gwasglue)
library(gwasvcf)
library(vcfR)
library(VariantAnnotation)
library(TwoSampleMR)
library(plyr)
library(data.table)

# Read in exposure data #
zheng_pQTLs <- read.csv("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/pQTL_instruments_Zheng.csv")

# Convert to exposure data for use in 2smr
zheng_pQTLS_formated <- format_data(zheng_pQTLs, type="exposure", phenotype_col = "Exposure", snp_col ="SNP", beta_col="Beta", se_col ="SE", 
                           eaf_col="Effect_allele_freq", effect_allele_col = "Effect_allele", other_allele_col = "Other_allele", pval_col = "P_value",
                           ncase_col = "Sample_size", chr_col = "CHR_SNP", pos_col = "POS_SNP")

# remove any rows where mr_keep.exposure=FALSE (means there is data missing) - think it is only FBLN1 that has data missing (doesnt have a value for "other_allele")
zheng_pQTLS_formated <- zheng_pQTLS_formated[!grepl("FALSE", zheng_pQTLS_formated$mr_keep.exposure), ]

# Split into list of individual proteins - this will keep the multiple instruments together for the time being
zheng_pQTLs_individual_proteins <- split(zheng_pQTLS_formated,zheng_pQTLS_formated$exposure, drop=TRUE)

# Read in outcome data 
prot_a_1 <- readVcf ("/Users/hw15842/Documents/PhD_work/protein_interactions/prot-a-1.vcf") 

# Convert outcome data for use in 2smr package 
prot_a_1_VCF <- gwasvcf_to_TwoSampleMR(prot_a_1, type="outcome")

#### Just read a subset of the vcf file 
vcfsubset <- query_gwas(prot_a_1, chrompos=c("1:1097291-1099437"))
vcfsubset2 <- gwasvcf_to_TwoSampleMR(vcfsubset, type="outcome")

## Need to set bcf tools
devtools::install_github('explodecomputer/genetics.binaRies')
set_plink()
set_bcftools()

## which snps do we want to extract - ones from zheng

snps_to_extract <- paste(zheng_pQTLS_formated$chr.exposure, zheng_pQTLS_formated$pos.exposure, sep=":")

prot_a_1_subset_extraction <- query_gwas("prot-a-1/prot-a-1.vcf.gz", chrompos=snps_to_extract)
prot_a_1_subset_extraction2 <- gwasvcf_to_TwoSampleMR(prot_a_1_subset_extraction, type="outcome")

## Run the MR for each exposure protein on the outcome protein


run_protein_MR <- function(exposure_protein){
  dat1 <- tryCatch(harmonise_data(exposure_protein, prot_a_1_subset_extraction2), error=function(e) NA)
  res1 <- tryCatch(mr_singlesnp(dat1, all_method = "mr_ivw"), error=function(e) NA)
  het <- tryCatch(mr_heterogeneity(dat1, method_list="mr_ivw"), error=function(e) NA)
  length <- tryCatch(nrow(res1), error=function(e) NA)
  a <- tryCatch(data.frame(matrix(NA, nrow=(length-1), ncol=3)), error=function(e) NA)
  tryCatch(colnames(a)<- c("Q", "Q_df", "Q_pval"), error=function(e) NA)
  b <- tryCatch(het[c(6:8)], error=function(e) NA)
  c <-tryCatch(rbind(a, b), error=function(e) NA)
  results <- tryCatch(cbind(res1, c), error=function(e) NA)
  return(results)
}


a <- zheng_pQTLs_individual_proteins[c(1:100)]
res_try <- lapply(a, run_protein_MR)
results_protein_interactions <- lapply(zheng_pQTLs_individual_proteins, run_protein_MR)
results_protein_interactions_table <- ldply(results_protein_interactions, data.table)




###########################################################################################################

## The above works for reading in "prot-a-1" need to now make it work to read all the "prot-a-"number"" in 

## See pQTL_MR_network_script.R 



#### Now need to do the colocalisation analysis ####

### need all the SNPs associated with the exposure and all snps that are associated with the outcome 
## got them for the outcome proteins - not sure where to get the exposure protein data from?
## can I just use the data on bc4 for the exposures as well?
## then need 500kb either side of the top SNP - which snp is top snp if more than one? 

library(coloc)


### need all the SNPs associated with the exposure and all snps that are associated with the outcome 
## got them for the outcome proteins - not sure where to get the exposure protein data from?
## can I just use the data on bc4 for the exposures as well?
## then need 500kb either side of the top SNP - which snp is top snp if more than one? 

library(coloc)


## Can use the inbuilt coloc function of gwasglue to run the coloc
## need to read in the subset 500kb either side of the SNP
## If a protein has multiple SNPs then do a coloc for each one

protein_linker_file <- read.table("prot-a_numbers_with_protein_names.txt")

## extract the vcf for each exposure (and each snp if more than one)

new_table <- merge(zheng_pQTLS_formated, protein_linker_file, by.x="exposure", by.y="V2")

outcome_protein <- paste("prot-a-2999/prot-a-2999.vcf.gz")

extract_from_VCF <- function(exposure_protein_number, exposure_protein_name, chr, rsid_position){
  bp_start <- (rsid_position - 500000)
  bp_finish <- (rsid_position + 500000)
  exposure_protein_VCF <- paste0(exposure_protein_number, "/", exposure_protein_number, ".vcf.gz")
  chr_position <- paste0(chr, ":", bp_start,"-", bp_finish)
  vcfsubset1 <- query_gwas(exposure_protein_VCF, chrompos=c(chr_position))
  vcfsubset2 <- query_gwas(outcome_protein, chrompos=c(chr_position))
  coloc_dataframes <- gwasvcf_to_coloc(vcfsubset1, vcfsubset2, chrompos = chr_position)
  coloc_results <- coloc.abf(coloc_dataframes[[1]], coloc_dataframes[[2]])
  coloc_results <- coloc_results[1]$summary
  coloc_results <- data.frame(coloc_results)
  coloc_results <- t(coloc_results)
  coloc_results <- data.frame(coloc_results)
  coloc_results$protein_name = as.character(exposure_protein_name)
  coloc_results$protein_number = as.character(exposure_protein_number)
  return(coloc_results)
}

prot_a_1_with_function <- extract_from_VCF(new_table$V1[1], new_table$exposure[1], new_table$chr.exposure[1], new_table$pos.exposure[1])

prot_a_1_with_function_formatted <- vcf_to_granges(prot_a_1_with_function) 

prot_mapply_try <- mapply(extract_from_VCF, new_table$V1, new_table$exposure, new_table$chr.exposure, new_table$pos.exposure)

new_table2 <- new_table[c(1,2),]
prot_mapply_try2 <- mapply(extract_from_VCF, new_table2$V1, new_table2$exposure, new_table2$chr.exposure, new_table2$pos.exposure)

prot_mapply_try2_table <- data.table(t(prot_mapply_try2))

df <- apply(prot_mapply_try2_table,2,as.character)

write.csv (df, file= "/Users/hw15842/Documents/PhD_work/protein_interactions/mapply_try.csv", quote=F, row.names=F, sep="\t")

## the protein name is in the file "prot-a-'number'_data.json
coloc_results <- gwasvcf_to_coloc(prot_mapply_try[[1]], prot_mapply_try[[2]], chrompos = as.character(3:185945324-186944705))




blah <- rbind.fill(lapply(results_protein_interactions, function(f) {
  as.data.frame(Filter(Negate(is.null), f))
}))
