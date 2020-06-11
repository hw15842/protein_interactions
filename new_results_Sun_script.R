#################################
#### new_results_Sun_script.R ###
#################################


args  <-  commandArgs(trailingOnly=TRUE)
outcome_protein   <-  toString(args[1])
my_data_location <- toString(args[2])
BC4_data_location <- toString(args[3])
results_location <- toString(args[4])

library(gwasglue)
library(gwasvcf)
library(vcfR)
library(VariantAnnotation)
library(TwoSampleMR)
library(plyr)
library(data.table)
library(coloc)

setwd(paste0(results_location))

load(paste0(my_data_location,"Zheng_pQTLs_with_file_names.rdata"))  ### This is the file with the correct pQTL to protein data and also has the names of the full summary stat files for the proteins


head(Zheng_pQTLs_with_file_names)

zheng_pQTLs_FOLKERSON_only <- Zheng_pQTLs_with_file_names[grepl("prot-b", Zheng_pQTLs_with_file_names$prot_file_name),]

head(zheng_pQTLs_FOLKERSON_only)

zheng_pQTLS_formated <- format_data(zheng_pQTLs_FOLKERSON_only, type="exposure", phenotype_col = "Platform_id", snp_col ="SNP", beta_col="Beta", se_col ="SE", 
                                    eaf_col="Effect_allele_freq", effect_allele_col = "Effect_allele", other_allele_col = "Other_allele", pval_col = "P_value",
                                    ncase_col = "Sample_size", chr_col = "CHR_SNP", pos_col = "POS_SNP")



# remove any rows where mr_keep.exposure=FALSE (means there is data missing) - think it is only FBLN1 that has data missing (doesnt have a value for "other_allele")
zheng_pQTLS_formated <- zheng_pQTLS_formated[!grepl("FALSE", zheng_pQTLS_formated$mr_keep.exposure), ]

# Split into list of individual proteins - this will keep the multiple instruments together for the time being
zheng_pQTLs_individual_proteins <- split(zheng_pQTLS_formated,zheng_pQTLS_formated$exposure, drop=TRUE)

## Need to set bcf tools


#devtools::install_github('explodecomputer/genetics.binaRies')
library(genetics.binaRies)
set_plink()
set_bcftools()


## which snps do we want to extract - ones from zheng

snps_to_extract <- paste(zheng_pQTLS_formated$chr.exposure, zheng_pQTLS_formated$pos.exposure, sep=":")

outcome_protein <- unlist(regmatches(outcome_protein, gregexpr('\\(?[0-9,.]+', outcome_protein))) ### This stops it pasting the full "--outcome_protein=1" and just keeps the "1" (or what ever number the array job is)

pastename1_outcome_protein_location <- paste0(BC4_data_location, "prot-a-", outcome_protein, "/", "prot-a-", outcome_protein, ".vcf.gz") # this is becuase need to go into the protein file and then get the vcf file 
paste(pastename1_outcome_protein_location)
outcome_protein_subset <- query_gwas(pastename1_outcome_protein_location, chrompos=snps_to_extract)
outcome_protein_subset <- gwasvcf_to_TwoSampleMR(outcome_protein_subset, type="outcome")


## Run the MR - using singlesnp MR as gives the wald ratio for each instrument and the IVW estimate if there are multiple instruments, also have added the heterogeneity values on the end 
## tryCatch just gives a null value if an error arrises when the snp is not in the outcome data 

run_protein_MR <- function(exposure_protein){
  dat1 <- tryCatch(harmonise_data(exposure_protein, outcome_protein_subset), error=function(e) NA)
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


Sun_results_protein_interactions_MR <- lapply(zheng_pQTLs_individual_proteins, run_protein_MR)
Sun_results_protein_interactions_MR_table <- ldply(Sun_results_protein_interactions_MR, data.table)

### Save the results 

pastefile1 <- paste0("Sun_prot-a-", outcome_protein, "_MR_results_Folkerson_only.rdata")
save(Sun_results_protein_interactions_MR_table, file=pastefile1)



########################
#### COLOCALISATION ####
########################

## Can use the inbuilt coloc function of gwasglue to run the coloc
## need to read in the subset 500kb either side of the SNP
## If a protein has multiple SNPs then do a coloc for each one

### Just need to do Folkerson on Sun as have already done Sun on Sun


##### Coloc of pQTLs from FOLKERSON proteins on SUN Proteins (which are being iterated through via the .sh file)

Coloc_FOLKERSON_exposure <- function(Folkerson_protein_file_name, chr, rsid_position, SNP, Platform_id){
  bp_start <- (rsid_position - 500000)
  bp_finish <- (rsid_position + 500000)
  bp_start[bp_start<0] <- 1
  Folkerson_protein_VCF_location <- paste0(BC4_data_location, Folkerson_protein_file_name, "/", Folkerson_protein_file_name, ".vcf.gz")
  chr_position <- paste0(chr, ":", bp_start,"-", bp_finish)
  vcfsubset1 <- tryCatch(query_gwas(Folkerson_protein_VCF_location, chrompos=c(chr_position)), error=function(e) NA)  ## VCF subset 1 is the Folkerson proteins
  vcfsubset2 <- tryCatch(query_gwas(pastename1_outcome_protein_location, chrompos=c(chr_position)), error=function(e) NA) ### VCF subset 2 is the SUN protein that the .sh file is iterating through
  coloc_dataframes <- tryCatch(gwasvcf_to_coloc(vcfsubset1, vcfsubset2, chrompos = chr_position), error=function(e) NA)
  coloc_results <- tryCatch(coloc.abf(coloc_dataframes[[1]], coloc_dataframes[[2]]), error=function(e) NA)
  coloc_results <- tryCatch(coloc_results[1]$summary, error=function(e) NA)
  coloc_results <- tryCatch(data.frame(coloc_results), error=function(e) NA)
  coloc_results <- tryCatch(t(coloc_results), error=function(e) NA)
  coloc_results <- tryCatch(data.frame(coloc_results), error=function(e) NA)
  coloc_results$exposure_protein_number = as.character(Folkerson_protein_file_name)
  coloc_results$exposure_platform_id = as.character(Platform_id)
  coloc_results$outcome_protein_number = as.character(paste0("prot-a-", outcome_protein))
  coloc_results$rsid = as.character(SNP)
  return(coloc_results)
}

Folkerson_on_Sun_coloc_results <- mapply(Coloc_FOLKERSON_exposure, zheng_pQTLs_FOLKERSON_only$prot_file_name, zheng_pQTLs_FOLKERSON_only$CHR_SNP, zheng_pQTLs_FOLKERSON_only$POS_SNP, zheng_pQTLs_FOLKERSON_only$SNP, zheng_pQTLs_FOLKERSON_only$Platform_id)
Folkerson_on_Sun_coloc_results_table <- data.table(t(Folkerson_on_Sun_coloc_results))

pastefile2 <- paste0("Folkerson_on_Sun_prot-a-", outcome_protein, "_coloc_results.rdata")
save(Folkerson_on_Sun_coloc_results_table, file=pastefile2)



