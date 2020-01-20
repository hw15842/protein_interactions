###########################
###### pQTL MR network ####
###########################


args  <-  commandArgs(trailingOnly=TRUE)
outcome_protein   <-  toString(args[1])
exposure_data_location <- toString(args[2])
outcome_data_location <- toString(args[3])
results_location <- toString(args[4])

library(gwasglue)
library(gwasvcf)
library(vcfR)
library(VariantAnnotation)
library(TwoSampleMR)
library(plyr)
library(data.table)
library(coloc)



# Read in exposure data #
zheng_pQTLs <- read.csv(paste0(exposure_data_location,"pQTL_instruments_Zheng.csv"))

head(zheng_pQTLs)

# Convert to exposure data for use in 2smr 
zheng_pQTLS_formated <- format_data(zheng_pQTLs, type="exposure", phenotype_col = "exposure_and_ID", snp_col ="SNP", beta_col="Beta", se_col ="SE", 
                                    eaf_col="Effect_allele_freq", effect_allele_col = "Effect_allele", other_allele_col = "Other_allele", pval_col = "P_value",
                                    ncase_col = "Sample_size", chr_col = "CHR_SNP", pos_col = "POS_SNP")


# remove any rows where mr_keep.exposure=FALSE (means there is data missing) - think it is only FBLN1 that has data missing (doesnt have a value for "other_allele")
zheng_pQTLS_formated <- zheng_pQTLS_formated[!grepl("FALSE", zheng_pQTLS_formated$mr_keep.exposure), ]

# Split into list of individual proteins - this will keep the multiple instruments together for the time being
zheng_pQTLs_individual_proteins <- split(zheng_pQTLS_formated,zheng_pQTLS_formated$exposure, drop=TRUE)


## Need to set bcf tools
devtools::install_github('explodecomputer/genetics.binaRies')
set_plink()
set_bcftools()


## which snps do we want to extract - ones from zheng

snps_to_extract <- paste(zheng_pQTLS_formated$chr.exposure, zheng_pQTLS_formated$pos.exposure, sep=":")

outcome_protein <- unlist(regmatches(outcome_protein, gregexpr('\\(?[0-9,.]+', outcome_protein))) ### This stops it pasting the full "--outcome_protein=1" and just keeps the "1" (or what ever number the array job is)

pastename1_outcome_protein_location <- paste0(outcome_data_location, "prot-a-", outcome_protein, "/", "prot-a-", outcome_protein, ".vcf.gz") # this is becuase need to go into the protein file and then get the vcf file 
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



results_protein_interactions <- lapply(zheng_pQTLs_individual_proteins, run_protein_MR)
results_protein_interactions_table <- ldply(results_protein_interactions, data.table)

### Save the results - saving list and table format, 

pastefile1 <- paste0(results_location, "prot-a-", outcome_protein, "_MR_results.txt")
write.table(results_protein_interactions_table, file=pastefile1, quote=F, row.names=F, sep="\t") 


########################
#### COLOCALISATION ####
########################

## Can use the inbuilt coloc function of gwasglue to run the coloc
## need to read in the subset 500kb either side of the SNP
## If a protein has multiple SNPs then do a coloc for each one



## read in the prot-a-'number' link protein name file
protein_linker_file <- read.table(paste0(exposure_data_location,"prot-a-numbers_protein_IDs.txt"), header=T)

head(protein_linker_file)

new_table <- merge(zheng_pQTLs, protein_linker_file, by.x="Platform_id", by.y="protein_ID")

head(new_table, n=20)
nrow(new_table)

## create function to extract 500kb either side of each pQTL instrument in Zheng from the BC4 full summary stats data 
## full summary stats data on BC4 is prot-a = Sun et al and prot-b = Folkersen et al
####### get the summary stats for the SNPs 500kb either side of the pQTL from Zheng data - this comes from the prot-a data on BC4 - thats why have the protein_linker_file, it uses the protein_IDs in zhengs data (platform_IDs)
####### then extract the same SNPs from the prot-a-1, prot-a-2, prot-a-3 etc files on BC4
####### then perform coloc using these two datasets


extract_from_VCF <- function(exposure_protein_number, exposure_protein_name, chr, rsid_position){
  bp_start <- tryCatch((rsid_position - 500000), error=function(e) NA)
  bp_finish <- tryCatch((rsid_position + 500000), error=function(e) NA)
  exposure_protein_VCF <- tryCatch(paste0(outcome_data_location, "/", exposure_protein_number, "/", exposure_protein_number, ".vcf.gz"), error=function(e) NA)
  chr_position <- tryCatch(paste0(chr, ":", bp_start,"-", bp_finish), error=function(e) NA)
  vcfsubset1 <- tryCatch(query_gwas(exposure_protein_VCF, chrompos=c(chr_position)), error=function(e) NA)
  vcfsubset2 <- tryCatch(query_gwas(pastename1_outcome_protein_location, chrompos=c(chr_position)), error=function(e) NA)
  coloc_dataframes <- tryCatch(gwasvcf_to_coloc(vcfsubset1, vcfsubset2, chrompos = chr_position), error=function(e) NA)
  coloc_results <- tryCatch(coloc.abf(coloc_dataframes[[1]], coloc_dataframes[[2]]), error=function(e) NA)
  coloc_results <- tryCatch(coloc_results[1]$summary, error=function(e) NA)
  coloc_results <- tryCatch(data.frame(coloc_results), error=function(e) NA)
  coloc_results <- tryCatch(t(coloc_results), error=function(e) NA)
  coloc_results <- tryCatch(data.frame(coloc_results), error=function(e) NA)
  coloc_results$protein_name = as.character(exposure_protein_name)
  coloc_results$protein_number = as.character(exposure_protein_number)
  return(coloc_results)
}

#a<- extract_from_VCF("prot-a-15","ACBD7",19, 54327869)
#a#

#b<- extract_from_VCF("prot-a-19", "ACP2", 14, 94844947)
#b#

#x<- tryCatch(extract_from_VCF("prot-a-18","ACP1",2,272203), error=function(e) NULL)
#x

## tester ## 

new_table2 <- new_table[c(1:100),]

results_coloc_tester <- mapply(extract_from_VCF, new_table2$protein_number, new_table2$Platform_id, new_table2$CHR_SNP, new_table2$POS_SNP)

results_coloc_tester_table <- data.table(t(results_coloc_tester))


results_coloc_tester_table <- apply(results_coloc_tester_table,2,as.character)

pastefile_test <- paste0(results_location, "/prot-a-", outcome_protein, "_coloc_results_tester.csv")
pastefile_test2 <- paste0(results_location, "/prot-a-", outcome_protein, "_coloc_results_tester.rdata")

save(results_coloc_tester_table, file=pastefile_test2)

write.csv (results_coloc_tester_table, file=pastefile_test, quote=F, row.names=F)


####


results_for_coloc_full <- mapply(extract_from_VCF, new_table$protein_number, new_table$Platform_id, new_table$CHR_SNP, new_table$POS_SNP)

results_for_coloc_full_table <- data.table(t(results_for_coloc_full))

results_for_coloc_full_table <- apply(results_for_coloc_full_table,2,as.character)

pastefile2 <- paste0(results_location, "/prot-a-", outcome_protein, "_coloc_results.csv")
write.csv (results_for_coloc_full_table, file=pastefile2, quote=F, row.names=F)


pastefile3 <- paste0(results_location, "/prot-a-", outcome_protein, "_coloc_results.csv")
save(results_for_coloc_full_table, file=pastefile3)



