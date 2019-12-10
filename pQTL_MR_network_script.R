###########################
###### pQTL MR network ####
###########################


args  <-  commandArgs(trailingOnly = TRUE)
outcome_protein   <-  toString(args[1]) #Phenotypic exposure or outcome of interest


library(gwasvcf)
library(vcfR)
library(VariantAnnotation)
library(TwoSampleMR)

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
outcome_protein <- readVcf ("/mnt/storage/private/mrcieu/data/IGD/data/public/prot-a-1.vcf") 

dir.create("new_folder")


