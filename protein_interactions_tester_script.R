
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
## 3282 proteins "a" on bc4 
## 83 proteins "b" on bc4



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


## Run the MR for each exposure protein on the outcome protein


run_protein_MR <- function(exposure_protein){
  dat1 <- tryCatch(harmonise_data(exposure_protein, prot_a_1_VCF), error=function(e) NULL)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

run_protein_MR <- function(exposure_protein){
  dat1 <- harmonise_data(exposure_protein, prot_a_1_VCF)
  res1 <- mr(dat1, method_list = "mr_ivw")
  return(res1)
}

a <- zheng_pQTLs_individual_proteins[c(645:655)]
res_try <- lapply(a, run_protein_MR)
results_protein_interactions <- lapply(zheng_pQTLs_individual_proteins, run_protein_MR)


###########################################################################################################

## The above works for reading in "prot-a-1" need to now make it work to read all the "prot-a-"number"" in 

## See pQTL_MR_network_script.R 

