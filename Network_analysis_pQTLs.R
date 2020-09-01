#####################################
#### Network analysis for pQTLs #####
####################################



setwd("/Users/hw15842/Documents/PhD_work/protein_interactions/new_results/")


load("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/Zheng_pQTLs_with_file_names.rdata")

load("unique_SNP_out_freq_table.rdata")   ### the SNPs that associate with one or more proteins 


library(ieugwasr)
ieugwasr::get_access_token()

### Get the genes for the SNPs that associate with 1 or more proteins ##

SNP_gene_query_1orMore <- variants_rsid(unique_SNP_out_freq_table$Var1)

SNP_gene_query_1orMore <- SNP_gene_query_1orMore[-c(grep("\\.", SNP_gene_query_1orMore$geneinfo)),]   ## remove ones that dont have gene info = 125 snps 

gene_freq_table_1orMore <- data.frame(table(SNP_gene_query_1orMore$geneinfo)) ## 341 genes for 490 SNPs


### Get the genes for all the pQTLs ###

SNP_gene_query_all <- variants_rsid(unique(Zheng_pQTLs_with_file_names$SNP))

SNP_gene_query_all <- SNP_gene_query_all[-c(grep("\\.", SNP_gene_query_all$geneinfo)),]   ## remove ones that dont have gene info = 236 snps 

gene_freq_table_all <- data.frame(table(SNP_gene_query_all$geneinfo)) ## 778 genes for 949 SNPs


### get genes for SNPs that associate with 1 or more but 50 or less proteins (1 to 50 inclusive) ###

rsids_1to50 <- subset(unique_SNP_out_freq_table, unique_SNP_out_freq_table$Freq <= 50)

SNP_gene_query_1to50 <- variants_rsid(rsids_1to50$Var1)

SNP_gene_query_1to50 <- SNP_gene_query_1to50[-c(grep("\\.", SNP_gene_query_1to50$geneinfo)),]   ## remove ones that dont have gene info = 122 snps 

gene_freq_table_1to50 <- data.frame(table(SNP_gene_query_1to50$geneinfo)) ## 330 genes for 439 SNPs






SNPs_not_in_gene_query <- subset(unique_SNP_out_freq_table_over50_assocs, !unique_SNP_out_freq_table_over50_assocs$SNP %in% SNP_gene_query$name) ### all on chr9... added in the SNPs that came out in the gene query but didnt have a gene 



