############################
#### Sorting the results ###
############################


args  <-  commandArgs(trailingOnly=TRUE)
results_location   <-  toString(args[1])
data_location <- toString(args[2])



library(plyr)
library(data.table)
library(ggplot2)

setwd(paste0(results_location))

load("all_MR_results_prot_a.rdata")
load("all_coloc_results_prot_a.rdata")



##############################
##### Remove "selfies" #######
##############################

############################
## Coloc "sefies" removed ##
############################

Coloc_results_prot_a_table$lower.outcome_prot_a_number <- tolower(Coloc_results_prot_a_table$prot_a_number)

Coloc_results_selfies_removed <- subset(Coloc_results_prot_a_table, Coloc_results_prot_a_table$protein_number != Coloc_results_prot_a_table$lower.outcome_prot_a_number)

save(Coloc_results_selfies_removed, file="Coloc_results_selfies_removed.rdata")


##########################
## MR "selfies" removed ##
##########################

protein_linker_file <- read.table(paste0(data_location,"/prot-a-numbers_protein_IDs.txt"), header=T)

protein_linker_file$merged_IDs <- paste0(protein_linker_file$protein_ID, "_", protein_linker_file$protein_number) 

MR_results_prot_a_table$lower_outcome <- tolower(MR_results_prot_a_table$outcome)

MR_results_prot_a_table_split <- data.frame(do.call('rbind', strsplit(as.character(MR_results_prot_a_table$exposure), ':', fixed=TRUE )))

MR_results_prot_a_table$exposure_protein_ID <- MR_results_prot_a_table_split$X2

MR_results_prot_a_table$joint_exp_out_IDs <- paste0(MR_results_prot_a_table$exposure_protein_ID , "_",  MR_results_prot_a_table$lower_outcome)

MR_results_selfies_removed <- MR_results_prot_a_table[!MR_results_prot_a_table$joint_exp_out_IDs %in% protein_linker_file$merged_IDs,]

save(MR_results_selfies_removed, file="MR_results_selfies_removed.rdata")

## Keep just the SNPs (remove the IVW estimates) ##

MR_selfies_removed_SNP_only <- MR_results_selfies_removed[!grepl("All - Inverse variance weighted", MR_results_selfies_removed$SNP),]
save(MR_selfies_removed_SNP_only, file="MR_selfies_removed_SNP_only.rdata") 


## keep just IVW ##
MR_selfies_removed_IVW_only <- MR_results_selfies_removed[grepl("All - Inverse variance weighted", MR_results_selfies_removed$SNP),]
save(MR_selfies_removed_IVW_only, file="MR_selfies_removed_IVW_only.rdata")

###################################################################################
###### Extract the coloc results where H4 is the largest (the two coloclalise) ####
###################################################################################


subset_all_coloc_results <- Coloc_results_selfies_removed[,2:6]
subset_all_coloc_results$max <- colnames(subset_all_coloc_results)[apply(subset_all_coloc_results,1,which.max)]
Coloc_results_selfies_removed$max <- subset_all_coloc_results$max
H4_largest <- Coloc_results_selfies_removed[Coloc_results_selfies_removed$max=="PP.H4.abf",]
H3_largest <- Coloc_results_selfies_removed[Coloc_results_selfies_removed$max=="PP.H3.abf",]
H2_largest <- Coloc_results_selfies_removed[Coloc_results_selfies_removed$max=="PP.H2.abf",]
H1_largest <- Coloc_results_selfies_removed[Coloc_results_selfies_removed$max=="PP.H1.abf",]
H0_largest <- Coloc_results_selfies_removed[Coloc_results_selfies_removed$max=="PP.H0.abf",]

save(H4_largest, file="H4_largest_coloc_prot_a.rdata")
save(H3_largest, file="H3_largest_coloc_prot_a.rdata")
save(H2_largest, file="H2_largest_coloc_prot_a.rdata")
save(H1_largest, file="H1_largest_coloc_prot_a.rdata")
save(H0_largest, file="H0_largest_coloc_prot_a.rdata")

## Select H4 results that are over 0.8 

H4_greater_0.8 <- subset(Coloc_results_selfies_removed, Coloc_results_selfies_removed$PP.H4.abf > 0.8)
save(H4_greater_0.8, file="H4_greater_0.8.rdata")

################################################
#### Extract all the significant MR results ####
################################################

sig_pval <- 0.05/(length(unique(MR_results_selfies_removed$exposure))) ## adjusting for multiple testing, just adjust for the proteins not all the SNPs 

MR_results_selfies_removed$p_num <- as.numeric(as.character(MR_results_selfies_removed$p))

significant_MRs <- subset(MR_results_selfies_removed, MR_results_selfies_removed$p_num<sig_pval) ### just the significant MR results

significant_MRs_plus <- subset(MR_results_selfies_removed, MR_results_selfies_removed$exposure %in% significant_MRs$exposure) ### for any significant MR results also get the other instruments and IVW for the protein the signifcant results is instrumenting 

save(significant_MRs, file="significant_MR_results_only.rdata")

save(significant_MRs_plus, file="significant_MR_results_and_corresponding_SNPs.rdata")

###############################
### Finding the new pQTLs #####
###############################


SNP_only_sig_MR_results <- significant_MRs[!grepl("All - Inverse variance weighted", significant_MRs$SNP),]
save(SNP_only_sig_MR_results, file="SNP_only_sig_MR_results.rdata")

### removing associations we already know about from zheng ###

zheng_pQTLs <- read.csv(paste0(data_location,  "/pQTL_instruments_Zheng.csv"), header=T)

zheng_pQTLs$SNP_and_protein <- paste0(zheng_pQTLs$SNP, "_", zheng_pQTLs$Platform_id)

sig_MRs_with_outcome_prot_IDs <- merge(SNP_only_sig_MR_results, protein_linker_file, by.x="lower_outcome", by.y="protein_number")

sig_MRs_with_outcome_prot_IDs$SNP_and_outcome <- paste0(sig_MRs_with_outcome_prot_IDs$SNP, "_", sig_MRs_with_outcome_prot_IDs$protein_ID)

sig_MRs_SNP_only_zheng_assoc_removed <- subset(sig_MRs_with_outcome_prot_IDs, !sig_MRs_with_outcome_prot_IDs$SNP_and_outcome %in% zheng_pQTLs$SNP_and_protein)

save(sig_MRs_SNP_only_zheng_assoc_removed, file="sig_MRs_SNP_only_zheng_assoc_removed.rdata")

##########################################################
#### Extract the coloc results for the significant MRs ###
##########################################################

SNP_only_sig_MR_results$exp_out_rsid_IDs <- paste0(SNP_only_sig_MR_results$joint_exp_out_IDs, "_", SNP_only_sig_MR_results$SNP)
Coloc_results_selfies_removed$exp_out_rsid_IDs <- paste0(Coloc_results_selfies_removed$protein_name, "_", Coloc_results_selfies_removed$lower.outcome_prot_a_number, "_", Coloc_results_selfies_removed$rsid)
sig_MRs_with_colocs <- merge(SNP_only_sig_MR_results, Coloc_results_selfies_removed, by="exp_out_rsid_IDs")
save(sig_MRs_with_colocs, file="significant_MRs_with_colocs.rdata")


## Find how many colocs are over 0.8 for H4 for sig MRs ##

sig_MRs_H4_greater_0.8 <- subset(sig_MRs_with_colocs, sig_MRs_with_colocs$PP.H4.abf > 0.8)
save(sig_MRs_H4_greater_0.8, file="sig_MRs_H4_greater_0.8.rdata")

sig_MRs_H3_largest <- subset(sig_MRs_with_colocs, sig_MRs_with_colocs$max=="PP.H3.abf")   ## 5109 sig MRs with H3 largest
sig_MRs_H1_largest <- subset(sig_MRs_with_colocs, sig_MRs_with_colocs$max=="PP.H1.abf")   ## 7 MRs with H1 largest
save(sig_MRs_H3_largest, file="sig_MRs_H3_largest.rdata")
save(sig_MRs_H1_largest, file="sig_MRs_H1_largest.rdata")

sig_MRs_H4_lower_0.8 <- subset(sig_MRs_with_colocs, sig_MRs_with_colocs$PP.H4.abf <= 0.8)
sig_MRs_H4_lower_0.8_freq_table <- data.frame(table(sig_MRs_H4_lower_0.8$max))
sig_MRs_H4_lower_0.8_freq_table

# Var1       Freq
# PP.H1.abf     7
# PP.H3.abf  5109
# PP.H4.abf  1546

### Sig MRs gretaer than 0.8 for H1 and H3
sig_MRs_H3_greater_0.8 <- subset(sig_MRs_H4_lower_0.8, sig_MRs_H4_lower_0.8$PP.H3.abf > 0.8)
nrow(sig_MRs_H3_greater_0.8)
# 4152
sig_MRs_H1_greater_0.8 <- subset(sig_MRs_H4_lower_0.8, sig_MRs_H4_lower_0.8$PP.H1.abf > 0.8)
nrow(sig_MRs_H1_greater_0.8)
# 0


###########################################################################
### Ploting the significant MR associations SNP against outcome protein ###
###########################################################################

## plot the associations between the SNPs and the outcome proteins that were significant ##
## Plotting outcome protein position against SNP position and different colours for cis vs trans associations ###

## Add exposure protein positions

protein_positions <- read.table(paste0(data_location, "protein_positions.txt"), header=T)

protein_positions_unique <- unique(protein_positions)

SNP_only_sig_MR_results_split <- data.frame(do.call('rbind', strsplit(as.character(SNP_only_sig_MR_results$exposure), ':', fixed=TRUE )))

SNP_only_sig_MR_results$exposure_protein_name <- SNP_only_sig_MR_results_split$X1

sig_MRs_with_prot_positions <- merge(SNP_only_sig_MR_results, protein_positions_unique, by.x="exposure_protein_name", by.y="Exposure")

SNP_positions <- read.table(paste0(data_location,"/SNP_positions.txt"), header=T)

SNP_positions_unique <- unique(SNP_positions) 

sig_MRs_with_prot_and_SNP_positions <- merge(sig_MRs_with_prot_positions, SNP_positions_unique, by="SNP") 

colnames(sig_MRs_with_prot_and_SNP_positions)[25:27] <- c("CHR.protein_exposure", "POS.start.protein_exposure", "POS.end.protein_exposure")

### now need to add the outcome protein positions 

load(paste0(data_location,"/Sun_prot_positions.rdata"))   ## see "finding_protein_positions.R" for how got these positions


sig_MRs_with_prot_and_SNP_positions_exp_and_out <- merge(sig_MRs_with_prot_and_SNP_positions, Sun_prot_positions, by.x="lower_outcome", by.y="protein_number")

colnames(sig_MRs_with_prot_and_SNP_positions_exp_and_out)[c(36,37,48)] <- c("POS.start.protein_outcome", "POS.end.protein_outcome", "CHR.protein_outcome")


## Add if the SNP is cis or trans from the OUTCOME protein gene position ##
## Cis variant within 1Mb of the gene(s) encoding the protein/protein complex 
## trans variant >1Mb from the gene(s) encoding the protein/protein complex
## also add column of trans on the same chromosome and trans on different chromosomes

sig_MRs_with_prot_and_SNP_positions_exp_and_out$trans.different.chr <- ifelse(sig_MRs_with_prot_and_SNP_positions_exp_and_out$CHR.protein_outcome == sig_MRs_with_prot_and_SNP_positions_exp_and_out$CHR_SNP, "same_chr", "trans_different_chr")

table(sig_MRs_with_prot_and_SNP_positions_exp_and_out$trans.different.chr)
#same_chr     trans_different_chr 
#26531              388486  

sig_MRs_with_prot_and_SNP_positions_exp_and_out$start.outcome.protein.minus1Mb <- as.numeric(sig_MRs_with_prot_and_SNP_positions_exp_and_out$POS.start.protein_outcome - 1000000) ## start position minus 1Mb

sig_MRs_with_prot_and_SNP_positions_exp_and_out$start.outcome.protein.minus1Mb[sig_MRs_with_prot_and_SNP_positions_exp_and_out$start.outcome.protein.minus1Mb<0] <- 1 ### makes sure no minu base pair positions 

sig_MRs_with_prot_and_SNP_positions_exp_and_out$end.outcome.protein.plus1Mb <- as.numeric(sig_MRs_with_prot_and_SNP_positions_exp_and_out$POS.end.protein_outcome + 1000000)

#sig_MRs_with_prot_and_SNP_positions_exp_and_out$CHR.protein_outcome <- toupper(sig_MRs_with_prot_and_SNP_positions$CHR.protein_outcome) ## has lower and uppercase X and Y for chromosomes




sig_MRs_with_prot_and_SNP_positions_exp_and_out$cis.or.trans <- ifelse(
															
															sig_MRs_with_prot_and_SNP_positions_exp_and_out$CHR.protein_outcome == sig_MRs_with_prot_and_SNP_positions_exp_and_out$CHR_SNP &
															sig_MRs_with_prot_and_SNP_positions_exp_and_out$POS_SNP > sig_MRs_with_prot_and_SNP_positions_exp_and_out$start.outcome.protein.minus1Mb &
															sig_MRs_with_prot_and_SNP_positions_exp_and_out$POS_SNP < sig_MRs_with_prot_and_SNP_positions_exp_and_out$end.outcome.protein.plus1Mb
															,
															"cis",
															"trans" 
																)

table(sig_MRs_with_prot_and_SNP_positions_exp_and_out$cis.or.trans)

#   cis  trans 
#  2630 412387  


save(sig_MRs_with_prot_and_SNP_positions_exp_and_out, file="sig_MRs_with_prot_and_SNP_positions_exp_and_out.rdata")

########################################################################
####### Plot the dotplot for the significant MRs as cis and trans ######
########################################################################

df1 <- sig_MRs_with_prot_and_SNP_positions_exp_and_out

#df$CHR_protein_levels <- factor(df$CHR.protein_outcome, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"))
#df$CHR_SNP_levels <- factor(df$CHR_SNP, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"))

plot1 <- ggplot(data = df1, aes(x=POS.start.protein_outcome, y =POS_SNP, color = as.factor(cis.or.trans))) +
  geom_point() + 
  facet_grid(rows=vars(CHR_SNP), cols=vars(CHR.protein_outcome),space="free", scales = 'free') +
  scale_x_discrete("CHR.protein_outcome") +
  scale_y_reverse("CHR_SNP", labels=NULL, breaks=NULL)


pdf("dotplot_sig_MRs.pdf")
print(plot1)    
dev.off() 

################################################################################
####### Plot the dotplot for the significant MRs with H3 vs H4 as colours ######
################################################################################


length(grep("PP.H3.abf", sig_MRs_with_colocs$max))
# 5109
length(grep("PP.H4.abf", sig_MRs_with_colocs$max))
# 408558
nrow(sig_MRs_with_colocs)
#413674
length(grep("PP.H1.abf", sig_MRs_with_colocs$max))
# 7

## If there is a significant MR the H1 is effectively an H3 

sig_MRs_with_colocs_and_positions <- merge(sig_MRs_with_colocs, sig_MRs_with_prot_and_SNP_positions_exp_and_out, by="exp_out_rsid_IDs")
sig_MRs_with_colocs_and_positions$H3.or.H4 <- ifelse(sig_MRs_with_colocs_and_positions$max == "PP.H4.abf", "H4", "H3")
save(sig_MRs_with_colocs_and_positions, file="sig_MRs_with_colocs_and_positions.rdata")
  

df2 <- sig_MRs_with_colocs_and_positions

plot2 <- ggplot(data = df2, aes(x=POS.start.protein_outcome, y =POS_SNP, color = as.factor(H3.or.H4))) +
  geom_point() + 
  facet_grid(rows=vars(CHR_SNP), cols=vars(CHR.protein_outcome),space="free", scales = 'free') +
  scale_x_discrete("CHR.protein_outcome") +
  scale_y_reverse("CHR_SNP", labels=NULL, breaks=NULL)

pdf("dotplot_sig_MRs_h3_vs_H4.pdf")
print(plot2)    
dev.off() 


##############################################################
#### How many pQTLs colocalise with more than one protein ####
##############################################################


out_prot_freq_table <- data.frame(table(H4_greater_0.8$lower.outcome_prot_a_number))
exp_prot_freq_table <- data.frame(table(unlist(H4_greater_0.8$protein_number)))
out_prot_freq_table <- out_prot_freq_table[order(out_prot_freq_table$Freq),]
exp_prot_freq_table <- exp_prot_freq_table[order(exp_prot_freq_table$Freq),]

save(out_prot_freq_table, file="H4_0.8_out_prot_freq_table.rdata")
save(exp_prot_freq_table, file="H4_0.8_exp_prot_freq_table.rdata")


SNP_freq_table_H4_0.8 <- data.frame(table(unlist(H4_greater_0.8$rsid)))
SNP_freq_table_H4_0.8 <- SNP_freq_table_H4_0.8[order(SNP_freq_table_H4_0.8$Freq),]
save(SNP_freq_table_H4_0.8, file="H4_0.8_SNP_freq_table.rdata")

## Find out how many SNPs colocalises with only 1 SNP ##

number_SNPs_colocalised_with <- data.frame(table(SNP_freq_table_H4_0.8$Freq))
save(number_SNPs_colocalised_with, file="number_SNPs_colocalised_with.rdata")








