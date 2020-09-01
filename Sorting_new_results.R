################################
#### Sorting the new results ###
################################



library(plyr)
library(data.table)
library(ggplot2)
library(dplyr)


setwd("/Users/hw15842/Documents/PhD_work/protein_interactions/new_results")


load("all_MR_results_prot_a.rdata")						    	      # Original MR results - all pQTLs on the Sun proteins - used wrong pQTLs file, need to remove the Folkerson SNPs and add the new ones 
load("all_coloc_results_prot_a.rdata")						        # Original coloc results of Sun pQTLs on Sun pQTLs - need to add the other pQTLs to it 

load("Folkerson_on_Folkerson_Coloc_results_table.rdata")	# Folkerson pQTLs on Folkerson proteins coloc analysis
load("Sun_on_Folkerson_Coloc_results_table.rdata")			  # Sun pQTLs on Folkerson proteins coloc analysis
load("Folkerson_MR_results_table.rdata")					        # All pQTLs on Folkerson proteins MR analysis 
load("SUN_COLOC_results_Folkerson_only_table.rdata")	  	# Folkerson only pQTLs on Sun proteins coloc analysis
load("SUN_MR_results_Folkerson_only_table.rdata")		    	# Folkerson only pQTLs on Sun proteins MR analysis

##########################################
##### Sort and Combine the MR results ####
##########################################

### remove the Folkerson pQTLs from the original MR results ###

wrong_Folkerson_IDs <- read.table("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/wrong_Folkerson_IDs.txt", header=T) 

MR_results_wrong_Folk_removed <- subset(MR_results_prot_a_table, !as.character(MR_results_prot_a_table$exposure) %in% as.character(wrong_Folkerson_IDs$wrong_Folkerson_IDs))

####### Add the new Folkerson pQTL MRs to the MR result file ###

# Add folkerson pQTLs only on Sun proteins

MR_results_wrong_Folk_removed[14:18] <- NULL
SUN_MR_results_Folkerson_only_table[14] <- NULL

MR_results_new_Folk_added <- rbind(MR_results_wrong_Folk_removed, SUN_MR_results_Folkerson_only_table)


# Add all pQTLs on FOlkerson proteins 

Folkerson_MR_results_table[14:18] <- NULL

MR_results_Sun_and_Folk <- rbind(MR_results_new_Folk_added, Folkerson_MR_results_table)

## remove MR selfies ##

load("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/Zheng_pQTLs_with_file_names.rdata")

MR_selfies_to_remove <- paste0(Zheng_pQTLs_with_file_names$Platform_id, "_", Zheng_pQTLs_with_file_names$prot_file_name)

MR_results_Sun_and_Folk_split_ID <- data.frame(do.call('rbind', strsplit(as.character(MR_results_Sun_and_Folk$exposure), ':', fixed=TRUE )))

MR_results_Sun_and_Folk$exposure_platform_id <- MR_results_Sun_and_Folk_split_ID$X2

MR_results_Sun_and_Folk$lower.outcome <- tolower(MR_results_Sun_and_Folk$outcome)

MR_results_Sun_and_Folk$exposure_on_outcome_IDs <- paste0(MR_results_Sun_and_Folk$exposure_platform_id, "_", MR_results_Sun_and_Folk$lower.outcome)

MR_results_Sun_and_Folk <- subset(MR_results_Sun_and_Folk, ! MR_results_Sun_and_Folk$exposure_on_outcome_IDs %in% MR_selfies_to_remove) 

save(MR_results_Sun_and_Folk, file="MR_results_Sun_and_Folk.rdata")

## Keep just the SNPs (remove the IVW estimates) ##

MR_results_SNP_only <- MR_results_Sun_and_Folk[!grepl("All - Inverse variance weighted", MR_results_Sun_and_Folk$SNP),]

## Have to remove results that have prot-b-12 and prot-b-24 as the outcome becuase something wierd going on in those dataset (figure this out lower down the script)
MR_results_SNP_only <- MR_results_SNP_only[!grepl("prot-b-12", MR_results_SNP_only$lower.outcome),]
MR_results_SNP_only <- MR_results_SNP_only[!grepl("prot-b-24", MR_results_SNP_only$lower.outcome),]
save(MR_results_SNP_only, file="MR_results_SNP_only.rdata") 

## keep just IVW ##
MR_results_IVW_only <- MR_results_Sun_and_Folk[grepl("All - Inverse variance weighted", MR_results_Sun_and_Folk$SNP),]
MR_results_IVW_only <- MR_results_IVW_only[!grepl("prot-b-12", MR_results_IVW_only$lower.outcome),]
MR_results_IVW_only <- MR_results_IVW_only[!grepl("prot-b-124", MR_results_IVW_only$lower.outcome),]

save(MR_results_IVW_only, file="MR_results_IVW_only.rdata") 



##########################################
####### Combine the Coloc results ########
##########################################

Coloc_results_prot_a_table$outcome_protein_number <- tolower(Coloc_results_prot_a_table$prot_a_number)
Coloc_results_prot_a_table$prot_a_number <- NULL
names(Coloc_results_prot_a_table)[7:8]<- c("exposure_platform_id", "exposure_protein_number")
Coloc_results_prot_a_table<- Coloc_results_prot_a_table[, c(1:6,8,7,10,9)]


names(Sun_on_Folkerson_Coloc_results_table)[7:8] <- c("exposure_protein_number", "exposure_platform_id")
Sun_on_Folkerson_Coloc_results_table$coloc_results <- NULL

SUN_COLOC_results_Folkerson_only_table <- data.frame(lapply(SUN_COLOC_results_Folkerson_only_table,unlist))


Coloc_results_Sun_and_Folk <- rbind(Coloc_results_prot_a_table, SUN_COLOC_results_Folkerson_only_table, Sun_on_Folkerson_Coloc_results_table)

## Remove Coloc selfies ##

Coloc_results_Sun_and_Folk <- subset(Coloc_results_Sun_and_Folk, !Coloc_results_Sun_and_Folk$exposure_protein_number == Coloc_results_Sun_and_Folk$outcome_protein_number)

row.names(Coloc_results_Sun_and_Folk) <- NULL

save(Coloc_results_Sun_and_Folk, file="Coloc_results_Sun_and_Folk.rdata")



#################################################################
####### Extract the new associations (significant pvals) ########
#################################################################

sig_pval <- 0.05/(length(unique(MR_results_SNP_only$outcome))) ## adjusting for multiple testing, just adjust for the proteins not all the SNPs  - COME BACK AND CHANGE TO OUTCOME! 

MR_results_SNP_only$p_val <- as.numeric(MR_results_SNP_only$p)

significant_MRs <- subset(MR_results_SNP_only, MR_results_SNP_only$p_val<sig_pval)

## Add columns for known/unknown results - i.e. associations already shown in zheng et al  ##

associations_already_known <- paste0(Zheng_pQTLs_with_file_names$prot_file_name, "_", Zheng_pQTLs_with_file_names$SNP)

significant_MRs$outcome_SNP_IDs <- paste0(significant_MRs$lower.outcome, "_", significant_MRs$SNP)

significant_MRs$known.unknown <- ifelse(significant_MRs$outcome_SNP_IDs %in% associations_already_known, "Known", "Unknown")

save(significant_MRs, file="significant_MRs.rdata")

significant_MRs_unique_SNP_out_associations <- significant_MRs[!duplicated(significant_MRs$outcome_SNP_IDs),]
nrow(significant_MRs_unique_SNP_out_associations)
# 20004
save(significant_MRs_unique_SNP_out_associations, file="significant_MRs_unique_SNP_out_associations.rdata")

###########################################################################
### Ploting the significant MR associations SNP against outcome protein ###
###########################################################################

## plot the associations between the SNPs and the outcome proteins that were significant ##
## Plotting outcome protein position against SNP position and different colours for cis vs trans associations ###

## Add the SNP positions

SNP_positions <- Zheng_pQTLs_with_file_names[12:14]
SNP_positions <- unique(SNP_positions)

significant_MRs_SNP_pos <- merge(significant_MRs, SNP_positions, by="SNP")
save(significant_MRs_SNP_pos, file="significant_MRs_SNP_pos.rdata")

### Add the outcome protein positions 

## Protein position data created in "Findin_protein_positions.R" ##

load("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/Outcome_protein_positions_SandF.rdata")

outcome_prot_position <- Outcome_protein_positions_SandF[,c(2:6,20)]
colnames(outcome_prot_position) <- paste("outcome", colnames(outcome_prot_position), sep="_")


significant_MRs_SNP_outProt_pos <- merge(significant_MRs_SNP_pos, outcome_prot_position, by.x="lower.outcome", by.y="outcome_protein_number")


## Add if the SNP is cis or trans from the OUTCOME protein gene position ##
## Cis variant within 1Mb of the gene(s) encoding the protein/protein complex 
## trans variant >1Mb from the gene(s) encoding the protein/protein complex

significant_MRs_SNP_outProt_pos$start.outcome.protein.minus1Mb <- as.numeric(as.numeric(significant_MRs_SNP_outProt_pos$outcome_start) - 1000000) ## start position minus 1Mb

significant_MRs_SNP_outProt_pos$start.outcome.protein.minus1Mb[significant_MRs_SNP_outProt_pos$start.outcome.protein.minus1Mb<0] <- 1 ### makes sure no minu base pair positions 

significant_MRs_SNP_outProt_pos$end.outcome.protein.plus1Mb <- as.numeric(as.numeric(significant_MRs_SNP_outProt_pos$outcome_end) + 1000000)

### Add in if the SNP is cis or trans from the outcome protein it is associated with ### 

significant_MRs_SNP_outProt_pos$cis.or.trans <- ifelse(
  
  significant_MRs_SNP_outProt_pos$outcome_CHR_protein_levels == significant_MRs_SNP_outProt_pos$CHR_SNP &
    significant_MRs_SNP_outProt_pos$POS_SNP > significant_MRs_SNP_outProt_pos$start.outcome.protein.minus1Mb &
    significant_MRs_SNP_outProt_pos$POS_SNP < significant_MRs_SNP_outProt_pos$end.outcome.protein.plus1Mb,
  "cis",
  "trans" 
)

table(significant_MRs_SNP_outProt_pos$cis.or.trans)

#   cis  trans 
#  2674 380828  


save(significant_MRs_SNP_outProt_pos, file="significant_MRs_SNP_outProt_pos.rdata")

significant_MRs_SNP_outProt_pos_unique <- significant_MRs_SNP_outProt_pos[!duplicated(significant_MRs_SNP_outProt_pos$outcome_SNP_IDs),]
nrow(significant_MRs_SNP_outProt_pos_unique)
# 19464
save(significant_MRs_SNP_outProt_pos_unique, file="significant_MRs_SNP_outProt_pos_unique.rdata")



################################# Dotplot of the new associations only ####################################
########## SNP position against outcome protein positions showing cis and trans associations #########

### Just the new associations and only the unique new SNP-outcome associations (there are SNP-expo-out associations that are different but the SNP-outcome associations are the same)

df1 <- subset(significant_MRs_SNP_outProt_pos_unique, significant_MRs_SNP_outProt_pos_unique$known.unknown == "Unknown")  ## we just want to plot the new associations in this one 

df1$CHR_SNP_levels <- factor(df1$CHR_SNP, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"))

plot1 <- ggplot(data = df1, aes(x=as.numeric(outcome_start), y =POS_SNP, color = as.factor(cis.or.trans))) +
  geom_point() + 
  facet_grid(rows=vars(CHR_SNP_levels), cols=vars(outcome_CHR_protein_levels),space="free", scales = 'free') +
  scale_x_discrete("Protein Position") +
  scale_y_reverse("SNP Position", labels=NULL, breaks=NULL) +
  scale_colour_discrete(name = "",
                        labels = c("Cis Associations", "Trans Associations")) 


pdf("dotplot_new_assoc_unique_only.pdf", width=10, height=10)
print(plot1)    
dev.off() 




###################################
#### Looking at the Coloc data ####
###################################

## Some of the pQTLs are not in the outcome protein files (prot-a and prot-b lightly different SNPs)

subset_all_coloc_results <- Coloc_results_Sun_and_Folk[,2:6]
subset_all_coloc_results$max <- colnames(subset_all_coloc_results)[max.col(subset_all_coloc_results,ties.method="first")]

num_coloc_highestPP <- table(subset_all_coloc_results$max)

# PP.H0.abf  PP.H1.abf  PP.H2.abf  PP.H3.abf  PP.H4.abf 
# 15241      4632418    70         117490     785287 

Coloc_results_Sun_and_Folk$max <- subset_all_coloc_results$max

H0_and_H2_colocs <- subset(Coloc_results_Sun_and_Folk, Coloc_results_Sun_and_Folk$max == "PP.H0.abf" |  Coloc_results_Sun_and_Folk$max == "PP.H2.abf" )
save(H0_and_H2_colocs, file="H0_and_H2_colocs.rdata")

exposure_proteins_H0andH2 <- subset(Zheng_pQTLs_with_file_names, Zheng_pQTLs_with_file_names$prot_file_name %in% H0_and_H2_colocs$exposure_protein_number & Zheng_pQTLs_with_file_names$SNP %in% H0_and_H2_colocs$rsid)
save(exposure_proteins_H0andH2, file="exposure_proteins_H0andH2.rdata")

exposure_proteins_H0andH2 <- exposure_proteins_H0andH2[c(12:14,24)]

H0_and_H2_colocs <- H0_and_H2_colocs[c(7,9,10)]

H0_and_H2_with_SNP_data <- merge(H0_and_H2_colocs, exposure_proteins_H0andH2, by.x=c("exposure_protein_number", "rsid"), by.y=c("prot_file_name", "SNP"))
save(H0_and_H2_with_SNP_data, file="H0_and_H2_with_SNP_data.rdata")

### Re run these exp prots on the out prots and double check that the reson they are H0 and H2 is becuase the pQTL not in the outcome protein data

#### All the H0 and H2s are either becuase the pQTL is not in the outcome protein or 
#### because the exposure protein is prot-b-12 or prot-b24 which are acting weirdly 
#### sowe can remove the H0 and H2 results, but also need to remove the outcome proteins that are protb-12 and -24 


Coloc_results_Sun_and_Folk_H0andH2removed <- subset(Coloc_results_Sun_and_Folk, Coloc_results_Sun_and_Folk$max == "PP.H1.abf" |  Coloc_results_Sun_and_Folk$max == "PP.H3.abf" |  Coloc_results_Sun_and_Folk$max == "PP.H4.abf")

freq_table_coloc <- table(Coloc_results_Sun_and_Folk_H0andH2removed$max)

#  PP.H1.abf  PP.H3.abf  PP.H4.abf 
#  4632418    117490     785287

save(Coloc_results_Sun_and_Folk_H0andH2removed, file="Coloc_results_Sun_and_Folk_H0andH2removed")

Coloc_results_for_analyses <- Coloc_results_Sun_and_Folk_H0andH2removed[!grepl("prot-b-12", Coloc_results_Sun_and_Folk_H0andH2removed$outcome_protein_number),]

Coloc_results_for_analyses <- Coloc_results_for_analyses[!grepl("prot-b-24", Coloc_results_for_analyses$outcome_protein_number),]

save(Coloc_results_for_analyses, file="Coloc_results_for_analyses.rdata") 
## this is with the H0 and H2 colocs removed becasue they didnt have the SNP in the outcome protein, 
## and the prot-b-12 and prot-b-24 proteins removed from the exposure and outcome becuase they were still giving H0 even when the SNP was present in both datasets
### The pvalues in the data set is different from what chris has so excluding prot-b-12 and prot-b-24 from all analysis

###############################################################
### How many significant MRs also have H4 greater than 0.8 ####
###############################################################


significant_MRs$expo_out_snp_ID <- paste0(significant_MRs$exposure_on_outcome_IDs, "_", significant_MRs$SNP)
Coloc_results_for_analyses$expo_out_snp_ID <- paste0(Coloc_results_for_analyses$exposure_protein_number,"_", Coloc_results_for_analyses$outcome_protein_number, "_",Coloc_results_for_analyses$rsid)


significant_MRs_with_coloc <- merge(significant_MRs, Coloc_results_for_analyses, by="expo_out_snp_ID")

table(significant_MRs_with_coloc$max)

#  PP.H1.abf   PP.H3.abf   PP.H4.abf 
#  4           4905        382855 
save(significant_MRs_with_coloc, file="significant_MRs_with_coloc.rdata")

significant_MRs_H4_greater0.8 <- subset(significant_MRs_with_coloc, significant_MRs_with_coloc$PP.H4.abf>0.8)
save(significant_MRs_H4_greater0.8, file="significant_MRs_H4_greater0.8.rdata")


################################################################################
####### Plot the dotplot for the significant MRs with H3 vs H4 as colours ######
################################################################################

## If there is a significant MR the H1 is effectively an H3 
significant_MRs_SNP_outProt_pos$expo_out_snp_ID <- paste0(significant_MRs_SNP_outProt_pos$exposure_on_outcome_IDs, "_", significant_MRs_SNP_outProt_pos$SNP)

significant_MRs_with_coloc_and_positions <- merge(significant_MRs_with_coloc, significant_MRs_SNP_outProt_pos, by="expo_out_snp_ID")
significant_MRs_with_coloc_and_positions$H3.or.H4 <- ifelse(significant_MRs_with_coloc_and_positions$max == "PP.H4.abf", "H4", "H3")

save(significant_MRs_with_coloc_and_positions, file="significant_MRs_with_coloc_and_positions.rdata")

## Just plot the unique and new associations between SNP and outcome 
significant_MRs_with_coloc_and_positions_unique <- significant_MRs_with_coloc_and_positions[!duplicated(significant_MRs_with_coloc_and_positions$outcome_SNP_IDs.x),]
save(significant_MRs_with_coloc_and_positions_unique, file="significant_MRs_with_coloc_and_positions_unique.rdata")

df2 <- subset(significant_MRs_with_coloc_and_positions_unique, significant_MRs_with_coloc_and_positions_unique$known.unknown.x == "Unknown")

df2$CHR_SNP_levels <- factor(df2$CHR_SNP, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"))


plot2 <- ggplot(data = df2, aes(x=as.numeric(outcome_start), y =POS_SNP, color = as.factor(H3.or.H4))) +
  geom_point() + 
  facet_grid(rows=vars(CHR_SNP_levels), cols=vars(outcome_CHR_protein_levels),space="free", scales = 'free') +
  scale_x_discrete("Protein Position") +
  scale_y_reverse("SNP Position", labels=NULL, breaks=NULL) +
  scale_colour_manual(name = "",
                      labels = c("Not Colocalised", "Colocalised"),
                      values = c("slateblue3", "springgreen3"))

pdf(file="dotplot_sig_MRs_H3_vs_H4_unique_and_new_assoc_only.pdf", width=10, height=10)
print(plot2)    
dev.off() 


####################################################################################################################################
####### Plot the dotplot for the significant MRs with H3 vs H4 as colours and diamond vs circle for new vs known associations ######
####################################################################################################################################

### Plot the unique SNP-outcome associations and distinguish between the known and unknown associations 

df3 <- significant_MRs_with_coloc_and_positions_unique

df3$CHR_SNP_levels <- factor(df3$CHR_SNP, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22", "X", "Y"))


plot3 <- ggplot(data = df3, aes(x=as.numeric(outcome_start), y =POS_SNP, color = as.factor(H3.or.H4), shape = as.factor(known.unknown.x))) +
  geom_point() + 
  facet_grid(rows=vars(CHR_SNP_levels), cols=vars(outcome_CHR_protein_levels),space="free", scales = 'free') +
  scale_x_discrete("Protein Position") +
  scale_y_reverse("SNP Position", labels=NULL, breaks=NULL) +
  scale_colour_discrete(name = "",
                        labels = c("Not Colocalised", "Colocalised")) +
  scale_shape_discrete(name = "",
                       labels = c("Known", "Unknown"))


pdf("dotplot_sig_MRs_unique_only_Coloc_known_and_unknown.pdf", width=12, height=8)
print(plot3)    
dev.off() 




################################################################
########## Reverse Causation -  reverse colocalisation #########
################################################################


H4_greater_0.8 <- subset(Coloc_results_for_analyses, Coloc_results_for_analyses$PP.H4.abf > 0.8)

H4_greater_0.8$out_on_expo_ID <- paste0(H4_greater_0.8$outcome_protein_number, "_", H4_greater_0.8$exposure_platform_id)

Coloc_results_for_analyses$expo_on_out_ID <- paste0(Coloc_results_for_analyses$exposure_platform_id, "_", Coloc_results_for_analyses$outcome_protein_number)

reverse_coloc_results <- subset(Coloc_results_for_analyses, Coloc_results_for_analyses$expo_on_out_ID %in% H4_greater_0.8$out_on_expo_ID)
save(reverse_coloc_results, file="reverse_coloc_results.rdata")
table(reverse_coloc_results$max)

# PP.H1.abf   PP.H3.abf   PP.H4.abf 
# 158587      2777    187728 

reverse_colocs_greater_H4 <- subset(reverse_coloc_results, reverse_coloc_results$PP.H4.abf> 0.8)

# 53% of the reverse colocs are H4 greater than 0.8

save(reverse_colocs_greater_H4, file="reverse_colocs_greater_H4.rdata")







#########################################################################################
########## Looking at SNPs assoc with num proteins vs SNP assoc with num traits #########
#########################################################################################


## Outcome traits that are in MRbase
outcome_traits_mrbase <- read.csv("/Users/hw15842/Documents/PhD_work/CCR5_paper/list_of_395_traits.csv", header=T)

## Remove NAs from mr.base.ids column in outcome_traits_mrbase ## - these NAs are the "other" traits that we need to get from elsewhere 

outcome_traits_mrbase <- outcome_traits_mrbase[!is.na(outcome_traits_mrbase$MR.Base.id),]

## extract the pQTLs from the MRbase traits ###
### Done on BC4 as times out otherwise ####
## See extract_SNPs_from_mrbase_traits.R  ##

## Load in the rdata file that was created from that with the SNPs extracted from the traits ## 

load("SNP_trait_results.rdata")


##### Calculate the r2 for the results to create the eigen vectors ####

### First Calculate F ###

SNP_trait_results$F_val <- ((SNP_trait_results$beta.outcome)/
                                     ((SNP_trait_results$se.outcome)))^2
 

## The the R2 ####
## This only holds when you have a linear regression with one exposure and a constant - fine for here but maybe not for other uses

SNP_trait_results$R_squ <- 1 / ((((SNP_trait_results$samplesize.outcome - 2)/SNP_trait_results$F_val)) +1)


df_to_turn_wide_matrix <- SNP_trait_results[c("SNP", "id.outcome", "R_squ")]

library(tidyr)

SNP_trait_wide_matrix <- pivot_wider(df_to_turn_wide_matrix, names_from = c("id.outcome"), values_from = c("R_squ"))

## this gives SNPs as rows, trait as columns and the r2 values for them 

### Now need to scale the data to have mean 0 and variance 1 ###
### We do this so that we are not giving a greater weighting to traits that are similar ###

SNP_trait_scaled <- SNP_trait_wide_matrix

SNP_trait_scaled[, -c(1)] <- scale(SNP_trait_scaled[, -c(1)], center = TRUE, scale = TRUE)

### Need to deal with the missingness in the dataset ##
## If the missingness is associated with pleiotropy then this will bias our results##
## However if the missingness is Random then should be fine ###
## I think the missingness is just due to which SNP-chip and imputation pannel is used ##
## BUT!!! need to check this ##

###  missMDA ##

library(missMDA)


imputed_results_trait <- imputePCA(SNP_trait_scaled[,-c(1)],method="Regularized",ncp=1)

imputed_results_trait <- data.frame(imputed_results_trait$completeObs)
PCA_results_trait <- prcomp(imputed_results_trait) ### Have to use prcomp because have to use prcomp for the proteins as there are more proteins than SNPs and princomp cant handle that
save(PCA_results_trait, file="PCA_results_trait.rdata")
trait_scores <- rowSums(PCA_results_trait$x)
trait_scores <- data.frame(trait_scores)
trait_scores$SNP <- SNP_trait_scaled$SNP
save(trait_scores, file="trait_scores.rdata")

### Need to do a regression of trait score on prot score with missingness as a covariant ###

### Need to do the same scaling for the protein scores ###
## Need to use the data extracted directly from the proteins - the beta ans se from the MR data set is for prot-prot not SNP-prot

load("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/SNP_only_extractionSNPs_from_FOLKERSON_table.rdata")
load("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/SNP_only_extractionSNPs_from_SUN_table.rdata")

SNP_prot_results <- rbind(SNPs_from_SUN_table, SNPs_from_FOLKERSON_table)

SNP_prot_results$F_val <- ((SNP_prot_results$beta.outcome)/
                              ((SNP_prot_results$se.outcome)))^2

SNP_prot_results$R_squ <- 1 / ((((SNP_prot_results$samplesize.outcome - 2)/SNP_prot_results$F_val)) +1)

SNP_prot_results$lower.outcome <- tolower(SNP_prot_results$outcome)

df_to_turn_wide_matrix_2 <- SNP_prot_results[c("SNP", "lower.outcome", "R_squ")]
SNP_prot_wide_matrix <- pivot_wider(df_to_turn_wide_matrix_2, names_from = c("lower.outcome"), values_from = c("R_squ"))

SNP_prot_scaled <- SNP_prot_wide_matrix
SNP_prot_scaled[, -c(1)] <- scale(SNP_prot_scaled[, -c(1)], center = TRUE, scale = TRUE)

imputed_results_prot <- imputePCA(SNP_prot_scaled[,-c(1)],method="Regularized",ncp=1)

imputed_results_prot <- data.frame(imputed_results_prot$completeObs)
PCA_results_prot <- prcomp(imputed_results_prot) #### used prcomp instead of princomp
save(PCA_results_prot, file="PCA_results_prot.rdata")
prot_scores <- rowSums(PCA_results_prot$x)

prot_scores <- data.frame(prot_scores)
prot_scores$SNP <- SNP_prot_scaled$SNP
save(prot_scores, file="prot_scores.rdata")


SNP_prot_scaled$na_count <- apply(SNP_prot_scaled, 1, function(x) sum(is.na(x)))
SNP_trait_scaled$na_count <- apply(SNP_trait_scaled, 1, function(x) sum(is.na(x)))

SNP_prot_NAs <- data.frame(SNP_prot_scaled[c("SNP", "na_count")])
SNP_trait_NAs <- data.frame(SNP_trait_scaled[c("SNP", "na_count")])
SNP_NAs <- merge(SNP_prot_NAs, SNP_trait_NAs, by="SNP")
names(SNP_NAs) <- c("SNP", "Prot_NAs", "Trait_NAs")

prot_trait_scores <- merge(prot_scores, trait_scores, by="SNP")

prot_trait_scores <- merge(prot_trait_scores, SNP_NAs, by="SNP")
save(prot_trait_scores, file="prot_trait_scores.rdata")

### Plotting the data ###

## ggplot function that will plot any linear regression ##
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


summary(lm(trait_scores ~ prot_scores + Prot_NAs + Trait_NAs, data=prot_trait_scores))

fit1 <- lm(trait_scores ~ prot_scores + Prot_NAs + Trait_NAs, data=prot_trait_scores)  ## missingness as covariate
regression_plot1 <- ggplotRegression(fit1)


pdf("linear_regression_scaled_imputed_trait_v_prot.pdf", width=12, height=8)
print(regression_plot1)    
dev.off() 

### see if MAF and LD score have an effect if added as covariates



fit2 <- 





##### Plot the raw data ####

head(SNP_prot_results)
head(SNP_trait_results)


SNP_trait_assoc_raw <- data.frame(table(SNP_trait_results$SNP))
SNP_prot_assoc_raw <- data.frame(table(SNP_prot_results$SNP))

SNP_assoc_raw <- merge(SNP_prot_assoc_raw, SNP_trait_assoc_raw, by="Var1")

names(SNP_assoc_raw)<- c("SNP", "Prot_assoc", "Trait_assoc")

SNP_assoc_raw <- merge(SNP_assoc_raw, SNP_NAs, by="SNP")

summary(lm(Prot_assoc ~ Trait_assoc + Prot_NAs + Trait_NAs, data=SNP_assoc_raw))

fit2 <- lm(Prot_assoc ~ Trait_assoc + Prot_NAs + Trait_NAs, data=SNP_assoc_raw)
ggplotRegression(fit2)

## Got 2 massive outliers, one in prot, a different one in trait ###
which.min(SNP_assoc_raw$Trait_assoc)
# 40
which.min(SNP_assoc_raw$Prot_assoc)
# 391

head(SNP_assoc_raw[order(SNP_assoc_raw$Prot_assoc),])
## And one SNP very low for both trait and prot row 384 

### plot again with these SNPs removed


SNP_assoc_raw_outlier_remove <- SNP_assoc_raw[-c(40,391,384),]

summary(lm(Prot_assoc ~ Trait_assoc + Prot_NAs + Trait_NAs, data=SNP_assoc_raw_outlier_remove))
 
fit3 <- lm(Prot_assoc ~ Trait_assoc + Prot_NAs + Trait_NAs, data=SNP_assoc_raw_outlier_remove)
ggplotRegression(fit3)

  






##################################
#### Gene enrichment analyses ####
##################################

unique_SNP_out_freq_table <- data.frame(table(significant_MRs_unique_SNP_out_associations$SNP))
unique_SNP_out_freq_table <- unique_SNP_out_freq_table[order(unique_SNP_out_freq_table$Freq),]
save(unique_SNP_out_freq_table, file="unique_SNP_out_freq_table.rdata")

unique_SNP_out_freq_table_over50_assocs <- subset(unique_SNP_out_freq_table, unique_SNP_out_freq_table$Freq > 50)
names(unique_SNP_out_freq_table_over50_assocs)[1] <- c("SNP")

SNP_positions <- Zheng_pQTLs_with_file_names[12:14]
SNP_positions <- unique(SNP_positions)

unique_SNP_out_freq_table_over50_assocs <- merge(unique_SNP_out_freq_table_over50_assocs, SNP_positions, by="SNP")
save(unique_SNP_out_freq_table_over50_assocs, file="unique_SNP_out_freq_table_over50_assocs.rdata")

library(ieugwasr)
ieugwasr::get_access_token()

SNP_gene_query <- variants_rsid(unique_SNP_out_freq_table_over50_assocs$SNP)
SNP_gene_query <- SNP_gene_query[-c(25,35,54),] ## remove ones that dont have gene info 

gene_freq_table <- data.frame(table(SNP_gene_query$geneinfo))

SNPs_not_in_gene_query <- subset(unique_SNP_out_freq_table_over50_assocs, !unique_SNP_out_freq_table_over50_assocs$SNP %in% SNP_gene_query$name) ### all on chr9... added in the SNPs that came out in the gene query but didnt have a gene 

snp_chrm_pos <- paste0(unique_SNP_out_freq_table_over50_assocs$CHR_SNP, ":", unique_SNP_out_freq_table_over50_assocs$POS_SNP)

SNP_gene_query_chrmpos <- variants_chrpos(snp_chrm_pos)

gene_IDs <- data.frame(do.call('rbind', strsplit(as.character(gene_freq_table$Var1), ':', fixed=TRUE )))

unique_SNP_out_freq_table_under50_assocs <- subset(unique_SNP_out_freq_table, unique_SNP_out_freq_table$Freq < 50)



## Are the genes in the NDEx query ##
## protein interactions ##
# SLC46A1 ABO and BCHE are not in the protein interactions network ##


## Add number of protein associations to the gene query table ##

SNP_gene_prot_assoc <- data.frame(SNP_gene_query[,2:5])

SNP_gene_prot_assoc <- merge(SNP_gene_prot_assoc, unique_SNP_out_freq_table_over50_assocs, by.x = "name", by.y = "SNP")
save(SNP_gene_prot_assoc, file="SNP_gene_prot_assoc.rdata")  ## the mass pleio SNPs and the genes they are associated with

genes_outside_prot_interactions <- grep("SLC46A1|ABO|BCHE", SNP_gene_prot_assoc$geneinfo)

SNP_gene_prot_assoc[genes_outside_prot_interactions,]


gene_interaction_network_1 <- grep("CFH|VTN|C7|SARM1|KLKB1|ZFPM2", SNP_gene_prot_assoc$geneinfo)
gin1 <- SNP_gene_prot_assoc[gene_interaction_network_1 ,]
gin1 <- gin1[order(gin1$geneinfo),]

gene_interaction_network_2 <- grep("APOE|BCHE|APOC1|SLC46A1", SNP_gene_prot_assoc$geneinfo)
gin2 <- SNP_gene_prot_assoc[gene_interaction_network_2 ,]
gin2 <- gin2[order(gin2$geneinfo),]


gene_interaction_network_3 <- grep("ABO|NLRP12", SNP_gene_prot_assoc$geneinfo)
gin3 <- SNP_gene_prot_assoc[gene_interaction_network_3 ,]
gin3 <- gin3[order(gin3$geneinfo),]



#### Add in the ABO SNPs that werent found in gene query but are ABO SNPs #####

SNPs_in_query <- SNP_gene_prot_assoc[c(1,5:7,4)]
names(SNPs_in_query)[1]<- c("SNP")
SNPs_manually_found_ABO <- subset(SNPs_not_in_gene_query, SNPs_not_in_gene_query$CHR_SNP == 9)
SNPs_manually_found_ABO$geneinfo <- c("ABO:28")

SNP_query_with_extra_ABO_snps <- rbind(SNPs_in_query, SNPs_manually_found_ABO)
save(SNP_query_with_extra_ABO_snps, file="SNP_query_with_extra_ABO_snps.rdata")


genes_outside_prot_interactions_ABO_snps_added <- grep("SLC46A1|ABO|BCHE", SNP_query_with_extra_ABO_snps$geneinfo)
SNP_query_with_extra_ABO_snps[genes_outside_prot_interactions_ABO_snps_added,]

gene_interaction_network_3_ABO_snps_added <- grep("ABO|NLRP12", SNP_query_with_extra_ABO_snps$geneinfo)
gin3_ABO_snps_added <- SNP_query_with_extra_ABO_snps[gene_interaction_network_3_ABO_snps_added ,]
gin3_ABO_snps_added <- gin3_ABO_snps_added[order(gin3_ABO_snps_added$geneinfo),]

### Number of gene associations from the gene association network NDEx ###

gene_assoc_network_assocs <- read.csv("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/num_assocs_gene_network.csv")

SNP_query_with_extra_ABO_snps_and_gene_assocs <- merge(SNP_query_with_extra_ABO_snps, gene_assoc_network_assocs, by.x = "geneinfo", by.y = "Gene")

fit4 <- lm(Assocs_gene_network ~ Freq, data=SNP_query_with_extra_ABO_snps_and_gene_assocs)  
regression_plot4 <- ggplotRegression(fit4)

## Remove NLRP12 and plot ###

SNP_query_with_extra_ABO_snps_and_gene_assocs_NLRP12_removed <- subset(SNP_query_with_extra_ABO_snps_and_gene_assocs, !SNP_query_with_extra_ABO_snps_and_gene_assocs$geneinfo=="NLRP12:91662")
fit5 <- lm(Assocs_gene_network ~ Freq, data=SNP_query_with_extra_ABO_snps_and_gene_assocs_NLRP12_removed)  
regression_plot5 <- ggplotRegression(fit5)


## Plot distribution of frequencies of SNPs ###

hist(unique_SNP_out_freq_table_over50_assocs$Freq)

unique_SNP_out_freq_table_under50_assocs <- subset(unique_SNP_out_freq_table, unique_SNP_out_freq_table$Freq < 50)
hist(unique_SNP_out_freq_table_under50_assocs$Freq)

hist(unique_SNP_out_freq_table$Freq)


over800_assocs <- subset(SNP_query_with_extra_ABO_snps, SNP_query_with_extra_ABO_snps$Freq > 800)

NLRP12_snps <- grep("NLRP12", SNP_query_with_extra_ABO_snps$geneinfo)
SNP_query_with_extra_ABO_snps[NLRP12_snps,]

genes_by_snp_freq_order <- SNP_query_with_extra_ABO_snps[order(SNP_query_with_extra_ABO_snps$Freq),]








#####################################################################
##################### SPECIFIC PLEIOTROPY ###########################
#####################################################################

sig_MR_with_coloc_less_than_50_assoc <- subset(significant_MRs_with_coloc_and_positions_unique, 
                                               significant_MRs_with_coloc_and_positions_unique$SNP.x %in% 
                                                 unique_SNP_out_freq_table_under50_assocs$Var1)

table(sig_MR_with_coloc_less_than_50_assoc$H3.or.H4)

# H3    H4 
# 1060  2302 
## So 31.5% do not colocalise and 68.4% do colocalise
## this is a smaller ratio than when looking at all the SNPs (H3=1617(9.3%) and H4=15648(90.6%))

sig_MR_with_coloc_more_than_50_assoc <- subset(significant_MRs_with_coloc_and_positions_unique, 
                                               significant_MRs_with_coloc_and_positions_unique$SNP.x %in% 
                                                 unique_SNP_out_freq_table_over50_assocs$SNP)


table(sig_MR_with_coloc_more_than_50_assoc$H3.or.H4)
# H3    H4 
# 517 13339
## much higher ratio of H4s in the mass pleiotropy SNPs





#############################################
########### COLIDER BIAS TESTING ############
#############################################

### Gib sent the eurpean 1000 genomes data via skype, I extracted the pQTls on bluecrystal. ####
thou_gen_mafs <- read.table("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/eur_1000_gen_pQTLs.txt", header=T, sep="\t")  ## this is build 37, gib sent link via skype http://fileserve.mrcieu.ac.uk/ld/data_maf0.01_rs_ref.tgz

Zheng_mafs <- Zheng_pQTLs_with_file_names[c("SNP", "Study", "CHR_SNP", "POS_SNP", "Effect_allele", "Other_allele", "Effect_allele_freq", "Sample_size")]


Zheng_thou_gen_mafs <- merge(Zheng_mafs, thou_gen_mafs, by.x="SNP", by.y="ID")  ### Keeps the duplicated Zheng SNPs from the different studies but not any duplicated thou gen SNPs

names(Zheng_thou_gen_mafs)[5:12] <- c("Effect_allele_Zheng", "Other_allele_Zheng", "Effect_allele_freq_Zheng", "Sample_size_Zheng", "Chrom_thou_gen", "Ref_allele_thou_gen", 
                                          "Alt_allele_thou_gen", "Alt_allele_freq_thou_gen")

## Need to make them all calculating the allele frequency for the effect allele as stated in effect allele for zheng data  ###

Zheng_thou_gen_mafs$Effect_allele_thou_gen <- ifelse(as.character(Zheng_thou_gen_mafs$Effect_allele_Zheng) == as.character(Zheng_thou_gen_mafs$Alt_allele_thou_gen), 
                                                     as.character(Zheng_thou_gen_mafs$Alt_allele_thou_gen),
                                                    as.character(Zheng_thou_gen_mafs$Ref_allele_thou_gen))

Zheng_thou_gen_mafs$Effect_allele_freq_thou_gen <- ifelse(as.character(Zheng_thou_gen_mafs$Effect_allele_Zheng) == as.character(Zheng_thou_gen_mafs$Alt_allele_thou_gen), 
                                                          Zheng_thou_gen_mafs$Alt_allele_freq_thou_gen,
                                                         (1 - Zheng_thou_gen_mafs$Alt_allele_freq_thou_gen))


Zheng_thou_gen_mafs$Sample_size_thou_gen <- 502  ## from plink2.log  on BC4

save(Zheng_thou_gen_mafs, file="Zheng_and_thou_gen_allele_frequencies.rdata")


# make sure that the effect allele is the same for the two populations
compare_af <- function(eaf_pop1, eaf_pop2, n_pop1, n_pop2)
{
  contingency <- matrix(
    c(2 * eaf_pop1 * n_pop1, 2 * (1-eaf_pop1) * n_pop1, 2 * eaf_pop2 * n_pop2, 2 * (1-eaf_pop2) * n_pop2), 2, 2)
  x <- fisher.test(contingency)
  x <- data.frame(t(unlist(x)))
  return(x)
}



allele_freq_comparison_1 <- mapply(compare_af, Zheng_thou_gen_mafs$Effect_allele_freq_Zheng, Zheng_thou_gen_mafs$Effect_allele_freq_thou_gen, Zheng_thou_gen_mafs$Sample_size_Zheng, Zheng_thou_gen_mafs$Sample_size_thou_gen)

allele_freq_comparison_2 <- data.table(t(allele_freq_comparison_1))
allele_freq_comparison_3 <- as.data.frame(matrix(unlist(allele_freq_comparison_2), nrow=nrow(allele_freq_comparison_2)))
names(allele_freq_comparison_3) <- names(allele_freq_comparison_2)

Zheng_thou_gen_maf_af_comparisons <- cbind(Zheng_thou_gen_mafs, allele_freq_comparison_3)


over_50_with_study <- merge(unique_SNP_out_freq_table_over50_assocs, unique(Zheng_pQTLs_with_file_names[c("SNP", "Study")]), by="SNP")
save(over_50_with_study, file="over_50_with_study.rdata")

Zheng_thou_gen_maf_af_comparisons$sig_pval <- ifelse(as.numeric(as.character(Zheng_thou_gen_maf_af_comparisons$p.value)) < (0.05/nrow(Zheng_thou_gen_maf_af_comparisons)),
                                                     "YES", "NO")

Zheng_thou_gen_maf_af_comparisons$Mass_pleio_SNP <- ifelse(Zheng_thou_gen_maf_af_comparisons$SNP %in% over_50_with_study$SNP,
                                                           "YES", "NO")

table(Zheng_thou_gen_maf_af_comparisons$Mass_pleio_SNP, Zheng_thou_gen_maf_af_comparisons$sig_pval)



#               Sig pval
# mass pleio    NO     YES
#         NO    1308   33
#         YES   757    0

##  None of the mass pleio ones are significant..... 

Zheng_thou_gen_maf_af_comparisons$num_assocs <- ifelse(Zheng_thou_gen_maf_af_comparisons$SNP %in% unique_SNP_out_freq_table$Var1,
                                                       unique_SNP_out_freq_table$Freq, "0")   
# Used the significant MRs as obviously those are the ones that are associated, but then we loose SNPs so need to add in 0 for snps that dont sig assoc with anything

save(Zheng_thou_gen_maf_af_comparisons, file="Zheng_thou_gen_maf_af_comparison.rdata")


### Make a regression with sig pval as a different shape and the mass pleio SNPs as a different colour

fit_MAFs_all_SNPs <- lm(Effect_allele_freq_thou_gen ~ Effect_allele_freq_Zheng, data=Zheng_thou_gen_maf_af_comparisons)
regression_plot_MAFs_all_SNPs <- ggplotRegression(fit_MAFs_all_SNPs) 
regression_plot_MAFs_all_SNPs 
                                  

plot_mafs <- ggplot(data = Zheng_thou_gen_maf_af_comparisons, aes(x=Effect_allele_freq_Zheng, y =Effect_allele_freq_thou_gen, 
                                                                  color = as.factor(Mass_pleio_SNP), shape = as.factor(sig_pval))) +
  geom_point() + 
  scale_colour_discrete(name = "", labels = c("Not mass pleio", "Mass pleio")) +
  scale_shape_discrete(name = "", labels = c("pvalue not sig", "pvalue sig"))

plot_mafs


## Make a plot of odds ratio (log odds?) against num prot associations per snp ###



plot_odds_vs_num_assocs <- ggplot(data = Zheng_thou_gen_maf_af_comparisons, aes(x=as.numeric(num_assocs), y =as.numeric(as.character(estimate.odds.ratio)), 
                                                                  color = as.factor(Mass_pleio_SNP), shape = as.factor(sig_pval))) +
  geom_point() + 
  scale_colour_discrete(name = "", labels = c("Not mass pleio", "Mass pleio")) +
  scale_shape_discrete(name = "", labels = c("pvalue not sig", "pvalue sig"))


plot_odds_vs_num_assocs



#### do with log odds as well

plot_log_odds_vs_num_assocs <- ggplot(data = Zheng_thou_gen_maf_af_comparisons, aes(x=as.numeric(num_assocs), y =abs(log(as.numeric(as.character(estimate.odds.ratio)))), 
                                                                                color = as.factor(Mass_pleio_SNP), shape = as.factor(sig_pval))) +
  geom_point() + 
  scale_colour_discrete(name = "", labels = c("Not mass pleio", "Mass pleio")) +
  scale_shape_discrete(name = "", labels = c("pvalue not sig", "pvalue sig"))


plot_log_odds_vs_num_assocs


sig_af_differences <- Zheng_thou_gen_maf_af_comparisons[as.numeric(as.character(Zheng_thou_gen_maf_af_comparisons$p.value)) < (0.05/nrow(Zheng_thou_gen_maf_af_comparisons)),]
## 33 SNPs that have a significant difference between the two allele frequencies 
plot_log_odds_vs_num_assocs_sig_only <- ggplot(data = sig_af_differences, aes(x=as.numeric(num_assocs), y =log(as.numeric(as.character(estimate.odds.ratio))), 
                                                                                    color = as.factor(Mass_pleio_SNP), shape = as.factor(sig_pval))) +
  geom_point() + 
  scale_colour_discrete(name = "", labels = c("Not mass pleio", "Mass pleio")) +
  scale_shape_discrete(name = "", labels = c("pvalue not sig", "pvalue sig"))


plot_log_odds_vs_num_assocs_sig_only

## Most of the SNPs that have sig diff allele frequenceis have zero protein associations 
# table(sig_af_differences$num_assocs)
# 0   1  2  3  78  8 
# 19  4  5  2  2   1 


library(ieugwasr)

#### 
batch_from_id <- function (id){
    sapply(strsplit(id, "-"), function(x) paste(x[1], x[2], sep = "-"))
  }

all_snps_func <- function(snp){
  o <- phewas(snp, pval=1, batch=c("prot-a", "prot-b", "prot-c"))
  o$batch <- batch_from_id(o$id)
  return(table(o$p < 0.05, o$batch))}

rsIDs <- over_50_with_study$SNP

all_snps <- lapply(rsIDs, all_snps_func)
names(all_snps) <- rsIDs

save(all_snps, file="over_50_snps_in_other_studies.rdata")

which_studies <- data.frame(rownames(summary(all_snps)), summary(all_snps)[,1])
names(which_studies) <- c("SNP", "Freq")

####
LD_proxies <- read.table("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/plink.tags.list_mass_pleio_LD_snps", header=T)
LD_proxies <- tibble::as_tibble(LD_proxies)
####

all_three_study_SNPs <- which_studies[which_studies$Freq == 6,]

LD_proxies_needed <- subset(LD_proxies, !LD_proxies$SNP %in% all_three_study_SNPs$SNP)

x <- strsplit(as.character(LD_proxies_needed$TAGS), "\\|")
names(x) <- LD_proxies_needed$SNP


SNPs_with_study_missing <- all_snps[LD_proxies_needed$SNP]
which_studies_missing <- data.frame(rownames(summary(SNPs_with_study_missing)), summary(SNPs_with_study_missing)[,1])
names(which_studies_missing) <- c("SNP", "Freq")


##############


prot_c_snp_list <- read.table("/Users/hw15842/Documents/PhD_work/protein_interactions/Data/protc_snplist.txt")$V1

protc_missing <- sapply(names(all_snps), function(x) {
  if(!"prot-c" %in% colnames(all_snps[[x]]))
  {
    return(x)
  }
}) %>% unlist()

LD_proxies_protc <- subset(LD_proxies, SNP %in% protc_missing)
x <- strsplit(as.character(LD_proxies_protc$TAGS), "\\|")
names(x) <- LD_proxies_protc$SNP



tags <- list()
for(i in names(x))
{
  tg <- x[[i]][x[[i]] %in% prot_c_snp_list]
  if(length(tg) >= 1)
  {
    tags[[i]] <- tibble(target=i, tag=tg[1])
  }
}

tags <- bind_rows(tags)

df <- lapply(tags$tag, all_snps_func)
names(df) <- tags$target

save(df, file="over_50_snps_in_other_studies_with_tagged_SNPs.rdata")

# make percentages

percent_func <- function(SNP){
  df_percent <- data.frame(((df[[SNP]][2]/(df[[SNP]][1]+df[[SNP]][2])))*100,
                           ((df[[SNP]][4]/(df[[SNP]][3]+df[[SNP]][4])))*100,
                           ((df[[SNP]][6]/(df[[SNP]][5]+df[[SNP]][6])))*100)
  names(df_percent) <- c("prot-a", "prot-b", "prot-c")
  return(df_percent)

}



percentages_taged_snps <- lapply(1:length(df), percent_func)
names(percentages_taged_snps) <- names(df)
percentages_taged_snps <- bind_rows(percentages_taged_snps, .id = "column_label")
percentages_taged_snps$tag_or_main <- "Tag"

other_snps <- all_snps[names(all_snps) %in% all_three_study_SNPs$SNP]

percent_func_other_snps <- function(SNP){
  other_snps_percent <- data.frame(((other_snps[[SNP]][2]/(other_snps[[SNP]][1]+other_snps[[SNP]][2])))*100,
                           ((other_snps[[SNP]][4]/(other_snps[[SNP]][3]+other_snps[[SNP]][4])))*100,
                           ((other_snps[[SNP]][6]/(other_snps[[SNP]][5]+other_snps[[SNP]][6])))*100)
  names(other_snps_percent) <- c("prot-a", "prot-b", "prot-c")
  return(other_snps_percent)
  
}

percentages_main_snps <- lapply(1:length(other_snps), percent_func_other_snps)
names(percentages_main_snps) <- names(other_snps)
percentages_main_snps <- bind_rows(percentages_main_snps, .id = "column_label")
percentages_main_snps$tag_or_main <- "Main"

percentages_SNPs <- rbind(percentages_main_snps, percentages_taged_snps)

save(percentages_SNPs, file = "percentages_SNPs_across_studies.rdata")


#### Looks like there are a lot more asociations in Sun compared to folkerson and suhre ####












################################
### Choose independent SNPs ####
################################

all_pQTLs <- Zheng_pQTLs_with_file_names[c("SNP", "CHR_SNP", "POS_SNP", "Study")]
all_pQTLs <- all_pQTLs[order(all_pQTLs$CHR_SNP, all_pQTLs$POS_SNP),]


over_50_with_study_sorted <- over_50_with_study[order(over_50_with_study$CHR_SNP, over_50_with_study$POS_SNP),]
rownames(over_50_with_study_sorted) <- NULL

### selecting the independent SNPs, if they are more than 100,000 base pairs apart from their closest SNP then independent ###
### the SNP selected is the one with the most protein associations ####
### if two have the same number of associations choose the Sun one ###
### if there are two Sun ones with the same amount of associations then pick the SNP that is more in the middle of all the SNPs location wise

independent_mass_pleio_snps <-  over_50_with_study_sorted[c(6,12,13,22,27,28,43,53,60,72),]   
save(independent_mass_pleio_snps, file = "independent_mass_pleio_snps.rdata")


#### Adding in the genes for each SNP ###

indep_SNPs_with_genes <- merge(independent_mass_pleio_snps, SNP_query_with_extra_ABO_snps, by="SNP", all.x=T)
indep_SNPs_with_genes <- subset(indep_SNPs_with_genes, select = -c(Freq.y, CHR_SNP.y, POS_SNP.y))
names(indep_SNPs_with_genes)[2:4] <- c("Freq", "CHR_SNP", "POS_SNP")
save(indep_SNPs_with_genes, file="indep_SNPs_with_genes.rdata")

indep_SNPs_with_genes_and_percent <- merge(indep_SNPs_with_genes, percentages_SNPs, by.x="SNP", by.y="column_label", all.x=T)   ## Tag here means that the tag SNP was used in SUhre t look up the SNP, not that the SNP is a tagged SNP
save(indep_SNPs_with_genes_and_percent, file="indep_SNPs_with_genes_and_percent.rdata")

### Add the gene info and percentages for all mass pleio SNPs ###
# doing this as not all the independent SNPs have percenatge data or gene info #

mass_pleio_SNPs_with_genes <- merge(over_50_with_study, SNP_query_with_extra_ABO_snps, by="SNP", all.x=T)
mass_pleio_SNPs_with_genes <- subset(mass_pleio_SNPs_with_genes, select = -c(Freq.y, CHR_SNP.y, POS_SNP.y))
names(mass_pleio_SNPs_with_genes)[2:4] <- c("Freq", "CHR_SNP", "POS_SNP")
save(mass_pleio_SNPs_with_genes, file="mass_pleio_SNPs_with_genes.rdata")

mass_pleio_SNPs_with_genes_and_percent <- merge(mass_pleio_SNPs_with_genes, percentages_SNPs, by.x="SNP", by.y="column_label", all.x=T)   ## Tag here means that the tag SNP was used in SUhre t look up the SNP, not that the SNP is a tagged SNP
save(mass_pleio_SNPs_with_genes_and_percent, file="mass_pleio_SNPs_with_genes_and_percent.rdata")

mass_pleio_SNPs_with_genes_and_percent_sorted <- mass_pleio_SNPs_with_genes_and_percent[order(mass_pleio_SNPs_with_genes_and_percent$CHR_SNP, mass_pleio_SNPs_with_genes_and_percent$POS_SNP),]



## The ABO snps seem to be the ones that are closest percentage wise between the three studies ###

## PLot the data ###
# three lines for each study percent, plot the independent SNPs gene next to it #
# or plot difference between SUn and Suhre and gene closest to... #


mass_pleio_SNPs_with_genes_and_percent_sorted$replicates <- ifelse(((mass_pleio_SNPs_with_genes_and_percent_sorted$`prot-a`) - (mass_pleio_SNPs_with_genes_and_percent_sorted$`prot-c`)) < 10,
                                                                   "Yes", "No")


indep_SNPs_with_genes_and_percent$replicates <- ifelse(((indep_SNPs_with_genes_and_percent$`prot-a`) - (indep_SNPs_with_genes_and_percent$`prot-c`)) < 10,
                                                                   "Yes", "No")

### might need to change indep snps based on closest genes as well ###


mass_pleio_SNPs_with_genes_and_percent_sorted$difference <- mass_pleio_SNPs_with_genes_and_percent_sorted$`prot-a` - mass_pleio_SNPs_with_genes_and_percent_sorted$`prot-c`

fit_percent_diff_vs_num_prot_assoc <- lm(difference ~ Freq, data=mass_pleio_SNPs_with_genes_and_percent_sorted)
regression_plot_fit_percent_diff_vs_num_prot_assoc <- ggplotRegression(fit_percent_diff_vs_num_prot_assoc) 
regression_plot_fit_percent_diff_vs_num_prot_assoc 


library(ggrepel)

row.names(mass_pleio_SNPs_with_genes_and_percent_sorted)<- NULL

label <- c(6,12,13,17,20,29,31,38,49,54,57,64,69)

mass_pleio_SNPs_with_genes_and_percent_sorted$label <- ifelse(row.names(mass_pleio_SNPs_with_genes_and_percent_sorted) %in% label, "YES", "NO")

plot_precentdiff_vs_prot_assoc <- ggplot(data = mass_pleio_SNPs_with_genes_and_percent_sorted, aes(x=as.numeric(Freq), y =difference, 
                                                                                                   color = as.factor(geneinfo),
                                                                                                   shape = as.factor(Study))) +
  geom_point(size = 3) + 
  geom_text_repel(aes(label=ifelse(label=="YES", geneinfo, ""),hjust=0,vjust=0)) 

plot_precentdiff_vs_prot_assoc



