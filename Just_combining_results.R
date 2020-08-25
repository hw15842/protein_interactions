######################################################################
### Reading all results files in and putting into one massive table ##
######################################################################

### Run interactively as issues with memory space ###

### For MR results for prot a ##

load_and_combine_MR <- function(prot_number){
  df<-get(load(paste0("prot-a-", prot_number, "_MR_results.rdata")))
  return(df)
}

prot_a_number_list <- 1:3282

MR_results_prot_a <- lapply(prot_a_number_list, load_and_combine_MR) 
MR_results_prot_a_table <- ldply(MR_results_prot_a, data.table)

### need to change all the sep-10 bits of this document to SEPT 10 (issue with the zheng pqtl.csv file - changed it to a date ## )
## Sep-10:SEPT10.12690.33.3 Sep-10:SEPT10.12690.33.3    --->   need to shange "Sep-10" to "SEPT10" in "ID" and "exposure" columns
MR_results_prot_a_table <- data.frame(lapply(MR_results_prot_a_table, function(x) {gsub("Sep-10", "SEPT10", x)}), stringsAsFactors=F)

save(MR_results_prot_a_table, file="all_MR_results_prot_a.rdata")
  
  
## For coloc results for prot a ##
  
load_and_combine_coloc <- function(prot_number){
  df<-get(load(paste0("prot-a-", prot_number, "_coloc_results.rdata")))
  df$prot_a_number <- paste0("PROT-a-", prot_number)
  return(df)
}

prot_a_number_list <- 1:3282

Coloc_results_prot_a <- lapply(prot_a_number_list, load_and_combine_coloc) 
Coloc_results_prot_a_table <- ldply(Coloc_results_prot_a, data.table)
Coloc_results_prot_a_table <- data.frame(lapply(Coloc_results_prot_a_table,unlist))
save(Coloc_results_prot_a_table, file="all_coloc_results_prot_a.rdata")
