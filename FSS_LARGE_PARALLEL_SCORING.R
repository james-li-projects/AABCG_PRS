library(data.table)
library(dplyr)
library(tidyr)

# importing passed argument for p-value threshold
args = commandArgs(trailingOnly=TRUE)
p_thresh = as.numeric(args[1])
p_thresh_string <- gsub("-","",toString(p_thresh))

# obtaining scoring file for Forward Stepwise Selection SNPs
for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  print(paste("GENERATING SCORING FILE FOR:",subtype,p_thresh_string))
  # loading in summary statistics
  sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/consistent_effect_allele/",subtype,".sumstats"),header=T) %>% arrange(`#CHROM`,POS) 
  
  # extracting lists of selected SNPs from the second window orientation 
  setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/Forward_Stepwise_Regression/",subtype,"/SELECTED_SNPS"))
  all_p_thresh_files <- list.files()[grep(p_thresh_string,list.files())]
  SELECTION_2_SNP_FILES <- all_p_thresh_files[grep("SELECTION_2_SNP_LIST_OFFSET_CHROM_WINDOW_SNP_LIST_pthresh_",all_p_thresh_files)]
  ENTIRE_SELECTION_2_SNP_LIST <- data.table(ID = as.character())
  for (snp_file in SELECTION_2_SNP_FILES) {
    load(snp_file)
    ENTIRE_SELECTION_2_SNP_LIST <- rbind(ENTIRE_SELECTION_2_SNP_LIST, SELECTION_2_SNP_LIST)
  }
  
  # filtering sumstats to only contain these selected SNPs from the second window orientation 
  sumstats <- sumstats %>% filter(ID %in% ENTIRE_SELECTION_2_SNP_LIST$ID)
  output_score_file <- sumstats %>% mutate(BETA=log(OR)) %>% select(ID,A1,BETA)
  write.table(output_score_file,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_FILES/",subtype,"_FSS_",p_thresh_string),row.names=F,col.names=F,quote=F,sep="\t")
  
}
