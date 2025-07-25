##########################
# Initializing libraries #
##########################
library(data.table)
library(tidyverse)
library(readxl)

# setting seed
set.seed(523)

####################################################################################
# Importing initial summary statistics and extracting dosages for significant CpGs #
####################################################################################
# importing passed argument for p-value threshold
args = commandArgs(trailingOnly=TRUE)
p_thresh = as.numeric(args[1])
p_thresh_string <- gsub("-","",toString(p_thresh))

# importing passed argument for subtype
subtype = args[2]

print("################################################")
print(paste0("STARTING PRS ANALYSIS FOR P-VALUE THRESHOLD OF: ",p_thresh_string))
print(Sys.time())
# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp/",subtype))

# reading initial iteration of summary statistics
print("IMPORTING TRAINING GWAS SUMMARY STATISTICS")
print(Sys.time())
sumstats <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/training_",subtype,"_PC_10.Status.glm.logistic.hybrid"),header=T)

# identifying all chromosomes that contain SNPs at a given p-value threshold
CHROM_INDEX_LIST <- unique((sumstats %>% filter(P<p_thresh))$`#CHROM`)

# writing out SNP lists for each chromosome window 
for (CHROM_INDEX in CHROM_INDEX_LIST) {
  # identifying SNPs under the p_threshold
  significant_SNP_df <- sumstats %>% filter(`#CHROM`==CHROM_INDEX) %>% filter(P<p_thresh) %>% arrange(`#CHROM`,POS)
  # creating a modified position value for each SNP that accounts for chromosomes being far apart prior to performing LD clumping
  significant_SNP_df_modified <- significant_SNP_df
  significant_SNP_df_modified <- significant_SNP_df_modified %>% mutate(index_SNP = seq(nrow(significant_SNP_df_modified)))
  significant_SNP_df_modified <- significant_SNP_df_modified %>% mutate(modified_POS = `#CHROM`*(1e14) + POS)
  # assigning each SNP to a window of 5 Mb
  significant_SNP_df_modified <- significant_SNP_df_modified %>% mutate(window = floor(modified_POS/5000000)) 
  # obtaining a list of SNPs for each window (code chunk directly below identifies unique windows)
  UNIQUE_WINDOW_LIST <- unique(significant_SNP_df_modified$window)
  TOTAL_WINDOW_NUM <- length(UNIQUE_WINDOW_LIST)
  # writing out chromosome window SNP lists
  for (index in seq(TOTAL_WINDOW_NUM)) {
    CHROM_WINDOW <- UNIQUE_WINDOW_LIST[index]
    CHROM_WINDOW_SNP_LIST <- significant_SNP_df_modified %>% filter(window == CHROM_WINDOW) %>% arrange(`#CHROM`,POS) %>% select(ID)
    write.table(CHROM_WINDOW_SNP_LIST,file = paste0("CHROM_WINDOW_SNP_LIST_pthresh_",p_thresh_string,"_chr_",CHROM_INDEX,"_window_",index,".txt"),quote=F,row.names=F,col.names=F)
    print(paste0("CHROM_WINDOW_SNP_LIST_pthresh_",p_thresh_string,"_chr_",CHROM_INDEX,"_window_",index))
  }
}
