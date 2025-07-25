library(data.table)
library(tidyverse)
library(readxl)
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,Status)
pheno_data$Status <- 3-pheno_data$Status
colnames(pheno_data)[1] <- "IID"

for (i in c(1:22)) {
  tmp_samp_list <- fread(paste0("chr",i,".samplist"),header=F)
  tmp_fam_df <- data.frame(
    FID = rep(0,nrow(tmp_samp_list)),
    IID = tmp_samp_list$V1,
    FATHER_ID = rep(0,nrow(tmp_samp_list)),
    MOTHER_ID = rep(0,nrow(tmp_samp_list)),
    SEX = rep(2,nrow(tmp_samp_list))
  )
  
  tmp_fam_pheno_df <- inner_join(tmp_fam_df, pheno_data, by = c("IID")) %>% select(FID,IID,FATHER_ID,MOTHER_ID,SEX,Status)
  write.table(tmp_fam_pheno_df, file = paste0("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr",i,".fam"), quote = F, row.names = F, col.names = F)
  print(paste0("Finished writing out chr",i,".fam"))
}
