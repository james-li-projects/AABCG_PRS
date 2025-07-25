library(dplyr)
library(data.table)
library(readxl)

# set working directory
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/ARISK_FAMH_SCORE_OUTPUT")

# importing all pheno data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,FamHistory,AFR_pro) %>% filter(FamHistory!=9)
colnames(pheno_data)[1] <- "#IID"
#pheno_data$Status <- 3-pheno_data$Status
pheno_data$FamHistory <- pheno_data$FamHistory-1

prs_score_list <- c(
  "OVERALL_PROPTUNE_ALLSUBTYPE.sscore",
  "ERPOS_ENSEMBLE_FSS.sscore",
  "ERNEG_ENSEMBLE_GLMNET.sscore",
  "TNBC_XANCESTRY_PRSICE2.sscore"
)

############################################
# loop to calculate fam history metrics
for (i in seq(length(prs_score_list))) {
  
  # print PRS name
  print(prs_score_list[i])
  
  # import current scoring file
  score_df <- fread(prs_score_list[i])
  
  # joining scoring file data with pheno data 
  JOINED_DF <- inner_join(score_df,pheno_data,by=c("#IID"))
  # obtaining SD standardized PRS scores
  JOINED_DF$SCORE_SD <- JOINED_DF$SCORE/sd(JOINED_DF$SCORE)
  
  # unadjusted for family history 
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ SCORE_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  logOR_sd <- fit["SCORE_SD","Estimate"]
  se_logOR_sd <- fit["SCORE_SD","Std. Error"]
  qnorm_val_sd <- qnorm(.975)
  logOR_sd_CI_95_lower_sd <- logOR_sd - qnorm_val_sd*se_logOR_sd
  logOR_sd_CI_95_upper_sd <- logOR_sd + qnorm_val_sd*se_logOR_sd
  OR <- sprintf("%.2f",exp(logOR_sd))
  LOWER <- sprintf("%.2f",exp(logOR_sd_CI_95_lower_sd))
  UPPER <- sprintf("%.2f",exp(logOR_sd_CI_95_upper_sd))
  print(paste0(OR, " (",LOWER,"-",UPPER,")"))
  
  # no family history 
  fit <- summary(glm(data = JOINED_DF %>% filter(FamHistory==0), formula = Status ~ SCORE_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  logOR_sd <- fit["SCORE_SD","Estimate"]
  se_logOR_sd <- fit["SCORE_SD","Std. Error"]
  qnorm_val_sd <- qnorm(.975)
  logOR_sd_CI_95_lower_sd <- logOR_sd - qnorm_val_sd*se_logOR_sd
  logOR_sd_CI_95_upper_sd <- logOR_sd + qnorm_val_sd*se_logOR_sd
  OR <- sprintf("%.2f",exp(logOR_sd))
  LOWER <- sprintf("%.2f",exp(logOR_sd_CI_95_lower_sd))
  UPPER <- sprintf("%.2f",exp(logOR_sd_CI_95_upper_sd))
  print(paste0(OR, " (",LOWER,"-",UPPER,")"))
  
  # has family history 
  fit <- summary(glm(data = JOINED_DF %>% filter(FamHistory==1), formula = Status ~ SCORE_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  logOR_sd <- fit["SCORE_SD","Estimate"]
  se_logOR_sd <- fit["SCORE_SD","Std. Error"]
  qnorm_val_sd <- qnorm(.975)
  logOR_sd_CI_95_lower_sd <- logOR_sd - qnorm_val_sd*se_logOR_sd
  logOR_sd_CI_95_upper_sd <- logOR_sd + qnorm_val_sd*se_logOR_sd
  OR <- sprintf("%.2f",exp(logOR_sd))
  LOWER <- sprintf("%.2f",exp(logOR_sd_CI_95_lower_sd))
  UPPER <- sprintf("%.2f",exp(logOR_sd_CI_95_upper_sd))
  print(paste0(OR, " (",LOWER,"-",UPPER,")"))
  
  # testing for interaction between PRS and family history
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ SCORE_SD*FamHistory + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  print(sprintf("%.2f",fit["SCORE_SD:FamHistory","Pr(>|z|)"]))
  
  # family history unadjusted for PRS
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ FamHistory + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  logOR_FamHistory <- fit["FamHistory","Estimate"]
  se_logOR_FamHistory <- fit["FamHistory","Std. Error"]
  qnorm_val_FamHistory <- qnorm(.975)
  logOR_FamHistory_CI_95_lower_FamHistory <- logOR_FamHistory - qnorm_val_FamHistory*se_logOR_FamHistory
  logOR_FamHistory_CI_95_upper_FamHistory <- logOR_FamHistory + qnorm_val_FamHistory*se_logOR_FamHistory
  OR <- sprintf("%.2f",exp(logOR_FamHistory))
  LOWER <- sprintf("%.2f",exp(logOR_FamHistory_CI_95_lower_FamHistory))
  UPPER <- sprintf("%.2f",exp(logOR_FamHistory_CI_95_upper_FamHistory))
  print(paste0(OR, " (",LOWER,"-",UPPER,")"))
  
  # family history adjusted for PRS
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ FamHistory + SCORE_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  logOR_FamHistory <- fit["FamHistory","Estimate"]
  se_logOR_FamHistory <- fit["FamHistory","Std. Error"]
  qnorm_val_FamHistory <- qnorm(.975)
  logOR_FamHistory_CI_95_lower_FamHistory <- logOR_FamHistory - qnorm_val_FamHistory*se_logOR_FamHistory
  logOR_FamHistory_CI_95_upper_FamHistory <- logOR_FamHistory + qnorm_val_FamHistory*se_logOR_FamHistory
  OR <- sprintf("%.2f",exp(logOR_FamHistory))
  LOWER <- sprintf("%.2f",exp(logOR_FamHistory_CI_95_lower_FamHistory))
  UPPER <- sprintf("%.2f",exp(logOR_FamHistory_CI_95_upper_FamHistory))
  print(paste0(OR, " (",LOWER,"-",UPPER,")"))
}


###############################
# loop to calculate ancestry metrics
# re-importing all pheno data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,AFR_pro)
colnames(pheno_data)[1] <- "#IID"

for (i in seq(length(prs_score_list))) {
  
  # print PRS name
  print(prs_score_list[i])
  
  # import current scoring file
  score_df <- fread(prs_score_list[i])
  
  # joining scoring file data with pheno data 
  JOINED_DF <- inner_join(score_df,pheno_data,by=c("#IID"))
  # obtaining SD standardized PRS scores
  JOINED_DF$SCORE_SD <- JOINED_DF$SCORE/sd(JOINED_DF$SCORE)
  # obtaining binary ancestry variable of 80%
  JOINED_DF <- JOINED_DF %>% mutate(BIN_ANCESTRY = ifelse(AFR_pro > 0.80,1,0)) 
  
  # interaction with ancestry
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ SCORE_SD*AFR_pro + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  print(fit["SCORE_SD:AFR_pro","Estimate"])
  print(fit["SCORE_SD:AFR_pro","Pr(>|z|)"])
}


###############################
# loop to calculate age metrics
for (i in seq(length(prs_score_list))) {
  
  # print PRS name
  print(prs_score_list[i])
  
  # import current scoring file
  score_df <- fread(prs_score_list[i])
  
  # joining scoring file data with pheno data 
  JOINED_DF <- inner_join(score_df,pheno_data,by=c("#IID"))
  # obtaining SD standardized PRS scores
  JOINED_DF$SCORE_SD <- JOINED_DF$SCORE/sd(JOINED_DF$SCORE)
  
  # interaction with age
  fit <- summary(glm(data = JOINED_DF, formula = Status ~ SCORE_SD*Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  print(sprintf("%.3f",fit["SCORE_SD:Age","Estimate"]))
  print(sprintf("%.3f",fit["SCORE_SD:Age","Pr(>|z|)"]))
}
