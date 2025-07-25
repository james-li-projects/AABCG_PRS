library(dplyr)
library(data.table)
library(readxl)
library(ROCnReg)
library(iCARE)
library(stringr)

SCORE_FILE_LIST <- c(
  "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT_NORM_SNP/OVERALL/OVERALL_ENSEMBLE_GLMNET.sscore",
  "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT_NORM_SNP/ERPOS/ERPOS_ENSEMBLE_GLMNET.sscore",
  "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT_NORM_SNP/ERNEG/ERNEG_XSUBTYPE_FSS_0.05.sscore",
  "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT_NORM_SNP/TNBC/TNBC_XSUBTYPE_FSS_0.05.sscore"
)

for (SCORE_FILE_NAME in SCORE_FILE_LIST) {
  # print score file name
  print(SCORE_FILE_NAME)
  SUBTYPE_STR <- str_split(SCORE_FILE_NAME,pattern="/")[[1]][10]
  
  # importing scoring df for the PRS of interest
  SCORE_FILE <- fread(SCORE_FILE_NAME) 
  NORM_FACTOR <- median(SCORE_FILE$ALLELE_CT)
  SCORE_FILE$SCORE <- SCORE_FILE$SCORE1_AVG*NORM_FACTOR
  SCORE_FILE <- SCORE_FILE %>% select(`#IID`,SCORE)
  
  # importing all pheno data
  pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
  pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,FamHistory) %>% filter(FamHistory!=9)
  colnames(pheno_data)[1] <- "#IID"
  pheno_data$FamHistory <- pheno_data$FamHistory-1
  
  # importing testing data
  testing.pheno_cov <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/testing.pheno_cov",header=T)
  testing.pheno_cov$Status <- testing.pheno_cov$Status-1
  
  # joining testing data with phenotype and PRS scoring data
  reg_tune_df <- inner_join(inner_join(testing.pheno_cov,pheno_data,by=c("#IID")),SCORE_FILE,by=c("#IID"))
  
  # performing the linear combination of PRS and fam history
  reg_tune_df <- reg_tune_df %>% dplyr::select(SCORE,FamHistory,Status,Age,Platform,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
  
  # storing coefficients
  tmp_coef_PRS <- summary(glm(data=reg_tune_df,formula=Status ~ .,family="binomial"))$coefficients["SCORE","Estimate"]
  tmp_coef_FamHistory <- summary(glm(data=reg_tune_df,formula=Status ~ .,family="binomial"))$coefficients["FamHistory","Estimate"]
  
  #########################################
  # obtaining combined scores in the validation set
  # importing validation data
  validation.pheno_cov <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation.pheno_cov",header=T)
  validation.pheno_cov$Status <- validation.pheno_cov$Status-1
  
  # joining testing data with phenotype and PRS scoring data
  reg_vad_df <- inner_join(inner_join(validation.pheno_cov,pheno_data,by=c("#IID")),SCORE_FILE,by=c("#IID"))
  reg_vad_df$COMBINED_PRS_FAMHIST <- tmp_coef_PRS*reg_vad_df$SCORE + tmp_coef_FamHistory*reg_vad_df$FamHistory
  
  # computing AUC
  output_AROC.sp <- AROC.sp(
    formula.h = COMBINED_PRS_FAMHIST~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
    group = "Status",
    data = reg_vad_df,
    tag.h = 0)
  print(output_AROC.sp$AUC)
  
  # outputting validation results
  write.table(reg_vad_df%>%select(`#IID`,COMBINED_PRS_FAMHIST),file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/COMBINED_PRS_FAMHIST/",SUBTYPE_STR),quote=F,row.names=F,col.names=F,sep="\t")
  
  
  # #########################################
  # # computing absolute risk
  # age_group_data_Black <- read_excel("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/DH_code/BreastCancerIncidence2021.xlsx")
  # colnames(age_group_data_Black) <- c("Age","OVERALL_INCIDENCE","Total_mortality","BC_mortality","Other_mortality","ERPOS_INCIDENCE","ERNEG_INCIDENCE","TNBC_INCIDENCE")
  # age_group_data_Black <- as.data.frame(age_group_data_Black)
  # age_black1 <- age_group_data_Black[rep(row.names(age_group_data_Black), each = 5), ]
  # age_black1[, 2:3] <- age_black1[, 2:3]/1e5
  # 
  # bc_model_log_or_PRS <- summary(glm(data=reg_vad_df,formula=Status ~ COMBINED_PRS_FAMHIST+Age+Platform+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,family="binomial"))$coefficients["COMBINED_PRS_FAMHIST","Estimate"]
  # 
  # new_cov_prof <- PRS1 %>% select(OVERALL_ENSEMBLE_GLMNET)
  # ref_cov_dat <- PRS1 %>% filter(case==0) %>% select(OVERALL_ENSEMBLE_GLMNET)
  # bc_model_formula_PRS <- case ~ OVERALL_ENSEMBLE_GLMNET 
  # 
  # sublist20<-list(name="OVERALL_ENSEMBLE_GLMNET",type="continuous")
  # bc_model_cov_info_PRS <-list(sublist20)
  # 
  # bc_model_log_or_PRS <- 8.9841  ## see logistic regression above 
  # names(bc_model_log_or_PRS) = c('OVERALL_ENSEMBLE_GLMNET')
  # bc_model_log_or_PRS
  # 
  # res_covs_snps_PRS = computeAbsoluteRisk(
  #   model.formula=bc_model_formula_PRS, 
  #   model.cov.info=bc_model_cov_info_PRS, 
  #   model.log.RR=bc_model_log_or_PRS, 
  #   model.ref.dataset=ref_cov_dat, 
  #   model.disease.incidence.rates=bc_inc,model.competing.incidence.rates=mort_inc,
  #   apply.age.start=20, 
  #   apply.age.interval.length=60,
  #   apply.cov.profile=new_cov_prof, 
  #   return.refs.risk=TRUE
  #   )
}



