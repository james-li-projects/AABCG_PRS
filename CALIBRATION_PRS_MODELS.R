# load libraries 
library(data.table)
library(dplyr)
library(devtools) # install.packages("devtools")
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)
library(ROCnReg)
library(ggplot2)
library(purrr)
library(tidyr)
library(ggplot2)

# set seed 
set.seed(1)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]

#for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  print(paste("OBTAINING CALIBRATED RESULTS FOR SUBTYPE:",subtype))
  # set working directory
  setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_SCORE_OUTPUT_TRAINING/",subtype))
  
  # reading in a list of all the PRS models
  score_file_list <- list.files()[grepl("sscore",list.files())]
 
  # importing validation set covariates
  pheno_cov <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/","validation_OVERALL.pheno_cov"))
  pheno_cov$Status <- pheno_cov$Status - 1
  pheno_cov <- pheno_cov %>% select(-paste0("PC",11:40))
  
  # initializing vectors to store variances
  NAME <- c()
  CALIBRATION_FACTOR <- c()
  
  # computing variances
  for (file in score_file_list) {
    # storing filename into a vector
    NAME <- c(NAME,file)
    
    # printing out the PRS name
    print(file)

    # importing testing and validation set PRS scores
    TESTING_VALIDATION_SCORE_FILE <- fread(paste0("../../FINAL_SCORE_OUTPUT/",subtype,"/",file)) 
    TESTING_VALIDATION_NORM_FACTOR <- median(TESTING_VALIDATION_SCORE_FILE$ALLELE_CT)
    
    # removing the plink normalization based on the number of SNPs
    TESTING_VALIDATION_SCORE_FILE$SCORE_RAW <- TESTING_VALIDATION_SCORE_FILE$SCORE1_AVG*TESTING_VALIDATION_NORM_FACTOR
    
    # obtaining calibration factor
    TMP_JOINED_DF <- inner_join(TESTING_VALIDATION_SCORE_FILE,pheno_cov,by=c("#IID"))
    CURRENT_CALIBRATION_FACTOR <- summary(glm(data=TMP_JOINED_DF,formula=Status ~ SCORE_RAW+Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,family="binomial"))$coefficients["SCORE_RAW","Estimate"]
    CALIBRATION_FACTOR <- c(CALIBRATION_FACTOR,CURRENT_CALIBRATION_FACTOR)
    
    # obtaining metrics
    METRICS_TMP_JOINED_DF <- TMP_JOINED_DF
    METRICS_TMP_JOINED_DF <- METRICS_TMP_JOINED_DF %>% mutate(SCORE_RAW = SCORE_RAW*CURRENT_CALIBRATION_FACTOR)
    print(summary(glm(data=TMP_JOINED_DF,formula=Status ~ SCORE_RAW+Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,family="binomial"))$coefficients["SCORE_RAW","Estimate"])
    print(summary(glm(data=METRICS_TMP_JOINED_DF,formula=Status ~ SCORE_RAW+Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,family="binomial"))$coefficients["SCORE_RAW","Estimate"])
    print(sqrt(mean((METRICS_TMP_JOINED_DF%>%filter(Status==1))$SCORE_RAW)-mean((METRICS_TMP_JOINED_DF%>%filter(Status==0))$SCORE_RAW)))
  }
  
  # assembling calibration DF for all PRS
  CALIBRATION_DF <- data.frame(NAME,CALIBRATION_FACTOR)
  
  # generate calibrated scoring files
  for (i in seq(nrow(CALIBRATION_DF))) {
    # importing scoring file
    SCORE_FILE_PREFIX <- gsub(".sscore","",CALIBRATION_DF$NAME[i])
    CURRENT_SCORE_FILE <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/FINAL_PRS_MODELS/",SCORE_FILE_PREFIX)) 
    # calibrating effect sizes
    CURRENT_CALIBRATION_FACTOR <- CALIBRATION_DF$CALIBRATION_FACTOR[i]
    CURRENT_SCORE_FILE$V3 <- CURRENT_SCORE_FILE$V3*CURRENT_CALIBRATION_FACTOR
    # outputting calibrated scoring file
    write.table(CURRENT_SCORE_FILE,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_PRS_MODELS/",SCORE_FILE_PREFIX),col.names=F,row.names=F,sep="\t",quote=F)
  }
#}


