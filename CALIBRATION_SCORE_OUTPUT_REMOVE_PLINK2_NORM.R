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

# iterating un-normalization of plink2 normalization on every score output file
for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  print(paste("OBTAINING CALIBRATED RESULTS WITH PLINK2 NORMALIZATION REMOVED FOR SUBTYPE:",subtype))
  
  # importing covariates for the validation set in order to identify scoring output for individuals in the validation set
  y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
  y_vad$Status <- y_vad$Status-1
  y_vad <- y_vad %>% select(-paste0("PC",11:40))
  
  # set working directory
  setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT_NORM_SNP/",subtype))
  
  # reading in a list of all the PRS models
  score_file_list <- list.files()[grepl("sscore",list.files())]
  
  # un-normalizing the score outputs
  for (file in score_file_list) {
    # reading in file
    SCORE_FILE <- fread(paste0(file)) 
    NORM_FACTOR <- median(SCORE_FILE$ALLELE_CT)
    SCORE_FILE$SCORE <- SCORE_FILE$SCORE1_AVG*NORM_FACTOR
    SCORE_FILE <- SCORE_FILE %>% select(`#IID`,SCORE)
    
    # joining with the validation set individual IDs and covariates
    JOINED_SCORE_FILE <- inner_join(SCORE_FILE,y_vad,by=c("#IID"))
    
    # writing out score file without plink normalization
    write.table(JOINED_SCORE_FILE,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT/",subtype,"/",file),row.names=F,col.names=T,sep="\t",quote=F)
  }
}


