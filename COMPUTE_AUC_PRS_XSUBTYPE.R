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
library(stringr)
library(tidyr)

# set seed 
set.seed(1)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
# subtype="OVERALL"

# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_OUTPUT/",subtype))

###################################################
# importing testing set covariates
y_tune <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/testing_",subtype,".pheno_cov"),header=T)
y_tune$Status <- y_tune$Status-1
# importing validation set covariates
y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
y_vad$Status <- y_vad$Status-1

###################################################
# importing calculated PRS scores
PRS_score_list <- list.files()[grep("sscore",list.files())]
imported_scores <- vector(mode="list",length=length(PRS_score_list))
for (i in 1:length(PRS_score_list)) {
  print(i)
  imported_scores[[i]] <- fread(PRS_score_list[i],header=T) %>% select(`#IID`,SCORE1_AVG)
  colnames(imported_scores[[i]])[2] <- PRS_score_list[i]
}
combined_PRS_scores <- Reduce(full_join,imported_scores)
combined_PRS_scores <- data.frame(combined_PRS_scores)
rownames(combined_PRS_scores) <- combined_PRS_scores$X.IID
colnames(combined_PRS_scores)[1] <- "#IID"
scores_tune <- combined_PRS_scores[y_tune$`#IID`,] 
scores_validate <- combined_PRS_scores[y_vad$`#IID`,]
scores_tune <- scores_tune %>% select(-`#IID`)
scores_validate <- scores_validate %>% select(-`#IID`)
# joining imported PRS scores with covariates
scores_tune_df <- scores_tune
scores_tune_df$`#IID` <- rownames(scores_tune_df)
reg_tune_df<-inner_join(scores_tune_df,y_tune,by=c("#IID"))
scores_validate_df <- scores_validate
scores_validate_df$`#IID` <- rownames(scores_validate_df)
reg_validate_df<-inner_join(scores_validate_df,y_vad,by=c("#IID"))

###################################################
# obtaining a list of unique PRS approaches
PRS_APPROACH_VEC <- c()
initial_PRS_APPROACH_VEC <- colnames(combined_PRS_scores)[-1]
initial_PRS_APPROACH_LIST <- str_split(initial_PRS_APPROACH_VEC,".sscore")
for (i in 1:length(initial_PRS_APPROACH_VEC)) {
  PRS_APPROACH_VEC[i] <- initial_PRS_APPROACH_LIST[[i]][1]
}
# only grep AFR PRS
PRS_APPROACH_VEC<-PRS_APPROACH_VEC[(!grepl("EUR",PRS_APPROACH_VEC))]
# stripping subtype names
PRS_APPROACH_VEC <- gsub("ERNEG_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("ERPOS_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("OVERALL_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("TNBC_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- unique(PRS_APPROACH_VEC)

#############################################
# Computing AUCs
XSUBTYPE_approach_name <- c()
XSUBTYPE_auc_vec_lower <- c()
XSUBTYPE_auc_vec_center <- c()
XSUBTYPE_auc_vec_upper <- c()
num_snp_list <- c()
for (current_approach in PRS_APPROACH_VEC) {
  print(paste("COMPUTING AUC OF XSUBTYPE PRS FOR:",current_approach))
  OVERALL_NAME <- paste0("OVERALL_",current_approach,".sscore")
  ERPOS_NAME <- paste0("ERPOS_",current_approach,".sscore")
  ERNEG_NAME <- paste0("ERNEG_",current_approach,".sscore")
  TNBC_NAME <- paste0("TNBC_",current_approach,".sscore")
  current_reg_tune_df <- reg_tune_df %>% select(
    Status,
    EUR_313_OVERALL.sscore,
    EUR_313_ERPOS.sscore,
    EUR_313_ERNEG.sscore,
    EUR_330_TNBC.sscore,
    all_of(OVERALL_NAME),
    all_of(ERPOS_NAME),
    all_of(ERNEG_NAME),
    all_of(TNBC_NAME),
    Age,Platform,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
  
  current_reg_tune_formula <- paste0(
  "Status ~ EUR_313_OVERALL.sscore + EUR_313_ERPOS.sscore + EUR_313_ERNEG.sscore + EUR_330_TNBC.sscore + ",
  OVERALL_NAME,
  " + ",
  ERPOS_NAME,
  " + ",
  ERNEG_NAME,
  " + ",
  TNBC_NAME,
  " + ",
  "Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
  
  coef_fit <- 
    summary(glm(
      data=current_reg_tune_df,
      formula=current_reg_tune_formula,
      family="binomial"))$coefficients
  
  current_reg_validate_df <- reg_validate_df
  current_reg_validate_df$SCORE <- 
    coef_fit["EUR_313_OVERALL.sscore","Estimate"]*current_reg_validate_df$EUR_313_OVERALL.sscore + 
    coef_fit["EUR_313_ERPOS.sscore","Estimate"]*current_reg_validate_df$EUR_313_ERPOS.sscore + 
    coef_fit["EUR_313_ERNEG.sscore","Estimate"]*current_reg_validate_df$EUR_313_ERNEG.sscore + 
    coef_fit["EUR_330_TNBC.sscore","Estimate"]*current_reg_validate_df$EUR_330_TNBC.sscore + 
    coef_fit[OVERALL_NAME,"Estimate"]*current_reg_validate_df[,OVERALL_NAME] + 
    coef_fit[ERPOS_NAME,"Estimate"]*current_reg_validate_df[,ERPOS_NAME] +
    coef_fit[ERNEG_NAME,"Estimate"]*current_reg_validate_df[,ERNEG_NAME] + 
    coef_fit[TNBC_NAME,"Estimate"]*current_reg_validate_df[,TNBC_NAME]
  
  
  #########################################
  # obtaining a final scoring file from input data.frame called cvfit_coef with column names of Name and Effect
  cvfit_coef_col1 <- c("EUR_313_OVERALL.sscore","EUR_313_ERPOS.sscore","EUR_313_ERNEG.sscore","EUR_330_TNBC.sscore",OVERALL_NAME,ERPOS_NAME,ERNEG_NAME,TNBC_NAME)
  cvfit_coef_col2 <- c()
  for (current_coef_name in cvfit_coef_col1) {
    cvfit_coef_col2 <- c(cvfit_coef_col2, coef_fit[current_coef_name,"Estimate"])
  }
  cvfit_coef <- data.frame(cvfit_coef_col1,cvfit_coef_col2)
  colnames(cvfit_coef) <- c("Name","Effect")
  
  # start of code to obtain scoring data.frame
  SCORE_FILE_DF <- data.frame(
    V1 = as.character(),
    V2 = as.character(),
    V3 = as.numeric(),
    V4 = as.numeric()
  )
  
  # obtaining weighted coefficients from each PRS model
  for (i in 1:nrow(cvfit_coef)) {
    # printing out index of PRS model
    print(i)
    # creating string to read in scoring file
    if(grepl("EUR_330",cvfit_coef$Name[i])) {
      TMP_SCORE_FILE_STR <- gsub("\\.","-",gsub(".sscore","",cvfit_coef$Name))[i]
    } else {
      TMP_SCORE_FILE_STR <- gsub(".sscore","",cvfit_coef$Name)[i]
    }
    # obtaining weight of the current PRS model and multiplying it for all coefficients in the scoring file
    TMP_SCORE_WEIGHT <- cvfit_coef$Effect[i]
    TMP_SCORE_FILE <- fread(paste0("../../SCORE_FILES_PROCESSED/",TMP_SCORE_FILE_STR)) 
    
    # account for PRS score normalization factor (allele counts)
    TMP_NORM_FACTOR <- median((fread(paste0("../../SCORE_OUTPUT/",subtype,"/",TMP_SCORE_FILE_STR,".sscore")))$ALLELE_CT)
    TMP_SCORE_FILE$V4 <- (TMP_SCORE_FILE$V3*TMP_SCORE_WEIGHT)/TMP_NORM_FACTOR

    # including the weighted scoring file to existing scoring file records 
    SCORE_FILE_DF <- rbind(SCORE_FILE_DF,TMP_SCORE_FILE)
  }
  colnames(SCORE_FILE_DF) <- c("V1","V2","V3","V4")
  
  # collating all the PRS weights
  COLLATED_SCORE_FILE_DF <- data.table(SCORE_FILE_DF %>% select(-V3) %>% dplyr::group_by(V1,V2) %>% summarise(Effect = sum(V4))) 
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% tidyr::separate(V1, into = c("chr","pos","a2","a1"),remove=F,sep=":")
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% mutate(AlignedEffect = ifelse(a1 == V2, Effect, -Effect))
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% select(V1,a1,AlignedEffect) 
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% dplyr::group_by(V1,a1) %>% summarise(FinalEffect = sum(AlignedEffect))
  COLLATED_SCORE_FILE_DF <- data.table(COLLATED_SCORE_FILE_DF)
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% filter(FinalEffect != 0)
  
  # writing out the final scoring file
  write.table(COLLATED_SCORE_FILE_DF,file=paste0("../../FINAL_PRS_MODELS/",subtype,"_XSUBTYPE_",current_approach),quote=F,row.names=F,col.names=F,sep="\t")
  #########################################
  
    # computing AUC of XSUBTYPE XSUBTYPE PRS
  output_AROC.sp <- AROC.sp(formula.h = SCORE~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                            group = "Status",
                            data = current_reg_validate_df,
                            tag.h = 0)
  print(paste("Covariate-adjusted AUC:", toString(sort(as.numeric(output_AROC.sp$AUC)))))
  
  # storing name of approach and AUC values prior to generating final output data.frame
  XSUBTYPE_approach_name <- c(XSUBTYPE_approach_name,current_approach)
  current_XSUBTYPE_auc_vec <- sort(output_AROC.sp$AUC)
  XSUBTYPE_auc_vec_lower <- c(XSUBTYPE_auc_vec_lower, as.numeric(current_XSUBTYPE_auc_vec[1]))
  XSUBTYPE_auc_vec_center <- c(XSUBTYPE_auc_vec_center, as.numeric(current_XSUBTYPE_auc_vec[2]))
  XSUBTYPE_auc_vec_upper <- c(XSUBTYPE_auc_vec_upper, as.numeric(current_XSUBTYPE_auc_vec[3]))
  
  # plotting histogram of PRS scores
  current_reg_validate_df$Legend <- NA
  current_reg_validate_df$Legend[current_reg_validate_df$Status==1] <- "Case"
  current_reg_validate_df$Legend[current_reg_validate_df$Status==0] <- "Control"
  png(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_AUC/",subtype,"/DENSITY_XSUBTYPE_PRS_",current_approach,".png"),units="in", width=5, height=5, res=2000)
  print(ggplot(current_reg_validate_df, aes(x=SCORE, fill=Legend)) + geom_density(alpha=0.4) + xlab("PRS Score") + ylab("Density"))
  dev.off()
}

# outputting final data.frame of AUCs 
XSUBTYPE_AUC_df <- data.frame(
  XSUBTYPE_approach_name,
  XSUBTYPE_auc_vec_lower,
  XSUBTYPE_auc_vec_center,
  XSUBTYPE_auc_vec_upper
)
colnames(XSUBTYPE_AUC_df) <- c("Approach","AUC_lower","AUC_center","AUC_upper")
print(XSUBTYPE_AUC_df %>% arrange(desc(AUC_center)))
save(XSUBTYPE_AUC_df,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_AUC/",subtype,"_XSUBTYPE_AUC_df.RData"))
