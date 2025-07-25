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

# set seed 
set.seed(1)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]
# subtype="OVERALL"

# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_OUTPUT/",subtype))

# initializing data.frame to store resulting AUCs
XENSEMBLE_AUC_df <- data.frame(Approach=as.character(),AUC_lower=as.numeric(),AUC_center=as.numeric(),AUC_upper=as.numeric())

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
# joining imported PRS scores with tuning covariates
scores_tune_df <- scores_tune
scores_tune_df$`#IID` <- rownames(scores_tune_df)
reg_tune_df<-inner_join(scores_tune_df,y_tune,by=c("#IID"))
# joining imported PRS scores with validate covariates
scores_validate_df <- scores_validate
scores_validate_df$`#IID` <- rownames(scores_validate_df)
reg_validate_df<-inner_join(scores_validate_df,y_vad,by=c("#IID"))


###########################################
# clumping PRS scores
first_clump_index = head(which(grepl("sscore",colnames(reg_tune_df))),1)
last_clump_index = tail(which(grepl("sscore",colnames(reg_tune_df))),1)
# generating p-values for every PRS model
clump_p_vec <- c()
remove_PRS_index_clump <- c()
for (clump_index in first_clump_index:last_clump_index) {
  print(clump_index)
  clump_PRS_name <- colnames(reg_tune_df)[clump_index]
  tmp_reg_tune_df <- reg_tune_df %>% dplyr::select(all_of(clump_index),Status,Age,Platform,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
  tmp_p_val <- summary(glm(data=tmp_reg_tune_df,formula=Status ~ .,family="binomial"))$coefficients[clump_PRS_name,4]
  clump_p_vec <- c(clump_p_vec,tmp_p_val)
}
# performing clumping
total_num_PRS <- length(c(first_clump_index:last_clump_index))
for (i in 1:(total_num_PRS-1)) {
  i_value <- reg_tune_df[,i]
  for (j in (i+1):total_num_PRS) {
    j_value <- reg_tune_df[,j]
    # making sure the SNPs that are compared are different SNPs
    if (i != j) {
      R2 = cor(i_value,j_value)^2
      i_moresig <- clump_p_vec[i] < clump_p_vec[j]
      if (((R2 > 0.95) == TRUE) & (i_moresig == FALSE)) {
        print(paste("Identified PRS",i,"for removal",R2,i_moresig)); 
        remove_PRS_index_clump <- c(remove_PRS_index_clump, i); 
        stop = TRUE; 
        break;
      } else {}  
    } else {}
  }
}
# removing clumped PRS models
reg_validate_df <- reg_validate_df %>% select(-all_of(remove_PRS_index_clump))
reg_tune_df <- reg_tune_df %>% select(-all_of(remove_PRS_index_clump))
###########################################


##########################################
####### FORWARD STEPWISE SELECTION #######
##########################################
data <- reg_tune_df %>% select(-(paste0("PC",c(11:40))))
platform_input<-data %>% select(`#IID`,Platform)
platform_output<-model.matrix(~ Platform - 1, data = platform_input)
platform_output<-platform_output[,-1]
data <- cbind(data,platform_output)
data <- data %>% select(-Platform,-`#IID`)

forced_variables<- colnames(data)
forced_variables<-forced_variables[(!grepl("sscore",forced_variables))]
forced_variables<-forced_variables[(!grepl("EUR",forced_variables))]
forced_variables<-forced_variables[(!grepl("Status",forced_variables))]
p_value_threshold<-0.05

###########################################
# listing out all the predictors
predictors <- setdiff(colnames(data), c("Status",forced_variables))  # Exclude forced variables
num_input_SNPs <- length(predictors)
selected_predictors <- character(0)  # Initialize an empty set of selected predictors
best_model <- NULL  # Initialize the best model
while (length(predictors) > 0) {
  min_p_value <- Inf
  best_predictor <- NULL
  
  for (i in seq_along(predictors)) {
    current_predictors <- c(forced_variables, selected_predictors, predictors[i])
    current_model <- glm(as.formula(paste("Status ~", paste(current_predictors, collapse = " + "))), 
                         data = data, 
                         family = binomial(link = "logit"))
    summary_info <- summary(current_model)
    
    # Extract p-value for the predictor being considered
    row_index_summary <- nrow(summary_info$coefficients)
    p_value <- summary_info$coefficients[row_index_summary, 4]
    
    if (p_value < min_p_value) {
      min_p_value <- p_value
      best_predictor <- predictors[i]
    }
  }
  
  if (min_p_value <= p_value_threshold) {
    selected_predictors <- c(selected_predictors, best_predictor)
    predictors <- setdiff(predictors, best_predictor)
    best_model <- glm(as.formula(paste("Status ~", paste(c(forced_variables, selected_predictors), collapse = " + "))), 
                      data = data, 
                      family = binomial(link = "logit"))
  } else {
    break
  }
}

# obtaining PRS model names and coefficients
best_model_coef_df <- data.frame(summary(best_model)$coefficients)
SELECTED_PRS_NAME_VEC <- rownames(best_model_coef_df)[grepl("sscore",rownames(best_model_coef_df))]
SELECTED_PRS_WEIGHT_VEC <- c()
for (i in 1:length(SELECTED_PRS_NAME_VEC)) {
  SELECTED_PRS_WEIGHT_VEC <- c(
    SELECTED_PRS_WEIGHT_VEC,
    best_model_coef_df[SELECTED_PRS_NAME_VEC[i],"Estimate"]
  )
}

# assembling data and covariates for the validation set
current_reg_validate_df <- reg_validate_df
SELECTED_PRS_LIST <- vector(mode="list",length=length(SELECTED_PRS_NAME_VEC))
for (i in 1:length(SELECTED_PRS_NAME_VEC)) {
  current_PRS_NAME <- SELECTED_PRS_NAME_VEC[i]
  SELECTED_PRS_LIST[[i]] <- current_reg_validate_df[,current_PRS_NAME] * SELECTED_PRS_WEIGHT_VEC[i]
}
SELECTED_PRS_VALIDATE_OUT_LIST <- transpose(SELECTED_PRS_LIST) %>% map(reduce, `+`)
SELECTED_PRS_VALIDATE_OUT_VEC <- unlist(SELECTED_PRS_VALIDATE_OUT_LIST)

current_reg_validate_df$SCORE <- 
  SELECTED_PRS_VALIDATE_OUT_VEC

# computing AUC of XSUBTYPE XSUBTYPE PRS
output_AROC.sp <- AROC.sp(formula.h = SCORE~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                          group = "Status",
                          data = current_reg_validate_df,
                          tag.h = 0)
print(paste("Covariate-adjusted AUC [ENSEMBLE | FSS]:", toString(sort(as.numeric(output_AROC.sp$AUC)))))

# storing AUC in final DF
XENSEMBLE_AUC_df[1,]<-data.frame(t(c("FSS",sort(as.numeric(output_AROC.sp$AUC)))))

#########################################
# assembling an input data.frame called cvfit_coef with column names of Name and Effect for each approach
cvfit_coef <- data.frame(SELECTED_PRS_NAME_VEC,SELECTED_PRS_WEIGHT_VEC)
colnames(cvfit_coef) <- c("Name","Effect")
print(cvfit_coef)

# obtaining a final scoring file from input data.frame called cvfit_coef with column names of Name and Effect
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
COLLATED_SCORE_FILE_DF <- data.frame(COLLATED_SCORE_FILE_DF)
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% filter(FinalEffect != 0)

# writing out the final scoring file
write.table(COLLATED_SCORE_FILE_DF,file=paste0("../../FINAL_PRS_MODELS/",subtype,"_ENSEMBLE_FSS"),quote=F,row.names=F,col.names=F,sep="\t")
#########################################


###################################################
###################################################
# performing 10-fold cross-validation with glmnet #
###################################################
###################################################
# setting up our glmnet inputs
joined_tuning_df <- reg_tune_df %>% select(-paste0("PC",c(1:40)),-"Age",-"Platform") %>% select(Status, everything())
joined_tuning_df_NO_SAMPLE_NAME <- joined_tuning_df %>% select(-`#IID`)
# running glmnet
Y=joined_tuning_df_NO_SAMPLE_NAME[,1]
X=as.matrix(joined_tuning_df_NO_SAMPLE_NAME[,-1])
cvfit<-cv.glmnet(x=X,y=Y,nfolds=10,type.measure="auc",family="binomial")
# obtaining coefficients of final model
cvfit$lambda.min
cvfit_coef <- data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
cvfit_coef$Name <- rownames(cvfit_coef)
cvfit_coef <- cvfit_coef %>% filter(s1 != 0) %>% mutate(Effect=s1) %>% select(Name,Effect)
cvfit_coef <- cvfit_coef[-1,]
print(cvfit_coef)

#########################################
# Evaluating model in the validation set
# assembling data and covariates for the validation set
current_reg_validate_df <- reg_validate_df
SELECTED_PRS_LIST <- vector(mode="list",length=nrow(cvfit_coef))
for (i in 1:nrow(cvfit_coef)) {
  current_PRS_NAME <- cvfit_coef$Name[i]
  SELECTED_PRS_LIST[[i]] <- current_reg_validate_df[,current_PRS_NAME] * cvfit_coef$Effect[i]
}
SELECTED_PRS_VALIDATE_OUT_LIST <- transpose(SELECTED_PRS_LIST) %>% map(reduce, `+`)
SELECTED_PRS_VALIDATE_OUT_VEC <- unlist(SELECTED_PRS_VALIDATE_OUT_LIST)

current_reg_validate_df$SCORE <- 
  SELECTED_PRS_VALIDATE_OUT_VEC

# computing AUC of ENSEMBLE PRS
output_AROC.sp <- AROC.sp(formula.h = SCORE~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                          group = "Status",
                          data = current_reg_validate_df,
                          tag.h = 0)
print(paste("Covariate-adjusted AUC [ENSEMBLE | GLMNET]:", toString(sort(as.numeric(output_AROC.sp$AUC)))))

# storing AUC in final DF
XENSEMBLE_AUC_df[2,]<-data.frame(t(c("GLMNET",sort(as.numeric(output_AROC.sp$AUC)))))
XENSEMBLE_AUC_df

#########################################


#########################################
# obtaining a final scoring file from input data.frame called cvfit_coef with column names of Name and Effect
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
COLLATED_SCORE_FILE_DF <- data.frame(COLLATED_SCORE_FILE_DF)
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% filter(FinalEffect != 0)

# writing out the final scoring file
write.table(COLLATED_SCORE_FILE_DF,file=paste0("../../FINAL_PRS_MODELS/",subtype,"_ENSEMBLE_GLMNET"),quote=F,row.names=F,col.names=F,sep="\t")
#########################################

save(XENSEMBLE_AUC_df,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_AUC/",subtype,"_XENSEMBLE_AUC_df.RData"))
