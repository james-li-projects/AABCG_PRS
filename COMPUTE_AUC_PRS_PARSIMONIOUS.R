# load libraries 
library(data.table)
library(dplyr)
library(devtools) # install.packages("devtools")
library(caret)
library(SuperLearner)
library(ranger)
library(glmnet)
library(ROCnReg)
library(stringr)
library(tidyr)

# set seed 
set.seed(1)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]

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
PRS_score_list <- PRS_score_list[grep(subtype,PRS_score_list)]
imported_scores <- vector(mode="list",length=length(PRS_score_list)+1)
num_regular_model <- length(PRS_score_list)
for (i in 1:num_regular_model) {
  print(i)
  imported_scores[[i]] <- fread(PRS_score_list[i],header=T) %>% select(`#IID`,SCORE1_AVG)
  colnames(imported_scores[[i]])[2] <- PRS_score_list[i]
}
# importing parsimonious scores
curr_index<-length(imported_scores)
imported_scores[[curr_index]] <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PARSIMONIOUS_330/AFR_SCORE_OUTPUT/",subtype,"_AFR_PARSIMONIOUS.sscore")) %>% select(`#IID`,SCORE1_AVG) 
colnames(imported_scores[[curr_index]])[2] <- paste0(subtype,"_AFR_PARSIMONIOUS")
PRS_score_list <- c(PRS_score_list,paste0(subtype,"_AFR_PARSIMONIOUS"))
# joining all scores
combined_PRS_scores <- Reduce(full_join,imported_scores)
combined_PRS_scores <- data.frame(combined_PRS_scores)
rownames(combined_PRS_scores) <- combined_PRS_scores$X.IID
colnames(combined_PRS_scores)[1] <- "#IID"

###################################################
# computing validation set AUCs
validation_auc_list <- vector(mode="list",length=length(PRS_score_list))
num_snp_list <- c()
for (i in 1:length(PRS_score_list)) {
  print(i)
  tmp_score_cov <- inner_join(imported_scores[[i]],y_vad,by=c("#IID"))
  colnames(tmp_score_cov)[2] <- "PRS"
  
  output_AROC.sp <- AROC.sp(formula.h = PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                            group = "Status",
                            data = tmp_score_cov,
                            tag.h = 0)
  validation_auc_list[[i]] <- sort(as.numeric(output_AROC.sp$AUC))
  
  # identifying number of SNPs in each model
  #name_approach <- colnames(imported_scores[[i]])[2]
  #name_approach <- str_split(name_approach,".sscore")[[1]][1]
  #approach_snp_list <- fread(paste0("../../SCORE_FILES_PROCESSED/",name_approach),header=F)
  #num_snp_list <- c(num_snp_list,length(unique(approach_snp_list$V1)))
}

# extracting validation set AUCs into a vector and ultimately into a data.frame
validation_auc_vec_lower <- c()
validation_auc_vec_center <- c()
validation_auc_vec_upper <- c()
for (i in 1:length(validation_auc_list)) {
  print(i)
  validation_auc_vec_lower <- c(validation_auc_vec_lower, validation_auc_list[[i]][1])
  validation_auc_vec_center <- c(validation_auc_vec_center, validation_auc_list[[i]][2])
  validation_auc_vec_upper <- c(validation_auc_vec_upper, validation_auc_list[[i]][3])
}
validation_AUC_df <- data.frame(colnames(combined_PRS_scores)[-1],validation_auc_vec_lower,validation_auc_vec_center,validation_auc_vec_upper)
colnames(validation_AUC_df) <- c("Approach","AUC_lower","AUC_center","AUC_upper")
print(validation_AUC_df %>% arrange(desc(AUC_center)))

###########################################################
# linearly combining genome-wide approaches with EUR PRS models in the tuning set
###########################################################
# extract subtype-specific PRS models
SUBTYPE_SPECIFIC_IID <- combined_PRS_scores %>% select(`#IID`)
SUBTYPE_SPECIFIC_AFR_SCORES <- combined_PRS_scores %>%
  select_if(grepl(subtype, names(.))) %>%
  select_if(!grepl("EUR", names(.))) %>%
  select_if(!grepl("PARSIMONIOUS", names(.))) %>%
  mutate(`#IID` = rownames(.))
SUBTYPE_SPECIFIC_EUR_SCORES <- combined_PRS_scores %>%
  select_if(grepl(subtype, names(.))) %>%
  select_if(grepl("PARSIMONIOUS", names(.))) %>% mutate(`#IID` = rownames(.))
# joining covariates of tuning set with PRS models
reg_df <- inner_join(SUBTYPE_SPECIFIC_IID,inner_join(SUBTYPE_SPECIFIC_EUR_SCORES,SUBTYPE_SPECIFIC_AFR_SCORES,by=c("#IID")),by=c("#IID"))
tune_reg_df <- inner_join(reg_df,y_tune,by=c("#IID"))
first_col_index_AFR_PRS <- head(which(grepl(subtype,colnames(tune_reg_df))),1)
last_col_index_AFR_PRS <- tail(which(grepl(subtype,colnames(tune_reg_df))),1) - first_col_index_AFR_PRS

# initializing vectors
name_vector_EUR <- c()
name_vector_AFR <- c()
coef_vector_EUR <- c()
coef_vector_AFR <- c()

# generating coefficients for linear combinations of XANCESTRY PRS
for(PRS_INDEX in c(first_col_index_AFR_PRS+1:last_col_index_AFR_PRS)) {
  # performing the linear combination of PRS models
  tmp_tune_reg_df <- tune_reg_df %>% dplyr::select(all_of(first_col_index_AFR_PRS),PRS_INDEX,Status,Age,Platform,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
  tmp_coef_EUR <- summary(glm(data=tmp_tune_reg_df,formula=Status ~ .,family="binomial"))$coefficients[2,1]
  tmp_coef_AFR <- summary(glm(data=tmp_tune_reg_df,formula=Status ~ .,family="binomial"))$coefficients[3,1]
  
  # storing coefficients
  coef_vector_EUR <- c(coef_vector_EUR,tmp_coef_EUR)
  coef_vector_AFR <- c(coef_vector_AFR,tmp_coef_AFR)
  
  # storing the name of the PRS models
  name_vector_EUR <- c(name_vector_EUR, colnames(tune_reg_df)[first_col_index_AFR_PRS])
  name_vector_AFR <- c(name_vector_AFR, colnames(tune_reg_df)[PRS_INDEX])
}

# assembling a data.frame of all the coefficients of PRS models
XANCESTRY_COEF_VEC <- data.frame(name_vector_EUR,coef_vector_EUR,name_vector_AFR,coef_vector_AFR)
colnames(XANCESTRY_COEF_VEC) <- c("NAME_EUR","COEF_EUR","NAME_AFR","COEF_AFR")

# joining covariates of validation set 
validation_reg_df <- inner_join(reg_df,y_vad,by=c("#IID"))
validation_prs_matrix <- matrix(ncol=nrow(XANCESTRY_COEF_VEC),nrow=nrow(validation_reg_df))

# computing linearly combined XANCESTRY PRS scores
for (i in 1:nrow(XANCESTRY_COEF_VEC)) {
  validation_prs_matrix[,i] <- validation_reg_df[,XANCESTRY_COEF_VEC$NAME_EUR[i]]*XANCESTRY_COEF_VEC$COEF_EUR[i] + validation_reg_df[,XANCESTRY_COEF_VEC$NAME_AFR[i]]*XANCESTRY_COEF_VEC$COEF_AFR[i]
}

# reformatting these XANCESTRY PRS scores into a data frame
validation_prs_df <- data.frame(validation_prs_matrix)
colnames(validation_prs_df) <- XANCESTRY_COEF_VEC$NAME_AFR
validation_prs_df$`#IID` <- validation_reg_df$`#IID`
rownames(validation_prs_df) <- validation_prs_df$`#IID`

# assessing AUCs of each linear combination
validation_prs_lincomb_df <- validation_reg_df %>%
  select_if(!grepl(subtype, names(.)))
validation_prs_lincomb_df <- inner_join(validation_prs_lincomb_df,validation_prs_df,by=c("#IID"))
first_col_index_XANCESTRY_PRS <- head(which(grepl(subtype,colnames(validation_prs_lincomb_df))),1)
last_col_index_XANCESTRY_PRS <- tail(which(grepl(subtype,colnames(validation_prs_lincomb_df))),1) 

# computing XANCESTRY set AUCs
XANCESTRY_auc_list <- vector(mode="list",length=length(c(first_col_index_XANCESTRY_PRS:last_col_index_XANCESTRY_PRS)))
num_snp_list <- c()
for (i in c(first_col_index_XANCESTRY_PRS:last_col_index_XANCESTRY_PRS)) {
  print(i - first_col_index_XANCESTRY_PRS + 1)
  tmp_score_cov <- validation_prs_lincomb_df %>% dplyr::select(Status,all_of(i),Age,Platform,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10)
  colnames(tmp_score_cov)[2] <- "PRS"
  output_AROC.sp <- AROC.sp(formula.h = PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                            group = "Status",
                            data = tmp_score_cov,
                            tag.h = 0)
  storing_index <- i - first_col_index_XANCESTRY_PRS + 1
  XANCESTRY_auc_list[[storing_index]] <- sort(as.numeric(output_AROC.sp$AUC))
  
  # determining number of snps in AFR model
  tmp_AFR_PRS_NAME <- colnames(validation_prs_lincomb_df)[i]
  tmp_AFR_PRS_NAME <- str_split(tmp_AFR_PRS_NAME,".sscore")[[1]][1]
  tmp_approach_snp_list_AFR <- fread(paste0("../../SCORE_FILES_PROCESSED/",tmp_AFR_PRS_NAME),header=F)
  # determining number of snps in EUR model
  tmp_EUR_PRS_NAME <- unique(name_vector_EUR)
  tmp_EUR_PRS_NAME <- str_split(tmp_EUR_PRS_NAME,".sscore")[[1]][1]
  tmp_approach_snp_list_EUR <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PARSIMONIOUS_330/AFR_SCORE_FILES/PARSIMONIOUS/",tmp_EUR_PRS_NAME),header=F)
  # determining number of snps in combined XANCESTRY model
  tmp_combined_snp_list <- unique(c(tmp_approach_snp_list_AFR$V1,tmp_approach_snp_list_EUR$V1))
  num_snp_list <- c(num_snp_list,length(tmp_combined_snp_list))
}

# extracting XANCESTRY set AUCs into a vector and ultimately into a data.frame
XANCESTRY_auc_vec_lower <- c()
XANCESTRY_auc_vec_center <- c()
XANCESTRY_auc_vec_upper <- c()
for (i in 1:length(XANCESTRY_auc_list)) {
  print(i)
  XANCESTRY_auc_vec_lower <- c(XANCESTRY_auc_vec_lower, XANCESTRY_auc_list[[i]][1])
  XANCESTRY_auc_vec_center <- c(XANCESTRY_auc_vec_center, XANCESTRY_auc_list[[i]][2])
  XANCESTRY_auc_vec_upper <- c(XANCESTRY_auc_vec_upper, XANCESTRY_auc_list[[i]][3])
}
XANCESTRY_AUC_df <- data.frame(colnames(validation_prs_lincomb_df[(first_col_index_XANCESTRY_PRS:last_col_index_XANCESTRY_PRS)]),XANCESTRY_auc_vec_lower,XANCESTRY_auc_vec_center,XANCESTRY_auc_vec_upper,num_snp_list)
colnames(XANCESTRY_AUC_df) <- c("Approach","AUC_lower","AUC_center","AUC_upper","Num_SNP")
print(XANCESTRY_AUC_df %>% arrange(desc(AUC_center)))


#########################################
for (row_index in 1:nrow(XANCESTRY_COEF_VEC)) {
  # assembling an input data.frame called cvfit_coef with column names of Name and Effect for each approach
  cvfit_coef_row1<-data.frame(XANCESTRY_COEF_VEC[row_index,c(1,2)])
  colnames(cvfit_coef_row1) <- c("Name","Effect")
  cvfit_coef_row2<-data.frame(XANCESTRY_COEF_VEC[row_index,c(3,4)])
  colnames(cvfit_coef_row2) <- c("Name","Effect")
  cvfit_coef <- rbind(cvfit_coef_row1,cvfit_coef_row2)
  
  #########################################
  # obtaining a final scoring file  
  SCORE_FILE_DF <- data.frame(
    V1 = as.character(),
    V2 = as.character(),
    V3 = as.numeric(),
    V4 = as.numeric()
  )
  
  ##########################################
  # obtaining weighted coefficients from the parsimonious model
  for (i in 1) {
    # printing out index of PRS model
    print(i)
    # creating string to read in scoring file
    TMP_SCORE_FILE_STR <- gsub(".sscore","",cvfit_coef$Name)[i]
    
    # obtaining weight of the current PRS model and multiplying it for all coefficients in the scoring file
    TMP_SCORE_WEIGHT <- cvfit_coef$Effect[i]
    TMP_SCORE_FILE <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PARSIMONIOUS_330/AFR_SCORE_FILES/PARSIMONIOUS/",TMP_SCORE_FILE_STR)) 
    
    # account for PRS score normalization factor (allele counts)
    TMP_NORM_FACTOR <- median((fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PARSIMONIOUS_330/AFR_SCORE_OUTPUT/",TMP_SCORE_FILE_STR,".sscore")))$ALLELE_CT)
    TMP_SCORE_FILE$V4 <- (TMP_SCORE_FILE$V3*TMP_SCORE_WEIGHT)/TMP_NORM_FACTOR
    
    # including the weighted scoring file to existing scoring file records 
    SCORE_FILE_DF <- rbind(SCORE_FILE_DF,TMP_SCORE_FILE)
  }
colnames(SCORE_FILE_DF) <- c("V1","V2","V3","V4")

  ##########################################
  # obtaining weighted coefficients from the AFR model
  for (i in 2) {
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
  
  ##########################################
  # collating all the PRS weights
  COLLATED_SCORE_FILE_DF <- data.table(SCORE_FILE_DF %>% select(-V3) %>% dplyr::group_by(V1,V2) %>% summarise(Effect = sum(V4))) 
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% tidyr::separate(V1, into = c("chr","pos","a2","a1"),remove=F,sep=":")
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% mutate(AlignedEffect = ifelse(a1 == V2, Effect, -Effect))
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% select(V1,a1,AlignedEffect) 
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% dplyr::group_by(V1,a1) %>% summarise(FinalEffect = sum(AlignedEffect))
  COLLATED_SCORE_FILE_DF <- data.frame(COLLATED_SCORE_FILE_DF)
  COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% filter(FinalEffect != 0)
  
  # writing out the final scoring file
  XANCESTRY_APPROACH_STRING <- gsub(".sscore","",cvfit_coef$Name[2])
  XANCESTRY_APPROACH_STRING <- gsub("TNBC_","",XANCESTRY_APPROACH_STRING)
  XANCESTRY_APPROACH_STRING <- gsub("ERNEG_","",XANCESTRY_APPROACH_STRING)
  XANCESTRY_APPROACH_STRING <- gsub("ERPOS_","",XANCESTRY_APPROACH_STRING)
  XANCESTRY_APPROACH_STRING <- gsub("OVERALL_","",XANCESTRY_APPROACH_STRING)
  write.table(COLLATED_SCORE_FILE_DF,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PARSIMONIOUS_330/FINAL_PARSIMONIOUS_MODELS/",subtype,"_XANCESTRY_",XANCESTRY_APPROACH_STRING),quote=F,row.names=F,col.names=F,sep="\t")
}
