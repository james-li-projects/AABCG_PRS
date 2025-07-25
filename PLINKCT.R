library(data.table)
library(dplyr)

subtype="OVERALL"

for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  print(subtype)
  setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PLINKCT/",subtype))
  sumstats<-fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/consistent_effect_allele/",subtype,".sumstats"))
  clump_list <- list.files()[grepl(".clumps",list.files())]
  
  for (clump_name in clump_list) {
    mod_clump_name <- paste0(subtype,"_","PLINKCT_",gsub(".clumps","",clump_name))
    print(mod_clump_name)
    clump_df <- fread(clump_name)
    clump_id_list <- clump_df$ID
    
    score_file <- sumstats %>% filter(ID %in% clump_id_list) %>% mutate(BETA=log(OR)) %>% select(ID,A1,BETA) %>% arrange(ID)
    
    write.table(score_file,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PLINKCT/SCORE_FILES/",mod_clump_name),row.names=F,col.names=F,quote=F,sep="\t")
  }
}

setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PLINKCT/SCORE_FILES"))

for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  score_file_list <- list.files()[grepl(subtype,list.files())]
  for (score_file_name in score_file_list) {
    print(score_file_name)
    pfile_path=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/split_validation_",subtype)
    output_path=paste0("../SCORE_OUTPUT/",score_file_name)
    system(paste0("plink2 --pfile ",pfile_path," --score ",score_file_name," --out ",output_path))
  }
}



library(ROCnReg)
library(pROC)
library(stringr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(data.table)

# set seed 
set.seed(1)

# importing subtype string
#!/usr/bin/env Rscript
subtype="TNBC"

# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PLINKCT/SCORE_OUTPUT"))

###################################################
# importing validation set covariates
y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
y_vad$Status <- y_vad$Status-1

###################################################
# importing calculated PRS scores
PRS_score_list <- list.files()[grep("sscore",list.files())]
PRS_score_list <- PRS_score_list[grepl(subtypem,PRS_score_list)]
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
scores_validate <- combined_PRS_scores[y_vad$`#IID`,]
scores_validate <- scores_validate %>% select(-`#IID`)
# joining imported PRS scores with covariates
scores_validate_df <- scores_validate
scores_validate_df$`#IID` <- rownames(scores_validate_df)
reg_validate_df<-inner_join(scores_validate_df,y_vad,by=c("#IID"))
colnames(reg_validate_df) <- gsub(paste0(subtype,"_"),"",colnames(reg_validate_df))
colnames(reg_validate_df) <- gsub(paste0(".sscore"),"",colnames(reg_validate_df))

###################################################
# obtaining a list of unique PRS approaches
PRS_APPROACH_VEC <- c()
initial_PRS_APPROACH_VEC <- colnames(combined_PRS_scores)[-1]
initial_PRS_APPROACH_LIST <- str_split(initial_PRS_APPROACH_VEC,".sscore")
for (i in 1:length(initial_PRS_APPROACH_VEC)) {
  PRS_APPROACH_VEC[i] <- initial_PRS_APPROACH_LIST[[i]][1]
}
# stripping subtype names
PRS_APPROACH_VEC <- gsub("ERNEG_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("ERPOS_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("OVERALL_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- gsub("TNBC_","",PRS_APPROACH_VEC)
PRS_APPROACH_VEC <- unique(PRS_APPROACH_VEC)

#############################################
# Computing AUCs
approach_name <- c()
auc_vec_lower <- c()
auc_vec_center <- c()
auc_vec_upper <- c()
num_snp_list <- c()
for (current_approach in PRS_APPROACH_VEC) {
  # current reg df
  current_reg_validate_df <- reg_validate_df
  current_reg_validate_df$SCORE1_AVG <- current_reg_validate_df[,current_approach]
  # computing AUC 
  current_reg_validate_df$Status <- as.factor(current_reg_validate_df$Status)
  output_AROC.sp <- AROC.sp(
    formula.h = SCORE1_AVG~Age+as.factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
    group = "Status",
    data = current_reg_validate_df,tag.h = 0)
  
  print("#################")
  print(current_approach)
  print(paste("Covariate-adjusted AUC:", toString(sort(as.numeric(output_AROC.sp$AUC)))))
  
  # storing name of approach and AUC values prior to generating final output data.frame
  approach_name <- c(approach_name,current_approach)
  current_auc_vec <- sort(output_AROC.sp$AUC)
  auc_vec_lower <- c(auc_vec_lower, as.numeric(current_auc_vec[1]))
  auc_vec_center <- c(auc_vec_center, as.numeric(current_auc_vec[2]))
  auc_vec_upper <- c(auc_vec_upper, as.numeric(current_auc_vec[3]))
}

# assembling final data.frame of AUCs 
AUC_DF <- data.frame(
  approach_name,
  auc_vec_lower,
  auc_vec_center,
  auc_vec_upper
)
colnames(AUC_DF) <- c("Approach","AUC_lower","AUC_center","AUC_upper")

# outputting the AUC DF
print(AUC_DF %>% arrange(desc(AUC_center)))










