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

# making a data frame to store all the subtype ENSEMBLE AUCs
COMBINED_FINAL_AUC_DF <- data.frame(
  Approach = as.character(),
  Method = as.character(),
  AUC_lower = as.numeric(),
  AUC_center = as.numeric(),
  AUC_upper = as.numeric(),
  logOR_sd_CI_95_lower_sd = as.numeric(),
  logOR_sd = as.numeric(),
  logOR_sd_CI_95_upper_sd = as.numeric(),
  Num_SNPs = as.numeric(),
  Subtype = as.character()
)

for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT/",subtype))

# initializing data.frame to store resulting AUCs
AUC_DF <- data.frame(Approach=as.character(),AUC_lower=as.numeric(),AUC_center=as.numeric(),AUC_upper=as.numeric(),Num_SNPs=as.numeric())

###################################################
# importing validation set covariates
y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
y_vad$Status <- y_vad$Status-1

###################################################
# importing calculated PRS scores
PRS_score_list <- list.files()[grep("sscore",list.files())]
PRS_score_list <- PRS_score_list[grep(subtype,PRS_score_list)]
PRS_score_list <- PRS_score_list[grep("ENSEMBLE",PRS_score_list)]

imported_scores <- vector(mode="list",length=length(PRS_score_list))
for (i in 1:length(PRS_score_list)) {
  print(i)
  imported_scores[[i]] <- fread(PRS_score_list[i],header=T) %>% select(`#IID`,SCORE)
  colnames(imported_scores[[i]])[2] <- PRS_score_list[i]
}
combined_PRS_scores <- Reduce(full_join,imported_scores)
combined_PRS_scores <- data.frame(combined_PRS_scores)
rownames(combined_PRS_scores) <- combined_PRS_scores$X.IID
colnames(combined_PRS_scores)[1] <- "#IID"
scores_validate <- combined_PRS_scores[y_vad$`#IID`,]
scores_validate <- scores_validate %>% select(-`#IID`)

# joining imported PRS scores with validate covariates
scores_validate_df <- scores_validate
scores_validate_df$`#IID` <- rownames(scores_validate_df)
reg_validate_df<-inner_join(scores_validate_df,y_vad,by=c("#IID"))

# assembling data and covariates for the validation set
current_reg_validate_df <- reg_validate_df
SELECTED_PRS_LIST <- colnames(current_reg_validate_df)[grep(paste0(subtype,"_"),colnames(current_reg_validate_df))]

# computing metrics of all PRS models for the specified subtype
for (m in c(1:length(SELECTED_PRS_LIST))) {
  print(paste("PROCESSING PRS",m,"of",length(SELECTED_PRS_LIST)))
  SELECTED_PRS <- SELECTED_PRS_LIST[m]
  
  # computing covariate adjusted AUC
  output_AROC.sp <- AROC.sp(
  formula.h = paste0(SELECTED_PRS,"~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"),
    group = "Status",
    data = current_reg_validate_df,
    tag.h = 0)
  print(paste("Covariate-adjusted AUC for",SELECTED_PRS,"method:", toString(signif(sort(as.numeric(output_AROC.sp$AUC)), digits = 4))))
  
  # identifying number of SNPs
  SELECTED_PRS_STR <- gsub(".sscore","",SELECTED_PRS)
  SELECTED_PRS_NUM_SNPS <- length(unique((fread(paste0("../../CALIBRATION_FINAL_PRS_MODELS/",SELECTED_PRS_STR)))$V1))
  
  
  ###############################
  # computing odds ratio per SD #
  ###############################
  or_sd_vad_reg_df <- reg_validate_df
  or_sd_vad_reg_df$SCORE <- or_sd_vad_reg_df[,SELECTED_PRS]
  control_SD <- sd((or_sd_vad_reg_df %>% filter(Status==0))$SCORE)
  or_sd_vad_reg_df$SCORE_SD <- or_sd_vad_reg_df$SCORE/control_SD
  or_sd_fit <- summary(glm(data = or_sd_vad_reg_df, formula = Status ~ SCORE_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  # preparing confidence interval of OR
  tmp_logOR_sd <- or_sd_fit["SCORE_SD","Estimate"]
  tmp_se_logOR_sd <- or_sd_fit["SCORE_SD","Std. Error"]
  qnorm_val_sd <- qnorm(.975)
  tmp_logOR_sd_CI_95_lower_sd <- tmp_logOR_sd - qnorm_val_sd*tmp_se_logOR_sd
  tmp_logOR_sd_CI_95_upper_sd <- tmp_logOR_sd + qnorm_val_sd*tmp_se_logOR_sd
  # storing results 
  logOR_sd <- c(logOR_sd,tmp_logOR_sd)
  logOR_sd_CI_95_lower_sd <- c(logOR_sd_CI_95_lower_sd,tmp_logOR_sd_CI_95_lower_sd)
  logOR_sd_CI_95_upper_sd <- c(logOR_sd_CI_95_upper_sd,tmp_logOR_sd_CI_95_upper_sd)
  
  
  TMP_AUC_DF<-data.frame(t(c(SELECTED_PRS,sort(as.numeric(output_AROC.sp$AUC)),logOR_sd_CI_95_lower_sd=tmp_logOR_sd_CI_95_lower_sd,logOR_sd=tmp_logOR_sd,logOR_sd_CI_95_upper_sd=tmp_logOR_sd_CI_95_upper_sd,SELECTED_PRS_NUM_SNPS)))
  colnames(TMP_AUC_DF) <- c("Approach","AUC_lower","AUC_center","AUC_upper","logOR_lower","logOR_center","logOR_upper","Num_SNPs")
  AUC_DF <- rbind(AUC_DF,TMP_AUC_DF)
}

# polishing the AUC DF
FINAL_AUC_DF <- AUC_DF
FINAL_AUC_DF$AUC_lower <- signif(as.numeric(FINAL_AUC_DF$AUC_lower), digits = 3)
FINAL_AUC_DF$AUC_center <- signif(as.numeric(FINAL_AUC_DF$AUC_center), digits = 3)
FINAL_AUC_DF$AUC_upper <- signif(as.numeric(FINAL_AUC_DF$AUC_upper), digits = 3)
FINAL_AUC_DF$Approach <- gsub("e","e-",gsub(paste0(subtype,"_"),"",gsub(".sscore","",FINAL_AUC_DF$Approach)))
# FINAL_AUC_DF %>% arrange(AUC_center)
FINAL_AUC_DF <- FINAL_AUC_DF %>% separate(Approach,sep="_",into=c("Class","Method"),remove=F) %>% select(-Class)
FINAL_AUC_DF$Subtype<-subtype

# storing the AUCs of the given subtype into the combined data frame
COMBINED_FINAL_AUC_DF <- rbind(COMBINED_FINAL_AUC_DF, FINAL_AUC_DF)
}

# make OR columns
COMBINED_FINAL_AUC_DF$OR_lower <- exp(as.numeric(COMBINED_FINAL_AUC_DF$logOR_lower))
COMBINED_FINAL_AUC_DF$OR_center <- exp(as.numeric(COMBINED_FINAL_AUC_DF$logOR_center))
COMBINED_FINAL_AUC_DF$OR_upper <- exp(as.numeric(COMBINED_FINAL_AUC_DF$logOR_upper))

# make OR columns exactly two decimal places
OR_cols <- c("OR_lower", "OR_center", "OR_upper")
COMBINED_FINAL_AUC_DF[OR_cols] <- lapply(COMBINED_FINAL_AUC_DF[OR_cols], function(x) ifelse(is.na(x), NA, sprintf("%.2f", round(x, 2))))

# making a combined OR column
COMBINED_FINAL_AUC_DF <- COMBINED_FINAL_AUC_DF %>% mutate(OR_combined=paste0(OR_center,"\n","(",OR_lower,", ",OR_upper,")"))

save(COMBINED_FINAL_AUC_DF,file=paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/ENSEMBLE_COMBINED_FINAL_AUC_DF.RData"))

######################################################
######################################################
######################################################
# parsing results
PARSED_DF <- COMBINED_FINAL_AUC_DF 
PARSED_DF$Approach <- paste(PARSED_DF$Subtype,PARSED_DF$Approach,sep="_")
PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="FSS","Forward Stepwise Selection","Elastic Net Regression")) %>% mutate(Subtype=ifelse(Subtype=="OVERALL","Overall",ifelse(Subtype=="ERPOS","ER-positive",ifelse(Subtype=="ERNEG","ER-negative","Triple-negative breast cancer"))))
PARSED_DF$Approach <- factor(PARSED_DF$Approach, levels = PARSED_DF$Approach)

# obtain labels for the different PRS methods
UNIQUE_SUBTYPE_LIST <- unique(PARSED_DF$Subtype)
SUBTYPE_PLOT_INDEX <- c()
for (UNIQUE_SUBTYPE in UNIQUE_SUBTYPE_LIST) {
  SUBTYPE_PLOT_INDEX <- c(SUBTYPE_PLOT_INDEX,median(which(PARSED_DF$Subtype==UNIQUE_SUBTYPE)))
}

png(paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/ENSEMBLE.png"),res=1200,units="in",height=4,width=11)
ggplot(data=PARSED_DF,aes(x=Approach,y=AUC_center,col=Subtype)) + geom_point() + geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper),width=0.4) + geom_text(aes(label=sprintf("%.3f",AUC_center),vjust = -5),size=3) + coord_cartesian(ylim = c(0.5, 0.7), xlim = c(0,nrow(PARSED_DF)+1), expand = FALSE, clip = "off") + theme_classic() + theme(plot.margin = unit(c(1, 8.5, 6.5, 1), "lines"), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") + 
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.48, label=prettyNum(PARSED_DF$Num_SNPs, big.mark = ",", scientific = FALSE),size=2) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+2, y=0.48, label="Number of SNPs in PRS",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.46, label=PARSED_DF$OR_combined,size=2) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+2, y=0.46, label="OR per SD",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.44, label=PARSED_DF$Method,size=2) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+2, y=0.44, label="PRS method",size=3) +
  
  annotate(geom = "text", x = SUBTYPE_PLOT_INDEX, y=0.42, label=UNIQUE_SUBTYPE_LIST,size=3) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+2, y=0.42, label="Breast Cancer Subtype",size=3) +
  
  ylab("Covariate-adjusted AUC")
dev.off()

