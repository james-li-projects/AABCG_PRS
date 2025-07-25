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

################################
# Identifying eligible samples #
# for overall BC no miss recep #
################################
# initializing libraries
library(readxl)
library(data.table)
library(dplyr)
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# converting breast cancer status to a binary variable
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 2-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
# outputting final covariates df
covariates <- pheno_data %>% select(Sample_Name,Status,ER,PR,HER2) %>% rename(case.control=Status)
# recoding missing receptor status in cases and controls
covariates <- covariates %>% mutate(ER=ifelse(ER==8,NA,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==8,NA,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==8,NA,HER2))
covariates <- covariates %>% mutate(ER=ifelse(ER==9,888,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==9,888,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==9,888,HER2))
# making receptor negative a 0 value
covariates <- covariates %>% mutate(ER=ifelse(ER==2,0,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==2,0,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==2,0,HER2))
# obtaining list of samples with no missing receptor data
eligible_samples <- (covariates %>% filter(ER %in% c(1,0,NA)) %>% filter(PR %in% c(1,0,NA))  %>% filter(HER2 %in% c(1,0,NA)))$Sample_Name
################################

# set seed 
set.seed(1)

subtype_vec = c(
  "OVERALL",
  "ERPOS",
  "ERNEG",
  "TNBC"
)
approach_name_vec = c(
  "OVERALL_PROPTUNE_ALLSUBTYPE.sscore",
  "ERPOS_ENSEMBLE_FSS.sscore",
  "ERNEG_ENSEMBLE_GLMNET.sscore",
  "TNBC_XANCESTRY_PRSICE2.sscore"
)

OR_SD_STR_OUTPUT_VEC <- c()

for (model_index in 1:length(subtype_vec)) {
subtype=subtype_vec[model_index]
approach_name=approach_name_vec[model_index]
print(paste(subtype,approach_name))

# setting working directory based on subtype
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT/",subtype))

###################################################
# importing validation set covariates
y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
y_vad$Status <- y_vad$Status-1

# filtering these for non-missing receptor samples if examining the overall BC risk model 
if (subtype=="OVERALL") {
  y_vad <- y_vad %>% filter(`#IID` %in% eligible_samples)
}

###################################################
# importing calculated PRS scores
PRS_score_list <- list.files()[grep(approach_name,list.files())]

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

###############################
# computing odds ratio per SD #
###############################
or_sd_vad_reg_df <- reg_validate_df
control_SD <- sd((or_sd_vad_reg_df %>% filter(Status==0))[,1])
or_sd_vad_reg_df$PRS_COMBINED <- or_sd_vad_reg_df[,1]
or_sd_vad_reg_df$PRS_COMBINED_SD <- or_sd_vad_reg_df$PRS_COMBINED/control_SD
or_sd_fit <- summary(glm(data = or_sd_vad_reg_df, formula = Status ~ PRS_COMBINED_SD + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
# preparing confidence interval of OR
logOR_sd <- or_sd_fit["PRS_COMBINED_SD","Estimate"]
se_logOR_sd <- or_sd_fit["PRS_COMBINED_SD","Std. Error"]
qnorm_val_sd <- qnorm(.975)
logOR_sd_CI_95_lower_sd <- logOR_sd - qnorm_val_sd*se_logOR_sd
logOR_sd_CI_95_upper_sd <- logOR_sd + qnorm_val_sd*se_logOR_sd
# outputting results
OR_SD_STR <- paste0("OR per SD: ", sprintf("%.3f",exp(logOR_sd)), " 95% CI: ", sprintf("%.3f",exp(logOR_sd_CI_95_lower_sd)), ", ", sprintf("%.3f",exp(logOR_sd_CI_95_upper_sd)))
print(OR_SD_STR)
print(exp(logOR_sd))
OR_SD_STR_OUTPUT_VEC <- c(OR_SD_STR_OUTPUT_VEC,OR_SD_STR)

######################################
# computing odds ratio by percentile #
######################################
# computing percentiles in control group participants
cov_validation3 <- reg_validate_df %>% filter(Status == 0) 
cov_validation3$PRS_COMBINED <- cov_validation3[,1]
cov_validation3 <- cov_validation3 %>% arrange(PRS_COMBINED)
cov_validation3$PRS_COMBINED_PERCENTILE <- c(1:nrow(cov_validation3)) / nrow(cov_validation3) * 100 
# obtaining PRS scores corresponding to each of our percentiles of interest
normal_lower <- 40
normal_upper <- 60
PRS_normal_lower <- cov_validation3$PRS_COMBINED[which.min(abs(cov_validation3$PRS_COMBINED_PERCENTILE - normal_lower))]
PRS_normal_upper <- cov_validation3$PRS_COMBINED[which.min(abs(cov_validation3$PRS_COMBINED_PERCENTILE - normal_upper))]

# listing out potential lower and upper bounds of percentiles we are interested in
vec_percent_lower <- c(0, 1, 5, 0, 10, 20, 60, 80, 90, 90, 95, 95, 99)
vec_percent_upper <- c(1, 5, 10, 10, 20, 40, 80, 90, 100, 95, 100, 99, 100)
vec_percent_range <- paste(vec_percent_lower,vec_percent_upper, sep = "-")

# obtaining ORs at each percentile
PCT_range <- c()
PCT_lower <- c()
PCT_upper <- c()
OR_lower <- c()
OR_center <- c()
OR_upper <- c()

for (i in 1:length(vec_percent_range)) {
  # defining current lower and upper bound of PRS percentile
  current_lower_pct <- vec_percent_lower[i]
  current_upper_pct <- vec_percent_upper[i]
  current_percent_range <- vec_percent_range[i]
  
  # extracting the equivalent PRS of the upper and lower bounds of PRS percentiles
  current_PRS_lower <- cov_validation3$PRS_COMBINED[which.min(abs(cov_validation3$PRS_COMBINED_PERCENTILE - current_lower_pct))]
  current_PRS_upper <- cov_validation3$PRS_COMBINED[which.min(abs(cov_validation3$PRS_COMBINED_PERCENTILE - current_upper_pct))]
  
  # creating a temporary annotation data frame that indicates whether a participant has a PRS in the specified percentile range or in the normal PRS range (between 40% and 60%)
  tmp_validation3 <- reg_validate_df
  tmp_validation3$PRS_COMBINED <- tmp_validation3[,1]
  tmp_validation3 <- tmp_validation3 %>% mutate(PRS_PERCENTILE_GROUP = ifelse((PRS_COMBINED >= current_PRS_lower) & (PRS_COMBINED <= current_PRS_upper), 1, ifelse((PRS_COMBINED > PRS_normal_lower) & (PRS_COMBINED < PRS_normal_upper), 0, NA))) %>% filter(!is.na(PRS_PERCENTILE_GROUP))
  # tmp_validation3$Platform <- factor(tmp_validation3$Platform)
  # setting WGS as the reference group in the sequencing platform variable
  #tmp_validation3 <- within(tmp_validation3, Platform <- relevel(Platform, ref = "WGS"))
  
  # obtaining OR and CI for the given PRS_PERCENTILE_GROUP
  fit_summary <- summary(glm(data = tmp_validation3, formula = Status ~ PRS_PERCENTILE_GROUP + Age + Platform + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family = "binomial"))$coefficients
  
  logOR <- fit_summary[2,1]
  se_logOR <- fit_summary[2,2]
  qnorm_val <- qnorm(.975)
  
  logOR_CI_95_lower <- logOR - qnorm_val*se_logOR
  logOR_CI_95_upper <- logOR + qnorm_val*se_logOR
  
  # outputting results
  print(paste0(current_percent_range, " OR: ", exp(logOR), " 95% CI: ", exp(logOR_CI_95_lower), ", ", exp(logOR_CI_95_upper)))
  PCT_range <- c(PCT_range, current_percent_range)
  OR_lower <- c(OR_lower,exp(logOR_CI_95_lower))
  OR_center <- c(OR_center,exp(logOR))
  OR_upper <- c(OR_upper,exp(logOR_CI_95_upper))
  PCT_lower <- c(PCT_lower,current_lower_pct)
  PCT_upper <- c(PCT_upper,current_upper_pct)
}

# resulting data frame of odds ratios
OR_DF <- data.frame(PCT_range,PCT_lower,PCT_upper,OR_lower,OR_center,OR_upper)
OR_DF <- OR_DF %>% mutate(PCT_center = (PCT_upper+PCT_lower)/2) %>% filter(!(PCT_range %in% c("0-10","90-100","95-100")))
OR_DF <- rbind(OR_DF,data.frame(PCT_range="40-60",PCT_lower=40,PCT_upper=60,OR_lower=1,OR_center=1,OR_upper=1,PCT_center=50))
OR_DF$PCT_range <- paste0(OR_DF$PCT_range,"%")
OR_DF$PCT_center <- OR_DF$PCT_center + 1.5
OR_DF$PCT_center[OR_DF$PCT_range=="1-5%"] <- OR_DF$PCT_center[OR_DF$PCT_range=="1-5%"] + 1
OR_DF$PCT_center[OR_DF$PCT_range=="99-100%"] <- OR_DF$PCT_center[OR_DF$PCT_range=="99-100%"] + 1

png(paste0("../../CALIBRATION_FINAL_SCORE_OR/",approach_name,".png"),res=1200,units="in",height=3.5,width=11.5)
p <- ggplot(data=OR_DF,aes(x=PCT_center,y=OR_center)) + geom_line() + geom_errorbar(aes(ymin = OR_lower, ymax = OR_upper),width=0.4) + geom_text(aes(label=sprintf("%.3f",OR_center),vjust = -5),size=3) + coord_cartesian(ylim = c(0, ceiling(max(OR_DF$OR_upper))), xlim = c(0,100), expand = FALSE, clip = "off") + theme_classic() + theme(plot.margin = unit(c(1, 3, 4, 1), "lines"), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") + 
  
  annotate(geom = "text", x = OR_DF$PCT_center, y=-1, label=OR_DF$PCT_range, size=2) + 
  annotate(geom = "text", x = 50, y=-1.5, label="PRS Percentile Range",size=4) +
  ylab("Covariate-adjusted Odds Ratio") +
  geom_line(aes(y = 1), col = "blue")
print(p)
dev.off()

}