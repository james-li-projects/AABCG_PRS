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

# set working directory
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT/",subtype))

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
# removing proportion-tuned PRS models from this list
PRS_score_list <- PRS_score_list[!grepl("PROP",PRS_score_list)]
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
logOR_sd <- c()
logOR_sd_CI_95_lower_sd <- c()
logOR_sd_CI_95_upper_sd <- c()

for (current_approach in PRS_APPROACH_VEC) {
  # current reg df
  current_reg_validate_df <- reg_validate_df
  current_reg_validate_df$SCORE <- current_reg_validate_df[,current_approach]
  # computing AUC 
  output_AROC.sp <- AROC.sp(formula.h = SCORE~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
                            group = "Status",
                            data = current_reg_validate_df,
                            tag.h = 0)
  current_auc_vec <- sort(output_AROC.sp$AUC)
  
  print(paste("Covariate-adjusted AUC:", toString(sort(as.numeric(output_AROC.sp$AUC)))))
  
  # storing name of approach and AUC values prior to generating final output data.frame
  approach_name <- c(approach_name,current_approach)
  current_auc_vec <- sort(output_AROC.sp$AUC)
  auc_vec_lower <- c(auc_vec_lower, as.numeric(current_auc_vec[1]))
  auc_vec_center <- c(auc_vec_center, as.numeric(current_auc_vec[2]))
  auc_vec_upper <- c(auc_vec_upper, as.numeric(current_auc_vec[3]))
  
  
  
  ###############################
  # computing odds ratio per SD #
  ###############################
  or_sd_vad_reg_df <- reg_validate_df
  or_sd_vad_reg_df$SCORE <- or_sd_vad_reg_df[,current_approach]
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

  # plotting histogram of PRS scores
  current_reg_validate_df$Legend <- NA
  current_reg_validate_df$Legend[current_reg_validate_df$Status==1] <- "Case"
  current_reg_validate_df$Legend[current_reg_validate_df$Status==0] <- "Control"
  png(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_AUC/DENSITY/",subtype,"/DENSITY_PRS_",current_approach,".png"),units="in", width=5, height=5, res=2000)
  print(ggplot(current_reg_validate_df, aes(x=SCORE, fill=Legend)) + geom_density(alpha=0.4) + xlab("PRS Score") + ylab("Density")) + theme_classic()
  dev.off()
}

# assembling final data.frame of AUCs 
AUC_DF <- data.frame(
  approach_name,
  auc_vec_lower,
  auc_vec_center,
  auc_vec_upper,
  logOR_sd_CI_95_lower_sd,
  logOR_sd,
  logOR_sd_CI_95_upper_sd
)
colnames(AUC_DF) <- c("Approach","AUC_lower","AUC_center","AUC_upper","logOR_lower","logOR_center","logOR_upper")

# make OR columns
AUC_DF$OR_lower <- exp(AUC_DF$logOR_lower)
AUC_DF$OR_center <- exp(AUC_DF$logOR_center)
AUC_DF$OR_upper <- exp(AUC_DF$logOR_upper)

# make OR columns exactly two decimal places
OR_cols <- c("OR_lower", "OR_center", "OR_upper")
AUC_DF[OR_cols] <- lapply(AUC_DF[OR_cols], function(x) ifelse(is.na(x), NA, sprintf("%.2f", round(x, 2))))

# making a combined OR column
AUC_DF <- AUC_DF %>% mutate(OR_combined=paste0(OR_center,"\n","(",OR_lower,", ",OR_upper,")"))

# adding a column for the number of SNPs in each PRS model
AUC_DF$Num_SNPs <- NA
for (k in 1:length(PRS_APPROACH_VEC)) {
  print(k)
  current_approach <- PRS_APPROACH_VEC[k]
  AUC_DF$Num_SNPs[k] <- nrow(fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_PRS_MODELS/",subtype,"_",current_approach),sep="\t",header=F))
}

# outputting the AUC DF
print(AUC_DF %>% arrange(desc(AUC_center)))
save(AUC_DF,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_AUC/AUC/",subtype,"/AUC_DF.RData"))

# polishing the AUC DF
FINAL_AUC_DF <- AUC_DF
FINAL_AUC_DF$AUC_lower <- signif(as.numeric(FINAL_AUC_DF$AUC_lower), digits = 3)
FINAL_AUC_DF$AUC_center <- signif(as.numeric(FINAL_AUC_DF$AUC_center), digits = 3)
FINAL_AUC_DF$AUC_upper <- signif(as.numeric(FINAL_AUC_DF$AUC_upper), digits = 3)
FINAL_AUC_DF$Approach <- gsub("e","e-",gsub(paste0(subtype,"_"),"",gsub(".sscore","",FINAL_AUC_DF$Approach)))
# FINAL_AUC_DF %>% arrange(AUC_center)
FINAL_AUC_DF <- FINAL_AUC_DF %>% separate(Approach,sep="_",into=c("Class","Method","Included SNPs"),remove=F) %>% mutate(Class = ifelse(Class=="SINGLEANCESTRY","Single-ancestry",ifelse(Class=="XANCESTRY","Cross-ancestry",ifelse(Class=="XSUBTYPE","Cross-subtype","Ensemble"))))
FINAL_AUC_DF$Class <- factor(FINAL_AUC_DF$Class,levels=c("Single-ancestry","Cross-ancestry","Cross-subtype"))
save(FINAL_AUC_DF,file=paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/",subtype,"_","FINAL_AUC_DF.RData"))

# setting the error bar label distance based on subtype since standard errors vary
vjust_subtype<-0
if (subtype=="OVERALL") {
  vjust_subtype <- -3.25
} else if (subtype=="ERPOS") {
  vjust_subtype <- -3.75
} else if (subtype=="ERNEG") {
  vjust_subtype <- -4.25
} else if (subtype=="TNBC") {
  vjust_subtype <- -5
} else {
  print("vjust_subtype NOT DEFINED")
}

######################################################
######################################################
######################################################
# parsing results for SINGLEANCESTRY PRS
PARSED_DF <- FINAL_AUC_DF %>% mutate(`Included SNPs`=ifelse(`Included SNPs`=="hm3",0,ifelse(`Included SNPs`=="CLUMPED",0.3,`Included SNPs`))) %>% mutate(`Included SNPs` = as.numeric(`Included SNPs`)) %>% filter(grepl("SINGLEANCESTRY",Approach),!grepl("ENSEMBLE",Approach),!grepl("CTSLEB",Approach),!grepl("PRScsx",Approach))
PARSED_DF <- PARSED_DF %>% arrange(Class,Method,`Included SNPs`) %>% mutate(`Included SNPs`=ifelse(`Included SNPs`==0,"HapMap3",ifelse(`Included SNPs`==1,"Clumped",`Included SNPs`)))
PARSED_DF$Approach <- factor(PARSED_DF$Approach, levels = PARSED_DF$Approach)
if (subtype %in% c("OVERALL","ERNEG","ERPOS")) {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
} else if (subtype=="TNBC") {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
}

# renaming PRSice-2
PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="PRSice-2","PLINK C+T",Method))

# finalizing the Included SNPs column
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`)] <- "All variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")] <- paste0("p<",PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")] )

# remove EUR PRS from the plot of single ancestry PRS models
PARSED_DF <- PARSED_DF %>% filter(!grepl("EUR",Approach))

# obtain labels for the different PRS methods
UNIQUE_METHOD_LIST <- unique(PARSED_DF$Method)
METHOD_PLOT_INDEX <- c()
for (UNIQUE_METHOD in UNIQUE_METHOD_LIST) {
  METHOD_PLOT_INDEX <- c(METHOD_PLOT_INDEX,median(which(PARSED_DF$Method==UNIQUE_METHOD)))
}

png(paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/",subtype,"_","SINGLEANCESTRY.png"),res=1200,units="in",height=4,width=11)
ggplot(data=PARSED_DF,aes(x=Approach,y=AUC_center,col=Method)) + geom_point() + geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper),width=0.8) + geom_text(aes(label=sprintf("%.3f",AUC_center),vjust = vjust_subtype),size=3) + coord_cartesian(ylim = c(0.5, 0.7), xlim = c(0,nrow(PARSED_DF)+1), expand = FALSE, clip = "off") + theme_classic() + theme(plot.margin = unit(c(1, 8.5, 6.5, 1), "lines"), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") + 
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.48, label=prettyNum(PARSED_DF$Num_SNPs, big.mark = ",", scientific = FALSE),size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.48, label="Number of SNPs in PRS",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.46, label=PARSED_DF$`Included SNPs`,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.46, label="SNPs Included",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.44, label=PARSED_DF$OR_combined,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.44, label="OR per SD",size=3) +
  
  annotate(
    geom = "text", 
    x = METHOD_PLOT_INDEX, 
    y = ifelse(UNIQUE_METHOD_LIST %in% c("Lassosum2", "LDpred2"), 0.412,0.415),
    label=UNIQUE_METHOD_LIST,size=3,
    angle = ifelse(UNIQUE_METHOD_LIST %in% c("Lassosum2", "LDpred2"), 30, 0)) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.415, label="PRS method",size=3) +
  
  ylab("Covariate-adjusted AUC")
dev.off()


######################################################
######################################################
######################################################
# parsing results for XANCESTRY PRS
PARSED_DF <- FINAL_AUC_DF %>% mutate(`Included SNPs`=ifelse(`Included SNPs`=="hm3",0,ifelse(`Included SNPs`=="CLUMPED",0.3,`Included SNPs`))) %>% mutate(`Included SNPs` = as.numeric(`Included SNPs`)) %>% filter((grepl("XANCESTRY",Approach) | grepl("CTSLEB",Approach) | grepl("PRScsx",Approach)),!grepl("ENSEMBLE",Approach),!grepl("XSUBTYPE",Approach))
# temporarily changing the name of the PRScsx method to sort correctly
PARSED_DF$Method[PARSED_DF$Method=="PRScsx"] <- "APRScsx"
PARSED_DF <- PARSED_DF %>% arrange(Method,Class,`Included SNPs`) %>% mutate(`Included SNPs`=ifelse(`Included SNPs`==0,"HapMap3",ifelse(`Included SNPs`==1,"Clumped",`Included SNPs`)))
# changing the name of the PRScsx method back after sorting above
PARSED_DF$Method[PARSED_DF$Method=="APRScsx"] <- "PRScsx"
PARSED_DF$Approach <- factor(PARSED_DF$Approach, levels = PARSED_DF$Approach)
if (subtype %in% c("OVERALL","ERNEG","ERPOS")) {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
} else if (subtype=="TNBC") {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
}

# modifying the Included SNPs column for PRSice-2
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & PARSED_DF$Method=="PRSice-2"] <- "All variants"

# renaming PRSice-2
PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="PRSice-2","PLINK C+T",Method))

# modifying the Included SNPs column (continued)
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & grepl("CTSLEB",PARSED_DF$Approach)] <- "0.3"
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & grepl("PRScsx",PARSED_DF$Approach)] <- "HapMap3"
# finalizing the Method column
PARSED_DF$Method[grepl("SINGLEANCESTRY",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "CT-SLEB"
PARSED_DF$Method[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "CT-SLEB"
PARSED_DF$Method[grepl("SINGLEANCESTRY",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "PRS-CSx"
PARSED_DF$Method[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "PRS-CSx"
# finalizing the Included SNPs column
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")] <- paste0("p<",PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")])

if (subtype %in% c("OVERALL","ERNEG","ERPOS")) {
  PARSED_DF$`Included SNPs`[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "p<0.3+EUR330"
  PARSED_DF$`Included SNPs`[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "HapMap3+EUR330"
} else if (subtype=="TNBC") {
  PARSED_DF$`Included SNPs`[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "p<0.3+EUR330"
  PARSED_DF$`Included SNPs`[grepl("XANCESTRY",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "HapMap3+EUR330"
}

# obtain labels for the different PRS methods
UNIQUE_METHOD_LIST <- unique(PARSED_DF$Method)
METHOD_PLOT_INDEX <- c()
for (UNIQUE_METHOD in UNIQUE_METHOD_LIST) {
  METHOD_PLOT_INDEX <- c(METHOD_PLOT_INDEX,median(which(PARSED_DF$Method==UNIQUE_METHOD)))
}

png(paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/",subtype,"_","XANCESTRY.png"),res=1200,units="in",height=4,width=11)
ggplot(data=PARSED_DF,aes(x=Approach,y=AUC_center,col=Method)) + geom_point() + geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper),width=0.8) + geom_text(aes(label=sprintf("%.3f",AUC_center),vjust = vjust_subtype),size=3) + coord_cartesian(ylim = c(0.5, 0.7), xlim = c(0,nrow(PARSED_DF)+1), expand = FALSE, clip = "off") + theme_classic() + theme(plot.margin = unit(c(1, 8.5, 6.5, 1), "lines"), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") + 
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.48, label=prettyNum(PARSED_DF$Num_SNPs, big.mark = ",", scientific = FALSE),size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.48, label="Number of SNPs in PRS",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.46, label=PARSED_DF$`Included SNPs`,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.46, label="SNPs Included",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.44, label=PARSED_DF$OR_combined,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.44, label="OR per SD",size=3) +
  
  annotate(geom = "text", x = METHOD_PLOT_INDEX, y = ifelse(UNIQUE_METHOD_LIST %in% c("Lassosum2", "LDpred2"), 0.412,0.415), label=UNIQUE_METHOD_LIST,size=3,angle = ifelse(UNIQUE_METHOD_LIST %in% c("Lassosum2", "LDpred2"), 30, 0)) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.415, label="PRS method",size=3) +
  
  ylab("Covariate-adjusted AUC")
dev.off()




######################################################
######################################################
######################################################
# parsing results for XSUBTYPE PRS
PARSED_DF <- FINAL_AUC_DF %>% mutate(`Included SNPs`=ifelse(`Included SNPs`=="hm3",0,ifelse(`Included SNPs`=="CLUMPED",0.3,`Included SNPs`))) %>% mutate(`Included SNPs` = as.numeric(`Included SNPs`)) %>% filter(grepl("XSUBTYPE",Approach),!grepl("ENSEMBLE",Approach))
# temporarily changing the name of the PRScsx method to sort correctly
PARSED_DF$Method[PARSED_DF$Method=="PRScsx"] <- "APRScsx"
PARSED_DF <- PARSED_DF %>% arrange(Method,Class,`Included SNPs`) %>% mutate(`Included SNPs`=ifelse(`Included SNPs`==0,"HapMap3",ifelse(`Included SNPs`==1,"Clumped",`Included SNPs`)))
# changing the name of the PRScsx method bak after sorting
PARSED_DF$Method[PARSED_DF$Method=="APRScsx"] <- "PRScsx"
PARSED_DF$Approach <- factor(PARSED_DF$Approach, levels = PARSED_DF$Approach)
if (subtype %in% c("OVERALL","ERNEG","ERPOS")) {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
} else if (subtype=="TNBC") {
  PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="EUR","EUR 330-variant",ifelse(Method=="FSS","Forward Stepwise Selection",ifelse(Method=="LASSOSUM2","Lassosum2",ifelse(Method=="LDPRED2AUTO","LDpred2",ifelse(Method=="PLINKCT","PLINK C+T",ifelse(Method=="PRSICE2","PRSice-2","")))))))
}

# finalizing the Included SNPs column for PRSice-2
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & PARSED_DF$Method=="PRSice-2"] <- "All variants"

# renaming PRSice-2
PARSED_DF <- PARSED_DF %>% mutate(Method=ifelse(Method=="PRSice-2","PLINK C+T",Method))

# finalizing the Included SNPs column (continued)
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & PARSED_DF$Method=="PRScsx"] <- "HapMap3"
PARSED_DF$`Included SNPs`[is.na(PARSED_DF$`Included SNPs`) & grepl("CTSLEB",PARSED_DF$Approach)] <- "0.3"
PARSED_DF$Method[grepl("XSUBTYPE",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "CT-SLEB"
PARSED_DF$Method[grepl("XSUBTYPE",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "PRS-CSx"


PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs`==330 | PARSED_DF$`Included SNPs`==313] <- "Selected variants"
PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")] <- paste0("p<",PARSED_DF$`Included SNPs`[PARSED_DF$`Included SNPs` %in% c("1e-08","1e-07","1e-06","1e-05","1e-04","0.001","0.01","0.05","0.1","0.2","0.3","0.4","0.5")])
PARSED_DF$`Included SNPs`[grepl("XSUBTYPE",PARSED_DF$Approach) & grepl("CTSLEB",PARSED_DF$Approach)] <- "p<0.3"
PARSED_DF$`Included SNPs`[grepl("XSUBTYPE",PARSED_DF$Approach) & grepl("PRScsx",PARSED_DF$Approach)] <- "HapMap3"


# obtain labels for the different PRS methods
UNIQUE_METHOD_LIST <- unique(PARSED_DF$Method)
METHOD_PLOT_INDEX <- c()
for (UNIQUE_METHOD in UNIQUE_METHOD_LIST) {
  METHOD_PLOT_INDEX <- c(METHOD_PLOT_INDEX,median(which(PARSED_DF$Method==UNIQUE_METHOD)))
}

png(paste0("../../CALIBRATION_FINAL_SCORE_AUC/AUC/",subtype,"_","XSUBTYPE.png"),res=1200,units="in",height=4,width=11)
ggplot(data=PARSED_DF,aes(x=Approach,y=AUC_center,col=Method)) + geom_point() + geom_errorbar(aes(ymin = AUC_lower, ymax = AUC_upper),width=0.8) + geom_text(aes(label=sprintf("%.3f",AUC_center),vjust = vjust_subtype),size=3) + coord_cartesian(ylim = c(0.5, 0.7), xlim = c(0,nrow(PARSED_DF)+1), expand = FALSE, clip = "off") + theme_classic() + theme(plot.margin = unit(c(1, 8.5, 6.5, 1), "lines"), axis.title.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), legend.position = "none") + 
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.48, label=prettyNum(PARSED_DF$Num_SNPs, big.mark = ",", scientific = FALSE),size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.48, label="Number of SNPs in PRS",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.46, label=PARSED_DF$`Included SNPs`,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.46, label="SNPs Included",size=3) +
  
  annotate(geom = "text", x = seq_len(nrow(PARSED_DF)), y=0.44, label=PARSED_DF$OR_combined,size=1.75,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.44, label="OR per SD",size=3) +
  
  annotate(geom = "text", x = METHOD_PLOT_INDEX[c(1:2)], y=0.415, label=UNIQUE_METHOD_LIST[c(1:2)],size=3,angle=30) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.415, label="PRS method",size=3) +
  annotate(geom = "text", x = METHOD_PLOT_INDEX[-c(1:2)], y=0.415, label=UNIQUE_METHOD_LIST[-c(1:2)],size=3) + 
  annotate(geom = "text", x = nrow(PARSED_DF)+3, y=0.415, label="PRS method",size=3) +
  
  ylab("Covariate-adjusted AUC")
dev.off()

