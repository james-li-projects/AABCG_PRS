library(dplyr)
library(data.table)
library(readxl)
library(ROCnReg)
library(ggplot2)

# set seed
set.seed(1)

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


#########################################
# importing scoring df for ERPOS 
SCORE_FILE_1 <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/ERPOS_ENSEMBLE_FSS.sscore") 
NORM_FACTOR <- median(SCORE_FILE_1$ALLELE_CT)
SCORE_FILE_1$SCORE <- SCORE_FILE_1$SCORE1_AVG*NORM_FACTOR
SCORE_FILE_1 <- SCORE_FILE_1 %>% select(`#IID`,SCORE)
colnames(SCORE_FILE_1) <- c("#IID","SCORE_ERPOS")
# importing scoring df for ERNEG 
SCORE_FILE_2 <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/ERNEG_ENSEMBLE_GLMNET.sscore") 
NORM_FACTOR <- median(SCORE_FILE_2$ALLELE_CT)
SCORE_FILE_2$SCORE <- SCORE_FILE_2$SCORE1_AVG*NORM_FACTOR
SCORE_FILE_2 <- SCORE_FILE_2 %>% select(`#IID`,SCORE)
colnames(SCORE_FILE_2) <- c("#IID","SCORE_ERNEG")
# importing scoring df for TNBC 
SCORE_FILE_3 <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/TNBC_XANCESTRY_PRSICE2.sscore") 
NORM_FACTOR <- median(SCORE_FILE_3$ALLELE_CT)
SCORE_FILE_3$SCORE <- SCORE_FILE_3$SCORE1_AVG*NORM_FACTOR
SCORE_FILE_3 <- SCORE_FILE_3 %>% select(`#IID`,SCORE)
colnames(SCORE_FILE_3) <- c("#IID","SCORE_TNBC")
# importing scoring df for OVERALL 
SCORE_FILE_4 <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/OVERALL_XANCESTRY_PRScsx.sscore") 
NORM_FACTOR <- median(SCORE_FILE_4$ALLELE_CT)
SCORE_FILE_4$SCORE <- SCORE_FILE_4$SCORE1_AVG*NORM_FACTOR
SCORE_FILE_4 <- SCORE_FILE_4 %>% select(`#IID`,SCORE)
colnames(SCORE_FILE_4) <- c("#IID","SCORE_OVERALL")

# combining score files
SCORE_FILE <- inner_join(SCORE_FILE_4,inner_join(SCORE_FILE_3,inner_join(SCORE_FILE_1,SCORE_FILE_2,by=c("#IID")),by=c("#IID")),by=c("#IID"))

# importing validation covariate data
validation.pheno_cov <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_OVERALL.pheno_cov",header=T)
validation.pheno_cov$Status <- validation.pheno_cov$Status-1

# importing testing covariate data
testing.pheno_cov <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/testing_OVERALL.pheno_cov",header=T)
testing.pheno_cov$Status <- testing.pheno_cov$Status-1

# importing all pheno data to identify proportions of subtypes in validation set
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx")) %>% filter(`AABCGS_ID...1` %in% validation.pheno_cov$`#IID`)
# identifying proportions of ER status in validation set
NUM_ERPOS <- nrow(pheno_data %>% filter(ER==1))
NUM_ERNEG <- nrow(pheno_data %>% filter(ER==2))
NUM_TNBC <- nrow(pheno_data %>% filter(ER==2) %>% filter(PR==2) %>% filter(HER2==2))
NUM_ERNEG <- NUM_ERNEG - NUM_TNBC
PROP_ERPOS <- NUM_ERPOS / (NUM_ERPOS + NUM_ERNEG + NUM_TNBC)
PROP_ERNEG <- NUM_ERNEG / (NUM_ERPOS + NUM_ERNEG + NUM_TNBC)
PROP_TNBC <- NUM_TNBC / (NUM_ERPOS + NUM_ERNEG + NUM_TNBC)
PROP_ERPOS
PROP_ERNEG
PROP_TNBC
# joining validation data with phenotype and PRS scoring data
reg_vad_df <- inner_join(validation.pheno_cov,SCORE_FILE,by=c("#IID"))
reg_vad_df$ER_TUNED_PRS <- 
  PROP_ERPOS*reg_vad_df$SCORE_ERPOS +
  PROP_ERNEG*reg_vad_df$SCORE_ERNEG +
  PROP_TNBC*reg_vad_df$SCORE_TNBC

# computing AUC for 3-subtype PRS
output_AROC.sp <- AROC.sp(
  formula.h = ER_TUNED_PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
  group = "Status",
  data = reg_vad_df,
  tag.h = 0)
print(output_AROC.sp$AUC)
three_subtype_AUC <- output_AROC.sp$AUC

# computing AUC for 3-subtype PRS in individuals with no missing receptor data
output_AROC.sp <- AROC.sp(
  formula.h = ER_TUNED_PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
  group = "Status",
  data = reg_vad_df %>% filter(`#IID` %in% eligible_samples),
  tag.h = 0)
print(output_AROC.sp$AUC)
three_subtype_AUC_noMiss <- output_AROC.sp$AUC


#########################################
# creating proportion calibration plot
#########################################
# joining testing data with phenotype and PRS scoring data
plot_df <- data.frame()
for (TMP_PROP in c(0,0.1,0.2,0.3,0.4,0.5,0.6,PROP_ERPOS,0.7,0.8,0.9,1)) {
  print(TMP_PROP)
  reg_vad_df <- inner_join(validation.pheno_cov,SCORE_FILE,by=c("#IID"))
  reg_vad_df$ER_TUNED_PRS <- 
    TMP_PROP*reg_vad_df$SCORE_ERPOS +
    (1-TMP_PROP)*reg_vad_df$SCORE_ERNEG
  # computing AUC
  output_AROC.sp <- AROC.sp(
    formula.h = ER_TUNED_PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
    group = "Status",
    data = reg_vad_df,
    tag.h = 0)
  print(output_AROC.sp$AUC)
  tmp_plot_df <- data.frame(cbind(Proportion_ERpos=TMP_PROP,t(output_AROC.sp$AUC)),group=NA)
  plot_df <- rbind(plot_df,tmp_plot_df)
}

# labeling internal proportion group differently
# plot_df <- plot_df %>% mutate(group=ifelse(Proportion_ERpos==PROP_ERPOS,"Two Subtype PRS (All Individuals)",group))

# adding in AUCs for 3-subtype PRS evaluations 
plot_df <- rbind(
  plot_df,
  cbind(Proportion_ERpos=PROP_ERPOS,t(three_subtype_AUC),group="Three Subtype PRS (All Individuals)")
)
plot_df <- rbind(
  plot_df,
  cbind(Proportion_ERpos=PROP_ERPOS,t(three_subtype_AUC_noMiss),group="Three Subtype PRS (Complete Subtype Data)")
)

# making the data numeric for the plot_df
plot_df <- plot_df %>% mutate(
  Proportion_ERpos=as.numeric(Proportion_ERpos),
  est=as.numeric(est),
  ql=as.numeric(ql),
  qh=as.numeric(qh)
)

# Load required libraries
library(ggplot2)
library(dplyr)

# Ensure 'group' is a factor with correct labels
plot_df <- plot_df %>%
  mutate(
    group = factor(group, levels = c(
      "Three Subtype PRS (All Individuals)",
      "Three Subtype PRS (Complete Subtype Data)"
    ))
  )

# Define custom colors for groups
group_colors <- c(
  "Three Subtype PRS (All Individuals)" = "#00008b",
  "Three Subtype PRS (Complete Subtype Data)" = "#40e0d0"
)

# Define proportion numerical labels every 0.1
x_labels <- seq(0, 1, by = 0.1) 

# Define the plot
p <- ggplot(plot_df, aes(x = Proportion_ERpos, y = est, ymin = ql, ymax = qh, color = group)) +
  # Error bars: colored by group, with NA values remaining black
  geom_errorbar(aes(color = group), width = 0.05, show.legend = FALSE) +  
  # Points: same colors as error bars
  geom_point(size = 3) +  
  scale_color_manual(
    values = c(group_colors, "gray"),  # Include black for NA values
    na.value = "gray",  # Assign black for NA values
    breaks = names(group_colors),  # Ensures only non-NA values appear in the legend
    labels = c("Three Subtype PRS (All Individuals)", "Three Subtype PRS (Complete Subtype Data)")  # Clean legend labels
  ) +
  labs(
    #title = "Performance of PRS tuned by proportion of cases by ER status in the internal validation set",
    x = "Assigned proportion of cases with ER-positive breast cancer",
    y = "Covariate-adjusted AUC"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),  # Removes legend title
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),  # **Removes all major grid lines**
    panel.grid.minor.x = element_blank(),  # **Removes vertical minor grid lines**
    panel.grid.minor.y = element_line(color = "#e6e6e6", linetype = "dashed", size = 0.3) 
  ) +
  scale_x_continuous(
    breaks = x_labels 
  ) 

# Save the plot
ggsave(
  filename = "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/PROP_ERROR_PLOT.png",
  plot = p,
  width = 10, height = 5, dpi = 1200
)


#########################################
# generating scoring file for each tuning proportion
#########################################
for (row_index in 1:nrow(plot_df)) {
  # assembling an input data.frame called cvfit_coef with column names of Name and Effect for each proportion
  current_prop_ERpos <- plot_df$Proportion_ERpos[row_index]
  cvfit_coef <- data.frame(
    Name=c("ERPOS_ENSEMBLE_FSS.sscore","ERNEG_ENSEMBLE_GLMNET.sscore"),
    Effect=c(current_prop_ERpos,1-current_prop_ERpos)
  )
  
  #########################################
  # obtaining a final scoring file  
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
    TMP_SCORE_FILE_STR <- gsub(".sscore","",cvfit_coef$Name)[i]
    
    # obtaining subtype string
    subtype <- sub("_.*", "", TMP_SCORE_FILE_STR)
    
    # obtaining weight of the current PRS model and multiplying it for all coefficients in the scoring file
    TMP_SCORE_WEIGHT <- cvfit_coef$Effect[i]
    TMP_SCORE_FILE <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/FINAL_PRS_MODELS/",TMP_SCORE_FILE_STR)) 
    
    # account for PRS score normalization factor (allele counts)
    TMP_NORM_FACTOR <- median((fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/",TMP_SCORE_FILE_STR,".sscore")))$ALLELE_CT)
    TMP_SCORE_FILE$V4 <- (TMP_SCORE_FILE$V3*TMP_SCORE_WEIGHT)#/TMP_NORM_FACTOR
    
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
  
  # obtaining string for final scoring file
  
  if (toString(current_prop_ERpos) == toString(PROP_ERPOS)) {
    current_prop_str = "InternalProp"
  } else {
    current_prop_str = toString(current_prop_ERpos)
  }
  
  # writing out the final scoring file
  write.table(COLLATED_SCORE_FILE_DF,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/FINAL_PRS_MODELS/OVERALL_PROPTUNE_",current_prop_str),quote=F,row.names=F,col.names=F,sep="\t")
}


#########################################
# generating scoring file incorporating all three subtypes
#########################################
# assembling an input data.frame called cvfit_coef with column names of Name and Effect for each proportion
cvfit_coef <- data.frame(
  Name=c("ERPOS_ENSEMBLE_FSS.sscore","ERNEG_ENSEMBLE_GLMNET.sscore","TNBC_XANCESTRY_PRSICE2.sscore"),
  Effect=c(PROP_ERPOS,PROP_ERNEG,PROP_TNBC)
)

#########################################
# obtaining a final scoring file  
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
  TMP_SCORE_FILE_STR <- gsub(".sscore","",cvfit_coef$Name)[i]
  
  # obtaining subtype string
  subtype <- sub("_.*", "", TMP_SCORE_FILE_STR)
  
  # obtaining weight of the current PRS model and multiplying it for all coefficients in the scoring file
  TMP_SCORE_WEIGHT <- cvfit_coef$Effect[i]
  TMP_SCORE_FILE <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/FINAL_PRS_MODELS/",TMP_SCORE_FILE_STR)) 
  print(nrow(TMP_SCORE_FILE))
  
  # account for PRS score normalization factor (allele counts)
  TMP_NORM_FACTOR <- median((fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PROPORTION_TUNED_PRS/",TMP_SCORE_FILE_STR,".sscore")))$ALLELE_CT)
  TMP_SCORE_FILE$V4 <- (TMP_SCORE_FILE$V3*TMP_SCORE_WEIGHT)#/TMP_NORM_FACTOR
  
  # including the weighted scoring file to existing scoring file records 
  SCORE_FILE_DF <- rbind(SCORE_FILE_DF,TMP_SCORE_FILE)
}
colnames(SCORE_FILE_DF) <- c("V1","V2","V3","V4")
print(nrow(SCORE_FILE_DF))

# collating all the PRS weights
COLLATED_SCORE_FILE_DF <- data.table(SCORE_FILE_DF %>% select(-V3) %>% dplyr::group_by(V1,V2) %>% summarise(Effect = sum(V4))) 
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% tidyr::separate(V1, into = c("chr","pos","a2","a1"),remove=F,sep=":")
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% mutate(AlignedEffect = ifelse(a1 == V2, Effect, -Effect))
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% select(V1,a1,AlignedEffect) 
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% dplyr::group_by(V1,a1) %>% summarise(FinalEffect = sum(AlignedEffect))
COLLATED_SCORE_FILE_DF <- data.frame(COLLATED_SCORE_FILE_DF)
COLLATED_SCORE_FILE_DF <- COLLATED_SCORE_FILE_DF %>% filter(FinalEffect != 0)

# writing out the final scoring file
write.table(COLLATED_SCORE_FILE_DF,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/FINAL_PRS_MODELS/OVERALL_PROPTUNE_ALLSUBTYPE"),quote=F,row.names=F,col.names=F,sep="\t")
