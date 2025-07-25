###########################
# loading packages and setting working directory
library(data.table)
library(tidyverse)
library(readxl)
setwd("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA")

# reading in seed number
args = commandArgs(trailingOnly=TRUE)
seed_num <- as.numeric(args[1])
print(paste("Seed Number:",seed_num))
set.seed(seed_num)

###########################
# assessing eigenvalues generated from PCA to determine the number of PCs to include 
eigenval_vec <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/AABCG_PCA_50PCs.eigenval",header=F)$V1
eigenval_df <- data.frame(eigenval_index = c(1:length(eigenval_vec)),eigenval = eigenval_vec)
p <- ggplot(data = eigenval_df, aes(x=eigenval_index,y=eigenval)) + geom_point() + theme_classic()
png("eigenval_plot.png")
p
dev.off()

###########################
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 3-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
colnames(pheno_data)[2] <- "Dataset"

# importing data of principal components
eigenvec <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/AABCG_PCA_50PCs.eigenvec",header=T)
colnames(eigenvec)[1] <- "Sample_Name"

# joining covariate and PC data.frames by sample IID
combined_cov <- inner_join(pheno_data,eigenvec, by = c("Sample_Name")) 

# recoding the dataset variable 
combined_cov$dataset_recode <- NA
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_01_WGS"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_02_WGS2"] <- "WGS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_03_MEGA_VANDY"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_04_MEGA_RP"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_05_MEGA_USC"] <- "MEGA" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_06_AMBER"] <- "AMBER" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_07_ROOT"] <- "ROOT" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_08_AABC"] <- "AABC" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_09_GBHS"] <- "GBHS" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_10_BCAC_OncoArray"] <- "BCAC_OncoArray" 
combined_cov$dataset_recode[combined_cov$Dataset=="Dataset_11_BCAC_iCOGS"] <- "BCAC_iCOGS" 

###########################
# assessing the relationship between global AFR ancestry proportion and the top 10 principal components
correlation_vec <- c(
    cor(combined_cov$AFR_pro, combined_cov$PC1),
    cor(combined_cov$AFR_pro, combined_cov$PC2),
    cor(combined_cov$AFR_pro, combined_cov$PC3),
    cor(combined_cov$AFR_pro, combined_cov$PC4),
    cor(combined_cov$AFR_pro, combined_cov$PC5),
    cor(combined_cov$AFR_pro, combined_cov$PC6),
    cor(combined_cov$AFR_pro, combined_cov$PC7),
    cor(combined_cov$AFR_pro, combined_cov$PC8),
    cor(combined_cov$AFR_pro, combined_cov$PC9),
    cor(combined_cov$AFR_pro, combined_cov$PC10)
)
cor_AFR_PC_df <- data.frame(PC_index = c(1:10), correlation = correlation_vec)

###########################
# assessing the relationship between PCs and datasets 
for (i in 1:10) {
    for (j in 1:10) {
        if (j > i) {
            PC_i <- paste0("PC",i)
            PC_j <- paste0("PC",j)
            x = combined_cov[,PC_i]
            y = combined_cov[,PC_j]
            Dataset = combined_cov$Dataset
            plot_df <- data.frame(x,y,Dataset)
            p <- ggplot(data = plot_df, aes(x=x,y=y,color=Dataset)) + geom_point() + xlab(PC_i) + ylab(PC_j) + theme_classic() 

            png(paste0("correlation_PC_dataset_",PC_i,"_",PC_j,".png"))
            print(p)
            dev.off()
            
        } else {
        # do nothing
        }
    } 
}

# regressing case control status on PCs and datasets 
summary(glm(family="binomial",data = combined_cov, formula = as.formula(gsub(",","",toString(c("as.factor(Status) ~ ",toString(paste(paste0(rep("PC",50),c(1:50)),rep("+",50))),"dataset_recode"))))))

# regressing PCs on datasets 
p_val_mat <- data.frame(summary(lm(data = combined_cov, formula = PC1 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC2 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC3 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC4 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC5 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC6 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC7 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC8 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC9 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC10 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC11 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC12 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC13 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC14 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC15 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC16 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC17 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC18 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC19 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC20 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC21 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC22 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC23 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC24 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC25 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC26 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC27 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC28 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC29 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC30 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC31 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC32 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC33 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC34 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC35 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC36 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC37 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC38 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC39 ~ Dataset))$coefficients[,4])
p_val_mat <- cbind(p_val_mat, summary(lm(data = combined_cov, formula = PC40 ~ Dataset))$coefficients[,4])


p_val_mat <- p_val_mat[-1,] 
colnames(p_val_mat) <- paste0(rep("PC",40),1:40)

################################
### creating covariate files ###
################################

# organizing the data to include relevant columns
organize_combined_cov <- combined_cov %>% select(
    Sample_Name,
    Status,
    Age_GWAS,
    dataset_recode,
    PC1,
    PC2,
    PC3,
    PC4,
    PC5,
    PC6,
    PC7,
    PC8,
    PC9,
    PC10,
    PC11,
    PC12,
    PC13,
    PC14,
    PC15,
    PC16,
    PC17,
    PC18,
    PC19,
    PC20,
    PC21,
    PC22,
    PC23,
    PC24,
    PC25,
    PC26,
    PC27,
    PC28,
    PC29,
    PC30,
    PC31,
    PC32,
    PC33,
    PC34,
    PC35,
    PC36,
    PC37,
    PC38,
    PC39,
    PC40
)
colnames(organize_combined_cov)[1:4] <- c("IID","Status","Age","Platform")

# transforming age into a z-score
stand_z_func <- function(x) {
  return( (x-mean(x)) / sd(x) )
}
organize_combined_cov$Age <- stand_z_func(organize_combined_cov$Age)


# randomly subsampling training/testing/validation sets indices separately from cases and controls
organize_cov_case <- organize_combined_cov %>% filter(Status == 2)
organize_cov_control <- organize_combined_cov %>% filter(Status == 1)
# obtaining sample indices for case and controls for each training/testing/validation sets
picked_case_training = sample(seq_len(nrow(organize_cov_case)),size = round(nrow(organize_cov_case)*.7))
picked_control_training = sample(seq_len(nrow(organize_cov_control)),size = round(nrow(organize_cov_control)*.7))
testing_validation_case = setdiff(seq_len(nrow(organize_cov_case)), picked_case_training)
testing_validation_control = setdiff(seq_len(nrow(organize_cov_control)), picked_control_training)
picked_case_testing = sample(testing_validation_case, size = round(nrow(organize_cov_case)*.1))
picked_control_testing = sample(testing_validation_control, size = round(nrow(organize_cov_control)*.1))
picked_case_validation = setdiff(testing_validation_case, picked_case_testing)
picked_control_validation = setdiff(testing_validation_control, picked_control_testing)

# generating subsampled training/testing/validation sets 
cov_case_training <- organize_cov_case[picked_case_training,]
cov_case_testing <- organize_cov_case[picked_case_testing,]
cov_case_validation <- organize_cov_case[picked_case_validation,]
cov_control_training <- organize_cov_control[picked_control_training,]
cov_control_testing <- organize_cov_control[picked_control_testing,]
cov_control_validation <- organize_cov_control[picked_control_validation,]
# assembling sets
training_combined <- rbind(cov_case_training,cov_control_training)
testing_combined <- rbind(cov_case_testing,cov_control_testing)
validation_combined <- rbind(cov_case_validation,cov_control_validation)
colnames(training_combined)[1] <- "#IID"
colnames(testing_combined)[1] <- "#IID"
colnames(validation_combined)[1] <- "#IID"

###################
# OVERALL BC RISK # 
###################
# outputting covariate and phenotype files
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined,file="training_OVERALL.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(testing_combined,file="testing_OVERALL.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(validation_combined,file="validation_OVERALL.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")

# outputting sample lists
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% select(`#IID`),file="training_OVERALL.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(testing_combined %>% select(`#IID`),file="testing_OVERALL.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(validation_combined %>% select(`#IID`),file="validation_OVERALL.samplist",quote=F,row.names=F,col.names=F,sep=" ")

##############################################
# GENERATING SAMPLE LISTS FOR OTHER SUBTYPES #
##############################################
# re-reading in the phenotype data to completely avoid any sort of confusion
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
iCOGS_SAMPLE_LIST <- (pheno_data %>% filter(Platform=="iCOGS"))$`AABCGS_ID...1`
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,Status,ER,PR,HER2)
colnames(pheno_data)[1] <- "#IID"

# identifying subtype sample lists
MISSING_UNKNOWN_ER <- (pheno_data %>% filter(ER==9,Status==1))$`#IID`
ERNEG_EXCLUDE_SAMPLE_LIST <- (pheno_data %>% filter(ER==1,Status==1))$`#IID`
ERNEG_EXCLUDE_SAMPLE_LIST <- c(ERNEG_EXCLUDE_SAMPLE_LIST, MISSING_UNKNOWN_ER)
ERPOS_EXCLUDE_SAMPLE_LIST <- (pheno_data %>% filter(ER==2,Status==1))$`#IID`
ERPOS_EXCLUDE_SAMPLE_LIST <- c(ERPOS_EXCLUDE_SAMPLE_LIST, MISSING_UNKNOWN_ER)
CASE_SAMPLE_LIST <- ((pheno_data %>% filter(Status==1)))$`#IID`
TNBC_SAMPLE_LIST <- (pheno_data %>% filter(ER==2,PR==2,HER2==2,Status==1))$`#IID`
TNBC_EXCLUDE_SAMPLE_LIST <- setdiff(CASE_SAMPLE_LIST,TNBC_SAMPLE_LIST)
# also excluding iCOGs from TNBC
TNBC_EXCLUDE_SAMPLE_LIST <- unique(c(TNBC_EXCLUDE_SAMPLE_LIST,iCOGS_SAMPLE_LIST))

#################
# ERPOS BC RISK # 
#################
# outputting covariate and phenotype files
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)) %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)),file="training_ERPOS.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)),file="testing_ERPOS.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)),file="validation_ERPOS.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")

# outputting sample lists
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="training_ERPOS.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="testing_ERPOS.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% ERPOS_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="validation_ERPOS.samplist",quote=F,row.names=F,col.names=F,sep=" ")

#################
# ERNEG BC RISK # 
#################
# outputting covariate and phenotype files
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)) %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)),file="training_ERNEG.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)),file="testing_ERNEG.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)),file="validation_ERNEG.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")

# outputting sample lists
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="training_ERNEG.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="testing_ERNEG.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% ERNEG_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="validation_ERNEG.samplist",quote=F,row.names=F,col.names=F,sep=" ")

#################
# TNBC BC RISK # 
#################
# outputting covariate and phenotype files
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)) %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)),file="training_TNBC.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)),file="testing_TNBC.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)),file="validation_TNBC.pheno_cov",quote=F,row.names=F,col.names=T,sep=" ")

# outputting sample lists
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input")
write.table(training_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="training_TNBC.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(testing_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="testing_TNBC.samplist",quote=F,row.names=F,col.names=F,sep=" ")
write.table(validation_combined %>% filter(!(`#IID` %in% TNBC_EXCLUDE_SAMPLE_LIST)) %>% select(`#IID`),file="validation_TNBC.samplist",quote=F,row.names=F,col.names=F,sep=" ")

##################################################
# subsetting pfiles using plink for each subtype #
##################################################
# identifying high missingness EUR 330 variants to exclude 
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/ALL_VARIANT/ALL_VARIANT --remove /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/exclude_bcac_icogs.samplist --extract /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/list_EUR660_variants_hg38 --geno 0.00028 dosage --make-pgen --out /scratch/jll1/tmp/eur660_included_training")
eur660_included_training <- fread("/scratch/jll1/tmp/eur660_included_training.pvar") %>% select(ID)
eur660_excluded_training <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/list_EUR660_variants_hg38",header=F) %>% rename(ID=V1)
eur660_excluded_training <- eur660_excluded_training %>% filter(!(ID%in%eur660_included_training$ID))
write.table(eur660_excluded_training,file="/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/EUR_330_DOSAGE/HIGH_MISS_EUR660_EXCLUDE_TRAINING.txt",quote=F,col.names=F,row.names=F)

# subsetting genetic data while excluding high missingness EUR 330 variants from the training set 
tmp_subtype_string_list <- c("OVERALL","ERPOS","ERNEG","TNBC")
tmp_set_string_list <- c("training","testing","validation")

for (tmp_set_string in tmp_set_string_list) {
  for (tmp_subtype_string in tmp_subtype_string_list) {
    if (tmp_set_string=="training") {
      plink_command_pfile <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr --keep /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/",tmp_set_string,"_",tmp_subtype_string,".samplist --exclude /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/EUR_330_DOSAGE/HIGH_MISS_EUR660_EXCLUDE_TRAINING.txt --make-pgen --out split_",tmp_set_string,"_",tmp_subtype_string)
      system(plink_command_pfile)
      plink_command_bfile <- paste0("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/bfile/FILTERED_all_combined_chr --keep /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/",tmp_set_string,"_",tmp_subtype_string,".samplist --exclude /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/EUR_330_DOSAGE/HIGH_MISS_EUR660_EXCLUDE_TRAINING.txt --make-bed --out split_",tmp_set_string,"_",tmp_subtype_string)
      system(plink_command_bfile)
    } else {
      plink_command_pfile <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr --keep /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/",tmp_set_string,"_",tmp_subtype_string,".samplist --make-pgen --out split_",tmp_set_string,"_",tmp_subtype_string)
      system(plink_command_pfile)
      plink_command_bfile <- paste0("plink2 -bfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/bfile/FILTERED_all_combined_chr --keep /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/",tmp_set_string,"_",tmp_subtype_string,".samplist --make-bed --out split_",tmp_set_string,"_",tmp_subtype_string)
      system(plink_command_bfile)
    }
  }
}

