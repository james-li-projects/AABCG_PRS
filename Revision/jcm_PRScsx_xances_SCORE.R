library(dplyr)
library(readr)
library(tidyr)
library(data.table)

# importing subtype and phi strings
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype = args[1]
phi = args[2]

# define output path
output_path <- paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_combined/PRScsx_",subtype,"_",phi,".tsv")

# exit early if output already exists
if (file.exists(output_path)) {
  cat("Output already exists at", output_path, "-- skipping.\n")
  quit(save = "no")
}

# setting seed 
set.seed(1)

# setting working directory 
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PRScsx/testing")

# loading in ancestry-specific PRS scores 
PRScsx_output_AFR <- fread(paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_scores/testing_AFR_",subtype,"_",phi,".sscore"), header = TRUE)
PRScsx_output_AFR$SCORE_AFR <- PRScsx_output_AFR$SCORE1_AVG
colnames(PRScsx_output_AFR)[1] <- "IID"
TMP_NORM_FACTOR_AFR <- median(PRScsx_output_AFR$ALLELE_CT)

PRScsx_output_EUR <- fread(paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_scores/testing_EUR_",subtype,"_",phi,".sscore"), header = TRUE)
PRScsx_output_EUR$SCORE_EUR <- PRScsx_output_EUR$SCORE1_AVG
colnames(PRScsx_output_EUR)[1] <- "IID"
TMP_NORM_FACTOR_EUR <- median(PRScsx_output_EUR$ALLELE_CT)

# importing testing set covariates
cov <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/testing_",subtype,".pheno_cov"),header=T)
cov$Status <- cov$Status-1
colnames(cov)[1] <- "IID"

# running logistic regression with covariates
PRScsx_output_joined <- inner_join(inner_join(PRScsx_output_AFR,PRScsx_output_EUR,by=c("IID")),cov,by=c("IID"))
fit <- summary(glm(data = PRScsx_output_joined, formula = Status ~ SCORE_AFR + SCORE_EUR + as.factor(Platform) + Age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10))
SCORE_AFR_weight <- fit$coefficients["SCORE_AFR",1]
SCORE_EUR_weight <- fit$coefficients["SCORE_EUR",1]

# read in scoring files
SCORE_AFR <- fread(paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx/AFR_",subtype,"_",phi,".txt"), header=FALSE) %>% mutate(EFFECT = V3 * SCORE_AFR_weight/TMP_NORM_FACTOR_AFR)
SCORE_EUR <- fread(paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx/EUR_",subtype,"_",phi,".txt"), header=FALSE) %>% mutate(EFFECT = V3 * SCORE_EUR_weight/TMP_NORM_FACTOR_EUR)

# combining scoring files with weights
SCORE_COMBINED <- rbind(SCORE_AFR,SCORE_EUR)
SCORE_COMBINED <- data.table(SCORE_COMBINED %>% dplyr::group_by(V1,V2) %>% summarise(COMBINED_EFFECT = sum(EFFECT)))
SCORE_COMBINED <- SCORE_COMBINED %>% filter(COMBINED_EFFECT != 0)

# writing out final scoring file
write.table(SCORE_COMBINED,
            file = output_path,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
