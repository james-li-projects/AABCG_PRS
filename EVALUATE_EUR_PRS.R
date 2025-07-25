library(data.table)
library(dplyr)
library(readxl)
library(pROC)
library(ROCnReg)

set.seed(1)

p_thresh_string <- "evaluateDataset"

AUC_vec <- c()
AUC_df <- data.frame()

for (current_set in c("training_OVERALL","testing_OVERALL","validation_OVERALL","training_ERPOS","testing_ERPOS","validation_ERPOS","training_ERNEG","testing_ERNEG","validation_ERNEG","training_TNBC","testing_TNBC","validation_TNBC")) {
# importing covariate data
cov <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/",current_set,".pheno_cov"),header=T)
colnames(cov)[1] <- "ID"
#cov <- cov %>% select(ID,Status)
cov$Status <- cov$Status-1
rownames(cov) <- cov$ID

#############################
# Computing EUR 313 SNP PRS #
#############################
print("[STEP 5A]: OBTAINING EFFECT SIZE ESTIMATES FOR THE EUR 313 SNP PRS")
print(Sys.time())
# set the working directory 
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp")

# importing the EUR 313 SNP PRS table with additional liftOver columns that were obtained using the liftOver website
#significant_SNP_EUR_313_SNP_PRS <- data.table(read_excel("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/EUR_313_SNP_PRS/mmc3_liftOver.xlsx")[c(1:313),])
#significant_SNP_EUR_313_SNP_PRS$effect_coef <- significant_SNP_EUR_313_SNP_PRS$`Overall Breast Cancerd`
#significant_SNP_EUR_313_SNP_PRS$A1 <- significant_SNP_EUR_313_SNP_PRS$`Effect Allele`
#write.table(significant_SNP_EUR_313_SNP_PRS %>% select(ID), file = paste0("significant_SNP_EUR_313_SNP_PRS_",p_thresh_string,".txt"), quote=F, row.names=F, col.names=F)

subtype=gsub("training_","",gsub("testing_","",gsub("validation_","",current_set)))
if (current_set %in% c("training_TNBC","testing_TNBC","validation_TNBC")) {
  significant_SNP_EUR_313_SNP_PRS <-fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_330_",subtype))
  
} else {
  significant_SNP_EUR_313_SNP_PRS <-fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_313_",subtype))
  
}
colnames(significant_SNP_EUR_313_SNP_PRS) <- c("ID","A1","effect_coef")
write.table(significant_SNP_EUR_313_SNP_PRS %>% select(ID), file = paste0("significant_SNP_EUR_313_SNP_PRS_",p_thresh_string,".txt"), quote=F, row.names=F, col.names=F)


print("[STEP 5A]: EXTRACTING AND IMPORTING DOSAGES OF THE EUR 313 SNP PRS FROM OUR evaluate SET")
print(Sys.time())
# run PLINK in R to extract dosages for evaluate set (214 of the 313 SNPs were successsfully extracted)
plink_command <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/split_",current_set," --extract ",paste0("significant_SNP_EUR_313_SNP_PRS_",p_thresh_string,".txt")," --export Av --out ./split_evaluate_",p_thresh_string)
system(plink_command)

# importing dosages for evaluate set into R
traw_file_EUR_313_SNP_PRS <- fread(paste0("./split_evaluate_",p_thresh_string,".traw"),header=T) %>% arrange(SNP)
print(nrow(traw_file_EUR_313_SNP_PRS))
# checking if sorting between the traw file and sumstats is identical
significant_SNP_EUR_313_SNP_PRS <- significant_SNP_EUR_313_SNP_PRS %>% filter(ID %in% traw_file_EUR_313_SNP_PRS$SNP) %>% arrange(ID)
identical(significant_SNP_EUR_313_SNP_PRS$ID, traw_file_EUR_313_SNP_PRS$SNP)

# swapping dosages depending on the counted allele in the traw file and the counted allele in the GWAS
index_different_counted_allele_EUR_313_SNP_PRS <- which(significant_SNP_EUR_313_SNP_PRS$A1 != traw_file_EUR_313_SNP_PRS$COUNTED)
traw_file_EUR_313_SNP_PRS[index_different_counted_allele_EUR_313_SNP_PRS, 7:ncol(traw_file_EUR_313_SNP_PRS)] <- 2 - traw_file_EUR_313_SNP_PRS[index_different_counted_allele_EUR_313_SNP_PRS, 7:ncol(traw_file_EUR_313_SNP_PRS)] 
dosage_mat_EUR_313_SNP_PRS <- data.matrix(traw_file_EUR_313_SNP_PRS[,7:ncol(traw_file_EUR_313_SNP_PRS)])

print("[STEP 5A]: COMPUTING PRS SCORES USING THE EUR 313 SNP PRS IN OUR evaluate SET")
print(Sys.time())
# computing EUR 313 SNP PRS in evaluate samples
PRS_weights_EUR_313_SNP_PRS <- significant_SNP_EUR_313_SNP_PRS$effect_coef



# computing PRS in evaluate samples
computed_EUR_313_SNP_PRS <- c()
PRS_weights <- significant_SNP_EUR_313_SNP_PRS$effect_coef
for (current_sample in cov$ID) {
  current_sample_string <- paste0("0_",current_sample)
  computed_EUR_313_SNP_PRS <- c(computed_EUR_313_SNP_PRS,
                               (sum(as.numeric(dosage_mat_EUR_313_SNP_PRS[,current_sample_string]) * PRS_weights))
  )
}

print("[STEP 5]: ASSESSING PREDICTIVE ABILITY OF PRS_AFR IN OUR evaluate SET")
print(Sys.time())
cov$PRS <- computed_EUR_313_SNP_PRS


# computing covariate adjusted AUC
output_AROC.sp <- AROC.sp(
  formula.h = paste0("PRS~Age+factor(Platform)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10"),
  group = "Status",
  data = cov,
  tag.h = 0)
print(paste("Covariate-adjusted AUC for",current_set,"method:", toString(signif(sort(as.numeric(output_AROC.sp$AUC)), digits = 4))))

UA_output_auc_evaluate_CI <- as.numeric(output_AROC.sp$AUC)
UA_output_auc_evaluate <- UA_output_auc_evaluate_CI[1]
UA_output_auc_evaluate_CI_lower <- UA_output_auc_evaluate_CI[2]
UA_output_auc_evaluate_CI_upper <- UA_output_auc_evaluate_CI[3]

print("##################################")
print(paste0("AUC METRICS OF THE SET: ",current_set))
print(UA_output_auc_evaluate)
print(UA_output_auc_evaluate_CI_lower)
print(UA_output_auc_evaluate_CI_upper)

AUC_vec <- c(AUC_vec,as.numeric(UA_output_auc_evaluate))
tmp_AUC_df <- cbind(current_set,t(data.frame(UA_output_auc_evaluate_CI)))
AUC_df <- rbind(AUC_df,tmp_AUC_df)
}

print("##################################")
print("EUR PRS AUCs of the training, testing, and validation sets:")
print(AUC_vec)
