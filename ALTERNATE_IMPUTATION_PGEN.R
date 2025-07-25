library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(ROCnReg)
library(pROC)

## setwd("/scratch/jll1/AABCG_PROCESSING/EUR330/traw")
setwd("/gpfs/data/huo-lab/BCAC/james.li/tmp/EUR330/traw")

# import traw
import_traw <- fread("/scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_geno_0.201.traw")

## calculate missing prop
miss_prop_vec <- rowSums(is.na(import_traw[,7:ncol(import_traw)]))/length(7:ncol(import_traw))
import_traw$prop_miss <- miss_prop_vec
# moving columns to beginning of file for cleaner subsequent processing
forward_cols <- c(head(colnames(import_traw),6),"prop_miss")
import_traw <- import_traw %>% relocate(all_of(forward_cols))
dim(import_traw)

# converting all dosages into numeric
import_traw[,8:ncol(import_traw)] <- lapply(import_traw[,8:ncol(import_traw)],as.numeric)

# harmonizing a SNP that is named differently from previous literature 
#### define conditions
condition1 <- "9:133271182:T:C" %in% import_traw$SNP
condition2 <- import_traw$ALT[import_traw$SNP=="9:133271182:T:C"] == "C"
condition3 <- import_traw$COUNTED[import_traw$SNP=="9:133271182:T:C"] == "T"
#### combine conditions 
boolean_result <- condition1 & condition2 & condition3
print(condition1)
print(condition2)
print(condition3)
print(boolean_result)
if (boolean_result == TRUE) {
  # swapping ALT and COUNTED alleles
  import_traw$ALT[import_traw$SNP=="9:133271182:T:C"] <- "T"
  import_traw$COUNTED[import_traw$SNP=="9:133271182:T:C"] <- "C"
  # reversing dosages
  print(summary(as.numeric(import_traw[which(import_traw$SNP=="9:133271182:T:C"),8:ncol(import_traw)])))
  import_traw[which(import_traw$SNP=="9:133271182:T:C"),8:ncol(import_traw)] <- 2 - import_traw[which(import_traw$SNP=="9:133271182:T:C"),8:ncol(import_traw)]
  print(summary(as.numeric(import_traw[which(import_traw$SNP=="9:133271182:T:C"),8:ncol(import_traw)])))
  # changing variant name
  import_traw$SNP[import_traw$SNP=="9:133271182:T:C"] <- "9:133271182:C:T"
  print("Harmonized 9:133271182:C:T")
}

############################################
# MEAN IMPUTATION LOW MISSINGNESS VARIANTS #
############################################
# mean imputing low missingness variants (less than 10/36191 individuals missing)
import_traw_1 <- import_traw %>% filter(prop_miss < 0.00028)
# saving a df of high missingness variants
import_traw_2 <- import_traw %>% filter(prop_miss > 0.00028)

# checking how many missing values there are in low missingness variants
## it looks like there is 0 missingness, but keeping this code in case another method of filtering requires it 
length(which(is.na(import_traw_1)))

# mean imputing values in import_traw_1
import_traw_1_backbone <- import_traw_1[,1:7]
import_traw_1_dosage <- import_traw_1[,8:ncol(import_traw_1)]
# replace NA values with row means
tmp_import_traw_1_dosage <- import_traw_1_dosage
tmp_import_traw_1_dosage <- t(apply(tmp_import_traw_1_dosage, 1, function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}))
# convert matrix back to data.frame
import_traw_1_dosage <- as.data.frame(tmp_import_traw_1_dosage)
# combining import_traw_1_backbone and import_traw_1_dosage that has been mean imputed
import_traw_1_imputed <- cbind(import_traw_1_backbone,import_traw_1_dosage)

# combining mean imputed dataframe for low missingness variants with the high missingness variants that have not been imputed yet
import_traw <- rbind(
  import_traw_1,
  import_traw_2
)

# subsetting traw to training set participants
import_traw <- import_traw %>% filter(prop_miss > 0.00028)

################################
# computing summary statistics #
################################
OVERALL_HIGHMISS <- data.frame()
ERPOS_HIGHMISS <- data.frame()
ERNEG_HIGHMISS <- data.frame()
TNBC_HIGHMISS <- data.frame()

# defining subtypes
for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  # importing training set covariates
  current.pheno_cov <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/training_",subtype,".pheno_cov"))
  current.pheno_cov <- current.pheno_cov %>% mutate(query_id=paste0("0_",`#IID`))
  
  # subsetting traw to backbone and dosages
  import_traw_BACKBONE <- import_traw[,1:6]
  select_ind<-current.pheno_cov$query_id
  import_traw_DOSAGE <- import_traw[,..select_ind]
  combined_import_traw <- cbind(import_traw_BACKBONE,import_traw_DOSAGE)
  combined_import_traw <- t(combined_import_traw)
  colnames(combined_import_traw) <- combined_import_traw[2,]
  combined_import_traw <- combined_import_traw[-(1:6),]
  combined_import_traw <- cbind(query_id = rownames(combined_import_traw), combined_import_traw)
  combined_import_traw <- data.frame(combined_import_traw)
  colnames(combined_import_traw) <- gsub("X","",colnames(combined_import_traw))
  colnames(combined_import_traw) <- gsub("\\.",":",colnames(combined_import_traw))

  # making a combined df for regression
  reg_df<-inner_join(current.pheno_cov,combined_import_traw,by=c("query_id")) %>% select(-query_id) %>% mutate(Status=Status - 1)
  
  # converting all values to numeric
  for (i in 45:ncol(reg_df)) {
    reg_df[,i]=as.numeric(unlist(reg_df[,..i]))
  }

  # running logistic regressions
  logistic_beta <- data.frame()
  for (i in 45:ncol(reg_df)) {
    print(i)
    # Define the outcome variable and predictor variables
    outcome_var <- "Status"
    predictor_vars <- c(paste0("`",colnames(reg_df)[i],"`"), "Age", "Platform", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")  # Predictors
    
    # Construct the formula
    formula <- as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + ")))
    
    # Fit the logistic regression model
    model <- glm(formula, data = reg_df, family = binomial)
    
    # Display the summary of the model
    fit<-summary(model)
    
    # storing values from regression output
    tmp_logistic_beta <- data.frame(SNP=colnames(reg_df)[i],beta=fit$coefficients[2,"Estimate"])
    logistic_beta <- rbind(logistic_beta,tmp_logistic_beta)
  }
  
  # join values with backbone
  joined_backbone_beta<-inner_join(import_traw_BACKBONE,logistic_beta,by=c("SNP")) 
  joined_backbone_beta <- joined_backbone_beta %>% mutate(BETA=ifelse(COUNTED != ALT, (-1)*beta, beta))
  
  joined_backbone_beta <- joined_backbone_beta %>% select(SNP,ALT,BETA) %>% rename(ID=SNP,A1=ALT)
  
  if (subtype == "OVERALL") {
    OVERALL_HIGHMISS <- joined_backbone_beta
    print(subtype)
  } else if (subtype == "ERPOS") {
    ERPOS_HIGHMISS <- joined_backbone_beta
    print(subtype)
  } else if (subtype == "ERNEG") {
    ERNEG_HIGHMISS <- joined_backbone_beta
    print(subtype)
  } else if (subtype == "TNBC") {
    TNBC_HIGHMISS <- joined_backbone_beta
    print(subtype)
  } else {
    stop("Unknown subtype value.")
  }
}

#######################################
# importing covariate data
library(readxl)
##pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- pheno_data %>% filter(`Dataset...3`!="Dataset_11_BCAC_iCOGS") %>% 
  select(`AABCGS_ID...1`,Status,ER,PR,HER2,Platform,`Dataset...3`)
colnames(pheno_data)[1] <- "#IID"
pheno_data <- pheno_data %>% rename(Dataset=`Dataset...3`)

# TNBC list
TNBC_SAMPLE_LIST <- (pheno_data %>% filter(ER==2,PR==2,HER2==2,Status==1))$`#IID`
TNBC_SAMPLE_LIST <- paste0("0_",TNBC_SAMPLE_LIST)

# ERPOS list
ERPOS_SAMPLE_LIST <- (pheno_data %>% filter(ER==1,Status==1))$`#IID`
ERPOS_SAMPLE_LIST <- paste0("0_",ERPOS_SAMPLE_LIST)
ERPOS_SAMPLE_LIST <- setdiff(ERPOS_SAMPLE_LIST,TNBC_SAMPLE_LIST)  

# ERNEG list
ERNEG_SAMPLE_LIST <- (pheno_data %>% filter(ER==2,Status==1))$`#IID`
ERNEG_SAMPLE_LIST <- paste0("0_",ERNEG_SAMPLE_LIST)
## ER- but not TNBC
ERNEG_SAMPLE_LIST <- setdiff(ERNEG_SAMPLE_LIST,TNBC_SAMPLE_LIST)

# Remaining overall breast cancer cases list 
OVERALL_SAMPLE_LIST <- (pheno_data %>% filter(Status==1))$`#IID`
OVERALL_SAMPLE_LIST <- paste0("0_",OVERALL_SAMPLE_LIST)
SUBTYPE_CASE_LIST <- c(TNBC_SAMPLE_LIST,ERPOS_SAMPLE_LIST,ERNEG_SAMPLE_LIST)
# All cases 
CASE_SAMPLE_LIST <- OVERALL_SAMPLE_LIST
OVERALL_SAMPLE_LIST <- setdiff(OVERALL_SAMPLE_LIST,SUBTYPE_CASE_LIST)

# Control sample list
CONTROL_SAMPLE_LIST <- (pheno_data %>% filter(Status==2))$`#IID`
CONTROL_SAMPLE_LIST <- paste0("0_",CONTROL_SAMPLE_LIST)

# subsetting traw into different partitions for each case type and controls
traw_BACKBONE<-import_traw[,..forward_cols]
traw_CONTROL<-import_traw[,..CONTROL_SAMPLE_LIST]
traw_OVERALL<-import_traw[,..OVERALL_SAMPLE_LIST]
traw_ERPOS<-import_traw[,..ERPOS_SAMPLE_LIST]
traw_ERNEG<-import_traw[,..ERNEG_SAMPLE_LIST]
traw_TNBC<-import_traw[,..TNBC_SAMPLE_LIST]

# obtaining mean imputed allele frequency in control samples
control_AF <- as.numeric(rowMeans(traw_CONTROL, na.rm=TRUE))/2

# retrieving list of variants
EUR_330_LIST <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_330_TNBC")$V1

# import score files (AFR for all subtypes except ERPOS)
OVERALL_SCORE_FILE <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/IMPUTE_330/SCORE_FILES/AFR_330_OVERALL")
colnames(OVERALL_SCORE_FILE) <- c("ID","A1","OVERALL_BETA")
ERPOS_SCORE_FILE <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/IMPUTE_330/SCORE_FILES/EUR_330_ERPOS")
colnames(ERPOS_SCORE_FILE) <- c("ID","A1","ERPOS_BETA")
ERNEG_SCORE_FILE <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/IMPUTE_330/SCORE_FILES/AFR_330_ERNEG")
colnames(ERNEG_SCORE_FILE) <- c("ID","A1","ERNEG_BETA")
TNBC_SCORE_FILE <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/IMPUTE_330/SCORE_FILES/AFR_330_TNBC")
colnames(TNBC_SCORE_FILE) <- c("ID","A1","TNBC_BETA")

# joining all beta coefficients from scoring files
ALL_BETA <- OVERALL_SCORE_FILE %>%
  inner_join(ERPOS_SCORE_FILE, by=c("ID","A1")) %>% 
  inner_join(ERNEG_SCORE_FILE, by=c("ID","A1")) %>%
  inner_join(TNBC_SCORE_FILE, by=c("ID","A1"))
ALL_BETA <- data.frame(ALL_BETA)
rownames(ALL_BETA) <- ALL_BETA$ID
ALL_BETA <- ALL_BETA %>% rename(V2=A1,V1=ID)

# check if risk alleles are the same as counted allele and reversing the sign of the beta coefficient if not the same
ALL_BETA_SELECT <- ALL_BETA[traw_BACKBONE$SNP,]
table(ALL_BETA_SELECT$V2 == traw_BACKBONE$COUNTED)
table(ALL_BETA_SELECT$V2 == traw_BACKBONE$ALT)

ALL_BETA_SELECT[, c(3:ncol(ALL_BETA_SELECT))] <- lapply(ALL_BETA_SELECT[, c(3:ncol(ALL_BETA_SELECT))],as.numeric)
ALL_BETA_SELECT <- ALL_BETA_SELECT %>% mutate( 
  across(names(ALL_BETA_SELECT[, c(3:ncol(ALL_BETA_SELECT))]), ~ case_when(
      ALL_BETA_SELECT$V2 == traw_BACKBONE$ALT ~ (-1) * .,
      TRUE ~ .) 
  ))

#####################################
# imputing dosages for each subtype #
#####################################
length(which(is.na(traw_CONTROL)))
control_rep <- matrix(rep(control_AF*2, times=ncol(traw_CONTROL)), nrow=nrow(traw_CONTROL), ncol=ncol(traw_CONTROL))
traw_CONTROL[which(is.na(traw_CONTROL), arr.ind=TRUE)] <- control_rep[which(is.na(traw_CONTROL), arr.ind=TRUE)]
rm(control_rep)

length(which(is.na(traw_ERPOS)))
ERpos_rep <- matrix(rep(
  2*(control_AF*exp(ALL_BETA_SELECT$ERPOS_BETA)/(1-control_AF+control_AF*exp(ALL_BETA_SELECT$ERPOS_BETA))), 
  times=ncol(traw_ERPOS)),  nrow=nrow(traw_ERPOS), ncol=ncol(traw_ERPOS))
traw_ERPOS[which(is.na(traw_ERPOS), arr.ind=TRUE)] <- ERpos_rep[which(is.na(traw_ERPOS), arr.ind=TRUE)]
rm(ERpos_rep)

length(which(is.na(traw_ERNEG)))
ERneg_rep <- matrix(rep(
  2*(control_AF*exp(ALL_BETA_SELECT$ERNEG_BETA)/(1-control_AF+control_AF*exp(ALL_BETA_SELECT$ERNEG_BETA))), 
  times=ncol(traw_ERNEG)),  nrow=nrow(traw_ERNEG), ncol=ncol(traw_ERNEG))
traw_ERNEG[which(is.na(traw_ERNEG), arr.ind=TRUE)] <- ERneg_rep[which(is.na(traw_ERNEG), arr.ind=TRUE)]
rm(ERneg_rep)

length(which(is.na(traw_TNBC)))
TNBC_rep <- matrix(rep(
  2*(control_AF*exp(ALL_BETA_SELECT$TNBC_BETA)/(1-control_AF+control_AF*exp(ALL_BETA_SELECT$TNBC_BETA))), 
  times=ncol(traw_TNBC)),  nrow=nrow(traw_TNBC), ncol=ncol(traw_TNBC))
traw_TNBC[which(is.na(traw_TNBC), arr.ind=TRUE)] <- TNBC_rep[which(is.na(traw_TNBC), arr.ind=TRUE)]
rm(TNBC_rep)

length(which(is.na(traw_OVERALL)))
OVERALL_rep <- matrix(rep(
  2*(control_AF*exp(ALL_BETA_SELECT$OVERALL_BETA)/(1-control_AF+control_AF*exp(ALL_BETA_SELECT$OVERALL_BETA))), 
  times=ncol(traw_OVERALL)),  nrow=nrow(traw_OVERALL), ncol=ncol(traw_OVERALL))
traw_OVERALL[which(is.na(traw_OVERALL), arr.ind=TRUE)] <- OVERALL_rep[which(is.na(traw_OVERALL), arr.ind=TRUE)]
rm(OVERALL_rep)

# combining all the imputed dosage matrices together with traw backbone
IMPUTED_330_DOSAGE <- cbind.data.frame(
  traw_BACKBONE,
  traw_CONTROL,
  traw_OVERALL,
  traw_ERPOS,
  traw_ERNEG,
  traw_TNBC
)
IMPUTED_330_DOSAGE <- IMPUTED_330_DOSAGE %>% select(-prop_miss)

# checking if all dosages got imputed
length(which(is.na(IMPUTED_330_DOSAGE)))
dim(IMPUTED_330_DOSAGE)

# writing out high missingness variant list
write.table(data.frame(IMPUTED_330_DOSAGE$SNP),file="/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/EUR_IMPUTE/highmiss_variant.list",quote=F,col.names=F,row.names=F)

# assembling final testing/validation pgen files for each subtype
for (subtype in c("OVERALL","ERPOS","ERNEG","TNBC")) {
  # exporting traw for testing/validation sets without high missingness variants
  system(paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/EUR_IMPUTE/combined_testing_validation_",subtype," --exclude /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/EUR_IMPUTE/highmiss_variant.list --export Av --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype))
  
  # identifying column/sample order for traw
  first_line <- readLines(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype,".traw"), n = 1)
  column_names <- strsplit(first_line, "\t")[[1]]
  
  # writing out imputed high missingness dosage for the current subtype
  CURRENT_IMPUTED_330_DOSAGE <- IMPUTED_330_DOSAGE %>% select(all_of(column_names))
  write.table(CURRENT_IMPUTED_330_DOSAGE,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/HIGHMISS_54_SNP_noheader_",subtype,".tsv"),quote=F,row.names=F,col.names=F,sep="\t")
  
  # catting this file with the larger traw
  system(paste0("cat /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/HIGHMISS_54_SNP_noheader_",subtype,".tsv >> /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype,".traw"))
  
  # creating fam file
  sample_order <- column_names[-c(1:6)]
  
  # creating fam files for these traw files
  example_fam <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr1.fam")
  example_fam <- data.frame(example_fam)
  example_fam <- example_fam %>% mutate(V7=paste0("0_",V2))
  rownames(example_fam) <- example_fam$V7
  example_fam<-example_fam[sample_order,]
  example_fam$V7 <- NULL
  write.table(example_fam,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype,".fam"),quote=F,row.names=F,col.names=F)
  
  # converting traw to pgen
  system(paste0("plink2 --import-dosage /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype,".traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --fam /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/FINAL_IMPUTE/combined_testing_validation_",subtype,".fam --sort-vars --make-pgen --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/combined_testing_validation_",subtype))
} 
