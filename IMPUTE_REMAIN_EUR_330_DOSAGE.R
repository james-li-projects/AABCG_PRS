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

######################################
# Intersecting 330 and BCAC sumstats
eur_330_snplist<-(fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_330/TNBC.txt") %>% mutate(ID=paste(chr_name,chr_position,other_allele,effect_allele,sep="_")))$ID
eur_sumstats <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/public_download/oncoarray_bcac_public_release_oct17.txt") 
eur_sumstats_330 <- eur_sumstats %>% filter(var_name %in% eur_330_snplist) %>% select(chr, position_b37, a0, a1,bcac_onco_icogs_gwas_beta,bcac_onco_icogs_gwas_erpos_beta,bcac_onco_icogs_gwas_erneg_beta)
nrow(eur_sumstats_330)

# lifting over coordinates for the 330
eur_sumstats_330 <- eur_sumstats_330 %>% mutate(col4=seq(nrow(eur_sumstats_330)))
tmp_liftOver <- eur_sumstats_330 %>% mutate(col1 = paste0(rep("chr",nrow(eur_sumstats_330)), chr)) %>% mutate(col2 = position_b37) %>% mutate(col2 = as.integer(col2)) %>% mutate(col3 = col2 + 1) %>% select(col1,col2,col3,col4) 
# writing temporary table out for liftOver
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp")
write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp/tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
colnames(liftOver_mapped) <- c("col1","col2","col3","col4")
# joining with original table 
liftOver_eur_sumstats_330 <- inner_join(eur_sumstats_330, liftOver_mapped, by = c("col4"))
liftOver_eur_sumstats_330 <- liftOver_eur_sumstats_330 %>% mutate(col1 = gsub("chr","",col1)) %>% mutate(V1=paste(col1,col2,a0,a1,sep=":"),V2=a1) %>% rename(OVERALL_BETA=bcac_onco_icogs_gwas_beta,ERPOS_BETA=bcac_onco_icogs_gwas_erpos_beta,ERNEG_BETA=bcac_onco_icogs_gwas_erneg_beta)
BCAC_330_BETA <- liftOver_eur_sumstats_330 %>% select(V1,V2,OVERALL_BETA,ERPOS_BETA,ERNEG_BETA)

# obtaining EUR-330 PRS for OVERALL, ERPOS, and ERNEG BC
## first retrieving EUR 313
EUR_313_OVERALL <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_313_OVERALL")
EUR_313_ERPOS <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_313_ERPOS")
EUR_313_ERNEG <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_313_ERNEG")
## removing redundant variants from BCAC
BCAC_330_BETA <- BCAC_330_BETA %>% filter(!(V1 %in% EUR_313_OVERALL$V1))
## obtaining EUR-330 for OVERALL, ERPOS, and ERNEG BC for imputation
### OVERALL
BCAC_330_BETA_OVERALL <- BCAC_330_BETA %>% select(V1,V2,OVERALL_BETA) %>% rename(V3=OVERALL_BETA)
OVERALL_330_BETA <- rbind(EUR_313_OVERALL,BCAC_330_BETA_OVERALL) %>% rename(OVERALL_BETA=V3)
### ERPOS
BCAC_330_BETA_ERPOS <- BCAC_330_BETA %>% select(V1,V2,ERPOS_BETA) %>% rename(V3=ERPOS_BETA)
ERPOS_330_BETA <- rbind(EUR_313_ERPOS,BCAC_330_BETA_ERPOS)
ERPOS_330_BETA <- rbind(EUR_313_ERPOS,BCAC_330_BETA_ERPOS) %>% rename(ERPOS_BETA=V3)
### ERNEG
BCAC_330_BETA_ERNEG <- BCAC_330_BETA %>% select(V1,V2,ERNEG_BETA) %>% rename(V3=ERNEG_BETA)
ERNEG_330_BETA <- rbind(EUR_313_ERNEG,BCAC_330_BETA_ERNEG)
ERNEG_330_BETA <- rbind(EUR_313_ERNEG,BCAC_330_BETA_ERNEG) %>% rename(ERNEG_BETA=V3)

# importing TNBC 330
EUR_330_BETA_TNBC <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/EUR_330_TNBC") %>% rename(TNBC_BETA=V3)

# joining all beta coefficients
ALL_BETA <- OVERALL_330_BETA %>%
  inner_join(ERPOS_330_BETA, by=c("V1","V2")) %>% 
  inner_join(ERNEG_330_BETA, by=c("V1","V2")) %>%
  inner_join(EUR_330_BETA_TNBC, by=c("V1","V2"))
ALL_BETA <- data.frame(ALL_BETA)
rownames(ALL_BETA) <- ALL_BETA$V1

# check if risk alleles are the same as counted allele and reversing the sign of the beta coefficient if not the same
ALL_BETA_SELECT <- ALL_BETA[traw_BACKBONE$SNP,]
table(ALL_BETA_SELECT$V2 == traw_BACKBONE$COUNTED)
table(ALL_BETA_SELECT$V2 == traw_BACKBONE$ALT)
ALL_BETA_SELECT[ALL_BETA_SELECT$V2 == traw_BACKBONE$COUNTED, ]

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

# changing the imputed traw to have the same column order as the traw with data from other genetic variants  
con_header <- file("/scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/OUTPUT/output_combined_traw/combined_mean_imputed.traw", "r")
lines_1 <- readLines(con_header, n = 1)
traw_column_order <- strsplit(lines_1, "\t")[[1]]
close(con_header)
IMPUTED_330_TRAW <- IMPUTED_330_DOSAGE[,traw_column_order]
write.table(IMPUTED_330_TRAW,file="/scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/IMPUTED_330_TRAW.traw",quote=F,row.names=F,col.names=T,sep="\t")

# combining this data with the other mean imputed variant data for PRS
system("cp /scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/OUTPUT/output_combined_traw/combined_mean_imputed.traw /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/PRS.traw")
system("tail -n +2 /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/IMPUTED_330_TRAW.traw >> /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/PRS.traw")

# creating fam file
fam <- read.table("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr1.fam")
fam <- fam %>% mutate(V7=paste0("0_",V2))
rownames(fam) <- fam$V7
fam_iid_order<-traw_column_order[-c(1:6)]
fam <- fam[fam_iid_order,]
fam <- fam %>% select(-V7)
write.table(fam,file="/scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/PRS.fam",quote=F,row.names=F,col.names=F)

# converting merged file into a pgen file
system("plink2 --import-dosage /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/PRS.traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --fam /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/PRS.fam --make-pgen --sort-vars --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/PRS/PRS")


