###########################
# loading packages and setting working directory
library(data.table)
library(tidyverse)
library(readxl)
library(tidyverse)
library(stringr)
library(tidyr)
set.seed(523)

###########################
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 3-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
colnames(pheno_data)[2] <- "Dataset"
combined_cov <- pheno_data

# writing out a list to exclude BCAC/iCOGs participants from the pgen
exclude_list <- combined_cov %>% filter(Dataset=="Dataset_11_BCAC_iCOGS") %>% select(Sample_Name)
write.table(exclude_list, file = paste0("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/exclude_bcac_icogs.samplist"),col.names=F,quote=F,row.names=F)
# excluding bcac/icogs from our AF calculation
combined_cov <- combined_cov %>% filter(Dataset!="Dataset_11_BCAC_iCOGS")

# recoding the dataset variable 
combined_cov$dataset_index <- NA
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_01_WGS"] <- 1 
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_02_WGS2"] <- 2
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_03_MEGA_VANDY"] <- 3
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_04_MEGA_RP"] <- 4
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_05_MEGA_USC"] <- 5
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_06_AMBER"] <- 6
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_07_ROOT"] <- 7
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_08_AABC"] <- 8
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_09_GBHS"] <- 9
combined_cov$dataset_index[combined_cov$Dataset=="Dataset_10_BCAC_OncoArray"] <- 10

# writing out control participant lists for each AABCG dataset
for (j in c(1,2,3,4,5,6,7,8,9,10)) {
  # writing out lists of participants
  tmp_control_list <- (combined_cov %>% filter(Status == 1,dataset_index == j) %>% select(Sample_Name)) 
  print(nrow(tmp_control_list))
  write.table(tmp_control_list, file = paste0("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/dataset_",j,".samplist"),col.names=F,quote=F,row.names=F)
  
  # plink command to extract participants and calculate allele frequencies for each dataset
  plink_command <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/ALL_VARIANT/ALL_VARIANT --keep /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/dataset_",j,".samplist --freq --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/dataset_",j)
  system(plink_command)
}

# reading in allele frequencies for each dataset
setwd("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP")

dataset_1_AF <- fread("dataset_1.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_2_AF <- fread("dataset_2.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_3_AF <- fread("dataset_3.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_4_AF <- fread("dataset_4.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_5_AF <- fread("dataset_5.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_6_AF <- fread("dataset_6.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_7_AF <- fread("dataset_7.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_8_AF <- fread("dataset_8.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_9_AF <- fread("dataset_9.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)
dataset_10_AF <- fread("dataset_10.afreq",header=T) %>% select(ID,ALT,ALT_FREQS)

# checking if all alt alleles are the same
all(sapply(list(dataset_2_AF$ALT,
                dataset_3_AF$ALT,
                dataset_4_AF$ALT,
                dataset_5_AF$ALT,
                dataset_6_AF$ALT,
                dataset_7_AF$ALT,
                dataset_8_AF$ALT,
                dataset_9_AF$ALT,
                dataset_10_AF$ALT
), FUN = identical, dataset_1_AF$ALT))

# checking if all SNP IDs are the same
all(sapply(list(dataset_2_AF$ID,
                dataset_3_AF$ID,
                dataset_4_AF$ID,
                dataset_5_AF$ID,
                dataset_6_AF$ID,
                dataset_7_AF$ID,
                dataset_8_AF$ID,
                dataset_9_AF$ID,
                dataset_10_AF$ID
), FUN = identical, dataset_1_AF$ID))

# combining all allele frequencies and identifying SNPs with wide allele frequencies
combined_AF_df <- dataset_1_AF <- fread("dataset_1.afreq",header=T) %>% select(ID,ALT_FREQS)
for (k in c(2,3,4,5,6,7,8,9,10)) {
  tmp_AF_vec <- fread(paste0("dataset_",k,".afreq"),header=T) %>% select(ALT_FREQS)
  combined_AF_df <- cbind(combined_AF_df, tmp_AF_vec)
}
colnames(combined_AF_df) <- c(
  "ID",paste0(rep("d",10),c(1,2,3,4,5,6,7,8,9,10))
)

combined_AF_df$maxAF <- pmax(combined_AF_df$d1,
                             combined_AF_df$d2,
                             combined_AF_df$d3,
                             combined_AF_df$d4,
                             combined_AF_df$d5,
                             combined_AF_df$d6,
                             combined_AF_df$d7,
                             combined_AF_df$d8,
                             combined_AF_df$d9,
                             combined_AF_df$d10,na.rm=TRUE)
combined_AF_df$minAF <- pmin(combined_AF_df$d1,
                             combined_AF_df$d2,
                             combined_AF_df$d3,
                             combined_AF_df$d4,
                             combined_AF_df$d5,
                             combined_AF_df$d6,
                             combined_AF_df$d7,
                             combined_AF_df$d8,
                             combined_AF_df$d9,
                             combined_AF_df$d10,na.rm=TRUE)
# compute difference between max and min allele frequencies
combined_AF_df$rangeAF <- combined_AF_df$maxAF - combined_AF_df$minAF
# count the number of non-NA values for each variant
columns_to_consider<-paste0("d",1:10)
combined_AF_df$non_na_count <- apply(combined_AF_df[,..columns_to_consider], 1, function(x) sum(!is.na(x)))

# identifying list of variants to exclude based on allele frequency
exclude_SNP_list_AF <- combined_AF_df %>% filter(rangeAF > 0.15) %>% select(ID)
print(nrow(exclude_SNP_list_AF))

# identifying list of variants to exclude that were present only in one dataset or only in WGS datasets
exclude_SNP_list_exclusive_1 <- combined_AF_df %>% filter(( (!is.na(d1) | !is.na(d2)) & rowSums(!is.na(select(., d3:d10))) == 0 )) %>% select(ID)
exclude_SNP_list_exclusive_2 <- combined_AF_df %>% filter(non_na_count == 1) %>% select(ID)
exclude_SNP_list_exclusive <- unique(rbind(
  exclude_SNP_list_exclusive_1,
  exclude_SNP_list_exclusive_2
))
print(nrow(exclude_SNP_list_exclusive))

# remove these exclusive variants prior to identifying duplicated/swapped allele variants 
combined_AF_df <- combined_AF_df %>% filter(non_na_count > 1) %>% filter(!( (!is.na(d1) | !is.na(d2)) & rowSums(!is.na(select(., d3:d10))) == 0 ))

# identifying list of variants to exclude that are either duplicated or have swapped alleles across all the chromosomes
modified_current_chr_snp_list <- combined_AF_df %>% separate(ID, sep=":",into=c("chr","pos","a1","a2"),remove=F) %>% mutate(ID_2 = paste(chr,pos,a2,a1,sep=":"))
all_ID <- c(modified_current_chr_snp_list$ID,modified_current_chr_snp_list$ID_2)
dup_index <- which(duplicated(all_ID))
dup_ID <- all_ID[dup_index]
exclude_SNP_list_dupSNP_df <- data.frame(dup_ID)
colnames(exclude_SNP_list_dupSNP_df)[1] <- "ID"
print(nrow(exclude_SNP_list_dupSNP_df))

# writing out list of variants to exclude
exclude_SNP_list <- rbind(exclude_SNP_list_AF, exclude_SNP_list_dupSNP_df, exclude_SNP_list_exclusive) %>% unique()
print(nrow(exclude_SNP_list))
write.table(exclude_SNP_list, file = "exclude_SNP_list.txt", col.names=F, row.names=F, quote=F)

###########################################
# INITIAL FILTERING: exclude BCAC iCOGS and problematic variants
plink_command_icogs_variants <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/ALL_VARIANT/ALL_VARIANT --remove /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/exclude_bcac_icogs.samplist --exclude /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/exclude_SNP_list.txt --make-pgen --out /scratch/jll1/AABCG_PROCESSING/PLINK2_FILTERING/initial_FILTERED_all_combined_chr")
system(plink_command_icogs_variants)
print("FINISHED FILTERING DUPLICATED, WIDE ALLELE FREQUENCY, AND REMOVED BCAC_iCOGS SNPS")

###########################################
# PRS FILTERING: filtering for variants that are missing in less than 10 individuals [10/36191=0.00028]
plink_command_filter_prs <- paste0("plink2 -pfile /scratch/jll1/AABCG_PROCESSING/PLINK2_FILTERING/initial_FILTERED_all_combined_chr --geno 0.00028 dosage --make-pgen --out /scratch/jll1/AABCG_PROCESSING/PLINK2_FILTERING/initial_FILTERED_all_combined_chr_geno0.00028")
system(plink_command_filter_prs)
print("FINISHED FILTERING ON GENOTYPE MISSINGNESS BASED ON DOSAGES")

###############################################
# Preparing variant lists for mean imputation #
###############################################
# making lists of variants for partitions
pvar <- fread("/scratch/jll1/AABCG_PROCESSING/PLINK2_FILTERING/initial_FILTERED_all_combined_chr_geno0.00028.pvar") %>% select(ID)
tot_files<-ceiling(nrow(pvar)/5000)
repeat_index<-rep(1:tot_files,each=5000)
repeat_index<-repeat_index[1:nrow(pvar)]
pvar <- pvar %>% mutate(write_index=repeat_index)
# writing out these variant lists
for (b in 1:tot_files) {
  tmp_pvar_df <- pvar %>% filter(write_index == b)
  write.table(tmp_pvar_df %>% select(ID), file=paste0("/scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/INPUT/input_snplist/snplist_",b,".txt"), quote=F, row.names=F, col.names=F)
}

# writing out partition lists for input into mean imputation
system("cd /scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/INPUT/input_snplist/; ls snplist_* | split -l 25 -d -a 3 - /scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/INPUT/input_partitionlist/partition_")

#########################################
# Preparing traw of EUR 330 variants    #
# that were filtered out for imputation #
#########################################
# make traw for 660 variants that were filtered out, predominantly from high missingness rates, but also wide AF 
## importing EUR 660
EUR660_list <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/list_EUR660_variants_hg38",header=F)
initial_retained_EUR660_list <- EUR660_list %>% filter(V1%in%pvar$ID)
filtered_out_EUR660_list <- EUR660_list %>% filter(!(V1%in%pvar$ID))
write.table(filtered_out_EUR660_list,file="/scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_list",quote=F,row.names=F,col.names=F)
## creating pgen
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/ALL_VARIANT/ALL_VARIANT --remove /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/wideAF_dupSNP/exclude_bcac_icogs.samplist --extract /scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_list --geno 0.201 dosage --make-pgen --out /scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_geno_0.201")
## creating traw
system("plink2 -pfile /scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_geno_0.201 --export Av --out /scratch/jll1/AABCG_PROCESSING/DU_IMPUTATION/input/filtered_out_EUR660_geno_0.201")

###########################################
# GWAS FILTERING: filtering for MAF and geno 
plink_command_filter_gwas <- paste0("plink2 -pfile /scratch/jll1/AABCG_PROCESSING/PLINK2_FILTERING/initial_FILTERED_all_combined_chr --maf 0.01 --geno 0.2 dosage --make-pgen --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/GWAS/GWAS")
system(plink_command_filter_gwas)
print("FINISHED FILTERING BASED ON MAF AND GENO")
