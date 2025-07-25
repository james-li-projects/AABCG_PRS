library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)

setwd("/scratch/jll1/AABCG_PROCESSING/EUR330/traw")

# importing dosages from AABCG WGS dataset
combined_traw_WGS <- data.frame()
for (j in 1:22) {
  print(j)
  
  # reading in file if it exists for indels
  if (file.exists(paste0("AABCG_WGS_indel_330_chr",j,".traw"))) {
    tmp_traw_WGS_1 <- fread(paste0("AABCG_WGS_indel_330_chr",j,".traw"))
  } else {
    tmp_traw_WGS_1 <- data.frame()
  }
  
  # reading in file if it exists for snps
  if (file.exists(paste0("AABCG_WGS_snp_330_chr",j,".traw"))) {
    tmp_traw_WGS_2 <- fread(paste0("AABCG_WGS_snp_330_chr",j,".traw"))
  } else {
    tmp_traw_WGS_2 <- data.frame()
  }
  
  tmp_traw_WGS <- rbind(tmp_traw_WGS_1,tmp_traw_WGS_2)
  combined_traw_WGS <- rbind(combined_traw_WGS,tmp_traw_WGS)
}

# list of datasets other AABCG datasets
aabcg_dataset_list <- c("AABC","AMBER","GBHS","MEGA_MEC","MEGA_RP","MEGA_Vandy","Onco","ROOT","WashU","iCOGS")

imported_traw <- vector(mode="list",length=length(aabcg_dataset_list))
for (i in 1:length(aabcg_dataset_list)) {
  study <- aabcg_dataset_list[i]
  print(study)
  combined_traw <- data.frame()
  for (j in 1:22) {
    print(j)
    tmp_traw <- fread(paste0("AABCG_",study,"_330_chr",j,".traw"))
    combined_traw <- rbind(combined_traw,tmp_traw)
  }
  imported_traw[[i]] <- combined_traw %>% select(-CHR,-`(C)M`,-POS)
}

combined_dataset_traw_noWGS <- Reduce(full_join,imported_traw)
combined_dataset_traw_noWGS <- data.frame(combined_dataset_traw_noWGS)

# combine the 330 variant data for WGS and all other datasets
combined_data_set_traw_ALL <- full_join(combined_traw_WGS, combined_dataset_traw_noWGS, by = c("SNP","COUNTED","ALT"))
combined_data_set_traw_ALL <- combined_data_set_traw_ALL %>% mutate(SNP=gsub("chr","",SNP))
# parsing this traw file to extract chromosome and position
combined_data_set_traw_ALL <- combined_data_set_traw_ALL %>% separate(SNP,into=c("F1","F2","F3","F4"),sep=":",remove=F) %>% mutate(CHR=F1,POS=F2,`(C)M`=0) %>% select(-F1,-F2,-F3,-F4)
# parsing individual IDs to match with previously processed data
colnames(combined_data_set_traw_ALL) <- str_extract(colnames(combined_data_set_traw_ALL), "[^_]+")
colnames(combined_data_set_traw_ALL) <- paste0("0_",colnames(combined_data_set_traw_ALL))
combined_data_set_traw_ALL <- combined_data_set_traw_ALL %>% rename(CHR=`0_CHR`,SNP=`0_SNP`,`(C)M`=`0_(C)M`,POS=`0_POS`,COUNTED=`0_COUNTED`,ALT=`0_ALT`)

# reading in the original dosages for 330 variants (277 of them are available)
original_660 <- read.table("/scratch/jll1/AABCG_PROCESSING/EUR330/original_traw/original_660.traw",header=T) 
original_660 <- original_660 %>% arrange(SNP) %>% rename(`(C)M`=`X.C.M`)
colnames(original_660) <- gsub("X0_","0_",colnames(original_660))

# parsing this new dosage matrix of 330 variants to match the original column ordering then dividing into a dataframe to harmonize with DS2 and a dataframe to combine later
combined_data_set_traw_ALL <- combined_data_set_traw_ALL %>% select(colnames(original_660))
# combined_traw_part_1 needs changed
# combined_traw_part_2 does not need changed
combined_traw_part_1 <- combined_data_set_traw_ALL %>% filter(SNP %in% original_660$SNP) %>% arrange(SNP)
combined_traw_part_2 <- combined_data_set_traw_ALL %>% filter(!(SNP %in% original_660$SNP)) %>% arrange(SNP)
original_660_part_1 <-  original_660 %>% filter(SNP %in% combined_data_set_traw_ALL$SNP) %>% arrange(SNP)
original_660_part_3 <-  original_660 %>% filter(!(SNP %in% combined_data_set_traw_ALL$SNP)) %>% arrange(SNP)
dim(combined_traw_part_1)
dim(combined_traw_part_2)
dim(original_660_part_1)
dim(original_660_part_3)

# obtaining original DS2 dosages if the current 330 variant import is NA
identical(original_660_part_1$SNP,combined_traw_part_1$SNP)
identical(original_660_part_1$COUNTED,combined_traw_part_1$COUNTED)
nrow(which(is.na(combined_traw_part_1), arr.ind=TRUE))
nrow(which(is.na(original_660_part_1), arr.ind=TRUE))

original_660_part_1[which(is.na(original_660_part_1), arr.ind=TRUE)] <- combined_traw_part_1[which(is.na(original_660_part_1), arr.ind=TRUE)]
nrow(which(is.na(original_660_part_1), arr.ind=TRUE))

length(which(is.na(original_660_part_3))) 
length(which(is.na(combined_traw_part_2))) 

# recombining 660 files: using original data as the master and using proper alignment to be read by plink
ALL_660_VARIANT <-rbind(
  original_660_part_1, combined_traw_part_2, original_660_part_3
)
dim(ALL_660_VARIANT)
ALL_660_VARIANT[7:ncol(ALL_660_VARIANT)] <- lapply(ALL_660_VARIANT[,7:ncol(ALL_660_VARIANT)],as.numeric)
ALL_660_VARIANT[7:ncol(ALL_660_VARIANT)] <- 2-ALL_660_VARIANT[7:ncol(ALL_660_VARIANT)]
ALL_660_VARIANT <- ALL_660_VARIANT %>% relocate(ALT, .after = POS)
ALL_660_VARIANT <- ALL_660_VARIANT %>% select(-`(C)M`) %>% rename(A1=ALT,A2=COUNTED)

# reading in psam file to obtain column order
#psam <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr.psam")

# reading in one of the fam files
fam <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr1.fam")
fam <- data.frame(fam)
fam <- fam %>% mutate(V7=paste0("0_",V2))
rownames(fam) <- fam$V7
fam$V7 <- NULL

# aligning order of fam files to that of the dosage df
individual_order <- (colnames(ALL_660_VARIANT))[-(1:5)]
fam <- fam[individual_order,]

# modified sample names to remove 0_ prefix
MODIFIED_ALL_660_VARIANT <- ALL_660_VARIANT
colnames(MODIFIED_ALL_660_VARIANT) <- gsub("0_","",colnames(MODIFIED_ALL_660_VARIANT))
MODIFIED_ALL_660_VARIANT[is.na(MODIFIED_ALL_660_VARIANT)] <- "."

# writing out dosage and fam files
write.table(MODIFIED_ALL_660_VARIANT,file="/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT.dosage",quote=F,row.names=F,col.names=F)
write.table(fam,file="/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT.fam",quote=F,row.names=F,col.names=F,sep="\t")

# convert dosage to pfile
system("
input_dosage_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT.dosage; 
input_fam_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT.fam; 
output_pfile=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT;
plink2 --import-dosage ${input_dosage_file} noheader chr-col-num=1 pos-col-num=3 format=1 skip0=1 skip1=1 --fam ${input_fam_file} --threads 1 --memory 100000 --sort-vars --make-pgen --out ${output_pfile} 
")

####################################################
# generating an entire master pfile by combining the 660 variant pgen with the rest of the dataset

## removing 660 variants 
remove_660_command <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr --exclude /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38/list_EUR660_variants_hg38 --sort-vars --make-pgen --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr_NO_660")
system(remove_660_command)

## converting to traw file
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr_NO_660 --export Av --out /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr_NO_660")
system("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/COMBINED_660_VARIANT/ALL_660_VARIANT --export Av --out /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/ALL_660_VARIANT")
## obtaining sample order from traw files
#### obtaining first line from first traw
con_1 <- file("/scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/ALL_660_VARIANT.traw", "r")
lines_1 <- readLines(con_1, n = 1)
data_1 <- strsplit(lines_1, "\t")[[1]]
data_1 <- data_1[7:length(data_1)]
close(con_1)
#### obtaining first line from second traw
con_2 <- file("/scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr_NO_660.traw", "r")
lines_2 <- readLines(con_2, n = 1)
data_2 <- strsplit(lines_2, "\t")[[1]]
data_2 <- data_2[7:length(data_2)]
close(con_2)
#### check if sample order is identical
if (identical(data_1,data_2) == TRUE) {
  sample_order <- data_1
} else {
  sample_order <- NULL
}

## obtaining a file removing header from 660 variant traw
system("tail -n+2 /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/ALL_660_VARIANT.traw > /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/ALL_660_VARIANT_noheader.traw")
## cat this file to the master traw
system("mv /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr_NO_660.traw /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr.traw")
system("cat /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/ALL_660_VARIANT_noheader.traw >> /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr.traw")

# creating fam files for these traw files
example_fam <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr1.fam")
example_fam <- data.frame(example_fam)
example_fam <- example_fam %>% mutate(V7=paste0("0_",V2))
rownames(example_fam) <- example_fam$V7
example_fam<-example_fam[sample_order,]
example_fam$V7 <- NULL
write.table(example_fam,file="/scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr.fam",quote=F,row.names=F,col.names=F)

# converting merged file into a pgen file
system("plink2 --import-dosage /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr.traw skip0=1 skip1=2 id-delim=_ chr-col-num=1 pos-col-num=4 ref-first --fam /scratch/jll1/AABCG_PROCESSING/PLINK2_ASSEMBLY/all_combined_chr.fam --make-pgen --sort-vars --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/MASTER_FILES/ALL_VARIANT/ALL_VARIANT")

