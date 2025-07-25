library(readxl)
library(bigsnpr)
library(dplyr)
library(readr)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(tidyr)
library(GenomicRanges)
library(data.table)
library(devtools)
library(tidyr)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
subtype=args[1]

# importing the entire list of rsids that can be included in the PRScsx algorithm
snpinfo_mult_1kg_hm3 <- fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/reference_panel/1kg/snpinfo_mult_1kg_hm3",header=T)

# import the rsid dictionary from Haoyu
info <- readRDS("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/rsid_dictionary_haoyu/SNP_GRCh37_38_match_update.rds")
info <- info %>% filter(rsid %in% snpinfo_mult_1kg_hm3$SNP)

# parsing this rsid dictionary 
info <- info %>% mutate(ID1 = paste(`chr`,`pos38`,`allele1_38`,`allele2_38`,sep=":")) %>% mutate(ID2 = paste(`chr`,`pos38`,`allele2_38`,`allele1_38`,sep=":"))
dict_orient_1 <- info %>% select(rsid,ID1)
colnames(dict_orient_1) <- c("rsid","SNP")
dict_orient_2 <- info %>% select(rsid,ID2)
colnames(dict_orient_2) <- c("rsid","SNP")

#################################
# setting working directory for reading and writing sumstats
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/parsed_xancestry_sumstats")

#################################
# reading in CT-SLEB parsed files and converting SNP IDs in these sumstats to rsIDs for AFR
CTSLEB_AFR <- fread(paste0("CTSLEB_hg38_sumstats_AFR_",subtype,".tsv"),header=T,sep=" ")
PRScsx_AFR_orient_1 <- inner_join(CTSLEB_AFR, dict_orient_1,by=c("SNP"))
PRScsx_AFR_orient_2 <- inner_join(CTSLEB_AFR, dict_orient_2,by=c("SNP"))
PRScsx_AFR <- rbind(PRScsx_AFR_orient_1,PRScsx_AFR_orient_2)
################################
# writing out SNP to rsID conversion dictionary
output_dictionary <- PRScsx_AFR %>% select(SNP,rsid)
write.table(output_dictionary,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/PRScsx/bim_dictionary_",subtype,".txt"),row.names=F,col.names=F,quote=F,sep="\t")
################################
PRScsx_AFR <- PRScsx_AFR %>% tidyr::separate(SNP,into=c("SNP_chr","SNP_pos","SNP_a2","SNP_a1"),sep=c(":"),remove=F)
# setting A2 alleles
PRScsx_AFR$A2[PRScsx_AFR$A1==PRScsx_AFR$SNP_a1] <- PRScsx_AFR$SNP_a2[PRScsx_AFR$A1==PRScsx_AFR$SNP_a1]
PRScsx_AFR$A2[PRScsx_AFR$A1!=PRScsx_AFR$SNP_a1] <- PRScsx_AFR$SNP_a1[PRScsx_AFR$A1!=PRScsx_AFR$SNP_a1]
# extracting relevant columns:
PRScsx_AFR <- PRScsx_AFR %>% select(rsid,A1,A2,BETA,SE)
colnames(PRScsx_AFR) <- c("SNP","A1","A2","BETA","SE")
# writing out sumstats
write.table(PRScsx_AFR, file=paste0("PRScsx_hg38_sumstats_AFR_",subtype,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")

#################################
# reading in CT-SLEB parsed files and converting SNP IDs in these sumstats to rsIDs for EUR
CTSLEB_EUR <- fread(paste0("CTSLEB_hg38_sumstats_EUR_",subtype,".tsv"),header=T,sep=" ")
PRScsx_EUR_orient_1 <- inner_join(CTSLEB_EUR, dict_orient_1,by=c("SNP"))
PRScsx_EUR_orient_2 <- inner_join(CTSLEB_EUR, dict_orient_2,by=c("SNP"))
PRScsx_EUR <- rbind(PRScsx_EUR_orient_1,PRScsx_EUR_orient_2)
PRScsx_EUR <- PRScsx_EUR %>% tidyr::separate(SNP,into=c("SNP_chr","SNP_pos","SNP_a2","SNP_a1"),sep=c(":"),remove=F)
# setting A2 alleles
PRScsx_EUR$A2[PRScsx_EUR$A1==PRScsx_EUR$SNP_a1] <- PRScsx_EUR$SNP_a2[PRScsx_EUR$A1==PRScsx_EUR$SNP_a1]
PRScsx_EUR$A2[PRScsx_EUR$A1!=PRScsx_EUR$SNP_a1] <- PRScsx_EUR$SNP_a1[PRScsx_EUR$A1!=PRScsx_EUR$SNP_a1]
# extracting relevant columns:
PRScsx_EUR <- PRScsx_EUR %>% select(rsid,A1,A2,BETA,SE)
colnames(PRScsx_EUR) <- c("SNP","A1","A2","BETA","SE")
# writing out sumstats
write.table(PRScsx_EUR, file=paste0("PRScsx_hg38_sumstats_EUR_",subtype,".tsv"),col.names=T,row.names=F,quote=F,sep="\t")
