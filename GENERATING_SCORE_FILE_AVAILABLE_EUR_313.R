library(dplyr)
library(data.table)
library(stringr)
library(pROC)
library(tidyr)

######################################
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_313")
subtype_list <- unlist(str_split(list.files(),pattern="\\."))
subtype_list <- subtype_list[!grepl("txt",subtype_list)]

######################################
# obtaining list of all variants that are available in our genetic dataset
available_variant_list <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr.pvar")$ID

######################################
# Intersecting 330 and BCAC sumstats
eur_330_snplist<-(fread("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_330/TNBC.txt") %>% mutate(ID=paste(chr_name,chr_position,other_allele,effect_allele,sep="_")))$ID
eur_sumstats <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/public_download/oncoarray_bcac_public_release_oct17.txt") 
eur_sumstats_330 <- eur_sumstats %>% filter(var_name %in% eur_330_snplist) %>% select(chr, position_b37, a0, a1,bcac_onco_icogs_gwas_beta,bcac_onco_icogs_gwas_erpos_beta,bcac_onco_icogs_gwas_erneg_beta)
nrow(eur_sumstats_330)

# obtaining liftOver of BCAC summary statistics
# lifting over coordinates for the 330
eur_sumstats_330 <- eur_sumstats_330 %>% mutate(col4=seq(nrow(eur_sumstats_330)))
tmp_liftOver_sumstats <- eur_sumstats_330 %>% mutate(col1 = paste0(rep("chr",nrow(eur_sumstats_330)), chr)) %>% mutate(col2 = position_b37) %>% mutate(col2 = as.integer(col2)) %>% mutate(col3 = col2 + 1) %>% select(col1,col2,col3,col4) 
# writing temporary table out for liftOver
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp")
write.table(tmp_liftOver_sumstats, file="tmp_liftOver_sumstats", quote=F,row.names=F,col.names=F,sep="\t")
# performing liftOver
system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp/tmp_liftOver_sumstats /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped_sumstats liftOver_unmapped")
# reading in liftOver mapped results
liftOver_mapped_sumstats <- fread("liftOver_mapped_sumstats",header=F,sep="\t")
colnames(liftOver_mapped_sumstats) <- c("col1","col2","col3","col4")
# joining with original table 
liftOver_eur_sumstats_330 <- inner_join(eur_sumstats_330, liftOver_mapped_sumstats, by = c("col4"))
liftOver_eur_sumstats_330 <- liftOver_eur_sumstats_330 %>% mutate(col1 = gsub("chr","",col1)) %>% mutate(V1=paste(col1,col2,a0,a1,sep=":"),V2=a1) %>% rename(OVERALL_BETA=bcac_onco_icogs_gwas_beta,ERPOS_BETA=bcac_onco_icogs_gwas_erpos_beta,ERNEG_BETA=bcac_onco_icogs_gwas_erneg_beta)
BCAC_330_BETA <- liftOver_eur_sumstats_330 %>% select(V1,V2,OVERALL_BETA,ERPOS_BETA,ERNEG_BETA)

######################################
# writing out subtype scoring files
for (subtype in subtype_list) {
  print(paste("PROCESSING SCORE FILE FOR:",subtype))
  setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_313")
  # importing scoring file and processing file for liftOver input
  PRS_313 <- fread(paste0(subtype,".txt"),header=F)
  PRS_313 <- PRS_313 %>% separate(V1,into=c("CHR","hg19","BASELINE_ALLELE","EFFECT_ALLELE"),sep="\\:") %>% rename(BETA=V3)
  print(identical(PRS_313$EFFECT_ALLELE,PRS_313$V2))
  PRS_313 <- PRS_313 %>% select(-V2)
  PRS_313 <- PRS_313 %>% mutate(col4 = seq(nrow(PRS_313))) %>% mutate(col4 = as.integer(col4))
  tmp_liftOver <- PRS_313 %>% mutate(col1 = paste0(rep("chr",nrow(PRS_313)), CHR)) %>% mutate(col2 = hg19) %>% mutate(col2 = as.integer(col2)) %>% mutate(col3 = col2 + 1) %>% select(col1,col2,col3,col4) 
  # writing temporary table out for liftOver
  setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp")
  write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
  # performing liftOver
  system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp/tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
  # reading in liftOver mapped results
  liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
  colnames(liftOver_mapped) <- c("col1","col2","col3","col4")
  # joining with original table 
  liftOver_PRS_313 <- inner_join(PRS_313, liftOver_mapped, by = c("col4"))
  
  # formatting the output score file
  liftOver_PRS_313$BETA <- as.numeric(liftOver_PRS_313$BETA)
  liftOver_PRS_313$A1 <- liftOver_PRS_313$EFFECT_ALLELE
  
  # formatting liftOver ID
  liftOver_PRS_313$ID <- paste(
    liftOver_PRS_313$CHR,
    liftOver_PRS_313$col2,
    liftOver_PRS_313$BASELINE_ALLELE,
    liftOver_PRS_313$EFFECT_ALLELE,
    sep=":"
  )
  
  # extract variants not in the 313
  current_additional <- BCAC_330_BETA %>% select(V1,V2,paste0(subtype,"_BETA"))
  colnames(current_additional) <- c("ID","A1","BETA")
  current_additional <- current_additional %>% filter(!(ID%in%liftOver_PRS_313$ID))
  
  # finalizing the output score file  
  liftOver_PRS_313 <- liftOver_PRS_313 %>% select(ID,A1,BETA)
  combined_330_score_file <- rbind(liftOver_PRS_313,current_additional)
  # filter for only variants that are available in our dataset
  output_score_file <- combined_330_score_file %>% filter(ID%in%available_variant_list)
  
  # outputting the final score file
  write.table(output_score_file,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/IMPUTE_330/SCORE_FILES/EUR_330_",subtype),quote=F,row.names=F,col.names=F,sep="\t")
}