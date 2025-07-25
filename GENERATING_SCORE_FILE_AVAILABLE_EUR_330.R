library(dplyr)
library(data.table)
library(stringr)
library(pROC)

setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_330")
intrinsic_subtype_list <- unlist(str_split(list.files(),pattern="\\."))
intrinsic_subtype_list <- intrinsic_subtype_list[!grepl("txt",intrinsic_subtype_list)]
intrinsic_subtype_list <- intrinsic_subtype_list[!grepl("ALL_SNP",intrinsic_subtype_list)]

# obtaining ID/A1 backbone of a consistent sumstats file
sumstats_backbone <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/consistent_effect_allele/OVERALL.sumstats"),header=T) %>% select(ID,A1)

# obtaining list of all variants that are available in our genetic dataset
available_variant_list <- fread("/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr.pvar")$ID

# writing out intrinsic subtype scoring files
for (intrinsic_subtype in intrinsic_subtype_list) {
  print(paste("PROCESSING SCORE FILE FOR:",intrinsic_subtype))
  setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/EUR_330")
  PRS_330 <- fread(paste0(intrinsic_subtype,".txt"),header=T)
  colnames(PRS_330) <- c("SNP","CHR","hg19","EFFECT_ALLELE","BASELINE_ALLELE","BETA","EFFECT_ALLELE_AF")
  PRS_330 <- PRS_330 %>% mutate(col4 = seq(nrow(PRS_330))) %>% mutate(col4 = as.integer(col4))
  tmp_liftOver <- PRS_330 %>% mutate(col1 = paste0(rep("chr",nrow(PRS_330)), CHR)) %>% mutate(col2 = hg19) %>% mutate(col2 = as.integer(col2)) %>% mutate(col3 = col2 + 1) %>% select(col1,col2,col3,col4) 
  # writing temporary table out for liftOver
  setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp")
  write.table(tmp_liftOver, file="tmp_liftOver", quote=F,row.names=F,col.names=F,sep="\t")
  # performing liftOver
  system("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/tmp/tmp_liftOver /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/liftOver/hg19ToHg38.over.chain.gz liftOver_mapped liftOver_unmapped")
  # reading in liftOver mapped results
  liftOver_mapped <- fread("liftOver_mapped",header=F,sep="\t")
  colnames(liftOver_mapped) <- c("col1","col2","col3","col4")
  # joining with original table 
  liftOver_PRS_330 <- inner_join(PRS_330, liftOver_mapped, by = c("col4"))
  
  # formatting the output score file
  liftOver_PRS_330$BETA <- as.numeric(liftOver_PRS_330$BETA)
  liftOver_PRS_330$A1 <- liftOver_PRS_330$EFFECT_ALLELE
  
  # formatting liftOver ID
  liftOver_PRS_330$ID <- paste(
    liftOver_PRS_330$CHR,
    liftOver_PRS_330$col2,
    liftOver_PRS_330$BASELINE_ALLELE,
    liftOver_PRS_330$EFFECT_ALLELE,
    sep=":"
  )
  
  # outputting the final score file
  liftOver_PRS_330 <- liftOver_PRS_330 %>% select(ID,A1,BETA)
  #output_score_file <- inner_join(liftOver_PRS_330,sumstats_backbone,by=c("ID","A1"))
  output_score_file <- liftOver_PRS_330 %>% filter(ID%in%available_variant_list)
  
  write.table(output_score_file,file=paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_FILES/EUR_330_",intrinsic_subtype),quote=F,row.names=F,col.names=F,sep="\t")
}