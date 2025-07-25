library(data.table)
library(dplyr)
library(tidyr)

setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/DAMON_VALIDATION/EUR_MODEL_LIFTOVER_HG38")

# importing hg38 330 variant coordinates
initial_330_SNP <- fread("EUR_330_Luminal-A-like")
initial_330_SNP <- initial_330_SNP %>% separate(V1,into=c("CHR","POS","A1","A2"),sep=":",remove=F)

# obtaining 660 variant IDs by flipping the A1 and A2 alleles
initial_330_SNP <- initial_330_SNP %>%
  mutate(ID1=paste(CHR,POS,A1,A2,sep=":")) %>%
  mutate(ID2=paste(CHR,POS,A2,A1,sep=":"))

# assembling 660 variant ID lists
id_list_1 <- initial_330_SNP %>% select(ID1)
colnames(id_list_1) <- c("ID")
id_list_2 <- initial_330_SNP %>% select(ID2)
colnames(id_list_2) <- c("ID")
id_list_all <- rbind(id_list_1,id_list_2)

# writing out the 660 variant ID list
write.table(id_list_all,file="list_EUR660_variants_hg38",quote=F,row.names=F,col.names=F)

# extracting traw for 660 variants
plink_traw_660_command <- paste0("plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr --extract list_EUR660_variants_hg38 --export Av --out /scratch/jll1/AABCG_PROCESSING/EUR330/original_traw/original_660")
system(plink_traw_660_command)

