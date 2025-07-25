# load libraries 
library(data.table)
library(dplyr)
library(tidyr)

# importing subtype string
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
input_filename=args[1]
output_filename=args[2]

# importing score file
SCORE_FILE_DF <- fread(input_filename,header=F)

# processing score file
colnames(SCORE_FILE_DF) <- c("V1","V2","V4")
PROCESSED_SCORE_FILE_DF <- data.table(SCORE_FILE_DF %>% dplyr::group_by(V1,V2) %>% summarise(Effect = sum(V4))) 
PROCESSED_SCORE_FILE_DF <- PROCESSED_SCORE_FILE_DF %>% tidyr::separate(V1, into = c("chr","pos","a2","a1"),remove=F,sep=":")
PROCESSED_SCORE_FILE_DF <- PROCESSED_SCORE_FILE_DF %>% mutate(AlignedEffect = ifelse(a1 == V2, Effect, -Effect))
PROCESSED_SCORE_FILE_DF <- PROCESSED_SCORE_FILE_DF %>% select(V1,a1,AlignedEffect) 
PROCESSED_SCORE_FILE_DF <- PROCESSED_SCORE_FILE_DF %>% dplyr::group_by(V1,a1) %>% summarise(FinalEffect = sum(AlignedEffect))
PROCESSED_SCORE_FILE_DF <- data.frame(PROCESSED_SCORE_FILE_DF)
PROCESSED_SCORE_FILE_DF <- PROCESSED_SCORE_FILE_DF %>% filter(FinalEffect != 0)

# output processed scoring file
write.table(PROCESSED_SCORE_FILE_DF,file=output_filename,quote=F,row.names=F,col.names=F,sep="\t")
print("FINISHED PROCESSING SCORE FILE")
