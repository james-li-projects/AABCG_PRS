library(dplyr)
library(data.table)
library(readxl)
library(ggplot2)
library(ROCnReg)

################################
# Identifying eligible samples #
# for overall BC no miss recep #
################################
# initializing libraries
library(readxl)
library(data.table)
library(dplyr)
# importing covariate data
pheno_data <- data.frame(read_excel("/gpfs/data/huo-lab/AABCG/data/AABCG_pheno_clean_Jan2023_share.xlsx"))
# converting breast cancer status to a binary variable
pheno_data <- pheno_data %>% select(`AABCGS_ID...1`,`Dataset...3`,Age_GWAS,Status,ER,PR,HER2,AFR_pro)
pheno_data$Status <- 2-pheno_data$Status
colnames(pheno_data)[1] <- "Sample_Name"
# outputting final covariates df
covariates <- pheno_data %>% select(Sample_Name,Status,ER,PR,HER2) %>% rename(case.control=Status)
# recoding missing receptor status in cases and controls
covariates <- covariates %>% mutate(ER=ifelse(ER==8,NA,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==8,NA,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==8,NA,HER2))
covariates <- covariates %>% mutate(ER=ifelse(ER==9,888,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==9,888,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==9,888,HER2))
# making receptor negative a 0 value
covariates <- covariates %>% mutate(ER=ifelse(ER==2,0,ER))
covariates <- covariates %>% mutate(PR=ifelse(PR==2,0,PR))
covariates <- covariates %>% mutate(HER2=ifelse(HER2==2,0,HER2))
# obtaining list of samples with no missing receptor data
eligible_samples <- (covariates %>% filter(ER %in% c(1,0,NA)) %>% filter(PR %in% c(1,0,NA))  %>% filter(HER2 %in% c(1,0,NA)))$Sample_Name
################################


# set working directory
setwd("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/ARISK_FAMH_SCORE_OUTPUT")

# define the absolute risk function
calculate_lifetime_risk <- function(sigma, ShortTime, rate_inc, rate_mor, Ystart, perc.u, perc.v, subtype) {
  
  # making a neat subtype string
  subtype_df <- data.frame(Neat=c("Overall","ER-positive","ER-negative","TNBC"),Crude=c("OVERALL","ERPOS","ERNEG","TNBC"))
  subtype_str <- (subtype_df %>% filter(Crude==subtype))$Neat
  
  prs_freq=perc.v - perc.u 
  
  OR <- (0.6-0.4)*(pnorm(qnorm(1-perc.u)+sigma) - pnorm(qnorm(1-perc.v)+sigma)) / 
    (perc.v-perc.u) / (pnorm(qnorm(0.6)+sigma) - pnorm(qnorm(0.4)+sigma)) 
  
  ORstring<- cbind.data.frame(g=c(1:length(OR)), 
                              prs_cat=paste0(perc.u*100, "-", perc.v*100, "%"), 
                              prs_freq, 
                              OR)
  print(paste0("Odds ratio per percentiles:"))
  print(ORstring[,2:4])
  
  Row = length(OR)
  Col = 81-Ystart
  
  S_g <- matrix(, nrow = Row, ncol = Col)  ## BC-free matrix, row for category, column for age
  S_g[,1] <- rep(1, Row) 
  AR_g_adj <- matrix(, nrow = Row, ncol = Col)  ## Cum risk matrix, row for category, column for age; 
  AR_g_adj[, 1] <- rep(0, Row)
  AR <- matrix(, nrow = Row, ncol = Col)
  AR[, 1] <- rep(0, Row)
  ## Breast cancer hazard of baseline PRS category for each age; 
  ## mu_0 record the hazard within time period t, while AR are for risk by the end of time period and the beginning of time period t+1). 
  mu_0 <- matrix(, nrow = 1, ncol = Col)  
  
  S_m <- matrix(, nrow = Row, ncol = Col)  ## matrix for survival from other diseases 
  mortality <- matrix(, nrow = Row, ncol = Col) ## dying from other diseases
  S_m[,1] <-  rep(1, Row) 
  
  for (t in 1:(Col-1)) {
    i_t <- rate_inc[t+Ystart-20]   ## BC incidence
    mu_0[t] <- solve(sum(prs_freq * OR * S_g[,t]) / sum(prs_freq * S_g[,t]), i_t)
    S_g[,t+1] <- S_g[,t] * (1 - mu_0[t] * OR)
    
    mortality[,t] <- rep(rate_mor[t+Ystart-20], Row) 
    S_m[,t+1] <- S_m[,t] * (1 - mortality[,t])
  }
  
  for (t in 1:(Col-1)){
    AR_g_adj[,t+1] <- mu_0[t] * OR * S_g[,t]*S_m[,t] + AR_g_adj[, t] 
    
    AR[,t+1] <-       mu_0[t] * OR * S_g[,t]           + AR[, t]  
  }
  
  rownames(AR) <- ORstring$prs_cat
  colnames(AR) <- Ystart:80
  AR_t <- t(AR*100)
  print(paste0("Cumulative risk from age ", Ystart, " to 80, without competing risk consideration: "))
  print(AR_t[nrow(AR_t),], digits=3)
  
  rownames(AR_g_adj) <- ORstring$prs_cat
  colnames(AR_g_adj) <- Ystart:80
  ## cumulative risk (%) adjusted for competing risk: by age ... 
  #print(t(AR_g_adj*100), digits=3)
  LT_risk <- t(AR_g_adj*100)
  print(paste0("Cumulative risk from age ", Ystart, " to 80, considering competing risk: "))
  print(LT_risk[nrow(LT_risk),], digits=3)
  
  # plotting function
  plot_LT_risk <- data.frame(t(data.frame(LT_risk)))
  plot_LT_risk <- plot_LT_risk %>% mutate(Percentile=rownames(plot_LT_risk))
  long <- melt(setDT(plot_LT_risk), id.vars = c("Percentile"), variable.name = "Age")
  long <- long %>% 
    mutate(Age = as.numeric(gsub("X","",Age))) %>%
    mutate(Percentile=gsub("X","",Percentile)) %>%
    mutate(Percentile=gsub("\\.","-",Percentile)) %>%
    mutate(Percentile=gsub("-$","%",Percentile))
  
  # defining palette
  c25 <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "darkturquoise",
    "#FF7F00", # orange
    "black",
    "#FB9A99", # lt pink
    "palegreen2",
    "#6A3D9A",
    "#FDBF6F", # lt orange
    "deeppink1","gray70", "khaki2",
    "maroon", "orchid1", "blue1", "steelblue4",
    "green1", "yellow4", "yellow3",
    "darkorange4", "brown"
  )
  
  # plotting risk
  png(paste0("plots/ARISK_",subtype,".png"),units="in",height=4.5,width=4.5,res=1200)
  # specifying palette
  cbPalette <- rev(c25[1:length(unique(long$Percentile))])
  long$Percentile = factor(long$Percentile, levels = c(
    "99-100%",
    "95-99%",
    "90-95%",
    "80-90%",
    "60-80%",
    "40-60%",
    "20-40%",
    "10-20%",
    "5-10%",
    "1-5%",
    "0-1%"
    ))
  p<-ggplot(long,aes(x=Age,y=value,col=Percentile)) + geom_line(linewidth=0.75) + theme(panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), panel.grid.major.y = element_line()) + scale_colour_manual(values = rainbow(11)) + ylab("Absolute risk estimate") + theme_classic() + theme(legend.title=element_blank(), legend.text = element_text(size=9), legend.position = c(0.1,0.60)) + theme(panel.grid.minor.y = element_line(colour="black", size=0.1)) + ggtitle(subtype_str) + theme(plot.title = element_text(hjust = 0.5,size=20),panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.position="none") # + ylim(c(0,1.75*max(long$value)))
  print(p)
  dev.off()
  
  ## Calculate the 5-year risk for women at age ... (rownames). 
  AR_ten <- matrix(, nrow = Row, ncol = (Col-ShortTime))
  for (t in 1:(Col-ShortTime)){
    AR_ten[,t] <- (AR_g_adj[, t+ShortTime] - AR_g_adj[, t]) / (S_g[,t] * S_m[,t])
  }
  rownames(AR_ten) <- ORstring$prs_cat
  colnames(AR_ten) <- Ystart:(80-ShortTime)
  ## For a woman at age (the row label) ..., the 10-year (%) risk is
  Risk_5y <- t(AR_ten*100)
  print(paste0(ShortTime, "-year risk for women aged ", Ystart, ": "))
  print(Risk_5y[1,], digits=3)
  
  # plotting function
  plot_Risk_5y <- data.frame(t(data.frame(Risk_5y)))
  plot_Risk_5y <- plot_Risk_5y %>% mutate(Percentile=rownames(plot_Risk_5y))
  long <- melt(setDT(plot_Risk_5y), id.vars = c("Percentile"), variable.name = "Age")
  long <- long %>% 
    mutate(Age = as.numeric(gsub("X","",Age))) %>%
    mutate(Percentile=gsub("X","",Percentile)) %>%
    mutate(Percentile=gsub("\\.","-",Percentile)) %>%
    mutate(Percentile=gsub("-$","%",Percentile))
  # plotting risk
  png(paste0("plots/TENRISK_",subtype,".png"),units="in",height=4.5,width=4.5,res=1200)
  # specifying palette
  cbPalette <- rev(c25[1:length(unique(long$Percentile))])
  long$Percentile = factor(long$Percentile, levels = c(
    "99-100%",
    "95-99%",
    "90-95%",
    "80-90%",
    "60-80%",
    "40-60%",
    "20-40%",
    "10-20%",
    "5-10%",
    "1-5%",
    "0-1%"
  ))
  if (subtype=="OVERALL") {
    p<-ggplot(long,aes(x=Age,y=value,col=Percentile)) + geom_line(linewidth=0.75) + theme(panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), panel.grid.major.y = element_line()) + scale_colour_manual(values = rainbow(11)) + ylab("10-year absolute risk estimate") + theme_classic() + theme(legend.title=element_blank(), legend.text = element_text(size=9), legend.position = c(0.1,0.60)) + theme(panel.grid.minor.y = element_line(colour="black", size=0.1)) + ggtitle(subtype_str) + theme(plot.title = element_text(hjust = 0.5,size=20),panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.position="none") + geom_hline(yintercept=2, linetype='dotted', col = 'red') # + ylim(c(0,1.75*max(long$value))) 
  } else {
    p<-ggplot(long,aes(x=Age,y=value,col=Percentile)) + geom_line(linewidth=0.75) + theme(panel.border = element_rect(colour = "black", fill=NA), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), panel.grid.major.y = element_line()) + scale_colour_manual(values = rainbow(11)) + ylab("10-year absolute risk estimate") + theme_classic() + theme(legend.title=element_blank(), legend.text = element_text(size=9), legend.position = c(0.1,0.60)) + theme(panel.grid.minor.y = element_line(colour="black", size=0.1)) + ggtitle(subtype_str) + theme(plot.title = element_text(hjust = 0.5,size=20),panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.position="none") # + ylim(c(0,1.75*max(long$value)))
  }
  
  print(p)
  dev.off()
  
  # printing out ages when PRS reaches above 2% 10-year risk
  if (subtype=="OVERALL") {
    print(long %>% filter(Percentile %in% c("99-100%"),value>=2) %>% head(1)) 
    print(long %>% filter(Percentile %in% c("95-99%"),value>=2) %>% head(1))
    print(long %>% filter(Percentile %in% c("90-95%"),value>=2) %>% head(1))
    } else {}
  # returning all summary metrics
  risk_combine <- list(LT_risk, Risk_5y) 
  return(risk_combine)  
}

# importing data for a given ancestry
age_group_data_Black <- read_excel("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/DH_code/BreastCancerIncidence2021.xlsx")
colnames(age_group_data_Black) <- c("Age","OVERALL_INCIDENCE","Total_mortality","BC_mortality","Other_mortality","ERPOS_INCIDENCE","ERNEG_INCIDENCE","TNBC_INCIDENCE")

# calculating life time risk
SUBTYPE_LIST <- c("OVERALL","ERPOS","ERNEG","TNBC")

####################################
# compute sigma for each PRS model #
####################################
# defining all imported PRS model scoring outputs
score_list <- c(
  "OVERALL_PROPTUNE_ALLSUBTYPE.sscore",
  "ERPOS_ENSEMBLE_FSS.sscore",
  "ERNEG_ENSEMBLE_GLMNET.sscore",
  "TNBC_XANCESTRY_PRSICE2.sscore"
)

# initializing vector to store sigma values
SIGMA_LIST <- c()
for (current_PRS_model in score_list) {
  # obtaining subtype string
  subtype <- sub("_.*", "", current_PRS_model)
  # importing score output
  score_validate_df <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CALIBRATION_FINAL_SCORE_OUTPUT/",subtype,"/",current_PRS_model)) %>% select(`#IID`,SCORE)
  # importing validation set covariates
  y_vad <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/validation_",subtype,".pheno_cov"),header=T)
  y_vad$Status <- y_vad$Status-1
  # filtering these for non-missing receptor samples if examining the overall BC risk model 
  if (subtype=="OVERALL") {
    y_vad <- y_vad %>% filter(`#IID` %in% eligible_samples)
  }
  # joining imported PRS scores with validate covariates
  reg_validate_df<-inner_join(score_validate_df,y_vad,by=c("#IID"))
  # computing AUC of PRS
  reg_validate_df$Status <- as.factor(reg_validate_df$Status)
  output_AROC.sp <- AROC.sp(
    formula.h = SCORE~Age+Platform+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,
    group = "Status",
    data = reg_validate_df,
    tag.h = 0)
  print(paste("Covariate-adjusted AUC:", output_AROC.sp$AUC[[1]]))
  # calculating sigma
  current_sigma <- qnorm(output_AROC.sp$AUC[[1]]) * sqrt(2)
  # obtaining sigma list 
  SIGMA_LIST <- c(SIGMA_LIST,current_sigma)
}

# print out sigma list for the subtypes
print(SIGMA_LIST)

# computing lifetime risk and 5/10-year risk 
for (i in 1:4) {
  subtype=SUBTYPE_LIST[i]
  SIGMA=SIGMA_LIST[i]
  print(paste("Evaluating Risk for subtype:",subtype))
  age_black <- cbind(age_group_data_Black[, c("Age",paste0(subtype,"_INCIDENCE"))],mort=age_group_data_Black$Total_mortality-age_group_data_Black$BC_mortality) 
  age_black1 <- age_black[rep(row.names(age_black), each = 5), ]
  age_black1[, 2:3] <- age_black1[, 2:3]/1e5
  perc.low <-c(0, 1, 5, 10, 20, 40, 60, 80, 90, 95, 99)/100
  perc.high <- c(perc.low[-1], 1)
  
  result_black <- calculate_lifetime_risk(sigma=SIGMA, Ystart=20, ShortTime=10,rate_inc= age_black1[,2], rate_mor= age_black1[,3],perc.u= perc.low, perc.v= perc.high,subtype=subtype)
}
