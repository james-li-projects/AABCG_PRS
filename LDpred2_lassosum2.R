library(remotes)
library(data.table)
library(magrittr)
library(bigsnpr)
library(dplyr)
library(pROC)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

set.seed(1)

args = commandArgs(trailingOnly=TRUE)
subtype = args[1]
snp_subset = args[2]

print(paste("PERFORMING ANALYSIS FOR SUBTYPE:",subtype))
print(paste("SNP FILTERING PROCEDURE:",snp_subset))

# set working directory
setwd(paste0("/scratch/jll1/PRS_AABCG/",subtype))

# Read in the phenotype and covariate files
phenotype <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pheno_cov/testing_",subtype,".pheno"))[,c("FID", "IID", "Status")]
covariate <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pheno_cov/testing_",subtype,".cov"),header=T)
pcs <- fread(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pheno_cov/testing_",subtype,".eigenvec"),header=T) %>%
  setnames(., colnames(.), c("FID","IID", paste0("PC",1:10)))
# rename columns
colnames(pcs) <- c("FID","IID", paste0("PC",1:10))
# generate required table
pheno <- merge(phenotype, covariate) %>%
  merge(., pcs)
pheno$Status <- pheno$Status - 1

# Load and transform the summary statistic file
# Read in the summary statistic file
sumstats <- bigreadr::fread2(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/sumstats/",snp_subset,"/",subtype,".sumstats")) 

# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = paste0("/scratch/jll1/PRS_AABCG/",subtype))
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)
snp_readBed(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/",snp_subset,"_split_testing_",subtype,".bed"))
# now attach the genotype object
obj.bigSNP <- snp_attach(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/",snp_subset,"_split_testing_",subtype,".rds"))
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
# calculate LD
for (chr in 1:22) {
  print(paste("CALCULATING LD FOR CHROMOSOME:",chr))
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

# perform LD score regression
print("PERFORMING LD SCORE REGRESSION")
df_beta <- info_snp #[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]

# calculate the null R2
print("CALCULATING NULL R2")
library(fmsb)
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- pheno[fam.order, on = c("FID", "IID")]
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste("PC", 1:10, sep = "", collapse = "+") %>% paste0("Status~Age+PlatformAABC+PlatformAMBER+PlatformBCAC_OncoArray+PlatformGBHS+PlatformMEGA+PlatformROOT+", .) %>% as.formula %>% glm(., data = y, family=binomial) %>% summary
#####null.r2 <- fmsb::NagelkerkeR2(null.model)
# Extract the null and residual deviance from the model summary
null_deviance <- null.model$null.deviance
residual_deviance <- null.model$deviance
# Get the number of observations
n <- nrow(y)
# Calculate Nagelkerke's R^2
R2_nagelkerke <- (1 - exp((residual_deviance - null_deviance) / n)) / (1 - exp(-null_deviance / n))
# Print the result
null.r2 <- R2_nagelkerke
null.r2


##############################################
##############################################
# Get adjusted beta from the auto model
print("OBTAINING ADJUSTED BETA FROM THE AUTO MODEL")
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2 = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
  ncores = NCORES#, 
  #use_MLE=TRUE,
  #allow_jump_sign = FALSE,
  #shrink_corr = coef_shrink
)
beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)

if(is.null(obj.bigSNP)){
  obj.bigSNP <- snp_attach(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/",snp_subset,"_split_testing_",subtype,".rds"))
}

#genotype <- obj.bigSNP$genotypes
print("PERFORMING MEAN IMPUTATION FOR MISSING VALUES")
genotype <- snp_fastImputeSimple(
  Gna = obj.bigSNP$genotypes,
  method = c("mean2"),
  ncores = NCORES
)
# calculate PRS for all samples
print("calculate PRS for all samples")
ind.test <- 1:nrow(genotype)
pred_auto <-
  big_prodMat(genotype,
              beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)
# scale the PRS generated from AUTO
print("scale the PRS generated from AUTO")
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
  rowMeans(beta_auto[,
                     abs(pred_scaled -
                           median(pred_scaled)) <
                       3 * mad(pred_scaled)])
print("COMPUTING PRED AUTO")
pred_auto <-
  big_prodVec(genotype,
              final_beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)
print("OUTPUTTING RELEVANT RESULTS")
# outputting relevant results
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/LDpred2_lassosum2/",subtype))
# append adjusted betas to snp info file and saving it
beta_info_snp <- cbind(info_snp,final_beta_auto)
save(beta_info_snp,file=paste0("beta_info_snp_",snp_subset,"_",subtype,"_auto.RData"))
# appending PRS score file and saving it
reg.dat <- y
reg.dat$PRS_pred_auto <- as.numeric(pred_auto)
save(reg.dat,file=paste0("reg.dat_",snp_subset,"_",subtype,"_auto.RData"))

print("OBTAINING AUC")
# Get the final performance of the LDpred models: auto model
roc_score=pROC::roc(reg.dat$Status, reg.dat$PRS_pred_auto) #AUC score
output_auc <- as.numeric(auc(roc_score))
output_auc_CI <- as.vector(ci.auc(roc_score))
output_auc_CI_lower <- output_auc_CI[1]
output_auc_CI_upper <- output_auc_CI[3]
print(paste("Unadjusted AUC (",subtype,"):", output_auc))
print(paste("Unadjusted AUC 95% CI (",subtype,"):", output_auc_CI_lower, output_auc_CI_upper))
save(output_auc_CI,file=paste0("output_auc_CI_",snp_subset,"_",subtype,"_auto.RData"))

##############################################
##############################################


##############################################
##############################################
# Running lassosum2
beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
(params2 <- attr(beta_lassosum2, "grid_param"))

pred_grid2 <- big_prodMat(genotype, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]])
params2$score <- apply(pred_grid2[ind.test, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(glm(y$Status ~ x, family = "binomial"))$coef["x", 3]
  # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

best_grid_lassosum2 <- params2 %>%
  mutate(id = row_number()) %>%
  arrange(desc(score)) %>%
  print() %>% 
  slice(1) %>%
  pull(id) %>% 
  beta_lassosum2[, .]

best_grid_lassosum2_index <- params2 %>%
  mutate(id = row_number()) %>%
  arrange(desc(score)) %>% 
  slice(1) %>%
  pull(id)

# Get the final performance of the lassosum2 model
reg.dat <- y
reg.dat$PRS_lassosum2 <- pred_grid2[,best_grid_lassosum2_index]
roc_score=pROC::roc(reg.dat$Status, reg.dat$PRS_lassosum2) #AUC score
output_auc <- as.numeric(auc(roc_score))
output_auc_CI <- as.vector(ci.auc(roc_score))
output_auc_CI_lower <- output_auc_CI[1]
output_auc_CI_upper <- output_auc_CI[3]
print(paste("Unadjusted AUC (",subtype,"):", output_auc))
print(paste("Unadjusted AUC 95% CI (",subtype,"):", output_auc_CI_lower, output_auc_CI_upper))
save(output_auc_CI,file=paste0("output_auc_CI_",snp_subset,"_",subtype,"_lassosum2.RData"))

# outputting relevant results
setwd(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/LDpred2_lassosum2/",subtype))
# append adjusted betas to snp info file and saving it
snp_betas<-beta_lassosum2[,best_grid_lassosum2_index]
beta_info_snp <- cbind(info_snp,snp_betas)
save(beta_info_snp,file=paste0("beta_info_snp_",snp_subset,"_",subtype,"_lassosum2.RData"))
save(reg.dat,file=paste0("reg.dat_",snp_subset,"_",subtype,"_lassosum2.RData"))
##############################################
##############################################

# remove temporary files generated for LDpred2 and lassosum2 
system_command <- paste0("rm /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/*.bk")
system(system_command)
system_command <- paste0("rm /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/*.rds")
system(system_command)
