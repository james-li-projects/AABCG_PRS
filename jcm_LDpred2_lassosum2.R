library(remotes)
library(data.table)
library(magrittr)
library(bigsnpr)
library(dplyr)
library(pROC)
library(random)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)

set.seed(1)

randString <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

args = commandArgs(trailingOnly=TRUE)
subtype = args[1]
snp_subset = args[2]
# subtype = "OVERALL"
# snp_subset = "hm3"

print(paste("PERFORMING ANALYSIS FOR SUBTYPE:",subtype))
print(paste("SNP FILTERING PROCEDURE:",snp_subset))

# set working directory
workdir <- paste0("/scratch/jmcclellan/james_ldpred2/", subtype)
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
setwd(workdir)

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
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
try({
  # preprocess the bed file (only need to do once for each data set)
  snp_readBed(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/", snp_subset, "_split_testing_", subtype, ".bed"))
}, silent = TRUE)
# now attach the genotype object
obj.bigSNP <- snp_attach(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/",snp_subset,"_split_testing_",subtype,".rds"))
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- runonce::save_run(file=glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/{subtype}/{snp_subset}_info_snp.rds"),
														 	{snp_match(sumstats, map)})
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
dir.create("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/intermediate_corr", recursive = TRUE, showWarnings = FALSE)

# calculate LD
corr <- runonce::save_run(
  file = glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/intermediate_corr/corr_{subtype}_{snp_subset}.rds"),
  {
    for (chr in 1:22) {
      print(paste("CALCULATING LD FOR CHROMOSOME:", chr))
      # Extract SNPs that are included in the chromosome
      ind.chr <- which(info_snp$chr == chr)
      ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]

      # Calculate the LD normally
      corr0 <- snp_cor(
        genotype,
        ind.col = ind.chr2,
        ncores = NCORES,
        infos.pos = POS2[ind.chr2],
        size = 3 / 1000
      )

      if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        tmp <- tempfile(pattern = randString(1), tmpdir = paste0("/scratch/jmcclellan/james_ldpred2/"))
        on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
        corr <- as_SFBM(corr0, tmp)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
    }
		ld <- runonce::save_run(file=glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/intermediate_corr/ld_{subtype}_{snp_subset}.rds"), {ld}
											)
    corr
  }
)

ld <- readRDS(glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/intermediate_corr/ld_{subtype}_{snp_subset}.rds"))

# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

# perform LD score regression
print("PERFORMING LD SCORE REGRESSION")
df_beta <- info_snp #[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]

h2_est <- runonce::save_run(file = glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/{subtype}/{snp_subset}_h2_est.rds"),
	{
	ldsc <- snp_ldsc(ld,
                   length(ld),
                   chi2 = (df_beta$beta / df_beta$beta_se)^2,
                   sample_size = df_beta$n_eff,
                   blocks = NULL)
  ldsc[["h2"]]
	}
)

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
# https://privefl.github.io/bigsnpr/articles/LDpred2.html#ldpred2-grid-grid-of-models
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

# Get adjusted beta from the grid model
print("OBTAINING ADJUSTED BETA FROM THE GRID MODEL")
beta_grid <- runonce::save_run(
  snp_ldpred2_grid(
    corr,
    df_beta,
    params,
    ncores = NCORES
    # use_MLE=TRUE,
    # allow_jump_sign = FALSE,
    # shrink_corr = coef_shrink
  ),
  file = glue::glue("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/beta_grid_{subtype}_{snp_subset}.rds")
)

if(is.null(obj.bigSNP)){
	try({
		# preprocess the bed file (only need to do once for each data set)
		snp_readBed(paste0("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/", snp_subset, "_split_testing_", subtype, ".bed"))
	}, silent = TRUE)
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
library(runonce)

outdir <- "/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pred_grid <- runonce::save_run(
  file = glue::glue("{outdir}/pred_grid_{subtype}_{snp_subset}.rds"),
  {
    big_prodMat(genotype,
                beta_grid,
                ind.row = ind.test,
                ind.col = info_snp$`_NUM_ID_`)
  }
)

NEED_THIS <- TRUE
if (isTRUE(NEED_THIS)){
  # outputting relevant results
  setwd(paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/", subtype))

  # append adjusted betas to snp info file and saving it
  # beta_info_snp <- cbind(info_snp, beta_grid)
  # save(beta_info_snp, file = paste0("beta_info_snp_", snp_subset, "_", subtype, "_grid.RData"))

  # initialize AUC storage columns in params
  params$AUC <- NA
  params$AUC_CI_lower <- NA
  params$AUC_CI_upper <- NA

  # loop over each model (column in pred_grid)
  for (i in seq_len(ncol(pred_grid))) {
    print(paste("PARAMS for model", i, ":"))
    print(params[i, ])

    reg.dat <- y
    reg.dat$PRS_pred_grid <- as.numeric(pred_grid[, i])
    save(reg.dat, file = paste0("reg.dat_", snp_subset, "_", subtype, "_grid_model", i, ".RData"))

		# append adjusted betas to snp info file and saving it
		beta_info_snp <- cbind(info_snp, beta_grid[,i])
		save(beta_info_snp, file = paste0("beta_info_snp_", snp_subset, "_", subtype, "_grid_model", i, ".RData"))


		# Write score files in the way james has liked to format them
    output_score_file <- beta_info_snp[,c(5,4,ncol(beta_info_snp))]
    colnames(output_score_file)[3] <- "beta"
    output_score_file <- output_score_file %>% filter(beta!=0)
    write.table(output_score_file,file=glue::glue("../SCORE_FILES/{subtype}_{snp_subset}_grid_model{i}.tsv"),row.names=F,col.names=F,quote=F,sep="\t")  
    print(paste("OBTAINING AUC for model", i))
    # Get the final performance of the LDpred models: grid model
    roc_score = pROC::roc(reg.dat$Status, reg.dat$PRS_pred_grid) # AUC score
    output_auc <- as.numeric(auc(roc_score))
    output_auc_CI <- as.vector(ci.auc(roc_score))
    output_auc_CI_lower <- output_auc_CI[1]
    output_auc_CI_upper <- output_auc_CI[3]
    print(paste("Unadjusted AUC (", subtype, ", model", i, "):", output_auc))
    print(paste("Unadjusted AUC 95% CI (", subtype, ", model", i, "):", output_auc_CI_lower, output_auc_CI_upper))

    # store in params
    params$AUC[i] <- output_auc
    params$AUC_CI_lower[i] <- output_auc_CI_lower
    params$AUC_CI_upper[i] <- output_auc_CI_upper

    save(output_auc_CI, file = paste0("output_auc_CI_", snp_subset, "_", subtype, "_grid_model", i, ".RData"))
  }

  # write updated params with AUC info to file
  out_path <- paste0("/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/", subtype, "_", snp_subset, "_params_and_auc.tsv")
  readr::write_tsv(params, out_path)
}

# remove temporary files generated for LDpred2 and lassosum2 
system_command <- paste0("rm /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/*.bk")
system(system_command)
system_command <- paste0("rm /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/*.rds")
system(system_command)
