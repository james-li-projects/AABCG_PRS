#!/bin/bash
#SBATCH --job-name=PREPARE_XANCESTRY_SUMSTATS
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PREPARE_XANCESTRY_SUMSTATS_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PREPARE_XANCESTRY_SUMSTATS_%A.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=460gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/2.0

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/PREPARE_XANCESTRY_SUMSTATS_EUR_SUBTYPES.R
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/PREPARE_XANCESTRY_SUMSTATS_EUR_TNBC.R
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/PREPARE_XANCESTRY_SUMSTATS_LIMIT_PVAL.R
