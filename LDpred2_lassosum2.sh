#!/bin/bash
#SBATCH --job-name=LDpred2_lassosum2
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/LDpred2_lassosum2_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/LDpred2_lassosum2_%A.err
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/1.9

# importing argument
subtype=${ARGS1}
snp_subset=${ARGS2}

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# running LDpred2_lassosum2
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/LDpred2_lassosum2.R ${subtype} ${snp_subset} >& /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/LDpred2_lassosum2/${subtype}/log_${snp_subset}.txt
