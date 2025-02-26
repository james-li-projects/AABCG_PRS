#!/bin/bash
#SBATCH --job-name=PREPROCESS_CLUMP_GENOME
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PREPROCESS_CLUMP_GENOME_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PREPROCESS_CLUMP_GENOME_%A.err
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/2.0

# importing argument
subtype=${ARGS1}
pthresh=0.3
r2thresh=0.8

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# performing LD clumping using plink2
gwas_sumstats=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/sumstats/training_${subtype}_PC_10.Status.glm.logistic.hybrid

plink2 -pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr --clump-p1 ${pthresh} --clump-p2 ${pthresh} --clump-r2 ${r2thresh} --clump-kb 5000 --clump ${gwas_sumstats} --clump-id-field SNP --clump-p-field P --clump-unphased --threads 1 --memory 100000 --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PLINKCT/PREPROCESS_CLUMP_GENOME/${subtype}/clump_p_${pthresh}_r2_${r2thresh}

