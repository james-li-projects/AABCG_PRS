#!/bin/bash
#SBATCH --job-name=SCORE_PRS
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/SCORE_PRS_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/SCORE_PRS_%A.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/2.0

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# importing arguments and setting paths
subtype=${ARGS1}
score_prefix=${ARGS2}
input_score_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_FILES_PROCESSED
output_score_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/SCORE_OUTPUT/${subtype}
pfile_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/pfile/combined_testing_validation_${subtype}

# generating PRS scores for the testing and validation sets
plink2 --pfile ${pfile_path} --score ${input_score_path}/${score_prefix} --out ${output_score_path}/${score_prefix}
