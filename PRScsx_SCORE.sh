#!/bin/bash
#SBATCH --job-name=PRScsx_SCORE
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PRScsx_SCORE_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PRScsx_SCORE_%A.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/2.0

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# importing arguments and setting paths
subtype=${ARGS1}
ancestry=${ARGS2}

pfile_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/split_testing_${subtype}
input_score_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PRScsx/ANCESTRY_SCORE_FILES
input_score_prefix=${subtype}_${ancestry}.txt
output_score_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/PRScsx/testing/
output_score_prefix=testing_${ancestry}_${subtype}

# generating PRS scores for the testing and validation sets
plink2 --pfile ${pfile_path} --score ${input_score_path}/${input_score_prefix} --threads 1 --memory 100000 --out ${output_score_path}/${output_score_prefix}
