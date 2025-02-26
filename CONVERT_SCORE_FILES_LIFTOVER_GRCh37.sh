#!/bin/bash
#SBATCH --job-name=CONVERT_SCORE_FILES_LIFTOVER_GRCh37
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CONVERT_SCORE_FILES_LIFTOVER_GRCh37_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CONVERT_SCORE_FILES_LIFTOVER_GRCh37_%A.err
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# importing arguments and setting paths
PRS_MODEL=${ARGS1}

# generating liftOver scoring files
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/CONVERT_SCORE_FILES_LIFTOVER_GRCh37.R ${PRS_MODEL}
