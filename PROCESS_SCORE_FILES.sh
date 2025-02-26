#!/bin/bash
#SBATCH --job-name=PROCESS_SCORE_FILES
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PROCESS_SCORE_FILES_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/PROCESS_SCORE_FILES_%A.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32gb
#SBATCH --partition=tier2q

i=${ARGS1}

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/2.0

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# process score files
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/PROCESS_SCORE_FILES.R ${i}
