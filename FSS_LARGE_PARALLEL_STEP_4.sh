#!/bin/bash
#SBATCH --job-name=FSS_LARGE_PARALLEL_STEP_4
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/FSS_LARGE_PARALLEL_STEP_4_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/FSS_LARGE_PARALLEL_STEP_4_%A.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb

# loading modules
module load gcc/12.1.0
module load R/4.2.1
module load plink/2.0

# importing argument that was passed
i=${ARGS1}
subtype=${ARGS2}
j=${ARGS3}

# running the FORWARD_STEPWISE_REGRESSION rscript
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/FSS_LARGE_PARALLEL_STEP_4.R ${i} ${subtype} ${j}
