#!/bin/bash
#SBATCH --job-name=CT-SLEB_RUN
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CT-SLEB_RUN_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CT-SLEB_RUN_%A.err
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100gb

# importing subtype argument
subtype=${ARGS1}

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/1.9

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# running CT-SLEB
Rscript /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/rscripts/CT-SLEB_RUN_FEWER_PARAM.R ${subtype} >& /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/CT-SLEB_FEWER_PARAM/${subtype}.log
