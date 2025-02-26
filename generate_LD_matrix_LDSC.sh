#!/bin/bash
#SBATCH --job-name=generate_LD_matrix_LDSC
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/generate_LD_matrix_LDSC_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/generate_LD_matrix_LDSC_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb

module load gcc/12.1.0
module load miniconda3
source activate ldsc

i=${ARGS1}

cd /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2

bfile=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/AFR_REF_${i}

python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/ldsc.py --bfile ${bfile} --l2 --ld-wind-kb 1000 --out ./ldsc/${i}

