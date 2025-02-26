#!/bin/bash
#SBATCH --job-name=PRS_MEAN_IMPUTATION
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=30gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/PRS_MEAN_IMPUTATION_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/PRS_MEAN_IMPUTATION_%A.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env
module load plink/2.0

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

cd /scratch/jll1/AABCG_PROCESSING/MEAN_IMPUTATION/INPUT/input_partitionlist

partition_file_list_name=${ARGS1}
for i in `cat ${partition_file_list_name}`
do
    Rscript /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/code/rscripts/PRS_MEAN_IMPUTATION.R ${i}
done

echo "FINISHED PARSING:"
echo ${partition_file_list_name}
