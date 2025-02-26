#!/bin/bash
#SBATCH --job-name=PARSE_DOSAGE
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/PARSE_DOSAGE_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/PARSE_DOSAGE_%A.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load miniconda3/23.1.0
source activate r_env 

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

cd /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/partition_filelist_forplink

partition_file_list_name=${ARGS1}
for i in `cat ${partition_file_list_name}`
do
    Rscript /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/code/rscripts/PARSE_DOSAGE.r ${i}
done

echo "FINISHED PARSING:"
echo ${partition_file_list_name}
