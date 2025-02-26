#!/bin/bash
#SBATCH --job-name=bcftools_query_variants
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_query_variants_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_query_variants_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=128gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load bcftools/1.17
module load samtools/1.17

# specifying tmp directory
export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

# setting file paths
i=${ARGS1}

cd /scratch/jll1/AABCG_PROCESSING/merge_chr_pass_001/chr_variant_ID_list
echo ${i}
bcftools view ../chr${i}.vcf.gz | bcftools query -f "%ID\n" > chr${i}_variant_ID_list.txt


