#!/bin/bash
#SBATCH --job-name=bcftools_concat_DATASET_01_SNP_indel
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_concat_DATASET_01_SNP_indel_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_concat_DATASET_01_SNP_indel_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=80gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load bcftools/1.17
module load samtools/1.17

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

i=${ARGS1}
in_path=/scratch/jll1/AABCG_PROCESSING/INTERMEDIATE_DATASET_01_FILES
out_path=/scratch/jll1/AABCG_PROCESSING/PREMERGE_VCF_FILES
fasta_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/hg38/hg38.fa

# concatenating the SNP and INDEL variants for DATASET_01
bcftools concat \
    ${in_path}/AABCG_dataset_01_WGS_SNP_chr${i}.multi.norm.vcf.gz \
    ${in_path}/AABCG_dataset_01_WGS_indel_chr${i}.multi.norm.vcf.gz | bcftools sort -Oz -o ${out_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz

# indexing the concatenated vcfs for each chromosome
bcftools index -t ${out_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz
