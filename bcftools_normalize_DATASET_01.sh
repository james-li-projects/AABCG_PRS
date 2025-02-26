#!/bin/bash
#SBATCH --job-name=bcftools_normalize_DATASET_01
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_normalize_DATASET_01_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_normalize_DATASET_01_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16gb

module load gcc/12.1.0
module load bcftools/1.17

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

first_path=${ARGS1}
prefix=${ARGS3}

output_path=/scratch/jll1/AABCG_PROCESSING/INTERMEDIATE_DATASET_01_FILES
input_vcf_file=${first_path}/${prefix}.vcf.gz
fasta_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/hg38/hg38.fa

# convert multiallelic to biallelic format, normalization to fix and swap REF & ALT, and export genotype in vcf format 
bcftools norm -m-both ${input_vcf_file} | bcftools norm -f ${fasta_file} -cs | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o ${output_path}/${prefix}.multi.norm.vcf.gz 
