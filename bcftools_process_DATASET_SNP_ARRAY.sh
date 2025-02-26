#!/bin/bash
#SBATCH --job-name=bcftools_process_DATASET_SNP_ARRAY
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_process_DATASET_SNP_ARRAY_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_process_DATASET_SNP_ARRAY_%A.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16gb

module load gcc/12.1.0
module load bcftools/1.17

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp

first_path=${ARGS1}
second_path=${ARGS2}
prefix=${ARGS3}

output_path=/scratch/jll1/AABCG_PROCESSING/PREMERGE_VCF_FILES
input_vcf_file=`ls ${first_path}/${second_path}/${prefix}.dose.vcf.gz`
fasta_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/data/hg38/hg38.fa

bcftools view ${input_vcf_file} | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' | bcftools sort -Oz -o ${output_path}/${prefix}.vcf.gz
bcftools index -t ${output_path}/${prefix}.vcf.gz
