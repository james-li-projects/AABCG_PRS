#!/bin/bash
#SBATCH --job-name=bcftools_merge_ALL_DATASET
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_merge_ALL_DATASET_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_merge_ALL_DATASET_%A.err
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=480gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load bcftools/1.17

# setting number of processes 
ulimit -n 10000
# printing out this set number of processes
ulimit -n 

i=${ARGS1}

export TMPDIR=/scratch/jll1/AABCG_PROCESSING/tmp/chr${i}
in_path=/scratch/jll1/AABCG_PROCESSING/PREMERGE_VCF_FILES
out_path=/scratch/jll1/AABCG_PROCESSING/merge_chr_pass_001

# merging vcfs, creating specific variant IDs, removing duplicate variants 
bcftools merge -m none \
    ${in_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_02_WGS2_to_share_extracted_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_03_MEGA_VANDY_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_04_MEGA_RP_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_05_MEGA_USC_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_06_AMBER_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_07_ROOT_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_08_AABC_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_09_GBHS_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_10_BCAC_OncoArray_chr${i}.vcf.gz \
    ${in_path}/AABCG_dataset_11_BCAC_iCOGS_chr${i}.vcf.gz -Oz -o ${out_path}/chr${i}.vcf.gz

# indexing this merged/processed vcf
bcftools index -t ${out_path}/chr${i}.vcf.gz
