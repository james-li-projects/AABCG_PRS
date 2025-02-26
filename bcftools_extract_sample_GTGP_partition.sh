#!/bin/bash
#SBATCH --job-name=bcftools_extract_sample_GTGP_partition
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_extract_sample_GTGP_partition_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_bcftools_extract_sample_GTGP_partition_%A.err
#SBATCH --time=36:00:00
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
input_dataset_chr_vcf_path=/scratch/jll1/AABCG_PROCESSING/PREMERGE_VCF_FILES
output_sample_path=/scratch/jll1/AABCG_PROCESSING/merge_chr_pass_001/sample_list/chr${i}
input_merge_vcf_path=/scratch/jll1/AABCG_PROCESSING/merge_chr_pass_001
output_partition_path=/scratch/jll1/AABCG_PROCESSING/merge_chr_pass_001/partition_raw

###########################
# Extracting sample lists #
###########################
# obtaining a list of samples with genotypes (WGS datasets)
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_01_WGS_chr${i}.vcf.gz > ${output_sample_path}/GTsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_02_WGS2_to_share_extracted_chr${i}.vcf.gz >> ${output_sample_path}/GTsamples.txt
# obtaining a list of samples with genotype probabilities (other datasets)
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_03_MEGA_VANDY_chr${i}.vcf.gz > ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_04_MEGA_RP_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_05_MEGA_USC_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_06_AMBER_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_07_ROOT_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_08_AABC_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_09_GBHS_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_10_BCAC_OncoArray_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt
bcftools query -l ${input_dataset_chr_vcf_path}/AABCG_dataset_11_BCAC_iCOGS_chr${i}.vcf.gz >> ${output_sample_path}/GPsamples.txt

####################################################################
# Generating GT and GP partition files accordingly for the samples #
####################################################################
# extracting GTs for samples with GTs into partitions
bcftools view -S ${output_sample_path}/GTsamples.txt ${input_merge_vcf_path}/chr${i}.vcf.gz | bcftools query -f "%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%GT]\n" | split -l 2500 -d -a 5 - ${output_partition_path}/chr${i}_GT_partition_
# extracting GPs for samples with GPs into partitions
bcftools view -S ${output_sample_path}/GPsamples.txt ${input_merge_vcf_path}/chr${i}.vcf.gz | bcftools query -f "%CHROM\t%ID\t%POS\t%REF\t%ALT[\t%GP]\n" | split -l 2500 -d -a 5 - ${output_partition_path}/chr${i}_GP_partition_ 
