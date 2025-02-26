#!/bin/bash
#SBATCH --job-name=convert_dosage_pfile
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_convert_dosage_pfile_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_convert_dosage_pfile_%A.err
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=150gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

i=${ARGS1}

input_dosage_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/combined_chr/dosage/chr${i}.dosage.gz
input_fam_file=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/sample_list/chr${i}.fam
output_pfile=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/combined_chr/pfile/chr${i}
plink2 --import-dosage ${input_dosage_file} noheader chr-col-num=1 pos-col-num=3 format=1 skip0=1 skip1=1 --fam ${input_fam_file} --threads 1 --memory 150000 --sort-vars --make-pgen --out ${output_pfile} 
