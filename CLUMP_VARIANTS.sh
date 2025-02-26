#!/bin/bash
#SBATCH --job-name=CLUMP_VARIANTS
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CLUMP_VARIANTS_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/CLUMP_VARIANTS_%A.err
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate r_env
module load plink/1.9

# importing argument
subtype=${ARGS1}

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jll1/tmp

# set working directory
cd /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile

# identify clumped variants
plink -bfile /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/split_training_${subtype} --clump-p1 0.3 --clump-p2 0.3 --clump-r2 0.8 --clump-kb 1000 --clump /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/sumstats/training_${subtype}_PC_10.Status.glm.logistic.hybrid --clump-snp-field SNP --clump-field P --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GENOME_WIDE_APPROACHES/bfile/split_training_${subtype}
awk '{print $3,$8}' split_training_${subtype}.clumped | cut -d" " -f1 | tail -n +2 > split_training_${subtype}.clumped.variant.list

# extracting clumped variants from each set and subtype plink bed file
plink -bfile split_training_${subtype} --extract split_training_${subtype}.clumped.variant.list --make-bed --out CLUMPED_split_training_${subtype}
plink -bfile split_testing_${subtype} --extract split_training_${subtype}.clumped.variant.list --make-bed --out CLUMPED_split_testing_${subtype}
plink -bfile split_validation_${subtype} --extract split_training_${subtype}.clumped.variant.list --make-bed --out CLUMPED_split_validation_${subtype}

