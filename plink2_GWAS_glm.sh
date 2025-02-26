#!/bin/bash
#SBATCH --job-name=plink2_GWAS_glm
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/plink2_GWAS_glm_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/plink2_GWAS_glm_%A.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

in_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input
pheno_cov_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input
out_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats

n=${ARGS1}
subtype=${ARGS2}
PC_index_list=`seq 1 ${n}`
PC_list=`for i in $PC_index_list; do echo PC${i}; done`

# running GWAS using plink2 [logistic regression] 
plink2 --pfile ${in_path}/split_training_${subtype} \
    --glm hide-covar \
    --pheno ${pheno_cov_path}/training_${subtype}.pheno_cov \
    --pheno-name Status \
    --covar ${pheno_cov_path}/training_${subtype}.pheno_cov \
    --covar-name Age Platform ${PC_list} \
    --threads 1 \
    --memory 256000 \
    --out ${out_path}/training_${subtype}_PC_${n}
