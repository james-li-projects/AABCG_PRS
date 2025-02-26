#!/bin/bash
#SBATCH --job-name=plink2_LD_prune
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_LD_prune.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_LD_prune.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

# identifying variants that have been pruned out 
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr \
    --indep-pairwise 500kb 0.2 \
    --threads 1 \
    --memory 256000 \
    --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/pruned

# obtaining a pfile that has been pruned
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/FILTERED_all_combined_chr \
    --extract /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/pruned.prune.in \
    --threads 1 \
    --memory 256000 \
    --make-pgen \
    --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/PRUNED
    
