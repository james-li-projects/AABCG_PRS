#!/bin/bash
#SBATCH --job-name=plink2_compute_pca
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_compute_pca.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_compute_pca.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

# generating PCs
plink2 --pfile /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/PRUNED \
    --pca 50 \
    --threads 1 \
    --memory 256000 \
    --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/PCA/AABCG_PCA_50PCs
