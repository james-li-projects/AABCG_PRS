#!/bin/bash
#SBATCH --job-name=plink2_merge_pfile
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=256gb
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_merge_pfile.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/logs/logs_plink2_merge_pfile.err
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load plink/2.0

# merging pfile into a single pfile 
cd /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/combined_chr/pfile
echo "/gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/combined_chr/pfile/chr1" > pmergelist.txt 
for i in `seq 2 22`
do
    echo /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/combined_chr/pfile/chr${i} >> pmergelist.txt
done
plink2 --pmerge-list pmergelist.txt pfile \
    --threads 1 \
    --memory 256000 \
    --make-pgen \
    --out /gpfs/data/huo-lab/BCAC/james.li/AABCG_PROCESSING/output/all_combined_chr/pfile/all_combined_chr
