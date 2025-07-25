#!/bin/bash
#SBATCH --job-name=LDpred2_lassosum2
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier1q,tier2q,tier3q

# loading modules 
source /etc/profile
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/huo-lab/jmcclellan/software/envs/r_mr

# importing argument
subtype=${ARGS1}
snp_subset=${ARGS2}

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jmcclellan/james_ldpred2/

# running LDpred2_lassosum2
Rscript /gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_LDpred2_lassosum2.R ${subtype} ${snp_subset}
