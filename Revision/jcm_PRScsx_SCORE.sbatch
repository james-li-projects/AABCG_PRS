#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --partition=express

# loading modules 
module load gcc/12.1.0
module load miniconda3
source activate /gpfs/data/huo-lab/jmcclellan/software/envs/r_mr
module load plink/2.0

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jmcclellan/tmp

pfile_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/split_testing_${subtype}
input_score_path=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx
input_score_prefix=${ancestry}_${subtype}_${phi}.txt
output_score_path=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_scores
output_score_prefix=testing_${ancestry}_${subtype}_${phi}

# generating PRS scores for the testing and validation sets
plink2 --pfile ${pfile_path} --score ${input_score_path}/${input_score_prefix} --threads 1 --memory 4000 --out ${output_score_path}/${output_score_prefix}

