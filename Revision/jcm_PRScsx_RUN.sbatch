#!/bin/bash
#SBATCH --job-name=PRScsx_RUN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier1q,tier2q,tier3q
#SBATCH --time=48:00:00

##SBATCH --mem-per-cpu=100gb
##SBATCH --partition=tier1q,tier2q,tier3q
##SBATCH --time=48:00:00

# loading modules 
source /etc/profile
cd $SLURM_SUBMIT_DIR
module load gcc/12.1.0
module load miniconda3/23.1.0
source activate /gpfs/data/huo-lab/BCAC/james.li/conda_env/pytorch_gnomix
module load plink/1.9

# importing arguments
i=${ARGS1}
subtype=${ARGS2}
phi=${ARGS3}

echo "chromosome ${i}"
echo "subtype: ${subtype}"
echo "phi: ${phi}"

# specifying the tmp directory to utilize in this job
export TMPDIR=/scratch/jmcclellan/tmp

# specifying number of threads to use 
N_THREADS=1
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

# Initializing paths for PRScsx
PRScsx_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/software/PRScsx/PRScsx.py
reference_path=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/reference_panel/1kg
sumstats_path_AFR=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/parsed_xancestry_sumstats/PRScsx_hg38_sumstats_AFR_${subtype}.tsv
sumstats_path_EUR=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/GWAS_sumstats/parsed_xancestry_sumstats/PRScsx_hg38_sumstats_EUR_${subtype}.tsv
VALIDATION_BIM_PREFIX=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/PRScsx/split_validation_${subtype}
out_path=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/${subtype}
mkdir -p ${out_path}

# identifying sample size for AFR GWAS summary statistics from training sets
subtype_sample_size_add1=`wc -l /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/split_training_${subtype}.psam | cut -d" " -f1`
subtype_sample_size="$((subtype_sample_size_add1-1))"

# Running PRScsx
python ${PRScsx_path} --ref_dir=${reference_path} --bim_prefix=${VALIDATION_BIM_PREFIX} --sst_file=${sumstats_path_AFR},${sumstats_path_EUR} --n_gwas=${subtype_sample_size},247173 --pop=AFR,EUR --phi=${phi} --chrom=${i} --out_dir=${out_path} --out_name=PRScsx --seed=1
