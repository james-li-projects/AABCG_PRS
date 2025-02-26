#!/bin/bash
#SBATCH --job-name=munge_sumstats_compute_h2
#SBATCH --output=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/munge_sumstats_compute_h2_%A.out
#SBATCH --error=/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/logs/munge_sumstats_compute_h2_%A.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=100gb
#SBATCH --partition=tier2q

module load gcc/12.1.0
module load miniconda3/23.1.0

source activate ldsc 

sumstats=${ARGS1}
sumstats_prefix=`echo $sumstats | rev | cut -d"." -f2- | rev`
echo ${sumstats}

if grep -q "GENOME_WIDE_SIG" "$sumstats"; then
  # munge sumstats
/gpfs/data/huo-lab/BCAC/james.li/conda_env/ldsc/bin/python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/munge_sumstats.py \
  --sumstats /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/sumstats/${sumstats} \
  --ignore beta \
  --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/munge_sumstats/${sumstats_prefix}
fi

if ! grep -q "GENOME_WIDE_SIG" "$sumstats"; then
  # munge sumstats
/gpfs/data/huo-lab/BCAC/james.li/conda_env/ldsc/bin/python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/munge_sumstats.py \
  --sumstats /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/sumstats/${sumstats} \
  --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/munge_sumstats/${sumstats_prefix}
fi


# computing heritability
# /gpfs/data/huo-lab/BCAC/james.li/conda_env/ldsc/bin/python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/ldsc.py --h2 /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/munge_sumstats/${sumstats_prefix}.sumstats.gz --ref-ld-chr /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/ldsc/ --w-ld-chr /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/ldsc/ --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/output/h2_${sumstats}

# computing heritability using AABCG LD reference panel
/gpfs/data/huo-lab/BCAC/james.li/conda_env/ldsc/bin/python /gpfs/data/huo-lab/BCAC/james.li/preconfluence/software/ldsc/ldsc.py --h2 /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/munge_sumstats/${sumstats_prefix}.sumstats.gz --ref-ld-chr /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/ --w-ld-chr /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/data/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/ --out /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/h2/output/h2_${sumstats}






