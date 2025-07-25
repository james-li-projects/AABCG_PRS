#!/bin/bash

subtypes=("OVERALL" "ERNEG" "ERPOS" "TNBC")
phis=("1e+00" "1e-04" "1e-06")
ancestries=("EUR" "AFR")

log_dir="/gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs"
sbatch_script="/gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_PRScsx_SCORE.sbatch"

mkdir -p "$log_dir"

for subtype in "${subtypes[@]}"; do
  for phi in "${phis[@]}"; do
    for ancestry in "${ancestries[@]}"; do

      jobname="${subtype}_${ancestry}_${phi}"
      sbatch --job-name="PRScsx_SCORE_$jobname" \
             --output="${log_dir}/${jobname}.log" \
             --export=ALL,subtype=${subtype},phi=${phi},ancestry=${ancestry} \
             "$sbatch_script"

      sleep 1

    done
  done
done
