#!/bin/bash

subtypes=("OVERALL" "ERPOS" "ERNEG" "TNBC")
phis=("1e+00" "1e-06" "1e-04")

for subtype in "${subtypes[@]}"; do
  for phi in "${phis[@]}"; do
    echo "Running for subtype=$subtype, phi=$phi"
    log_file="/gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs/PRScsx_xances_${subtype}_${phi}.log"
    Rscript /gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_PRScsx_xances_SCORE.R "$subtype" "$phi" &> "$log_file"
  done
done
