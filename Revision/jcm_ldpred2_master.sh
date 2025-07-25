mkdir -p /gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs

for snp_subset in hm3; do
  prev_jobid=""
  # for subtype in OVERALL ERPOS ERNEG TNBC; do
  for subtype in OVERALL ERPOS ERNEG TNBC; do
    outfile="/gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/${subtype}_${snp_subset}_params_and_auc.tsv"
    if [ ! -f "$outfile" ]; then
      if [ -z "$prev_jobid" ]; then
        jobid=$(sbatch --parsable \
               --job-name=LDpred2_${subtype}_${snp_subset} \
               --output=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs/LDpred2_${subtype}_${snp_subset}_%j.log \
               --export=ARGS1=${subtype},ARGS2=${snp_subset} \
               /gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_ldpred2.sbatch)
      else
        jobid=$(sbatch --parsable \
               --dependency=afterok:${prev_jobid} \
               --job-name=LDpred2_${subtype}_${snp_subset} \
               --output=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs/LDpred2_${subtype}_${snp_subset}_%j.log \
               --export=ARGS1=${subtype},ARGS2=${snp_subset} \
               /gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_ldpred2.sbatch)
      fi
      echo "Submitted ${subtype} ${snp_subset} as job ${jobid}"
      prev_jobid="$jobid"
    else
      echo "Skipping ${subtype} ${snp_subset} - output already exists."
    fi
  done
done
