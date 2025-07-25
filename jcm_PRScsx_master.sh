for phi in 1e-6 1e-4 1 # 1e-2 already run by James
do
  for subtype in OVERALL ERPOS ERNEG TNBC
  do
    for i in $(seq 1 22)
    do
      echo "Running PRScsx for ${subtype}, chromosome ${i}, phi ${phi}"
      sbatch --job-name=PRScsx_${subtype}_${i}_${phi} \
             --output=/gpfs/data/huo-lab/jmcclellan/james_paper_additional/logs/PRScsx_${subtype}_chr${i}_phi${phi}_%j.log \
             --export=ARGS1=${i},ARGS2=${subtype},ARGS3=${phi} \
             /gpfs/data/huo-lab/jmcclellan/james_paper_additional/jcm_PRScsx_RUN.sbatch
    done
  done

done
