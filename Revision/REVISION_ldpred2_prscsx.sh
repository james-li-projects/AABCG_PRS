#########################################################
# code to process revisions to PRS-CSx and LDpred2-grid #
#########################################################
cd /gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_combined
cd /gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/SCORE_FILES

# copying prscsx files into scoring file directory
for file in /gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/PRScsx_combined/PRScsx_*.tsv; do
basename=$(basename "$file")
subtype=$(echo "$basename" | cut -d'_' -f2)
suffix=$(echo "$basename" | cut -d'_' -f3)
newname="${subtype}_PRScsx_${suffix}"
cp "$file" "/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/REVISIONS/ldpred2_prscsx/SCORE_FILES/$newname"
done

# copying ldpred2-grid files into scoring file directory
cp /gpfs/data/huo-lab/jmcclellan/james_paper_additional/output/ldpred2/SCORE_FILES/* /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/REVISIONS/ldpred2_prscsx/SCORE_FILES/
  
# renaming all + to - in file strings
cd /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/REVISIONS/ldpred2_prscsx/SCORE_FILES/
for f in *+*; do mv "$f" "${f//+/-}"; done
# removing .tsv
for f in *.tsv; do mv "$f" "${f%.tsv}"; done

# scoring files
cd /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/output/REVISIONS/ldpred2_prscsx/SCORE_FILES/
for subtype in OVERALL ERPOS ERNEG TNBC
do
for score_file in ${subtype}_*
  do
sbatch --export=ARGS1=${subtype},ARGS2=${score_file} /gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/code/sbatch/REVISION_SCORE_PRS_ldpred2_prscsx.sh
done
done
