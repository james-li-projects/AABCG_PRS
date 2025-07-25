library(readr)
library(dplyr)
library(tidyr)
library(glue)

subtypes <- c("ERNEG", "ERPOS", "OVERALL", "TNBC")
ancestries <- c("EUR", "AFR")
phis <- c("1e+00", "1e-04", "1e-06")

for (subtype in subtypes) {
  # Read and transform dictionary
  dict_path <- glue("/gpfs/data/huo-lab/BCAC/james.li/PRS_AABCG/input/PRScsx/bim_dictionary_{subtype}.txt")
  dict <- read_table(dict_path, col_names = FALSE)
  dict <- dict %>% separate(X1, into = c("chr", "pos", "a2", "a1"), remove = FALSE, sep = ":")

  for (ancestry in ancestries) {
    for (phi in phis) {
      
      # Combine per-chromosome files
      file_paths <- glue("output/{subtype}/PRScsx_{ancestry}_pst_eff_a1_b0.5_phi{phi}_chr{1:22}.txt")
      score_file <- bind_rows(lapply(file_paths, read_table, col_names = FALSE, show_col_types = FALSE)) %>%
				select(X2,X4,X6)
      
      # Join with dictionary and align effect sizes
      joined_score <- inner_join(score_file, dict, by = c("X2" = "X2"))
      joined_score <- joined_score %>%
        mutate(AlignedEffect = ifelse(a1 == X2, X6, -X6)) %>%
        select(X1, a1, AlignedEffect) %>%
        filter(AlignedEffect != 0) %>%
        arrange(X1)
      
      # Output file
      out_dir <- "output/PRScsx"
      if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
      out_file <- glue("{out_dir}/{ancestry}_{subtype}_{phi}.txt")
      write_tsv(joined_score, out_file, col_names = FALSE)
    }
  }
}
