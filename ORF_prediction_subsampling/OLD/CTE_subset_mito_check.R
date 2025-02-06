library(ORFik)
out_dir <- "~/livemount/shared_results/ribocrypt_paper/mito_check/"
all_types <- list(c("dORF", "doORF"), c("uORF", "uoORF"), "NTE", "CTE")
min_lengths <- c(15*3, 15*3, 15, 15)
species <- c("all_merged-Homo_sapiens", "all_merged-Saccharomyces_cerevisiae")
s <- species[1]
types <- all_types[3]
for (s in species) {
  message(s)
  df <- read.experiment(s)
  ref_dir <- dirname(df@fafile)
  translon_dir <- file.path(ref_dir, "predicted_translons")
  if (s == "all_merged-Homo_sapiens") translon_dir <- file.path(translon_dir, "mane_isoforms_only")

  table <- setDT(fst::read_fst(file.path(translon_dir, "predicted_translons_with_sequence.fst")))
  grl <- readRDS(file.path(translon_dir, "predicted_translons_with_sequence_ranges.rds"))
  stopifnot(nrow(table) == length(grl))
  # Define filter
  for (types in all_types) {
    message(types[1])
    predicted <- table$predicted == TRUE
    type_hits <- table$type %in% types
    length_valid <- table$length > 15*3
    # Subset
    subset <- predicted & type_hits & length_valid
    if (sum(subset) == 0) {
      message("No hits, next")
      next
    }
    length(grl)
    grl_subset <- grl[subset]
    length(grl_subset)
    dt <- table[subset,]
    dir.create(file.path(out_dir, s), showWarnings = FALSE, recursive = TRUE)
    fst::write_fst(dt, file.path(out_dir, s, paste0(types[1], "_for_mito_check.fst")))
    saveRDS(grl_subset, file.path(out_dir, s, paste0(types[1], "_for_mito_check_ranges.rds")))
  }
}

