library(ORFik)
out_dir <- "~/livemount/shared_results/ribocrypt_paper/mito_check/"
all_types <- list(c("dORF", "doORF"), c("uORF", "uoORF"), "NTE", "CTE")
min_lengths <- c(15*3, 15*3, 15, 15)
# species <- c("all_merged-Homo_sapiens", "all_merged-Saccharomyces_cerevisiae")
ref_dir <- ORFik::config()["ref"]
all_exp <- list.experiments(validate = FALSE, pattern = "all_merged-", libtypeExclusive = "RFP", BPPARAM = SerialParam())
all_exp <- all_exp[libtypes == "RFP"]
organisms <- all_exp$name
organisms <- organisms[grep("_modalities$|HEK293|_04|Homo_sapiens", organisms, invert = TRUE)]
org <- "all_merged-Homo_sapiens"
organisms <-  c(org, organisms)

s <- organisms[1]
types <- all_types[3]
for (s in organisms) {
  message(s)
  df <- read.experiment(s)
  ref_dir <- dirname(df@fafile)
  translon_dir_sORFs <- translon_dir <- file.path(ref_dir, "predicted_translons")
  if (s == "all_merged-Homo_sapiens") translon_dir_sORFs <- file.path(translon_dir, "mane_isoforms_only")
  table_sORF <- setDT(fst::read_fst(file.path(translon_dir_sORFs, "predicted_translons_with_sequence.fst")))
  grl_sORF <- readRDS(file.path(translon_dir_sORFs, "predicted_translons_with_sequence_ranges.rds"))

  table_cte <- setDT(fst::read_fst(file.path(translon_dir, "ORFik", "CTE_candidates", "CTE_candidates.fst")))
  table_nte <- setDT(fst::read_fst(file.path(translon_dir, "ORFik", "NTE_candidates", "NTE_candidates.fst")))
  table_cte[, type := "CTE"]
  table_nte[, type := "NTE"]
  table_cte[, length := ORF_length_nt]
  table_nte[, length := ORF_length_nt]
  grl_cte <- readRDS(file.path(translon_dir, "ORFik", "CTE_candidates", "CDS_CTE_candidates_ranges.rds"))
  grl_nte <- readRDS(file.path(translon_dir, "ORFik", "NTE_candidates", "CDS_NTE_candidates_ranges.rds"))

  stopifnot(nrow(table_sORF) == length(grl_sORF))
  stopifnot(nrow(table_cte) == length(grl_cte))
  stopifnot(nrow(table_nte) == length(grl_nte))
  # Define filter
  for (types in all_types) {
    message("- ", types[1])
    if (all(types %in% c("dORF", "doORF", "uORF", "uoORF"))) {
      table <- table_sORF
      grl <- grl_sORF
    } else if (types %in% "CTE") {
      table <- table_cte
      grl <- grl_cte
    } else if (types %in% "NTE") {
      table <- table_nte
      grl <- grl_nte
    } else stop("Type not found")

    predicted <- table$predicted == TRUE
    type_hits <- table$type %in% types
    length_valid <- table$length > 15*3
    # Subset
    subset <- predicted & type_hits & length_valid
    message("-- total_predicted (Subset size): ", sum(predicted), "(", sum(subset), ")")
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
list.files(dirname(result_dir), "[NTE|CTE|INTRON]_predicted_ranges\\.rds", recursive = TRUE, full.names = T)
