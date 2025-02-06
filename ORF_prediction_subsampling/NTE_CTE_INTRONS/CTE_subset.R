source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")
github_forks_repo <- "~/livemount/forks"
conservation_dir <- file.path(github_forks_repo, "translon-conservation/")
conservation_raw_dir <- file.path(conservation_dir, "data", "raw")

ref_dir <- ORFik::config()["ref"]
all_exp <- list.experiments(validate = FALSE, pattern = "all_merged-", libtypeExclusive = "RFP")
all_exp <- all_exp[libtypes == "RFP"]
organisms <- all_exp$name
organisms <- organisms[grep("_modalities$|HEK293|_04|Homo_sapiens", organisms, invert = TRUE)]
remake <- F
org <- "all_merged-Homo_sapiens_04_oct_2024_all"
organisms <- c(org, organisms)
organisms <- organisms[-1]
for (org in organisms) {
  df <- read.experiment(org, validate = FALSE)
  org_short <- gsub(" ", "_", tolower(organism(df)))
  message("- ", org_short)
  org_dir <- result_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik", "CTE_candidates")
  if (file.exists(file.path(org_dir, "CDS_CTE_candidates_ranges.rds")) & !remake) next
  if (!dir.exists(org_dir)) dir.create(org_dir, recursive = TRUE)

  RFP <- fimport(filepath(df, "cov"))
  txdb <- loadTxdb(df)
  all_cds <- cds <- loadRegion(txdb, "cds")
  mrna <- loadRegion(txdb, "mrna")
  stopifnot(length(cds) == length(mrna)); length(cds)
  subset <- integer()
  if (org_short == "homo_sapiens") {
    mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
    subset <- mane
  } else try(subset <- filterTranscripts(txdb, NULL, 30, 31, longestPerGene = TRUE), silent = TRUE)
  if (length(subset) < length(mrna)*0.2) {
    if (length(subset) == 0) {
      warning("Organism has no trailers, extending 150nt on mrnas!")
    } else warning("Organism has little trailers, extending 150nt on mrnas!")
    mrna <- extendTrailers(mrna, 150)
    subset <- filterTranscripts(txdb, NULL, 30, NULL, longestPerGene = TRUE)
  }

  cds <- cds[names(cds) %in% subset]
  mrna <- mrna[names(mrna) %in% subset]
  trailers <- loadRegion(df, "trailers")
  cds_stops <- stopSites(cds, TRUE, TRUE, TRUE)
  # cds_stops <- split(cds_stops, names(cds_stops))
  windows_all <- windowPerGroup(cds_stops, tx = mrna, upstream = -1, downstream = 300)
  summary(widthPerGroup(windows_all))

  orfs_gr <- CTE_orfs(windows_all, df)
  length(orfs_gr)
  stopifnot(length(unique(txNames(orfs_gr))) == length(txNames(orfs_gr)))
  names(orfs_gr) <- txNames(orfs_gr)
  message("Coverage statistics and filtering...")
  orfScores <- coverage_statistics(orfs_gr, RFP)

  overlap_filter <- orfs_gr %over% all_cds
  filter <- orfScores$ORFScores > 5 & orfScores$frame_zero_RP > 300 &
    orfScores$ORF_F0_codons_covered > 3 & orfScores$`ORF_F0_codons_covered(%)` > 10 &
    orfScores$frame_bias_relative > 0.5 & !overlap_filter
  orfScores[, predicted := filter]

  final <- append_gene_ids(orfs_gr, df, orfScores)
  final <- orf_to_cds_coverage_statistcs(cds[final$tx_ids], RFP, final)

  message("-- Final candidates: ", nrow(final))
  final_predicted <- final[filter]
  message("-- Final Predicted: ", nrow(final_predicted))

  predicted_gr <- orfs_gr[filter]
  unl_gr <- c(unlistGrl(cds[txNames(orfs_gr)]), unlistGrl(orfs_gr))
  cds_cte <- split(unl_gr, names(unl_gr))
  length(cds_cte)
  cds_cte <- reduce(cds_cte)
  length(cds_cte)

  fwrite(final, file = file.path(org_dir, "CTE_candidates.csv"))
  fst::write_fst(final, file.path(org_dir, "CTE_candidates.fst"))
  saveRDS(orfs_gr, file = file.path(org_dir, "CTE_candidates_ranges.rds"))
  saveRDS(cds_cte, file = file.path(org_dir, "CDS_CTE_candidates_ranges.rds"))
  fwrite(final_predicted, file = file.path(org_dir, "CTE_predicted.csv"))
  fst::write_fst(final_predicted, file.path(org_dir, "CTE_predicted.fst"))
  saveRDS(cds_cte[txNames(predicted_gr)], file = file.path(org_dir, "CDS_predicted_ranges.rds"))
  saveRDS(predicted_gr, file = file.path(org_dir, "CTE_predicted_ranges.rds"))

  if (org_short == "homo_sapiens") { # The sets for phylo analysis ->
    not_predicted <- orfs_gr[!filter & !overlap_filter & !(names(orfs_gr) %in% names(predicted_gr))]
    not_predicted_gr <-
      sample_to_matching_lengths_quantiles_and_size(predicted_gr, not_predicted)
    background_from_predicted_subset <- trailers[unique(names(predicted_gr))]
    background_from_predicted_subset_equal_size <-
      window_next_to_window(extendTrailers(background_from_predicted_subset, 200), predicted_gr)
    stopifnot(length(background_from_predicted_subset_equal_size) == length(predicted_gr))
    region <- "CTEs"
    background <- "trailer"
    file_names <- c("_predicted", "_not_predicted",
                    paste0("_predicted_whole_", background),
                    paste0("_predicted_whole_", background, "_equal_size"))
    file_names <- paste0(region, file_names, ".bed12")
    files_to_send <- file.path(result_dir, file_names)
    ORFik::export.bed12(predicted_gr, files_to_send[1])
    ORFik::export.bed12(not_predicted_gr, files_to_send[2])
    ORFik::export.bed12(background_from_predicted_subset, files_to_send[3])
    ORFik::export.bed12(background_from_predicted_subset_equal_size, files_to_send[4])

    stopifnot(all(file.exists(files_to_send)))
    file.copy(files_to_send, file.path(conservation_raw_dir, file_names))
  }
}
