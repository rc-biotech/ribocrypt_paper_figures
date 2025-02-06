# Setup
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")

github_forks_repo <- "~/livemount/forks"
conservation_dir <- file.path(github_forks_repo, "translon-conservation/")
conservation_raw_dir <- file.path(conservation_dir, "data", "raw")

ref_dir <- ORFik::config()["ref"]
all_exp <- list.experiments(validate = FALSE, pattern = "all_merged-", libtypeExclusive = "RFP", BPPARAM = SerialParam())
all_exp <- all_exp[libtypes == "RFP"]

organisms <- all_exp$name
organisms <- organisms[grep("_modalities$|HEK293|_04|Homo_sapiens", organisms, invert = TRUE)]
remake <- F
org <- "all_merged-Homo_sapiens_04_oct_2024_all"
organisms <- c(org, organisms)
# organisms <- organisms[1]
for (org in organisms) {
  # Load data
  df <- read.experiment(org)
  org_short <- gsub(" ", "_", tolower(organism(df)))
  message("- ", org_short)
  org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik")
  org_dir_out <- result_dir <- file.path(org_dir, "INTRON_candidates")
  dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

  already_done <- file.exists(file.path(org_dir_out, "intron_candidates_whole_intron_equal_size_1_to_7_ranges.qs"))
  has_no_introns <- file.exists(file.path(org_dir_out, "no_introns.rds"))
  if ((already_done & !remake) | has_no_introns) next

  introns <- loadRegion(df, "introns")
  introns <- introns[lengths(introns) > 0]
  if (length(introns) == 0) {
    message("- Species has no introns: ", org_short)
    saveRDS(TRUE, file.path(org_dir_out, "no_introns.rds"))
    next
  }

  subset <- integer()
  if (org_short == "homo_sapiens") {
    mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
    subset <- mane
  } else try(subset <- filterTranscripts(df, NULL, 150, NULL), silent = TRUE)

  RFP <- fimport(filepath(df, "cov")) #if (!exists("RFP", mode = "S4"))
  leaders <- loadRegion(df, "leaders")
  leaders <- leaders[lengths(leaders) == 1] # First intron must be in CDS


  all_cds <- cds <- loadRegion(df, "cds")
  mrna <- loadRegion(df, "mrna")

  # Subset to intronic
  coding_tx <- names(cds)
  name_subset <- coding_tx[coding_tx %in% subset & coding_tx %in% names(introns) & coding_tx %in% names(leaders)]
  intronic_genes <- length(name_subset)
  cds <- cds[names(cds) %in% name_subset]
  mrna <- mrna[names(mrna) %in% name_subset]
  introns <- introns[names(introns) %in% name_subset]

  stopifnot(length(introns) > 0)
  stopifnot(length(cds) == length(mrna)); length(cds)
  stopifnot(length(cds) == length(introns)); length(introns)

  cds_starts_all <- startSites(cds, keep.names = TRUE, is.sorted = TRUE)
  cds_stops_all <- stopSites(cds, keep.names = TRUE, is.sorted = TRUE)

  # Run per intron index
  intron_indices <- seq(7)
  intron_index <- intron_indices[1]
  all_indices_dt <- data.table()
  candidates_gr_all <- GRangesList()
  all_introns_gr_all <- GRangesList()
  all_orfs_gr_all <- GRangesList()
  for (intron_index in intron_indices) {
    message(intron_index)
    introns_first_intron <- subset_to_group_index(introns, index = intron_index, cds = cds)
    if (length(introns_first_intron) == 0) next

    phase <- mcols(introns_first_intron)$phase
    cds_exon_length <- mcols(introns_first_intron)$cds_exon_length
    intron_starts <- startSites(introns_first_intron, TRUE, TRUE, TRUE)

    windows_all <- windowPerGroup(intron_starts, tx = introns_first_intron, upstream = -(phase), downstream = 300)
    summary(widthPerGroup(windows_all))

    orfs_gr <- CTE_orfs(windows_all, df)
    names(orfs_gr) <- paste0(sub("_[0-9]$", "", names(orfs_gr)), "_", intron_index)
    orfs_gr@unlistData$intron_index <- intron_index
    length(orfs_gr)
    tx_names <- names(unlistGrl(orfs_gr))

    introns_first_intron <- introns_first_intron[tx_names]
    introns_first_intron@unlistData$intron_index <- intron_index
    names(introns_first_intron) <- paste0(names(introns_first_intron), "_", intron_index)
    # Detect coverage
    message("-- Coverage on stop codon extensions")
    orfScores <- coverage_statistics(orfs_gr, RFP)

    overlap_filter <- orfs_gr %over% all_cds
    filter <- orfScores$ORFScores > 5 & orfScores$frame_zero_RP > 300 &
      orfScores$ORF_F0_codons_covered > 3 & orfScores$`ORF_F0_codons_covered(%)` > 10 & !overlap_filter

    # TODO check if!
    if (sum(orfScores[filter]$ORF_sum) == 0) final <- data.table()
    final <- append_gene_ids(orfs_gr, df, orfScores, tx_names)
    final <- orf_to_cds_coverage_statistcs(cds[final$tx_ids], RFP, final)

    N_introns_filtered <- sum(filter)
    N_introns_gt30nt <- mcols(orfs_gr)[1,]$width_gt_min
    N_introns_has_stop <- length(orfs_gr)
    N_introns_codon_1_is_stop <- mcols(orfs_gr)[1,]$codon_1_is_stop
    N_introns <- length(introns_first_intron)
    intron_specific_stats <- cbind(intron_index, N_introns, N_introns_gt30nt, N_introns_codon_1_is_stop, N_introns_has_stop,
                                   N_introns_filtered, intron_length = widthPerGroup(introns_first_intron[names(orfs_gr)]))
    final <- cbind(intron_specific_stats, final, predicted = filter)

    message("-- Final candidates: ", nrow(final))
    final_predicted <- final[filter]
    message("-- Final Predicted: ", nrow(final_predicted))

    all_indices_dt <- rbindlist(list(all_indices_dt, final))
    candidates_gr_all <- c(candidates_gr_all, orfs_gr[filter])
    all_orfs_gr_all <- c(all_orfs_gr_all, orfs_gr)
    all_introns_gr_all <- c(all_introns_gr_all, introns_first_intron)
  }
  # Save
  stopifnot(nrow(all_indices_dt) == length(all_orfs_gr_all))
  stopifnot(nrow(all_indices_dt) == length(all_introns_gr_all))
  duplicate_filter <- countOverlaps(all_orfs_gr_all, all_orfs_gr_all, type = "equal") > 1
  all_orfs_gr_all <- all_orfs_gr_all[!duplicate_filter]
  all_introns_gr_all <- all_introns_gr_all[!duplicate_filter]
  all_indices_dt <- all_indices_dt[!duplicate_filter, ]
  stopifnot(nrow(all_indices_dt) == length(all_orfs_gr_all))
  predicted_gr <- candidates_gr_all[!(countOverlaps(candidates_gr_all, candidates_gr_all, type = "equal") > 1)]
  not_predicted_gr_all <-  all_orfs_gr_all[countOverlaps(all_orfs_gr_all, predicted_gr, type = "equal") == 0]
  stopifnot(length(all_orfs_gr_all) == length(predicted_gr) + length(not_predicted_gr_all))

  # The quantile and size correct nonpredicted set
  not_predicted_gr <-
    try(sample_to_matching_lengths_quantiles_and_size(predicted_gr, not_predicted_gr_all), silent = TRUE)
  if (is(not_predicted_gr, "try-error")) {
    message("Could not sample negative set, use set of size 0!")
    not_predicted_gr <- GRangesList()
  }
  if (org_short == "homo_sapiens") {
    # Make sure we use same for human
    not_predicted_gr <-
      sample_to_matching_lengths_quantiles_and_size(predicted_gr, not_predicted_gr_all)
    not_predicted_gr <- not_predicted_gr_all[fread("~/livemount/forks/translon-conservation/data/raw/introns_not_predicted.bed12")$V4]
  }


  # The whole intron and whole intron size correct set
  background_from_predicted_subset <- all_introns_gr_all[names(predicted_gr)]
  background_from_predicted_subset_equal_size <- window_next_to_window(background_from_predicted_subset,
                                                                       predicted_gr, removeEmpty = FALSE)
  summary(widthPerGroup(background_from_predicted_subset_equal_size))

  fst::write_fst(all_indices_dt,  file.path(result_dir, "intron_candidates_intron_1_to_7.fst"))
  saveRDS(predicted_gr, file.path(result_dir, "INTRON_predicted_ranges.rds"))
  qs::qsave(not_predicted_gr_all, file.path(result_dir, "intron_all_not_predicted_1_to_7_ranges.qs"))
  qs::qsave(all_orfs_gr_all, file.path(result_dir, "intron_all_orfs_intron_1_to_7_ranges.qs"))

  file.path(result_dir, c("intron_candidates_intron_1_to_7_ranges.qs", "intron_all_not_predicted_equal_size_1_to_7_ranges.qs",
                          "intron_all_whole_intron_1_to_7_ranges.qs", "intron_candidates_whole_intron_equal_size_1_to_7_ranges.qs"))
  qs::qsave(predicted_gr, file.path(result_dir, "intron_candidates_intron_1_to_7_ranges.qs"))
  qs::qsave(not_predicted_gr, file.path(result_dir, "intron_all_not_predicted_equal_size_1_to_7_ranges.qs"))
  qs::qsave(all_introns_gr_all, file.path(result_dir, "intron_all_whole_intron_1_to_7_ranges.qs"))
  qs::qsave(background_from_predicted_subset_equal_size, file.path(result_dir, "intron_candidates_whole_intron_equal_size_1_to_7_ranges.qs"))

  if (org_short == "homo_sapiens") { # The sets for phylo analysis ->
    message("- Making human bed12 files for phylo scores")
    region <- "introns"
    background <- "intron"
    file_names <- c("_predicted", "_not_predicted",
                    paste0("_predicted_whole_", background),
                    paste0("_predicted_whole_", background, "_equal_size"))
    file_names <- paste0(region, file_names, ".bed12")
    files_to_send <- file.path(result_dir, file_names)
    ORFik::export.bed12(predicted_gr, files_to_send[1])
    ORFik::export.bed12(noncandidates_gr_all_subsample, files_to_send[2])
    ORFik::export.bed12(candidates_whole_introns_gr, files_to_send[3])
    ORFik::export.bed12(candidates_whole_introns_equal_size_gr, files_to_send[4])

    stopifnot(all(file.exists(files_to_send)))
    file.copy(files_to_send, file.path(conservation_raw_dir, file_names))
  }
}

list.files(dirname(result_dir), "CDS_[NTE|CTE|INTRON]_predicted_ranges\\.rds", recursive = TRUE, full.names = T)

# # Intron retention index check
# dt <- data.table(tx_id = txNames(predicted_gr), intron_id = gsub(".*_", "", names(predicted_gr)))
# table(dt$intron_id)
# dt[, id_string := paste(intron_id, sep = "_", collapse = "&"), by = tx_id]
# sort(table(dt$id_string), decreasing = T)[1:15]

# RNA-seq check
# df_modalities <- read.experiment("all_merged-Homo_sapiens_modalities")
# RNA <- fimport(filepath(df_modalities, "cov"))
# all_introns_gr_all <- qs::qread(file.path(result_dir, "intron_candidates_intron_1_to_7_ranges.qs"))
# table(all_introns_gr_all@unlistData$intron_index)
# gr <-all_introns_gr_all@unlistData
# counts <- countOverlapsW(gr, RFP)
# dt <- data.table(index = gr$intron_index, counts, width = as.integer(width(gr)))
# dt[, normalized_counts := (counts / width)*1000]
# plot(dt[, mean(normalized_counts), by = index]$V1)
#
# all_introns_gr_all <- lapply(seq(7), function(intron_index) subset_to_group_index(introns, index = intron_index, cds = cds))
# all_introns_gr_all <- lapply(seq(7), function(intron_index) {res <- all_introns_gr_all[[intron_index]]; res@unlistData$intron_index <- intron_index; res})
# res <- do.call(c, all_introns_gr_all[1:7])
# all_introns_gr_all <- res
# table(all_introns_gr_all@unlistData$intron_index)
# gr <-all_introns_gr_all@unlistData
# counts <- countOverlapsW(gr, RFP)
# dt2 <- data.table(index = gr$intron_index, counts, width = as.integer(width(gr)))
# dt2[, normalized_counts := (counts / width)*1000]
# plot(dt2[, mean(normalized_counts), by = index]$V1)
#
# dt_final <- rbindlist(list(dt, dt2), idcol = T)
# dt_final[, percentage := 100*(normalized_counts / sum(normalized_counts)), by = .(.id, index)]
# dt_final[, zscore := (normalized_counts - mean(normalized_counts)), by = .(.id, index)]
# dt_final[, Set := as.factor(ifelse(.id == 1, "Predicted retianed intron", "All Introns"))]
#
# ggplot(dt_final, aes(as.factor(index), zscore, fill = Set)) +
#   geom_boxplot(position = position_dodge(width = 0.75), outliers = FALSE) +  # Dodge boxes
#   scale_fill_manual(values = c("blue", "red")) +  # Custom colors
#   theme_minimal() + ggtitle("RNA-seq coverage (zscore normalized) by intron Index", "Coverage normalized by length before zscore")
