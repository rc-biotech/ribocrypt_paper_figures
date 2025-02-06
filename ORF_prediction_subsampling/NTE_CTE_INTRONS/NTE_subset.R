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
  df <- read.experiment(org, validate = FALSE)
  org_short <- gsub(" ", "_", tolower(organism(df)))
  message("- ", org_short)
  org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik")
  org_dir_in <- file.path(org_dir, "longest_predicted")
  org_dir_out <- result_dir <- file.path(org_dir, "NTE_candidates")
  if (file.exists(file.path(org_dir_out, "CDS_NTE_candidates_ranges.rds")) & !remake) next
  if (!dir.exists(org_dir_out)) dir.create(org_dir_out, recursive = TRUE)

  RFP <- fimport(filepath(df, "cov"))
  txdb <- loadTxdb(df)
  all_cds <- cds <- loadRegion(txdb, "cds")
  mrna <- loadRegion(txdb, "mrna")
  leaders <- loadRegion(df, "leaders")
  stopifnot(length(cds) == length(mrna)); length(cds)
  subset <- integer()
  if (org_short == "homo_sapiens") {
    mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
    subset <- mane
  } else try(subset <- filterTranscripts(txdb, 31, 30, NULL, longestPerGene = TRUE), silent = TRUE)
  if (length(subset) < length(mrna)*0.2) {
    if (length(subset) == 0) {
      warning("Organism has no leaders, extending 150nt on mrnas!")
    } else warning("Organism has too few leader defined mrnas, extending 150nt on mrnas!")
    mrna <- extendLeaders(mrna, 150)
    subset <- filterTranscripts(txdb, NULL, 30, NULL, longestPerGene = TRUE)
  }

  cds <- cds[names(cds) %in% subset]
  mrna <- mrna[names(mrna) %in% subset]

  cds_starts_p_vec_all <- pmapToTranscriptF(startSites(cds, TRUE, TRUE, TRUE), mrna[names(cds)])
  # Detect all and subset
  orfs_table_files <- list.files(org_dir_in, "prediction_table.rds", full.names = TRUE)
  if (length(orfs_table_files) == 0) {
    longestORF <- FALSE
    ORF_categories_to_keep = "NTE"
    startCodons <- paste("ATG", "CTG", "TTG", "GTG", sep = "|")
    orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df, TRUE),
                                    longestORF = longestORF, startCodon = startCodons)
    all_NTE_candidates <- ORFik:::categorize_and_filter_ORFs(orf_candidate_ranges, ORF_categories_to_keep, cds, mrna)
    saveRDS(all_NTE_candidates, file = file.path(org_dir_out, "NTE_all_candidates_ranges.rds"))

    detect_ribo_orfs(df, org_dir_in, ORF_categories_to_keep = ORF_categories_to_keep,
                     longestORF = longestORF, cds = cds, mrna = mrna,
                     orf_candidate_ranges = orf_candidate_ranges, orfs_gr = all_NTE_candidates,
                     startCodon = startCodons, libraries = list(RFP))
    table <- riboORFs(df[1,], type = "table", org_dir_in)
  } else all_NTE_candidates <- readRDS(file.path(org_dir_out, "NTE_all_candidates_ranges.rds"))

  orfs_table_files <- list.files(org_dir_in, "prediction_table.rds", full.names = TRUE)
  stopifnot(length(orfs_table_files) == 1)
  orfs_table <- readRDS(orfs_table_files)
  if (colnames(orfs_table)[1] == "gene") colnames(orfs_table)[1:2] <- c("ensembl_gene_id", "ensembl_tx_name")
  ids_pred <- orfs_table[predicted == TRUE,][!duplicated(ensembl_tx_name),]$id
  ids_non_pred <- orfs_table[predicted == FALSE][!(ensembl_tx_name %in% orfs_table[id %in% ids_pred]$ensembl_gene_id)][!duplicated(ensembl_tx_name),]$id
  stopifnot(!any(ids_non_pred %in% ids_pred))
  stopifnot(all(orfs_table$type == "NTE"))
  gr_files <- list.files(org_dir_in, "_candidates.rds", full.names = TRUE)
  stopifnot(length(gr_files) == 1)
  orfs_gr_raw <- readRDS(gr_files)
  stopifnot(nrow(orfs_table) == length(orfs_gr_raw))
  orfs_gr_raw <- orfs_gr_raw[c(ids_pred, ids_non_pred)]

  length(orfs_gr_raw)
  orfs_gr_raw <- orfs_gr_raw[names(orfs_gr_raw) %in% names(mrna)]
  length(orfs_gr_raw)

  orfs_starts_p_vec <- pmapToTranscriptF(startSites(orfs_gr_raw, TRUE, TRUE, TRUE), mrna[names(orfs_gr_raw)])

  cds_starts_p_vec <- cds_starts_p_vec_all[names(orfs_starts_p_vec)]

  NTE_section <- GRanges(seqnames(orfs_starts_p_vec), IRanges(start(orfs_starts_p_vec), start(cds_starts_p_vec) - 1), strand(orfs_starts_p_vec))
  summary(width(NTE_section))
  orfs_gr <- pmapFromTranscriptF(NTE_section, mrna[names(orfs_gr_raw)], removeEmpty = TRUE)

  gt_min_nt_CTE <- widthPerGroup(orfs_gr) > 30
  print(paste(round(sum(100 * gt_min_nt_CTE) / length(gt_min_nt_CTE), 2), "% (gt30)"))
  orfs_gr <- orfs_gr[gt_min_nt_CTE]

  message("Coverage statistics and filtering...")
  orfScores <- coverage_statistics(orfs_gr, RFP)

  overlap_filter <- orfs_gr %over% all_cds
  filter <- orfScores$ORFScores > 5 & orfScores$frame_zero_RP > 300 &
    orfScores$ORF_F0_codons_covered > 3 & orfScores$`ORF_F0_codons_covered(%)` > 10 &
    orfScores$frame_bias_relative > 0.4 & !overlap_filter
  final <- append_gene_ids(orfs_gr, df, orfScores)
  final[, longest_active_isoform := FALSE]
  final[filter, longest_active_isoform := seq(.N) == which.max(ORF_length_nt), by = tx_ids]
  filter <- filter & final$longest_active_isoform

  final[, predicted := filter]
  final[, overlaps_other_cds := overlap_filter]
  final <- orf_to_cds_coverage_statistcs(cds[final$tx_ids], RFP, final)
  message("-- Final candidates: ", nrow(final))
  final_predicted <- final[filter]
  message("-- Final Predicted: ", nrow(final_predicted))

  mcols(orfs_gr) <- DataFrame(pred = filter)
  predicted_gr <- orfs_gr[filter]
  cds_nte <- orfs_gr_raw[names(orfs_gr)]
  predicted_cds_nte <- cds_nte[names(predicted_gr)]

  fwrite(final, file = file.path(org_dir_out, "NTE_candidates.csv"))
  fst::write_fst(final, file.path(org_dir_out, "NTE_candidates.fst"))
  saveRDS(orfs_gr, file = file.path(org_dir_out, "NTE_candidates_ranges.rds"))
  saveRDS(cds_nte, file = file.path(org_dir_out, "CDS_NTE_candidates_ranges.rds"))
  fwrite(final_predicted, file = file.path(org_dir_out, "NTE_predicted.csv"))
  fst::write_fst(final_predicted, file.path(org_dir_out, "NTE_predicted.fst"))
  saveRDS(predicted_cds_nte, file = file.path(org_dir_out, "CDS_predicted_ranges.rds"))
  saveRDS(predicted_gr, file = file.path(org_dir_out, "NTE_predicted_ranges.rds"))

  if (org_short == "homo_sapiens") { # The sets for phylo analysis ->
    not_predicted <- orfs_gr[!filter & !overlap_filter & !(txNames(orfs_gr) %in% txNames(predicted_gr))]
    not_predicted_gr <-
      sample_to_matching_lengths_quantiles_and_size(predicted_gr, not_predicted)
    background_from_predicted_subset <- leaders[unique(names(predicted_gr))]
    background_from_predicted_subset_equal_size <-
      window_next_to_window(extendLeaders(background_from_predicted_subset, 200), predicted_gr,
                            start_direction = "3'")
    stopifnot(length(background_from_predicted_subset_equal_size) == length(predicted_gr))
    region <- "NTEs"
    background <- "leader"
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

