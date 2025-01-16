# Setup
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")

result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
df <- read.experiment("all_merged-Homo_sapiens_04_oct_2024_all")


RFP <- fimport(filepath(df, "cov"))
mane <- data.table::fread("~/livemount/Bio_data/references/homo_sapiens/canonical_isoforms.txt")[[1]]
# mane <- filterTranscripts(df, NULL, 150, NULL)
leaders <- loadRegion(df, "leaders")
leaders <- leaders[lengths(leaders) == 1] # First intron must be in CDS
introns <- loadRegion(df, "introns")
introns <- introns[lengths(introns) > 0]

all_cds <- cds <- loadRegion(df, "cds")
mrna <- loadRegion(df, "mrna")

# Subset to intronic
coding_tx <- names(cds)
name_subset <- coding_tx[coding_tx %in% mane & coding_tx %in% names(introns) & coding_tx %in% names(leaders)]
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
  phase <- mcols(introns_first_intron)$phase
  cds_exon_length <- mcols(introns_first_intron)$cds_exon_length
  intron_starts <- startSites(introns_first_intron, TRUE, TRUE, TRUE)

  windows_all <- windowPerGroup(intron_starts, tx = introns_first_intron, upstream = -(phase), downstream = 300)
  summary(widthPerGroup(windows_all))

  orfs_gr <- CTE_orfs(windows_all, df)

  orfs_gr@unlistData$intron_index <- intron_index
  length(orfs_gr)
  introns_first_intron <- introns_first_intron[txNames(orfs_gr)]
  # Detect coverage
  message("-- Coverage on stop codon extensions")
  orfScores <- coverage_statistics(orfs_gr, RFP)

  overlap_filter <- orfs_gr %over% all_cds
  filter <- orfScores$ORFScores > 5 & orfScores$frame_zero_RP > 300 &
    orfScores$ORF_F0_codons_covered > 3 & orfScores$`ORF_F0_codons_covered(%)` > 10 & !overlap_filter

  # TODO check if!
  if (sum(orfScores[filter]$ORF_sum) == 0) final <- data.table()
  final <- append_gene_ids(orfs_gr, df, orfScores)
  final <- orf_to_cds_coverage_statistcs(cds[final$tx_ids], RFP, final)

  N_introns_filtered <- sum(filter)
  N_introns_gt30nt <- mcols(orfs_gr)[1,]$width_gt_min
  N_introns_has_stop <- length(orfs_gr)
  N_introns_codon_1_is_stop <- mcols(orfs_gr)[1,]$codon_1_is_stop
  N_introns <- length(introns_first_intron)
  intron_specific_stats <- cbind(intron_index, N_introns, N_introns_gt30nt, N_introns_codon_1_is_stop, N_introns_has_stop,
                                 N_introns_filtered, intron_length = widthPerGroup(introns_first_intron[txNames(orfs_gr)]))
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
candidates_gr_all <- candidates_gr_all[!(countOverlaps(candidates_gr_all, candidates_gr_all, type = "equal") > 1)]
noncandidates_gr_all <-  all_orfs_gr_all[countOverlaps(all_orfs_gr_all, candidates_gr_all, type = "equal") == 0]
stopifnot(length(all_orfs_gr_all) == length(noncandidates_gr_all) + length(candidates_gr_all))

fst::write_fst(all_indices_dt,  file.path(result_dir, "intron_candidates_intron_1_to_7.fst"))
qs::qsave(candidates_gr_all, file.path(result_dir, "intron_candidates_intron_1_to_7_ranges.qs"))
qs::qsave(all_orfs_gr_all, file.path(result_dir, "intron_all_orfs_intron_1_to_7_ranges.qs"))
qs::qsave(all_introns_gr_all, file.path(result_dir, "intron_all_whole_intron_1_to_7_ranges.qs"))
qs::qsave(noncandidates_gr_all, file.path(result_dir, "intron_all_not_predicted_1_to_7_ranges.qs"))

ORFik::export.bed12(candidates_gr_all, file.path(result_dir, "intron_candidates.bed12"))
ORFik::export.bed12(noncandidates_gr_all, file.path(result_dir, "intron_not_predicted.bed12"))
ORFik::export.bed12(candidates_gr_all, file.path(result_dir, "intron_candidates_whole_intron_1_to_7.bed12"))
# candidates_gr_all_unlist <- unlistGrl(candidates_gr_all)
# noncandidates_gr_all_unlist <- unlistGrl(noncandidates_gr_all)
# rtracklayer::export.bed(candidates_gr_all_unlist, file.path(result_dir, "intron_candidates_intron_1_to_7.bed"))


