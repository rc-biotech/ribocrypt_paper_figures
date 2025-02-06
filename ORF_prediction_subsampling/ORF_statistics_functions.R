library(ORFik)
library(data.table)
library(massiveNGSpipe)
library(ggplot2)


# Function definitions
subset_to_group_index <- function(introns, index, up_to_index = FALSE, cds = NULL) {
  introns_second_intron <- introns
  strand_bool <- strandBool(introns)
  introns_second_intron[strand_bool] <- heads(introns[strand_bool], index)
  introns_second_intron[!strand_bool] <- tails(introns[!strand_bool], index)
  introns_second_intron <- introns_second_intron[lengths(introns_second_intron) == index]
  stopifnot(all(lengths(introns_second_intron) == index))

  if (!up_to_index & index != 1) {
    strand_bool <- strandBool(introns_second_intron)
    introns_second_intron[strand_bool] <- tails(introns_second_intron[strand_bool], 1)
    introns_second_intron[!strand_bool] <- heads(introns_second_intron[!strand_bool], 1)
    stopifnot(all(lengths(introns_second_intron) == 1))
    length(introns_second_intron) / length(introns_first_intron)
  }

  if (!is.null(cds)) {
    if (is.null(cds_starts_all) | is.null(cds_stops_all)) {
      cds_starts_all <- startSites(cds, keep.names = TRUE, is.sorted = TRUE)
      cds_stops_all <- stopSites(cds, keep.names = TRUE, is.sorted = TRUE)
    }
    introns_second_intron <- intron_in_cds_only(introns_second_intron, cds,
                                                cds_starts_all, cds_stops_all)
  }
  return(introns_second_intron)
}

sample_to_matching_lengths_quantiles_and_size <- function(candidates_gr_all,
                                                          noncandidates_gr_all,
                                                          max_diff = 4,
                                                          size = length(candidates_gr_all),
                                                          sampling_attempts = 100) {
  if (size == 0) return(GRangesList())
  A <- widthPerGroup(candidates_gr_all, FALSE)
  B <- widthPerGroup(noncandidates_gr_all, FALSE)
  Q_A <- quantile(A, probs = c(0.25, 0.50, 0.75))

  for (sampling in seq(sampling_attempts)) {
    sampling_indices <- match_quantiles_indices(A, B, size)
    B_resampled <- B[sampling_indices]

    Q_B <- quantile(B_resampled, probs = c(0.25, 0.50, 0.75))
    invalid_sampling <- mean(A) > mean(B_resampled) + max_diff | mean(A) < mean(B_resampled) - max_diff
    if (!invalid_sampling) break
  }
  stopifnot(!invalid_sampling)
  message("Completed on sampling attempt: ", sampling)
  message("New quantiles:")
  print(Q_A)
  print(Q_B)
  return(noncandidates_gr_all[sampling_indices])
}

#' Match quantiles from 2 sets
#'
#' Given vectors A and B, subset to length of A and match quantiles in B
#' to quantiles in A
#' @return The indices in B to subset on to get matching quantiles and size
match_quantiles_indices <- function(A, B, size = length(A)) {
  # Compute quantiles of A
  A_quants <- quantile(A, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

  # Assign each value in A to a quantile bin
  A_bins <- cut(A, breaks = c(-Inf, A_quants, Inf), labels = FALSE, include.lowest = TRUE)

  # Get proportion of A in each bin
  A_bin_counts <- table(A_bins) / length(A)

  # Compute corresponding bins for B
  B_bins <- cut(B, breaks = c(-Inf, A_quants, Inf), labels = FALSE, include.lowest = TRUE)

  # Sample indices from B based on A's bin proportions
  selected_indices <- c()

  for (bin in seq_along(A_bin_counts)) {
    bin_indices <- which(B_bins == bin)

    # Number to sample from this bin
    sample_size <- round(A_bin_counts[bin] * size)

    # Ensure at least one sample if available
    if (length(bin_indices) > 0) {
      selected_indices <- c(selected_indices, sample(bin_indices, min(sample_size, length(bin_indices)), replace = FALSE))
    }
  }

  # Ensure the final sample size matches the requested size
  if (length(selected_indices) > size) {
    selected_indices <- sample(selected_indices, size, replace = FALSE)
  } else if (length(selected_indices) < size) {
    extra_indices <- sample(setdiff(seq_along(B), selected_indices), size - length(selected_indices), replace = FALSE)
    selected_indices <- c(selected_indices, extra_indices)
  }

  return(selected_indices)
}

intron_in_cds_only <- function(introns_first_intron, cds, cds_starts_all,
                               cds_stops_all) {
  stopifnot(all(lengths(introns_first_intron) == 1))
  cds <- cds[names(introns_first_intron)]
  stopifnot(length(cds) == length(introns_first_intron))
  intron_starts <- startSites(introns_first_intron, TRUE, TRUE, TRUE)
  intron_stops <- stopSites(introns_first_intron, TRUE, TRUE, TRUE)
  cds_starts <- cds_starts_all[names(intron_starts)]
  cds_stops <- cds_stops_all[names(intron_starts)]

  diff <- rep(as.integer(NA), length(intron_starts))
  strand_bool <- strandBool(intron_starts)
  diff[strand_bool] <-  start(intron_starts[strand_bool]) - cds_starts[strand_bool]
  diff[!strand_bool] <- cds_starts[!strand_bool] - end(intron_starts[!strand_bool])

  diff_stop <- rep(as.integer(NA), length(intron_stops))
  diff_stop[strand_bool] <-  end(intron_stops[strand_bool]) - cds_stops[strand_bool]
  diff_stop[!strand_bool] <- cds_stops[!strand_bool] - start(intron_stops[!strand_bool])

  stopifnot(!anyNA(diff))
  stopifnot(!anyNA(diff_stop))
  # stopifnot(all(diff > 0))
  summary(diff)
  summary(diff[strand_bool])
  # Update
  valid <- diff > 0 & diff_stop < 0
  strand_bool <- strand_bool[valid]
  intron_starts <- intron_starts[valid]
  introns_first_intron <- introns_first_intron[valid]
  introns <- introns[valid]
  cds <- cds[valid]
  mrna <- mrna[valid]
  diff <- diff[valid]


  phase <- abs(diff %% -3)
  mcols(introns_first_intron)$phase <- phase
  mcols(introns_first_intron)$cds_exon_length <- widthPerGroup(subset_to_group_index(cds, 1))
  return(introns_first_intron)
}

append_gene_ids <- function(candidates_gr, df, candidates, tx_ids = txNames(candidates_gr)) {
  tx_ids <- data.table(CTE_id = names(candidates_gr), tx_ids)
  symbols <- symbols(df)
  if (nrow(symbols) > 0) {
    if (!is.null(symbols$uniprot_id) )symbols[, uniprot_id := NULL]
    merge <- data.table::merge.data.table(tx_ids, symbols, by.x = "tx_ids", by.y = "ensembl_tx_name",
                                          all.x = TRUE, all.y = FALSE, sort = FALSE)
  } else merge <- tx_ids

  final <- cbind(merge, candidates)
  final[]
  return(final)
}

CTE_orfs <- function(windows, faFile, min_width = 30) {
  gt_min_nt_CTE <- widthPerGroup(windows) > min_width
  windows <- windows[gt_min_nt_CTE]
  summary(widthPerGroup(windows))
  # find stop codon
  message("- Find stop codons")
  seqs <- txSeqsFromFa(windows, faFile, TRUE, TRUE)
  codon_1_is_stop <- grepl(stopDefinition(1), subseq(seqs, 1, 3))
  windows <- windows[!codon_1_is_stop]
  seqs <- seqs[!codon_1_is_stop]
  subseq(seqs, 1, 3) <- "ATG"
  orfs <- findORFs(seqs, startCodon = "ATG")
  orfs <- orfs[start(orfs) == 1] # Subset back to only start
  length(orfs)
  orfs_gr <- ORFik:::mapToGRanges(windows, orfs, groupByTx = FALSE, grl_is_sorted = TRUE)
  mcols(orfs_gr) <- DataFrame(codon_1_is_stop = sum(codon_1_is_stop),
                              width_gt_min = sum(gt_min_nt_CTE))[1,]
  return(orfs_gr)
}

coverage_statistics <- function(orfs_gr, RFP) {
  cov <- coveragePerTiling(orfs_gr, RFP, as.data.table = TRUE, withFrames = TRUE)
  cov[, codon := ceiling(position / 3) ]
  cov[, codon_sum := sum(count), by = .(genes, codon)]
  cov[, count_relative_to_codon := count / codon_sum, by = .(genes, codon)]
  cov[!is.finite(count_relative_to_codon), count_relative_to_codon := 0]

  cov_rel <- cov[, .(sum_rel_frame = sum(count_relative_to_codon),
                     non_zero_codons_frame = sum(count_relative_to_codon > 0)), by = .(genes, frame)]
  cov_rel[, non_zero_codons_total := as.integer(rep(cov[frame == 0, sum(codon_sum > 0), by = genes]$V1, each = 3))]
  cov_rel <- cov_rel[frame == 0,]
  cov_rel[, frame_bias_relative := round(sum_rel_frame / non_zero_codons_total, 2)]
  cov_rel[!is.finite(frame_bias_relative), frame_bias_relative := 0]

  orfScores <- orfScore(orfs_gr, RFP, is.sorted = TRUE, coverage = cov)
  orfScores[, frame_bias := round(frame_zero_RP / (frame_zero_RP + frame_one_RP + frame_two_RP), 2)]
  orfScores[!is.finite(frame_bias), frame_bias := 0]
  orfScores[, frame_bias_relative := cov_rel$frame_bias_relative]
  orfScores[, ORF_length_nt := widthPerGroup(orfs_gr)]
  orfScores[, ORF_length_codons := ORF_length_nt / 3]

  cov_hits <- cov[, sum(count > 0), by = .(genes, frame)]
  orfScores[, `ORF_F0_codons_covered` := cov_hits[frame == 0]$V1]
  orfScores[, `ORF_F1_codons_covered` := cov_hits[frame == 1]$V1]
  orfScores[, `ORF_F2_codons_covered` := cov_hits[frame == 2]$V1]
  orfScores[, `ORF_F0_codons_covered(%)` := 100*(ORF_F0_codons_covered / ORF_length_codons)]
  orfScores[, `ORF_F1_codons_covered(%)` := 100*(ORF_F1_codons_covered / ORF_length_codons)]
  orfScores[, `ORF_F2_codons_covered(%)` := 100*(ORF_F2_codons_covered / ORF_length_codons)]

  orfScores[, ORF_sum := rowSums(orfScores[, .(frame_zero_RP, frame_one_RP, frame_two_RP)])]
  orfScores[, ORF_count_length_ratio := ORF_sum / ORF_length_nt]
  orfScores[, ORF_F0_count_length_ratio :=  frame_zero_RP/ ORF_length_nt]
  orfScores[]
  return(orfScores)
}

orf_to_cds_coverage_statistcs <- function(cds, RFP, dt) {
  dt[, cds_sum := countOverlaps(cds, RFP)]
  dt[, cds_length_nt := widthPerGroup(cds)]
  dt[, cds_count_length_ratio := cds_sum / cds_length_nt]
  dt[, `ORF_cds_count_length_ratio(%)` := (ORF_count_length_ratio / cds_count_length_ratio)*100]
  dt[]
  return(dt)
}

window_next_to_window <- function(whole_window,
                                  first_window,
                                  start_direction = "5'",
                                  width_multiplier = 1L,
                                  removeEmpty = TRUE) {
  stopifnot(length(whole_window) == length(first_window))
  if (length(whole_window) == 0) return(GRangesList())
  whole_window_grl <- whole_window
  whole_window <- pmapToTranscriptF(whole_window, whole_window, FALSE, TRUE, TRUE)
  window_width <- widthPerGroup(first_window, FALSE)
  if (start_direction == "5'") {
    start(whole_window) <- start(whole_window) + window_width
    end(whole_window) <- start(whole_window) + window_width*width_multiplier + 1L
  } else {
    end(whole_window) <- end(whole_window) - window_width
    start(whole_window) <- end(whole_window) - window_width*width_multiplier - 1L
  }

  whole_window <- trim(whole_window)
  to_map_back <- unlistGrl(whole_window)
  names(to_map_back) <- NULL
  second_window <- pmapFromTranscriptF(to_map_back, whole_window_grl, removeEmpty = removeEmpty)
  return(second_window)
}
