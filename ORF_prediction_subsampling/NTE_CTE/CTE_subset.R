source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")

ref_dir <- ORFik::config()["ref"]
all_exp <- list.experiments(validate = FALSE, pattern = "all_merged-", libtypeExclusive = "RFP")
all_exp <- all_exp[libtypes == "RFP"]
organisms <- all_exp$name
organisms <- organisms[grep("_modalities$|HEK293|_04|Homo_sapiens", organisms, invert = TRUE)]
org <- "all_merged-Homo_sapiens_04_oct_2024_all"
organisms <- c(org, organisms)[1]
for (org in organisms) {
  df <- read.experiment(org, validate = FALSE)
  org_short <- gsub(" ", "_", tolower(organism(df)))
  message("- ", org_short)
  org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik", "CTE_candidates")
  if (file.exists(file.path(org_dir, "CDS_CTE_candidates_ranges.rds"))) next
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

  cds_stops <- stopSites(cds, TRUE, TRUE, TRUE)
  # cds_stops <- split(cds_stops, names(cds_stops))
  windows_all <- windowPerGroup(cds_stops, tx = mrna, upstream = -1, downstream = 300)
  summary(widthPerGroup(windows_all))

  orfs_gr <- CTE_orfs(windows_all, df)
  length(orfs_gr)

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
}

org <- "all_merged-Homo_sapiens_04_oct_2024_all"
df <- read.experiment(org, validate = FALSE)

org_short <- gsub(" ", "_", tolower(organism(df)))
message("- ", org_short)
org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik")
org_dir_out <- file.path(org_dir, "CTE_candidates")
# Annotation and track
RFP <- fimport(filepath(df, "cov"))
txdb <- loadTxdb(df)
all_cds <- cds <- loadRegion(txdb, "cds")
mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
subset <- mane
mrna <- loadRegion(txdb, "mrna")
trailers <- loadRegion(df, "leaders")
cds <- all_cds[names(all_cds) %in% subset]
mrna <- mrna[names(mrna) %in% subset]
# all_CTE_candidates <- readRDS(file.path(org_dir_out, "CTE_all_candidates_ranges.rds"))
stopifnot(length(cds) == length(mrna)); length(cds)
# Results
final <- setDT(fst::read_fst(file.path(org_dir_out, "CTE_candidates.fst")))
orfs_gr <- readRDS(file = file.path(org_dir_out, "CTE_candidates_ranges.rds"))

all_indices_dt_merged <- copy(final)
all_indices_dt_merged[, predicted := as.factor(predicted)]
predicted <- all_indices_dt_merged$predicted
all_indices_dt_merged$predicted <- NULL
melt_all <- melt(all_indices_dt_merged[, -c("tx_ids", "CTE_id", "ensembl_gene_id", "external_gene_name")])
melt_all[, predicted := rep(predicted, length.out = nrow(melt_all))]
melt_all[ variable == "ORF_cds_count_length_ratio(%)", variable := "Relative expression(%)"]
melt_all <- rbindlist(list(melt_all, data.table(variable = "N_retained_CTE(%)",
                                                value = round(100*sum(final$predicted) / nrow(final), 2),
                                                              predicted = FALSE)))
melt_all[ variable == "cds_sum", variable := "CDS expression"]
melt_all[, variable := factor(variable, levels = unique(variable))]
melt_all[, predicted := factor(predicted, levels = c(T, F))]

melt_all_subset <- melt_all[variable %in% c("N_retained_CTE(%)", "ORF_F0_codons_covered(%)",
                                            "frame_bias_relative", "ORF_length_codons", "CDS expression",
                                            "Relative expression(%)")]
melt_all_subset[variable == "ORF_length_codons", variable := "Distance to first inframe stop"]
melt_all_subset[variable == "ORF_F0_codons_covered(%)", variable := "Codons covered (%)"]
melt_all_subset[variable == "frame_bias_relative", variable := "Frame bias"]
melt_all_subset[, variable := factor(variable, levels = unique(variable)[c(6,2,1,3,4,5)])]
stats_plot_CTE <- ggplot(melt_all_subset, aes(predicted, value, fill = predicted)) + geom_boxplot(outliers = FALSE) +
  facet_wrap(~ variable, scales = "free", ncol = 2) + theme_classic() + xlab("Prediction status") +
  ggtitle("Human CTE analysis", subtitle = paste0("On all mane isoforms with 3' UTRs (total tx:", nrow(final), ")")) +
  theme(legend.position = "top"); stats_plot_CTE

# MetaWindow plot
extend <- 18
candidates_gr_45 <- orfs_gr[widthPerGroup(orfs_gr) > 20]
gr_to_plot <- windowPerGroup(candidates_gr_45, candidates_gr_45, 0, 20)
gr_to_plot <- extendLeaders(gr_to_plot, extend)
mcols(gr_to_plot)$predicted <- final[chmatch(CTE_id, names(gr_to_plot)),][!is.na(tx_ids),]$predicted
cov_final_all <- rbindlist(lapply(c(TRUE, FALSE), function(prediction_status) {
  dt <- coveragePerTiling(gr_to_plot[mcols(gr_to_plot)$predicted == prediction_status], RFP, as.data.table = TRUE, withFrames = TRUE)
  dt[, fraction := prediction_status]
}))
cov_final <- cov_final_all[position <= 50 + extend,]
cov_final <- coverageScorings(cov_final, "mean")
cov_final[, frame := (position-1) %% 3]
cov_final[, position := position - extend]
cov_final[, prediction_status := fraction]

plot_with_cds_line <- ggplot(cov_final, aes(x = position, y = score, col = prediction_status)) +
  geom_line() + theme_classic() + xlab("Position relative to CDS stop site (TES) (by nt)") +
  ylab("Z-score  counts") + scale_color_brewer(palette="Dark2"); plot_with_cds_line
metacoverage_plot_CTE <- plot_with_cds_line


# cov_final <- copy(cov_final_all)
# cov_final[, position := position - extend]
# cov_final <- cov_final[position <= 50,]
# cov_final <- coverageScorings(cov_final, "transcriptNormalized")
# cov_final[, frame := (position-1) %% 3]
#
# plot_with_cds <- pSitePlot(cov_final, facet = TRUE, title = "Coverage by prediction status (with CDS flank)")
# plot_with_cds <- plot_with_cds + theme_classic() + xlab("Position of Stop codon (from first in frame nt relative to CDS)") + ylab("ZScore normalized")
# plot_with_cds
#
cov_final <- copy(cov_final_all)
cov_final[, position := position - extend]
cov_final <- cov_final[position <= 50,][position >= 0,]
cov_final <- coverageScorings(cov_final, "transcriptNormalized")
cov_final[, frame := (position-1) %% 3]

plot_with_cds <- pSitePlot(cov_final, facet = TRUE, title = "Coverage by prediction status (with CDS flank)")
plot_with_cds <- plot_with_cds + theme_classic() + xlab("Position of Stop codon (from first in frame nt relative to CDS)") + ylab("ZScore normalized")
plot_with_cds
