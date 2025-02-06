source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")
ref_dir <- ORFik::config()["ref"]

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
candidates_gr_45 <- orfs_gr#[widthPerGroup(orfs_gr) > 20]
gr_to_plot <- windowPerGroup(startSites(candidates_gr_45, TRUE, TRUE, TRUE), mrna[names(candidates_gr_45)], extend, 44)
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
cov_final <- coverageScorings(cov_final, "sum")
cov_final[, frame := (position-1) %% 3]


plot_with_cds <- pSitePlot(cov_final, facet = TRUE, title = "Coverage by prediction status (with CDS flank)", frameSum = T)
plot_with_cds <- plot_with_cds + theme_classic() + xlab("Position of Stop codon (from first in frame nt relative to CDS)") + ylab("ZScore normalized")
plot_with_cds
