source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")

result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
intronic_genes <- 11575


# Load
all_orfs_gr_all <- qs::qread(file.path(result_dir, "intron_all_orfs_intron_1_to_7_ranges.qs"))
candidates_gr_all <- qs::qread(file.path(result_dir, "intron_candidates_intron_1_to_7_ranges.qs"))
noncandidates_gr_all <- qs::qread(file.path(result_dir, "intron_all_not_predicted_1_to_7_ranges.qs"))
all_introns_gr_all <- qs::qread(file.path(result_dir, "intron_all_whole_intron_1_to_7_ranges.qs"))
all_indices_dt <- setDT(fst::read_fst(file.path(result_dir, "intron_candidates_intron_1_to_7.fst")))


# Compute tables
all_indices_dt_merged <- copy(all_indices_dt)
all_indices_dt_merged[, predicted := as.factor(predicted)]
predicted <- all_indices_dt_merged$predicted
all_indices_dt_merged$predicted <- NULL
melt_all <- melt(all_indices_dt_merged[, -c("tx_ids", "CTE_id", "ensembl_gene_id", "external_gene_name")])
melt_all[, intron_index := as.factor(rep(melt_all[variable == "intron_index"]$value, length.out = nrow(melt_all)))]
melt_all <- melt_all[variable != "intron_index"]
melt_all[, predicted := rep(predicted, length.out = nrow(melt_all))]
melt_all[ variable == "ORF_cds_count_length_ratio(%)", variable := "Relative expression(%)"]
melt_all[ variable == "N_introns_filtered", variable := "N_retained_introns(%)"]
melt_all[ variable == "cds_sum", variable := "CDS expression"]
melt_all[variable == "N_retained_introns(%)", value := round(100*value / melt_all[variable == "N_introns"]$value, 2)]
melt_all <- melt_all[!(predicted == "FALSE" & variable %in% c("N_introns", "N_retained_introns(%)"))]
melt_all[, variable := factor(variable, levels = unique(variable))]
melt_all[, predicted := factor(predicted, levels = c(T, F))]

melt_all_subset <- melt_all[variable %in% c("N_retained_introns(%)", "ORF_F0_codons_covered(%)",
                                            "frame_bias_relative", "ORF_length_codons", "CDS expression",
                                            "Relative expression(%)")]
melt_all_subset[variable == "ORF_length_codons", variable := "Distance to first inframe stop"]
melt_all_subset[variable == "ORF_F0_codons_covered(%)", variable := "Codons covered (%)"]
melt_all_subset[variable == "frame_bias_relative", variable := "Frame bias"]
melt_all_subset[, variable := factor(variable, levels = unique(variable)[c(1,3,2,4,5,6)])]
melt_all_subset[, intron_index := as.numeric(as.character(intron_index))]
melt_all_subset[, orf_index := seq.int(.N), by = .(variable, intron_index, predicted)]
melt_all_subset[, intron_index := as.factor(ifelse(intron_index > 2, "3+", intron_index))]
melt_all_subset <- melt_all_subset[,.(value = mean(value)),by = .(variable, intron_index, predicted, orf_index)]
intronic_stats_plot <- ggplot(melt_all_subset, aes(intron_index, value, fill = predicted)) + geom_boxplot(outliers = FALSE) +
  facet_wrap(~ variable, scales = "free", ncol = 2) + theme_classic() + xlab("Intron index") +
  ggtitle("Human intron analysis", subtitle = paste0("On all CDS spliced mane isoforms (total tx:", intronic_genes, ")")) +
  theme(legend.position = "top"); intronic_stats_plot


# MetaWindow plot
df <- read.experiment("all_merged-Homo_sapiens_04_oct_2024_all")
RFP <- fimport(filepath(df, "cov"))

extend <- 18
intron_indices <- seq(max(as.numeric(as.character(melt_all$intron_index))))
candidates_gr_all_45 <- candidates_gr_all[widthPerGroup(candidates_gr_all) > 44]
names(candidates_gr_all_45) <- make.unique(names(candidates_gr_all_45))
gr_to_plot <- windowPerGroup(candidates_gr_all_45, candidates_gr_all_45, 0, 40)
gr_to_plot <- extendLeaders(gr_to_plot, extend)
gr_to_plot@unlistData$intron_index <- candidates_gr_all_45@unlistData$intron_index
cov_final_all <- rbindlist(lapply(intron_indices, function(intron_index) {
  dt <- coveragePerTiling(gr_to_plot[gr_to_plot@unlistData$intron_index == intron_index], RFP, as.data.table = TRUE, withFrames = TRUE)
  dt[, fraction := intron_index]
}))
cov_final <- cov_final_all[position <= 50 + extend,]
cov_final <- coverageScorings(cov_final, "fracPos")
cov_final[, frame := (position-1) %% 3]
cov_final[, position := position - extend]

# Barplot per intron index
plot_with_cds <- pSitePlot(cov_final, facet = TRUE, title = "Coverage per Intron index (with CDS flank)")
plot_with_cds <- plot_with_cds + theme_classic() + xlab("Position of Intron (from first in frame nt relative to CDS)") + ylab("Transcript Normalized counts")
plot_with_cds

# Lineplot with CDS extension
cov_final[,Intron := ifelse(fraction > 2, "3+", fraction)]
cov_final_filtered_intron <- cov_final[,.(score = mean(score)),by = .(Intron, position)]
plot_with_cds_line <- ggplot(cov_final_filtered_intron, aes(x = position, y = score, col = Intron)) +
  geom_line() + theme_classic() + xlab("Position relative to exon/intron junction (by codon)") +
  ylab("Z-score normalized counts") + scale_color_brewer(palette="Dark2"); plot_with_cds_line
metacoverage_plot <- plot_with_cds_line

# Lineplot without CDS extension
cov_final_filtered_intron <- copy(cov_final_filtered_intron)
cov_final_filtered_intron[, log_score := log10(score*1000)]
plot_with_cds_line <- ggplot(cov_final_filtered_intron, aes(x = position, y = log_score, col = Intron)) +
  geom_line() + theme_classic() + xlab("Position relative to exon/intron junction (by codon)") +
  ylab("Z-score normalized counts") + scale_color_brewer(palette="Dark2"); plot_with_cds_line
metacoverage_plot <- plot_with_cds_line




cov_final_filtered <- cov_final[position <= 45 & position > (extend + 3),]
cov_final_filtered[, position := position - (extend + 3)]
plot <- pSitePlot(cov_final_filtered, facet = TRUE, title = "Coverage per Intron index")
plot <- plot + theme_classic() +
  xlab("Position of Intron (from first in frame nt relative to CDS)") + ylab("Transcript Normalized counts")
plot
metacoverage_plot_old <- pSitePlot(cov_final_filtered, facet = FALSE, title = "Coverage of exon/intron junction")
metacoverage_plot_old <- metacoverage_plot_old + theme_classic() +
  xlab("Position of Intron (from first in frame nt relative to CDS)") + ylab("Transcript Normalized counts")
metacoverage_plot_old

cov_final_filtered_log <- copy(cov_final_filtered)
# cov_final_filtered_log[, score := log2(score + 1)]
coverageHeatMap(cov_final_filtered_log, title = "Coverage per Intron index (with CDS flank)", addFracPlot = TRUE,
                ylab = "Intron Index", xlab = "Position relative to exon/intron junction", scoring = "Transcript Normalized (log2)")

plots <- cowplot::plot_grid(plot_with_cds, plot_with_cds_line, ncol = 1)
# massiveNGSpipe:::plot_all_versions(plots, file.path(result_dir, "intron_indices_coverage"), send_to_discord = TRUE)

# OLD code
# ggplot(melt_all, aes(intron_index, value, fill = predicted)) + geom_boxplot(outliers = FALSE) +
#   facet_wrap(~ variable, scales = "free", ncol = 5) + theme_classic() + xlab("Intron index") +
#   ggtitle("Human intron analysis", subtitle = paste0("On all CDS spliced mane isoforms (total tx:", intronic_genes, ")"))
#
# melt_all_subset <- melt_all[variable %in% c("N_introns", "N_retained_introns(%)", "ORFScores", "ORF_length_codons", "CDS expression", "Relative expression(%)")]
# intronic_stats_plot_old <- ggplot(melt_all_subset, aes(intron_index, value, fill = predicted)) + geom_boxplot(outliers = FALSE) +
#   facet_wrap(~ variable, scales = "free", ncol = 2) + theme_classic() + xlab("Intron index") +
#   ggtitle("Human intron analysis", subtitle = paste0("On all CDS spliced mane isoforms (total tx:", intronic_genes, ")")) +
#   theme(legend.position = "top"); intronic_stats_plot_old

# Inspect
# cov_final_all_sub <- cov_final_all[position > (extend + 3),]
# for (g in seq(unique(cov_final_all_sub$genes))) {
#   print(symbols[ensembl_tx_name == txNames(candidates_gr_all)[g]])
#   plot(pSitePlot(cov_final_all_sub[genes == g & fraction == 1,]))
#   readline(prompt="Press [enter] to continue")
# }

# mcols(candidates_gr_all)$intron_index <- all_indices_dt$intron_index
#
# seqs_all <- txSeqsFromFa(split(unlistGrl(introns), seq_along(introns@unlistData)), df, TRUE, TRUE)
# seqs_all_sub <- subseq(seqs_all, 1,1)
# names(seqs_all_sub) <- rep(txNames(introns), lengths(introns))
#
# seqs_sub <- seqs_all_sub[paste0(txNames(candidates_gr_all), )]
#
#
# tab_seqs <- table(as.data.table(seqs_sub))
# tab_seqs_all <- table(as.data.table(seqs_all_sub))
# plot(tab_seqs[tab_seqs > 5])
