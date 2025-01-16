source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")
ref_dir <- ORFik::config()["ref"]


org <- "all_merged-Homo_sapiens_04_oct_2024_all"
df <- read.experiment(org, validate = FALSE)
org_short <- gsub(" ", "_", tolower(organism(df)))
message("- ", org_short)
org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik")
org_dir_in <- file.path(org_dir, "longest_predicted")
org_dir_out <- file.path(org_dir, "NTE_candidates")
# Annotation and track
RFP <- fimport(filepath(df, "cov"))
txdb <- loadTxdb(df)
all_cds <- cds <- loadRegion(txdb, "cds")
mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
subset <- mane
mrna <- loadRegion(txdb, "mrna")
leaders <- loadRegion(df, "leaders")
cds <- cds[names(cds) %in% subset]
mrna <- mrna[names(mrna) %in% subset]
all_NTE_candidates <- readRDS(file.path(org_dir_out, "NTE_all_candidates_ranges.rds"))
stopifnot(length(cds) == length(mrna)); length(cds)
# Results
final <- setDT(fst::read_fst(file.path(org_dir_out, "NTE_candidates.fst")))
orfs_gr <- readRDS(file = file.path(org_dir_out, "NTE_candidates_ranges.rds"))

all_indices_dt_merged <- copy(final)
all_indices_dt_merged[, predicted := as.factor(predicted)]
predicted <- all_indices_dt_merged$predicted
all_indices_dt_merged$predicted <- NULL
melt_all <- melt(all_indices_dt_merged[, -c("tx_ids", "CTE_id", "ensembl_gene_id", "external_gene_name")])
melt_all[, predicted := rep(predicted, length.out = nrow(melt_all))]
melt_all[ variable == "ORF_cds_count_length_ratio(%)", variable := "Relative expression(%)"]
melt_all <- rbindlist(list(melt_all, data.table(variable = "N_retained_CTE(%)",
                                                value = round(100*sum(final$predicted) / length(unique(final[predicted == FALSE]$tx_ids)), 2),
                                                predicted = FALSE)))
melt_all[ variable == "cds_sum", variable := "CDS expression"]
melt_all[, variable := factor(variable, levels = unique(variable))]
melt_all[, predicted := factor(predicted, levels = c(T, F))]

melt_all_subset <- melt_all[variable %in% c("N_retained_CTE(%)", "ORF_F0_codons_covered(%)",
                                            "frame_bias_relative", "ORF_length_codons", "CDS expression",
                                            "Relative expression(%)")]
melt_all_subset[variable == "ORF_length_codons", variable := "Distance to CDS_TIS"]
melt_all_subset[variable == "ORF_F0_codons_covered(%)", variable := "Codons covered (%)"]
melt_all_subset[variable == "frame_bias_relative", variable := "Frame bias"]
melt_all_subset[, variable := gsub("CTE", "NTE", variable)]
melt_all_subset[, variable := factor(variable, levels = unique(variable)[c(6,2,1,3,4,5)])]
stats_plot_NTE <- ggplot(melt_all_subset, aes(predicted, value, fill = predicted)) + geom_boxplot(outliers = FALSE) +
  facet_wrap(~ variable, scales = "free", ncol = 2) + theme_classic() + xlab("Prediction status") +
  ggtitle("Human NTE analysis", subtitle = paste0("On all mane/canonical isoforms with coverage on 5' UTRs\n",
                                                  "Total tx: ", length(mrna),
                                                  "\nTotal tx with NTE:", length(unique(names(all_NTE_candidates))),
                                                  "\nTotal tx with candidate NTE:", length(unique(final$tx_ids)), " (> 30nt & 10 reads on NTE)")) +
  theme(legend.position = "top"); stats_plot_NTE

# MetaWindow plot
extend <- 18
extend_up <- 30
candidates_gr_45 <- orfs_gr[widthPerGroup(orfs_gr) > extend_up]
gr_to_plot <- windowPerGroup(stopSites(candidates_gr_45, TRUE, TRUE, TRUE), mrna, extend_up, extend)
mcols(gr_to_plot) <- mcols(candidates_gr_45)
gr_to_plot <- gr_to_plot[widthPerGroup(gr_to_plot) == 49]
pred_status <- mcols(gr_to_plot)$pred
cov_final_all <- rbindlist(lapply(c(TRUE, FALSE), function(prediction_status) {
  dt <- coveragePerTiling(gr_to_plot[pred_status == prediction_status], RFP, as.data.table = TRUE, withFrames = TRUE)
  dt[, fraction := prediction_status]
}))

cov_final <- cov_final_all[position <= 50 + extend,]
cov_final <- coverageScorings(cov_final, "mean")
cov_final[, frame := (position-1) %% 3]
cov_final[, position := position - extend_up]
cov_final[, prediction_status := fraction]

plot_with_cds_line <- ggplot(cov_final, aes(x = position, y = score, col = prediction_status)) +
  geom_line() + theme_classic() + xlab("Position relative to CDS TIS (by nt)") +
  ylab("Mean counts") + scale_color_brewer(palette="Dark2"); plot_with_cds_line
metacoverage_plot_NTE <- plot_with_cds_line

# Inspected:
# Case where next start would fail
#ribocrypt.org/?&dff=all_merged-Homo_sapiens_04_oct_2024_all&gene=SPSB1-ENSG00000171621&tx=ENST00000328089&library=RFP&frames_type=columns&kmer=1&log_scale=FALSE&extendLeaders=0&extendTrailers=0&viewMode=FALSE&other_tx=FALSE&add_uorfs=FALSE&genomic_region=&customSequence=NTG&phyloP=FALSE&summary_track=FALSE&go=TRUE
# Case where next start would work
#ribocrypt.org/?&dff=all_merged-Homo_sapiens_04_oct_2024_all&gene=PCSK9-ENSG00000169174&tx=ENST00000302118&library=RFP&frames_type=columns&kmer=1&log_scale=FALSE&extendLeaders=50&extendTrailers=0&viewMode=FALSE&other_tx=FALSE&add_uorfs=FALSE&genomic_region=&customSequence=NTG&phyloP=FALSE&summary_track=FALSE&go=TRUE
# Whcih are found by this model, but not in full model ?
# final_predicted[!(tx_ids %in% orfs_table[predicted == TRUE & type == "NTE" & ensembl_tx_name %in% names(mrna)]$ensembl_tx_name)]
# Whcih are found by full model, but not in this model ?
# orfs_table[predicted == TRUE & type == "NTE" & ensembl_tx_name %in% names(mrna)][ensembl_tx_name %in% final[overlaps_other_cds == FALSE]$tx_ids][!(ensembl_tx_name %in% final_predicted$tx_ids),]
