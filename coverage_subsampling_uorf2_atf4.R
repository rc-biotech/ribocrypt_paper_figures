source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ORFik)

# Setup experiment and parameters
all_predictions_folder <- "~/livemount/shared_results/predicted_orfs"
df_merged_h <- read.experiment("all_merged-Homo_sapiens")
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
rds_dir <- file.path(plot_dir, "rds_temp_objects")
atf4_uorf2_dt.rds <- file.path(rds_dir, "atf4_uorf2_dt.rds")
metadata_used <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]
org_orf_out <- ORFik:::riboORFsFolder(df_merged_h, all_predictions_folder)
orf_pred_files <- list.files(org_orf_out, full.names = TRUE)
prediction_table_files <- grep("prediction_table\\.rds$", orf_pred_files, value = T)


if (!file.exists(atf4_uorf2_dt.rds)) {
  # Merged table

  table_all <- readRDS(grep(pattern = "all_merged", prediction_table_files, value = TRUE))
  table_all[, `:=`(study = "all_merged", predicted_on = "all_merged")]

  candidates_gr <- readRDS(file.path(all_predictions_folder, "Ribo_orfs_Homo_sapiens/uORF_uoORF_annotated_Homo_sapiens_all_merged-Homo_sapiens_RFP_candidates.rds"))

  is_predicted_uorf <- table_all$type %in% c("uorf", "uoORF") & table_all$predicted == TRUE

  all_predicted_uorfs <- candidates_gr[is_predicted_uorf]

  #Pick isoform from gene you want
  gene <- "ATF4"
  fst_isoform <- "ENST00000337304"
  table_all[(external_gene_name %in% gene), ]
  isoform <- "ENST00000674920"
  table_all[predicted == TRUE & ensembl_tx_name %in% isoform,]

  atf4_uorfs_predicted <- candidates_gr[table_all$predicted == TRUE & table_all$ensembl_tx_name %in% isoform,]
  atf4_uorfs_predicted04 <- atf4_uorfs_predicted[2]
  names(atf4_uorfs_predicted04) <- fst_isoform
  RiboCrypt::multiOmicsPlot_ORFikExp(loadRegion(df_merged_h, "mrna", names.keep = fst_isoform),
                                     df_merged_h,
                                     annotation = atf4_uorfs_predicted,
                                     reads = outputLibs(df_merged_h, type = "bigwig", output.mode = "envirlist", force = TRUE,
                                                        naming = "full", BPPARAM = BiocParallel::SerialParam()))

  # Now check fst

  gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
  table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                          normalization = "tpm", value.var = "count",
                                                          kmer = 1, metadata = metadata_used,
                                                          min_count = 0, as_list = TRUE)$table
  table_cds_nt_uorf <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst),
                                            subset = loadRegion(df_all, "leaders")[names(gene_path_fst)])
  table_cds_nt_uorf <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst),
                                            subset = atf4_uorfs_predicted04)
  table_uorf_flanks <- table_all_coord[seq(80, 100),]
  # Animation
  subset_sizes <- c(1, 3, 5, 10, 20, 50, 100, 500, 1000)
  subsampling_fst_animate(table_cds_nt_uorf, subset_sizes = subset_sizes, normalize = T)
  all_subsets <- rbindlist(lapply(seq(200), function(x) cbind(rbindlist(sample_libs_from_fst(table_uorf_flanks, subset_sizes)), sampling = x)))
  all_subsets[,position := position + 79]
  all_subsets <- all_subsets[, .(counts = min(counts)), by = .(position, size)]

  all_subsets[, tx_norm_counts_subset := counts / sum(counts + 1) , by = size]

  saveRDS(all_subsets, atf4_uorf2_dt.rds)
} else all_subsets <- readRDS(atf4_uorf2_dt.rds)


ggplot(all_subsets) + geom_raster(aes(x = position, y = as.factor(size), fill = counts))
gg_uorf <- ggplot(all_subsets[!(size %in% c(20, 50, 100))]) + geom_raster(aes(x = position, y = as.factor(size), fill = tx_norm_counts_subset)) + theme_classic() +
  scale_fill_gradientn(colours=(c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", 3)))) + ylab("Libraries sampled") + xlab("Position") +
  guides(fill=guide_legend(title="Normalized counts\n(Worst case sampling)")) + theme(legend.position="top"); gg_uorf
in_frame <- ranges(ORFik::pmapToTranscriptF(atf4_uorfs_predicted04, loadRegion(df_all)["ENST00000337304"]))
in_orf <- seq.int(from = as.integer(start(in_frame)), to = as.integer(end(in_frame)))
in_frame <- seq.int(from = as.integer(start(in_frame)), to = as.integer(end(in_frame)), 3)

# Summary with frame Frame
summary_track_type <- "columns"
summary_profile <- data.table(count = all_subsets[size == 1000,]$counts)
summary_profile[, `:=`(position = seq.int(.N)) ]
summary_profile[, `:=`(frame = factor((position-1) %% 3)) ]
summary_plot <- RiboCrypt:::createSinglePlot(summary_profile, TRUE, 1, "",
                                             FALSE, lines = NULL,
                                             type = summary_track_type,
                                             flip_ylabel = FALSE, as_plotly = FALSE) +
  theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank()) ; summary_plot
summary_plot <- summary_plot + theme_classic() + theme(legend.position = "none", axis.text=element_blank(),  axis.ticks=element_blank())
# Annotation
lines <- NULL
withFrames <- TRUE
colors <- TRUE
tx_width <- 21
start <- 1
end <- tx_width
start <- start + 6
end <- start  + 12 - 3
grl <- GRangesList(GRanges("1", IRanges(start, end)))
names(grl) <- "ATF4: uORF2"

ranges <- unlistGrl(grl)
ranges <- c(GRanges("1", IRanges(1, tx_width)), ranges)
dt <- RiboCrypt:::geneBoxFromRanges(ranges, tx_width,
                        cols = c("#FFFFFF", c("#F8766D","#00BA38","#619CFF")[start(ranges[-1]) %% 3 + 1]))[[1]]
gene_model_panel <- RiboCrypt:::geneModelPanelPlot(dt)

grob <- gg_uorf
to_use_logicals <- c(T, TRUE, T)
final_plot <- cowplot::plot_grid(plotlist = list(summary_plot, grob, gene_model_panel)[to_use_logicals],
                                 ncol = 1, rel_heights = c(0.2, 0.75, 0.05)[to_use_logicals]); plot(final_plot)

ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4.png"), gg_uorf, width = 6, height = 6)
ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4.svg"), gg_uorf, width = 6, height = 6)

ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4_all.png"), final_plot, width = 6, height = 6)
ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4_all.svg"), final_plot, width = 6, height = 6)

ratio <- all_subsets[position %in% in_frame, sum(tx_norm_counts_subset), by = size]$V1 / all_subsets[!(position %in% in_frame), sum(tx_norm_counts_subset), by = size]$V1 + 1
names(ratio) <- as.factor(unique(all_subsets$size)); ratio
plot(ratio)

ratio <- all_subsets[position %in% in_frame, sum(tx_norm_counts_subset), by = size]$V1 / all_subsets[!(position %in% in_orf), sum(tx_norm_counts_subset), by = size]$V1 + 1
names(ratio) <- as.factor(unique(all_subsets$size)); ratio
plot(ratio)

ratio <- all_subsets[position %in% in_orf, sum(tx_norm_counts_subset), by = size]$V1 / all_subsets[!(position %in% in_orf), sum(tx_norm_counts_subset), by = size]$V1 + 1
names(ratio) <- as.factor(unique(all_subsets$size)); ratio
plot(ratio)


