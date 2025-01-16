library(ORFik); library(ggplot2)
source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
all_tables <- list.files("~/livemount/shared_results/predicted_orfs/Ribo_orfs_Homo_sapiens/", "_prediction_table.rds")
table <- readRDS("~/livemount/shared_results/predicted_orfs/Ribo_orfs_Homo_sapiens/uORF_uoORF_annotated_Homo_sapiens_all_merged-Homo_sapiens_RFP_prediction_table.rds")

df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
lib_sizes <- readRDS(file.path(QCfolder(df_all), "totalCounts_mrna.rds"))
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")

collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]

# PTPN7: 1047:1047
gene_path_fst <- load_allsamples_fst("PTPN7", df_all, collection_genes, isoform_index = 2)
point <- 1047
TIS <- as.numeric(start(pmapToTranscriptF(loadRegion(df_all, "cds")[names(gene_path_fst)], loadRegion(df_all, "mrna")[names(gene_path_fst)])))
point_iORF <- TIS + 13

gene_total <- RiboCrypt:::load_collection(gene_path_fst)
matchings <- RiboCrypt:::match_collection_to_exp(m, df_all)
term <- "CELL_LINE" #TISSUE
meta_sub <- m[matchings, term, with = FALSE][[1]]
region_type <- "leader+cds"
subset <- RiboCrypt:::subset_fst_coord_by_region(df_all, names(gene_path_fst), region_type)
gene <- RiboCrypt:::compute_collection_table(gene_path_fst, lib_sizes, df_all, term,
                                             normalization = "maxNormalized", value.var = "logscore", subset = subset,
                                             kmer = 9, metadata = m, min_count = 500, as_list = TRUE, ratio_interval = c(point, point))

table <- gene[[1]]
tail(colnames(table), 20)

meta_sub[match(colnames(table), levels(gene_total$library))]

ratios <- unlist(table[point, ] / colMeans(table[(point+100):nrow(table),]))
length(ratios)
plot(log2(ratios + 1))
plot(log2(gene[[2]] + 1))


order <- order(ratios)
color <- c("default (White-Blue)", "Matrix (black,green,red)")
color_index <- 2

cut_off_high <- (quantile(log2(ratios), 0.80)); cut_off_high
col_subset <- names(which(log2(ratios + 1) > cut_off_high))
sub_table_high <- table[, ..col_subset] #[,..order, which = FALSE]
plot_high <- RiboCrypt:::get_meta_browser_plot(sub_table_high, color[color_index], color_mult = 2)
plot_high

cut_off_low <- (quantile(log2(ratios), 0.6)); cut_off_low
cut_off_low_min <- (quantile(log2(ratios), 0.31)); cut_off_low_min
col_subset <- names(which(log2(ratios + 1) < cut_off_low & log2(ratios + 1) > cut_off_low_min))
sub_table_low <- table[, ..col_subset]
plot_low <- RiboCrypt:::get_meta_browser_plot(sub_table_low, color[color_index], color_mult = 2)
plot_low

grob_high <- grid::grid.grabExpr(ComplexHeatmap::draw(plot_high))
grob_low <- grid::grid.grabExpr(ComplexHeatmap::draw(plot_low))

# plot <- RiboCrypt:::get_meta_browser_plot(table[, order, with = F], color[color_index], color_mult = 2)
# plot

# Enrichment
metadata_value <- meta_sub[match(colnames(table), levels(gene_total$library))]
score <- log2(ratios + 1) > cut_off_high
meta_dt <- data.table(grouping = metadata_value, order = seq_along(metadata_value), cluster = ifelse(score, "uORF up", "uORF down"))
enrich_dt_temp <- RiboCrypt:::allsamples_meta_stats(meta_dt)
enrich <- copy(enrich_dt_temp)

{
enrich_dt <- as.data.table(enrich, keep.rownames = T)
enrich_dt <- suppressWarnings(melt(enrich_dt))
enrich_dt[, `:=`(variable, factor(as.character(variable)))]
enrich_dt <- enrich_dt[variable != "uORF down",]
enrich_dt <- enrich_dt[abs(value) > 2,]
enrich_dt <- enrich_dt[rn != "", ]
text_size <- 13
text_size_title <- text_size + 10
enrichment_plot <- ggplot(enrich_dt) + geom_bar(aes(x = rn,
                                                    y = value, fill = variable), stat = "identity", position = position_dodge()) +
  theme_minimal() + labs(fill = "Cluster") + xlab("Cell line") +
  ylab("Enrichment") + geom_hline(yintercept = c(3, -3),
                                  linetype = "dashed", color = "red", linewidth = 1) +
  theme(axis.title = element_text(size = 32), axis.text.x = element_text(size = (text_size),
                                                                         angle = 45), axis.text.y = element_text(size = (text_size)),
        legend.text = element_text(size = text_size), legend.title = element_text(size = 32))
enrichment_plot
}

geneModelPanelPlotFst <- function(df, id, region_type, render_text = id, interval = NULL, names_interval = "T", shift_cds_start = 0) {
  lengths <- ORFik:::optimizedTranscriptLengths(df)
  length <- lengths[tx_name == id]
  # Custom UTR lengths for Sac cer (yeast)
  if (organism(df) == "Saccharomyces cerevisiae")
    length$utr5_len <- 650
  tx_width <- ncol(heatmap) # length$tx_len

  start <- 1
  end <- tx_width

  if (length$cds_len > 0 & region_type %in% c("mrna", "leader+cds")) {
    start <- start + length$utr5_len
    end <- start + length$cds_len
    if (shift_cds_start != 0) start <- start + shift_cds_start
  }
  grl <- GRangesList(GRanges("1", IRanges(start, end)))
  names(grl) <- render_text

  ranges <- unlistGrl(grl)
  ranges <- c(GRanges("1", IRanges(1, tx_width)), ranges)
  if (!is.null(interval)) {
    target <- GRanges("1", IRanges(interval))
    names(target) <- names_interval
    ranges <- c(ranges, target)
  }
  dt <- RiboCrypt:::geneBoxFromRanges(ranges, tx_width,
                          cols = c("#FFFFFF", c("#F8766D","#00BA38","#619CFF")[start(ranges[-1]) %% 3 + 1]))[[1]]
  gene_model_panel <- RiboCrypt:::geneModelPanelPlot(dt)
}

summary_track_allsamples <- function(m, summary_track_type = "columns", as_plotly = FALSE) {

  summary_profile <- data.table(count = rowSums(m))
  summary_profile[, `:=`(position = seq.int(.N)) ]
  summary_profile[, `:=`(frame = factor((position-1) %% 3)) ]
  summary_plot <- RiboCrypt:::createSinglePlot(summary_profile, TRUE, 1, "",
                                               FALSE, lines = NULL,
                                               type = summary_track_type,
                                               flip_ylabel = FALSE, as_plotly = as_plotly)
  if (!as_plotly) {
    summary_plot +
      theme(axis.text.y=element_blank(),  axis.ticks.y=element_blank())
  }
  return(summary_plot)
}

shift <- 27
gene_model_panel <- geneModelPanelPlotFst(df_all, names(gene_path_fst), region_type, render_text = "PTPN7 (CDS)", c((point:(point+12))-shift, ((point_iORF-shift):(point_iORF+32))), c("uORF", "iORF"), -shift)

final_plot <- cowplot::plot_grid(plotlist = list(NULL, grob_low, NULL, grob_high, gene_model_panel, enrichment_plot), labels = c("", "uORF down", "","uORF up", NULL, NULL),
                                 ncol = 1, rel_heights = c(0.015, 0.335, 0.015, 0.335, 0.05, 0.25), greedy = FALSE,  vjust = 0)
final_plot
# # final_plot
# ggsave(final_plot, filename = "~/livemount/shared_results/ribocrypt_paper/Figure3_PTPN7.svg", width = 10, height = 10)
# browseURL("~/livemount/shared_results/ribocrypt_paper/Figure3_PTPN7.jpg")
#
# plot_all <- RiboCrypt:::get_meta_browser_plot(table, color[color_index], color_mult = 2)
# grob_all <- grid::grid.grabExpr(ComplexHeatmap::draw(plot_all))
# final_plot_all <- cowplot::plot_grid(plotlist = list(grob_all, gene_model_panel, enrichment_plot),
#                                  ncol = 1, rel_heights = c(0.7, 0.05, 0.25), greedy = FALSE,  vjust = 0)
# # final_plot_all
# ggsave(final_plot_all, filename = "~/livemount/shared_results/ribocrypt_paper/Figure3_PTPN7_all.svg", width = 10, height = 10)
# # browseURL("~/livemount/shared_results/ribocrypt_paper/Figure3_PTPN7_all.jpg")
