source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ggplot2)
library(gridExtra)
library(plotly)
# devtools::install_github("m-swirski/RiboCrypt")
# all_exp <- list.experiments(validate = FALSE)

m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
shared_res_dir <- "~/livemount/shared_results"; stopifnot(dir.exists(shared_res_dir))
all_predictions_folder <- file.path(shared_res_dir, "predicted_orfs")
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
rds_dir <- file.path(plot_dir, "rds_temp_objects")
dir.create(rds_dir, FALSE, TRUE); dir.create(plot_dir, FALSE, TRUE);
translon_file <- file.path(rds_dir, "ORF_prediction_final.rds")

# Similarity to all merged (From Figure 4)


# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- seq(20)
indices_max <- seq(20001) # Number of ORFs to use (max per species, yeast has ~ 5k max)
ORF_categories_keep <- c("annotated", "uORF", "dORF")
longestORF = FALSE
startCodon = startDefinition(1)
stopCodon = stopDefinition(1)
minimumLength = 0



## Figure 4D: Ribo-seq ORF prediction downsampling
species <- c("Saccharomyces_cerevisiae", "Homo_sapiens")
if (file.exists(translon_file)) {
  dt_translon_all <- readRDS(translon_file)
} else {
  dt_translon_all <- data.table()
  for (specie in species) {
    message("Species: ", specie)

    df_merged <- read.experiment(paste0("all_merged-", specie), validate = FALSE)

    if (specie == "Homo_sapiens") {
      txnames <- filterTranscripts(df_merged, longestPerGene = TRUE)

      mrna <- loadRegion(df_merged, "mrna", names.keep = txnames)
      cds <- loadRegion(df_merged, "cds", names.keep = txnames)
    } else {
      mrna <- loadRegion(df_merged, "mrna")
      mrna <- extendLeaders(mrna, 600)
      cds <- loadRegion(df_merged, "cds")
    }

    result_folder <- ORFik:::riboORFsFolder(df_merged, all_predictions_folder)
    res_files <- list.files(result_folder, full.names = TRUE)
    all_done_list <- grep(pattern = "prediction_table\\.rds", basename(res_files), value = T)
    merged_not_done <- !any(grepl("_all_merged-", all_done_list))
    if (merged_not_done) {
      orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df_merged, TRUE),
                                      longestORF = longestORF, startCodon = startCodon,
                                      stopCodon = stopCodon, minimumLength = minimumLength)

      prefix_result <- paste(c(ORF_categories_keep, organism(df_merged), "all_merged"), collapse = "_")
      res <- ORFik:::detect_ribo_orfs(df_merged, out_folder = result_folder, prefix_result,
                                      ORF_categories_to_keep = ORF_categories_keep,
                                      orf_candidate_ranges = orf_candidate_ranges)
    }
    res_files <- list.files(result_folder, full.names = TRUE)
    table_all <- readRDS(grep(pattern = "all_merged", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE))
    grl_all <- readRDS(grep(pattern = "all_merged", grep("candidates\\.rds$", res_files, value = T), value = TRUE))

    if (specie == "Homo_sapiens") {
      table_all <- table_all[ensembl_tx_name %in% filterTranscripts(df_merged, longestPerGene = TRUE),]
      grl_all <- grl_all[names(grl_all) %in% filterTranscripts(df_merged, longestPerGene = TRUE)]
    }
    stopifnot((nrow(table_all) == length(grl_all)) & length(grl_all) > 0)
    table_all[, `:=`(study = "all_merged", predicted_on = "all_merged")]
    subset <- (table_all$type %in% ORF_categories_keep) & table_all$predicted
    table_all_pred <- table_all[subset,]
    message("- Running for prediction groups:")
    print(table(table_all_pred$type))
    message("- Unique genes:")
    print(length(unique(table_all_pred$ensembl_tx_name)))


    # FST sampling
    df_all <- read.experiment(paste0("all_samples-", specie), validate = FALSE)
    collection_folder <- file.path(resFolder(df_all), "collection_tables/")
    collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
    collection_genes[, gene_symbol := gsub("-.*", "", label)]
    subset_fst <- subset & (table_all$ensembl_tx_name %in% collection_genes$value)
    table_all_pred_fst <- table_all[subset_fst,]
    message("- Unique genes with fst:")
    unique_genes <- length(unique(table_all_pred_fst$ensembl_tx_name))
    unique_ORFs <- length(unique(table_all_pred_fst$id))
    print(unique_genes)
    message("- Unique ORFs with fst:")
    print(unique_ORFs)
    if (length(indices_max) > unique_ORFs) {
      message("Note: Subsetting ORFs to max ORFs")
    }
    indices <- indices_max[indices_max <= unique_ORFs]
    genes <- table_all_pred_fst$ensembl_tx_name[indices]

    mrna_fst <- mrna[genes]
    grl_fst <- grl_all[subset_fst]
    grl_fst <- grl_fst[indices]
    stopifnot(length(unique(names(mrna_fst))) == length(unique(names(grl_fst))))
    orf_start_gr <- startSites(grl_fst, TRUE, TRUE, TRUE)
    upstream_gr <- windowPerGroup(orf_start_gr, mrna, 20, -2)
    valid_upstream <- as.logical(widthPerGroup(upstream_gr) > 0)
    message("- Loading gene fsts")
    bins <- cut_width(seq_along(genes), 100, center = 51)

    all_bins <- lapply(levels(bins), function(bin) {
      message("-- Bin:", bin)
      file_grl <- file.path(rds_dir, paste0("ORF_merge_grl", specie, bin, "_subsamp_", length(samplings), ".fst"))
      file_upstream <- file.path(rds_dir, paste0("ORF_merge_upstream", specie, bin, "_subsamp_", length(samplings), ".fst"))
      if (all(file.exists(file_grl, file_upstream))) {
        message("--- Preloaded")
        return(list(setDT(fst::read_fst(file_grl)), setDT(fst::read_fst(file_upstream))))
      } else {
        index_bin <- indices[bins == bin]
        list <- lapply(index_bin, function(index) {
          id <- genes[index]
          names(id) <- index
          message("--- ID: ", id, if (!is.null(names(id))) " (", names(id), ")")
          gene_path_fst <- paste0(collection_folder, id, ".fst")
          table_long <- RiboCrypt:::load_collection(gene_path_fst)

          table_grl <- subset_fst_by_region(df_all, table_long, id = id,
                                            gene_mrna = mrna_fst[index], subset = grl_fst[index])
          table_grl[, ORF := index]
          if (valid_upstream[index]) {
            table_up <- subset_fst_by_region(df_all, table_long, id = id,
                                             gene_mrna = mrna_fst[index], subset = upstream_gr[index])
            table_up <- table_up[, .(mean = mean(count), ORF = index), by = library]
          } else {
            table_up <- data.table(library = unique(table_grl$library), mean = 0, ORF = index)
          }

          list(table_grl, table_up)
        })
        message("---- Summary statistics for bin: ", bin)
        table_grls <- rbindlist(lapply(list, function(x) x[[1]]))
        table_ups <- rbindlist(lapply(list, function(x) x[[2]]))
        stopifnot(length(index_bin) == length(unique(table_ups$ORF)))
        orfs_cov_stats <- table_grls[, .(mean = mean(count), median = median(count), sum = sum(count)), by = .(ORF, library)]
        orfs_cov_stats[, reads_start := table_grls[table_grls[ , .I[position == min(position)], by = .(ORF, library)]$V1, count]]
        table_grls[, frame := rep(seq.int(0, 2), length.out = .N), by = .(ORF, library)]
        ORFscores <- dcast(table_grls[, .(frame_sum = sum(count)), by = .(ORF, library, frame)],
                           formula = ORF + library ~ frame, fun.aggregate = sum, value.var = "frame_sum")[,-c(1,2), with = FALSE]
        colnames(ORFscores) <- paste0("frame", colnames(ORFscores))
        orfs_cov_stats <- data.table(orfs_cov_stats, ORFscores)

        fst::write.fst(orfs_cov_stats, file_grl)
        fst::write.fst(table_ups, file_upstream)
        return(list(orfs_cov_stats, table_ups))
      }
    })
    orfs_cov_stats <- rbindlist(lapply(all_bins, function(x) x[[1]]))
    upstream_cov_stats <- rbindlist(lapply(all_bins, function(x) x[[2]]))
    stopifnot(nrow(orfs_cov_stats) == nrow(upstream_cov_stats) & nrow(upstream_cov_stats) != 0)


    # Merged run

    subset_sizes_use <- c(subset_sizes[subset_sizes < nrow(df_all)], nrow(df_all))
    message("Starting simulation")
    sample_list <- lapply(subset_sizes_use, function(sub) {
      message("\nSubset size: ", sub)
      cat("-- Sampling nr. ")
      dt <- rbindlist(lapply(samplings, function(samp) {
        cat(paste(samp, ", "))
        sample_lib <- sample(nrow(df_all), size = sub, replace = FALSE)
        sample_cov <- orfs_cov_stats[as.integer(library) %in% sample_lib, ]
        sample_cov_merged <- sample_cov[, .(reads_start = sum(reads_start), mean = sum(mean), median = sum(median)),
                                        by = ORF]
        sample_cov_merged[, upstream_mean := upstream_cov_stats[as.integer(library) %in% sample_lib,][, .(mean = sum(mean)), by = ORF]$mean]
        ORFscores <- sample_cov[, .(frame0 = sum(frame0), frame1 = sum(frame1), frame2 = sum(frame2), sum = sum(sum)), by = ORF]
        ORFscores[, `:=`(frame0_wins = (frame0 > frame1) & (frame0 > frame2))]
        ORFscores[, ORFscore := log2(rowSums((ORFscores[, .(frame0, frame1, frame2)] - sum)^2 / sum) + 1)]
        ORFscores[frame0_wins == FALSE, ORFscore := -1*ORFscore]
        ORFscores[is.nan(ORFscore), ORFscore := 0]

        sample_cov_merged[, `:=`(ORFScores = ORFscores$ORFscore, sum = ORFscores$sum)]

        sample_cov_merged[, predicted := (reads_start > 3) & (sum > 10) &
                            (mean > (upstream_mean * 1.3)) &
                            (ORFScores > 2.5) &
                            ((reads_start + 3) > median)]

        data.table(predicted = sum(sample_cov_merged$predicted), sampling_nr = samp)
      }
      ))
      dt[, subset_size := sub]
      return(dt)
    })

    dt <- rbindlist(sample_list)
    dt[, subset_size := as.factor(subset_size)]
    dt[, species := specie]
    dt[, prediction_orfs := round((predicted / length(indices))*100, 2)]

    # dt_translon_all <- dt_translon_all[species != "Homo_sapiens",]
    dt_translon_all <- rbindlist(list(dt_translon_all, dt))
  }
  saveRDS(dt_translon_all, translon_file)
}


translon_plot <- ggplot(data = dt_translon_all[as.numeric(as.character(subset_size)) <= 1000,], aes(x = subset_size, y = prediction_orfs, fill = species)) +
  geom_boxplot(fatten=0.7, outlier.size = 0.3) +
  theme_minimal() + ylab("Predicted ORFs (%)") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
  theme(legend.position="top") + guides(color="none") +
  guides(fill=guide_legend(title="Species")); plot(translon_plot)

translon_plot_orfs <- ggplot(data = dt_translon_all, aes(x = subset_size, y = predicted, fill = species, color = species)) +
  geom_boxplot() +
  theme_minimal() + ylab("# Predicted ORFs") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
  theme(legend.position="top") + guides(color="none") +
  guides(fill=guide_legend(title="Species")); plot(translon_plot_orfs)

# Percentage output
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred.png"),
       translon_plot, width = 5, height = 3, dpi = 400)
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred.jpg"),
       translon_plot, width = 5, height = 3, dpi = 400)
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred.svg"),
       translon_plot, width = 6, height = 6, dpi = 400)

# Numeric output
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred_numeric.png"),
       translon_plot_orfs, width = 5, height = 3, dpi = 400)
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred_numeric.jpg"),
       translon_plot_orfs, width = 5, height = 3, dpi = 400)
ggsave(file.path(plot_dir, "Figure_4D_subsampling_orf_pred_numeric.svg"),
       translon_plot_orfs, width = 6, height = 6, dpi = 400)

# Split by species
species_finished <- unique(dt_translon_all$species)
for (s in species_finished) {
  translon_plot_orfs <- ggplot(data = dt_translon_all[species == s,], aes(x = subset_size, y = predicted, fill = species, color = species)) +
    geom_boxplot() +
    theme_minimal() + ylab("# Predicted ORFs") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
    theme(legend.position="top") + guides(color="none") +
    guides(fill=guide_legend(title="Species")); plot(translon_plot_orfs)
  ggsave(file.path(plot_dir, paste0("Figure_4D_subsampling_orf_pred_numeric_", s, ".png")),
         translon_plot_orfs, width = 5, height = 3, dpi = 400)
  ggsave(file.path(plot_dir, paste0("Figure_4D_subsampling_orf_pred_numeric_", s, ".jpg")),
         translon_plot_orfs, width = 5, height = 3, dpi = 400)
  ggsave(file.path(plot_dir, paste0("Figure_4D_subsampling_orf_pred_numeric_", s, ".pdf")),
         translon_plot_orfs, width = 5, height = 3, dpi = 400)

}

# BY CDS
dt <- readRDS("~/livemount/shared_results/predicted_orfs/Ribo_orfs_Homo_sapiens/uORF_uoORF_annotated_Homo_sapiens_all_merged-Homo_sapiens_RFP_prediction_table.rds")
df <- read.experiment("all_merged-Homo_sapiens")
values <- c(length(filterTranscripts(df, 0,1,0,longestPerGene = TRUE)),
            length(unique(dt[type == "annotated",]$ensembl_gene_id)),
            length(unique(dt[type == "annotated" & predicted == TRUE,]$ensembl_gene_id)))
names(values) <- c("all genes", "10 reads genes", "predicted genes")
values_rel_human <- round((values / values[1])*100, 2)
values_rel_human # As percentage
values # Raw numbers

dt <- readRDS("~/livemount/shared_results/predicted_orfs/Ribo_orfs_Saccharomyces_cerevisiae/annotated_uORF_dORF_Saccharomyces cerevisiae_all_merged_all_merged-Saccharomyces_cerevisiae_RFP_WT_prediction_table.rds")
df <- read.experiment("all_merged-Saccharomyces_cerevisiae")
values <- c(length(filterTranscripts(df, 0,1,0,longestPerGene = TRUE)),
            length(unique(dt[type == "annotated",]$ensembl_gene_id)),
            length(unique(dt[type == "annotated" & predicted == TRUE,]$ensembl_gene_id)))
names(values) <- c("all genes", "10 reads genes", "predicted genes")
values_rel_yeast <- round((values / values[1])*100, 2)
values_rel_yeast # As percentage
values # Raw numbers

dt_cds_normalized <- dt_translon_all[as.numeric(as.character(subset_size)) <= 1000,]
dt_cds_normalized[species == "Saccharomyces_cerevisiae", prediction_orfs := prediction_orfs * (values_rel_yeast["predicted genes"] / 100)]
dt_cds_normalized[species == "Homo_sapiens", prediction_orfs := prediction_orfs * (values_rel_human["predicted genes"] / 100)]
translon_plot_cds <- ggplot(data = dt_cds_normalized, aes(x = subset_size, y = prediction_orfs, fill = species)) +
  geom_boxplot(fatten=0.7, outlier.size = 0.3) +
  theme_minimal() + ylab("Predicted CDSs (%)") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
  theme(legend.position="top") + guides(color="none") +
  geom_hline(yintercept = c(values_rel_human["10 reads genes"], values_rel_yeast["10 reads genes"]), color = c("red", "#00BFC4"), linetype = 'dotted') +
  scale_y_continuous(limits = c(0, 96)) +
  guides(fill=guide_legend(title="Species")); plot(translon_plot_cds)

save_image_formats <- function(plot, dir, filename, formats = c(".svg", ".jpg"),
                               google_dir = NULL,
                               width = 6, height = 6, dpi = 400) {
  file_names <- paste0(filename, formats)
  for (file in file_names) {
    message(file)
    ggsave(file.path(dir, file),
           plot, width = width, height = height, dpi = dpi)
    if (!is.null(google_dir)) {
      googledrive::drive_put(file.path(dir, file),
                             google_dir,
                             name = file)
    }
  }
}
save_image_formats(translon_plot_cds, plot_dir, "Figure_4F_subsampling_cds_pred",
                   google_dir = "https://drive.google.com/drive/folders/1SYsf0cf-gwU9xFw0Kht7FJaMSOhWzHGY")





# glm_smooth <- function(...) {
#   geom_smooth(method = "glm", method.args = list(family = quasi(link = "log", variance = "mu")), ...)
# }

# translon_plot_glm <- ggplot(data = dt_translon_all[as.numeric(as.character(subset_size)) <= 1000,], aes(x = as.numeric(as.character(subset_size)), y = prediction_orfs, color = species)) +
#   glm_smooth() + theme_minimal() + ylab("Predicted ORFs (%)") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
#   theme(legend.position="top") + guides(color="none") +
#   guides(fill=guide_legend(title="Species")); plot(translon_plot_glm)
#
# p<-p+geom_ribbon(aes(ymin=data$lower, ymax=data$upper), linetype=2, alpha=0.1)
#   theme_minimal() + ylab("Predicted ORFs (%)") +  xlab(paste("# Libraries (of", length(samplings), "subsamplings)")) +
#   theme(legend.position="top") + guides(color="none") +
#   guides(fill=guide_legend(title="Species")); plot(translon_plot)

# Inspect
# table_all_pred_fst
# index <- 6
# table_all_pred_fst[index,]
# gene <- genes[index]
# gene_path_fst <- paste0(collection_folder, gene, ".fst")
# table_long <- RiboCrypt:::load_collection(gene_path_fst)
# all_subsets <- table_long[, .(count = sum(count)), by = position]
# all_subsets[, frame := factor(rep(seq.int(0,2), length.out = .N))]
# table_grl <- subset_fst_by_region(df_all, table_long, id = id,
#                                   gene_mrna = mrna_fst[index], subset = grl_fst[index])
#
#
# ggplot(data = all_subsets, mapping = aes(x = position, y = count)) +
#   geom_line(aes(color = frame)) + ggtitle(gene) + geom_vline(xintercept = (c(min(table_grl$position), max(table_grl$position))))









