# Single gene metabrowser viewer

source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ggplot2)
library(gridExtra)
library(plotly)
# devtools::install_github("m-swirski/RiboCrypt")
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"

collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]



## Samplings correlation to all_merged

# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- 100

# RPL7
gene <- c("RPL7")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "count",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst))

# Animation
subsampling_fst_animate(table_cds_nt_tex101, subset_sizes = subset_sizes, normalize = FALSE)%>% RiboCrypt:::lineDeSimplify()


# Correlation tests (nt level)
cor_test_nt <- correlation_repeat_samplings_fst(table_cds_nt, subset_sizes, samplings)

# Correlation tests (codon level)

cor_test_codons <- correlation_repeat_samplings_fst(nt_to_codons_fst(df_all, gene_path_fst), subset_sizes, samplings)


# TEX101
gene <- c("TEX101")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "logscore",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt_tex101 <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst))

# Animation
subsampling_fst_animate(table_cds_nt_tex101, subset_sizes = subset_sizes, normalize = FALSE)

# Correlation tests (nt level)
cor_test_nt_tex101 <- correlation_repeat_samplings_fst(table_cds_nt_tex101, subset_sizes, samplings)

# Correlation tests (codon level)
table_long <- RiboCrypt:::load_collection(gene_path_fst)
table_cds <- subset_fst_by_region(df_all, table_long, id = names(gene_path_fst))
table_cds[, position := position - min(position) + 1]
table_cds[, codon := floor((position-1) / 3)]
table_cds_codons <- table_cds[, .(count = sum(count)), by = .(library, codon)]
colnames(table_cds_codons)[2] <- "position"
cor_test_codons_tex101 <- correlation_repeat_samplings_fst(table_cds_codons, subset_sizes, samplings)

# ATF4
gene <- c("ATF4")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "logscore",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt_ATF4 <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst))

# Animation
subsampling_fst_animate(table_cds_nt_ATF4, subset_sizes = subset_sizes, normalize = T)

# Correlation tests (nt level)
cor_test_nt_ATF4 <- correlation_repeat_samplings_fst(table_cds_nt_ATF4, subset_sizes, samplings)

# Correlation tests (codon level)

cor_test_codons_ATF4 <- correlation_repeat_samplings_fst(table_cds_codons, subset_sizes, samplings)



# Final figure
cor_all_list <- list(RPL7_nt = cor_test_nt, RPL7_codon = cor_test_codons,
                     TEX101_nt = cor_test_nt_tex101, TEX101_codon = cor_test_codons_tex101,
                     ATF4_nt = cor_test_nt_ATF4, ATF4_codon = cor_test_codons_ATF4)
cor_merged <- merged_cor_subsamplings(cor_all_list)
cor_plot <- ggplot(data = cor_merged) + theme_minimal() +
  geom_boxplot(aes(y = value, x = subsampling, fill = variable)) +
  ylab("Correlation (to all merged)") + xlab("# Libraries (of 100 subsamplings)"); cor_plot

# Frame
table_100 <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                  normalization = "tpm", value.var = "count",
                                                  kmer = 1, metadata = m, min_count = 100,
                                                  as_list = FALSE, format = "long")
table_100[, `:=`(score_tpm = NULL, score = NULL, logscore = NULL)]
table_cds <- subset_fst_by_region(df_all, table_100, id = names(gene_path_fst))
table_cds[, frame := (position-1) %% 3]
table_cds[, codon := floor((position-1) / 3)]
table_cds_frame <- copy(table_cds)
table_cds_frame[, position := frame]
length(unique(table_100$library))

samplings_res_frame <- lapply(seq(samplings), function(x) {
  message("Sampling: ", x)
  subsets <- sample_libs_from_fst(table_cds_frame, subsets = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                                               1500, 2000, 2492))

  as_mat <- matrix(unlist(lapply(subsets, function(s) s$tx_norm_counts)),
                   ncol = length(subsets))
  colnames(as_mat) <- subset_sizes
  rownames(as_mat) <- c("1","2","3")
  as_mat
})


best_frame <- sapply(samplings_res_frame, function(x) which.max(x[1,]))
best_frame <- factor(names(best_frame), levels = subset_sizes)
best_frame <- as.factor(sort(table((best_frame)), decreasing = T))
plot(x = factor(names(best_frame), levels = names(best_frame), ordered = T), y = as.numeric(as.character(best_frame)))

# frame: all samples check
table_100_wide <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                       normalization = "tpm", value.var = "count",
                                                       kmer = 1, metadata = m, min_count = 100,
                                                       as_list = FALSE, format = "wide")
table_cds_wide <- subset_fst_by_region(df_all, table_100_wide, id = names(gene_path_fst))
frame_usage <- coverageSim::frame_usage(table_cds_wide)
a[, frame_ratio := V1 / sum(V1), by = .(library)]

# Samples as melted
table <- RiboCrypt:::load_collection(gene_path_fst)
aggregate <- table[, .(agg_sum = sum(count)), by = position]
plot(aggregate$position, aggregate$agg_sum)


# Ribo-seq ORF prediction downsampling

df_merged <- read.experiment("all_merged-Saccharomyces_cerevisiae")
remove.experiments(df_merged)
all_predictions_folder <- "~/livemount/shared_results/predicted_orfs"
result_folder <- ORFik:::riboORFsFolder(df_merged, all_predictions_folder)
ORF_categories_keep <- c("uORF", "uoORF", "annotated", "NTE", "NTT", "internal", "doORF", "dORF")
ORF_categories_keep <- c("annotated")
prefix_result <- paste(c(ORF_categories_keep, organism(df_merged), "all_merged"), collapse = "_")
res <- ORFik:::detect_ribo_orfs(df_merged, out_folder = result_folder, prefix_result,
                                ORF_categories_to_keep = ORF_categories_keep)
res_files <- list.files(result_folder, full.names = TRUE)
table_all <- readRDS(grep(pattern = "all_merged", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE))
table_all[, `:=`(study = "all_merged", predicted_on = "all_merged")]

# Study and Single lib
yeast_studies <- all_exp[organism == "Saccharomyces cerevisiae",]
yeast_studies <- yeast_studies[!grepl("all_merged|all_samples|RNA-seq|_modalities", name),]
aligned_folders <- file.path(ORFik::config()["bam"], yeast_studies$name, "aligned")
bigwig_folders <- file.path(aligned_folders, "bigwig")
yeast_studies <- yeast_studies[dir.exists(bigwig_folders),]

mrna <- loadRegion(df_merged, "mrna")
cds <- loadRegion(df_merged, "cds")
longestORF = FALSE
startCodon = startDefinition(1)
stopCodon = stopDefinition(1)
minimumLength = 0
orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df_merged, TRUE),
                         longestORF = longestORF, startCodon = startCodon,
                         stopCodon = stopCodon, minimumLength = minimumLength)

for (exp in yeast_studies$name[1:10]) {
  message(exp)
  message("- Study merged")
  df_study <- read.experiment(exp, validate = FALSE)
  study_path <- file.path(libFolder(df_study), "pshifted_merged", "RFP_merged.ofst")
  if (!file.exists(study_path)) {
    study_path <- file.path(libFolder(df_study), "pshifted_merged", "RFP.ofst")
    if (!file.exists(study_path)) {
      stop("Could not find merged study track for: ", exp)
    }
  }

  prefix_result <- paste(c(ORF_categories_keep, organism(df_merged), "by_study"), collapse = "_")
  res <- ORFik:::detect_ribo_orfs(df_study[1,], out_folder = result_folder, prefix_result,
                                  ORF_categories_to_keep = ORF_categories_keep,
                                  mrna = mrna, cds = cds, orf_sequences = orf_candidate_ranges,
                                  libraries = list(fimport(study_path)))
  message("- Single Library")
  df_sample <- df_study[1,]
  remove.experiments(df_sample)
  prefix_result <- paste(c(ORF_categories_keep, organism(df_sample), "by_sample"), collapse = "_")
  res <- ORFik:::detect_ribo_orfs(df_sample, out_folder = result_folder, prefix_result,
                                  ORF_categories_to_keep = ORF_categories_keep,
                                  mrna = mrna, cds = cds, orf_sequences = orf_candidate_ranges,
                                  libraries = list(fimport(filepath(df_sample[1,], type = "pshifted"))))

}

res_files <- list.files(result_folder, full.names = TRUE)
pred_studies_paths <- grep(pattern = "by_study", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE)
pred_studies <- rbindlist(lapply(pred_studies_paths, function(x) cbind(readRDS(x), study = x, predicted_on = "study")))

res_files <- list.files(result_folder, full.names = TRUE)
pred_sample_paths <- grep(pattern = "by_sample", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE)
pred_samples <- rbindlist(lapply(pred_sample_paths, function(x) cbind(readRDS(x), study = x, predicted_on = "single_library")))

stopifnot(length(pred_studies_paths) == length(pred_sample_paths))
dt_translon <- rbindlist(list(table_all, pred_studies, pred_samples))
dt_translon <- dt_translon[, .(prediction_orfs = sum(predicted)), by = .(study, predicted_on)]
dt_translon[, predicted_on := factor(predicted_on, levels = rev(c("all_merged", "study", "single_library")), ordered = TRUE)]
dt_translon[, predicted_on]

translon_plot <- ggplot(data = dt_translon, aes(x = predicted_on, y = prediction_orfs)) +
  geom_boxplot() +
  theme_minimal() + ylab("# Predicted ORFs") +  xlab(paste("# Libraries (of", length(pred_studies_paths), "subsamplings)")); translon_plot

Figure_subsampling <- do.call("grid.arrange", c(list(cor_plot, translon_plot), ncol = 1))

ggsave("~/livemount/shared_results/ribocrypt_paper/Figure_subsampling.png",
       Figure_subsampling, width = 6, height = 6)

# ATF4 uORF

gene <- c("ATF4")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "count",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt_uorf <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst),
                                            subset = loadRegion(df_all, "leaders")[names(gene_path_fst)])

# Animation
subsampling_fst_animate(table_cds_nt_uorf, subset_sizes = subset_sizes, normalize = T)
all_subsets <- rbindlist(sample_libs_from_fst(table_all_coord, subset_sizes))
all_subsets <- all_subsets[position %in% seq(80, 100),]
ggplot(all_subsets) + geom_raster(aes(x = position, y = as.factor(size), fill = counts))
gg_uorf <- ggplot(all_subsets) + geom_raster(aes(x = position, y = as.factor(size), fill = tx_norm_counts)) +
  scale_fill_gradientn(colours=(c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", 3))))
ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4.png"), gg_uorf, width = 6, height = 6)

# SAMD11 uORF

gene <- c("SAMD11")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "count",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt_uorf <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst),
                                          subset = loadRegion(df_all, "leaders")[names(gene_path_fst)])

# Animation
uorf_subset <- table_cds_nt_uorf[220:370,]
subsampling_fst_animate(uorf_subset, subset_sizes = c(subset_sizes,200,300,500,1000), normalize = F)
all_subsets <- rbindlist(sample_libs_from_fst(table_all_coord, subset_sizes))
all_subsets <- all_subsets[position %in% seq(250, 380),]
ggplot(all_subsets) + geom_raster(aes(x = position, y = as.factor(size), fill = counts))
gg_uorf <- ggplot(all_subsets) + geom_raster(aes(x = position, y = as.factor(size), fill = as.numeric(tx_norm_counts))) +
  scale_fill_gradientn(colours=(c("#000000", "#2CFA1F", "yellow2", rep("#FF2400", 3))))
ggsave(file.path(plot_dir, "coverage_subsampling_uorf2_atf4.png"), gg_uorf, width = 6, height = 6)

subsampling_fst_animate_col <- function(table, subset_sizes = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                 1500, 2000, 2500), normalize = TRUE) {
  all_subsets <- rbindlist(sample_libs_from_fst(table, subset_sizes))
  all_subsets[,frame := as.factor((position-1) %% 3)]
  browser()
  y_counts <- if(normalize) {
    all_subsets$tx_norm_counts
  } else all_subsets$counts
  plot <- ggplot(data = all_subsets, mapping = aes(x = position, y = y_counts, frame = size, color = frame)) + theme_bw() +
    geom_line(size = 1) + ggtitle(gene)
  anime <- ggplotly(plot) %>% animation_slider(currentvalue = list(prefix = "Samples ", font = list(color="red"))) %>%
    RiboCrypt:::lineDeSimplify()
  return(anime)
}

subsampling_fst_animate_col(uorf_subset, normalize = F)

ggplot() + theme_bw() + geom_col(aes(x = 1:151, y = apply(uorf_subset,1,sum), fill = as.factor(0:150 %% 3)))
ggplot() + geom_rect(mapping = aes(xmin = 1, ymin = 1, xmax = 2, ymax = 2, fill = 1))
p <- ggplot() + geom_rect(mapping = aes(xmin = rep(1:3,3) - 0.5, ymax = c(10,3,2,15,3,3,25,4,5), xmax = rep(1:3,3) + 0.5, ymin = rep(0,9), frame = c(1,1,1,2,2,2,3,3,3), fill = as.factor(rep(1:3,3))))
