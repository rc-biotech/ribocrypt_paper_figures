source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ggplot2)
library(gridExtra)
library(plotly)
# devtools::install_github("m-swirski/RiboCrypt")
all_exp <- list.experiments(validate = FALSE)
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
shared_res_dir <- "~/livemount/shared_results"; stopifnot(dir.exists(shared_res_dir))
all_predictions_folder <- file.path(shared_res_dir, "predicted_orfs")
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
rds_dir <- file.path(plot_dir, "rds_temp_objects")
dir.create(rds_dir, FALSE, TRUE); dir.create(plot_dir, FALSE, TRUE);
collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]

# Similarity to all merged (From Figure 4)

## Figure 4C correlation of libraries merged

# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- 100
genes <- c("RPL7", "TEX101", "ATF4")

cor_to_allmerged_rds <- file.path(rds_dir, "cor_to_allmerged.rds")

if (!file.exists(cor_to_allmerged_rds)) {
  cor_all_list <- lapply(genes, function(gene) {
    gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
    cov_test_codons <- correlation_repeat_samplings_fst(nt_to_codons_fst(df_all, gene_path_fst),
                                                        subset_sizes, samplings)
  })
  names(cor_all_list) <- genes

  cor_merged <- merged_cor_subsamplings(cor_all_list)
  saveRDS(cor_merged, cor_to_allmerged_rds)
}


cor_plot <- ggplot(data = cor_merged) + theme_minimal() +
  geom_boxplot(aes(y = value, x = subsampling, fill = variable), outlier.shape=NA) +
  ylab("Correlation (to all merged)") + xlab(paste("# Libraries (of", samplings, "subsamplings)")) +
  theme(legend.position="top") + guides(fill=guide_legend(title="Gene")); cor_plot

## Figure 4D Ribo-seq ORF prediction downsampling
species <- c("Saccharomyces_cerevisiae", "Homo_sapiens")
dt_translon_all <- data.table()

specie <- species[2]
df_merged <- read.experiment(paste0("all_merged-", specie), validate = FALSE)
suppressWarnings(remove.experiments(df_merged))
result_folder <- ORFik:::riboORFsFolder(df_merged, all_predictions_folder)
res_files <- list.files(result_folder, full.names = TRUE)

all_done_list <- grep(pattern = "prediction_table\\.rds", basename(res_files), value = T)
all_done_length <- length(all_done_list) == 21

if (!all_done_length) {
  # Merged run
  ORF_categories_keep <- c("annotated")
  if (specie == "Homo_sapiens") {
    txnames <- filterTranscripts(df_merged, longestPerGene = TRUE)

    mrna <- loadRegion(df_merged, "mrna", names.keep = txnames)
    cds <- loadRegion(df_merged, "cds", names.keep = txnames)
  } else {
    mrna <- loadRegion(df_merged, "mrna")
    cds <- loadRegion(df_merged, "cds")
  }

  longestORF = FALSE
  startCodon = startDefinition(1)
  stopCodon = stopDefinition(1)
  minimumLength = 0
  orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df_merged, TRUE),
                                  longestORF = longestORF, startCodon = startCodon,
                                  stopCodon = stopCodon, minimumLength = minimumLength)

  merged_not_done <- !any(grepl("_all_merged-", all_done_list))
  if (merged_not_done) {
    prefix_result <- paste(c(ORF_categories_keep, organism(df_merged), "all_merged"), collapse = "_")
    res <- ORFik:::detect_ribo_orfs(df_merged, out_folder = result_folder, prefix_result,
                                    ORF_categories_to_keep = ORF_categories_keep,
                                    orf_candidate_ranges = orf_candidate_ranges)
  }

  # Study and Single lib
  org_studies <- all_exp[organism == organism(df_merged),]
  org_studies <- org_studies[!grepl("all_merged|all_samples|RNA-seq|_modalities", name),]
  aligned_folders <- file.path(ORFik::config()["bam"], org_studies$name, "aligned")
  bigwig_folders <- file.path(aligned_folders, "bigwig")
  org_studies <- org_studies[dir.exists(bigwig_folders),]

  for (exp in org_studies$name[1:10]) {
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
                                    mrna = mrna, cds = cds, orf_candidate_ranges = orf_candidate_ranges,
                                    libraries = list(fimport(study_path)))
    message("- Single Library")
    df_sample <- df_study[1,]
    prefix_result <- paste(c(ORF_categories_keep, organism(df_sample), "by_sample"), collapse = "_")
    res <- ORFik:::detect_ribo_orfs(df_sample, out_folder = result_folder, prefix_result,
                                    ORF_categories_to_keep = ORF_categories_keep,
                                    mrna = mrna, cds = cds, orf_candidate_ranges = orf_candidate_ranges,
                                    libraries = list(fimport(filepath(df_sample[1,], type = "pshifted"))))
  }
}


res_files <- list.files(result_folder, full.names = TRUE)
table_all <- readRDS(grep(pattern = "all_merged", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE))[type %in% ORF_categories_keep,]
if (specie == "Homo_sapiens") {
  table_all <- table_all[ensembl_tx_name %in% filterTranscripts(df_merged, longestPerGene = TRUE),]
}
table_all[, `:=`(study = "all_merged", predicted_on = "all_merged")]

pred_studies_paths <- grep(pattern = "by_study", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE)
pred_studies <- rbindlist(lapply(pred_studies_paths, function(x) cbind(readRDS(x), study = x, predicted_on = "study")))

pred_sample_paths <- grep(pattern = "by_sample", grep("prediction_table\\.rds$", res_files, value = T), value = TRUE)
pred_samples <- rbindlist(lapply(pred_sample_paths, function(x) cbind(readRDS(x), study = x, predicted_on = "single_library")))

stopifnot(length(pred_studies_paths) == length(pred_sample_paths) & length(pred_studies_paths) != 0)
dt_translon <- rbindlist(list(table_all, pred_studies, pred_samples))
dt_translon <- dt_translon[, .(prediction_orfs = sum(predicted)), by = .(study, predicted_on)]
dt_translon[, predicted_on := factor(predicted_on, levels = rev(c("all_merged", "study", "single_library")), ordered = TRUE)]
dt_translon[, predicted_on]
all_cds <- length(cds)
dt_translon[, prediction_orfs := round((prediction_orfs / all_cds)*100, 2)]
dt_translon[, species := organism(df_merged)]

dt_translon_all <- rbindlist(list(dt_translon_all, dt_translon))


dt_translon_all[, coloring := species]
dt_translon_all[predicted_on != "all_merged", coloring := NA]
translon_plot <- ggplot(data = dt_translon_all, aes(x = predicted_on, y = prediction_orfs, fill = species, color = coloring)) +
  geom_boxplot() +
  theme_minimal() + ylab("Predicted CDS' (%)") +  xlab(paste("# Libraries (of", length(pred_studies_paths), "subsamplings)")) +
  theme(legend.position="top") + guides(color="none") +
  guides(fill=guide_legend(title="Species")); plot(translon_plot)




## Merged plots

Figure_subsampling <- do.call("grid.arrange", c(list(cor_plot,
                                                     translon_plot),
                                                ncol = 1))

ggsave(file.path(plot_dir, "Figure_4C&D_subsampling.png"),
       Figure_subsampling, width = 6, height = 6, dpi = 400)
ggsave(file.path(plot_dir, "Figure_4C&D_subsampling.svg"),
       Figure_subsampling, width = 6, height = 6, dpi = 400)


