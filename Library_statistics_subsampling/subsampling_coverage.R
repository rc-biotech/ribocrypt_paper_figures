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

# Similarity to all merged

## coverage % number of libraries merged

# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- 20
n_genes <- 20
min_covs <- c(0, 2, 10, 100, 1000) # Minimum count to filter

genes <- c("RPL7", "TEX101", "ATF4")
genes <- unique(c(genes, sample(collection_genes$gene_symbol,
                                size = n_genes - length(genes), replace = FALSE)))

cov_all_list <- lapply(genes, function(gene) {
  gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
  cov_test_codons <- coverage_repeat_samplings_fst(nt_to_codons_fst(df_all, gene_path_fst),
                                                   subset_sizes, samplings, min_covs)
})
names(cov_all_list) <- genes
# cov figure
cov_merged <- merged_cor_subsamplings(cov_all_list)

cov_merged[, min_count := rep(rep(rep(min_covs, each = length(subset_sizes)), samplings), length(cov_all_list))]
cov_merged[min_count %in% c(0, 2), min_count := min_count + 1]
cov_merged[, min_count := as.factor(min_count)]
fwrite(cov_merged, file.path(plot_dir, "coverage_subsampling_table_genes.csv"))
cov_merged_mean <- cov_merged[, .(value = mean(value)), by = .(min_count, subsampling)]

cov_plot_line <- ggplot(data = cov_merged_mean) + theme_minimal() +
  geom_line(aes(y = value, x = subsampling, group = min_count, color = min_count)) + labs("Min. Count") +
  ylab("Coverage (%)") + xlab(paste("# Libraries (of", samplings, "subsamplings)")); cov_plot_line

ggsave(file.path(plot_dir, "coverage_subsampling.png"), cov_plot_line, width = 6, height = 6)
ggsave(file.path(plot_dir, "coverage_subsampling.pdf"), cov_plot_line, width = 6, height = 6)

# Replicates vs Size
ct_all_h <- countTable(df_all, "mrna", type = "summarized")[filterTranscripts(df_all, longestPerGene = T)]
ct_all_h <- as.data.table(assay(ct_all_h))
libsizes <- colSums(ct_all_h)
libnames <- bamVarName(df_all)
names(libsizes) <- libnames

libsizes <- sort(libsizes)
summary(libsizes) / 1e6

lib_size_total <- c(10, 20, 100, 500)
lib_size_mean <- c(2, 6, 20, 50)
# x = lib_size_total, y = coverage (%), group = lib_size_mean
x <- data.table(id = seq_along(libsizes), value = libsizes, libnames = libnames)

sum_game <- function(x, attempts = 50, wanted = 100000, fudge = .1, draws = 45) {
  x <- x[value <= wanted,]
  stopifnot(nrow(x) >= draws)
  for (i in seq(attempts)) {
    e <- sample(x = x$id, size = draws)
    d <- sum(x$value[x$id %in% e])
    if (d < (wanted + (wanted * fudge)) & d > (wanted - (wanted * fudge))) {
      print(paste("I found a sum of", d, "on attempt", i, "(", draws, "draws)"))
      d <- x$value[x$id %in% e]
      names(d) <- e
      return(d)
    } else {
      d <- NULL
      e <- NULL
    }
  }
  print("Abject failure")
  return(e)
}

subset_sizes <- c(1, 3, 10)
samplings <- 10
#solution <- sum_game(x, wanted = 2e6, draws = subset_sizes[1], attempts = 20000)
total_counts <- c(2, 6, 20, 100) * 1e6
genes_3 <- c("RPL7", "TEX101", "ATF4")
genes_sel <- unique(c(genes_3, sample(collection_genes$gene_symbol,
                                size = 400 - length(genes_3), replace = FALSE)))

cov_merged_mean <- lapply(total_counts,  function(total_count) {
  message("- Total counts: ", total_count)
cov_all_list_lib <- lapply(genes_sel, function(gene) {
  gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
  cov_test_codons <- coverage_repeat_samplings_fst(nt_to_codons_fst(df_all, gene_path_fst),
                                                   subset_sizes,
                                                   samplings, min_covs, x, total_count)
})

names(cov_all_list_lib) <- genes_sel
# cov figure
cov_merged <- merged_cor_subsamplings(cov_all_list_lib)
cov_merged[, min_count := rep(rep(rep(min_covs, each = length(subset_sizes)), samplings), length(genes_sel))]
cov_merged[min_count %in% c(0, 2), min_count := min_count + 1]
cov_merged[, min_count := as.factor(min_count)]
cov_merged_mean <- cov_merged[, .(value = mean(value)), by = .(min_count, subsampling)]
cov_merged_mean[, total_count := total_count]
return(cov_merged_mean)
})

cov_merged_mean_dt <- rbindlist(cov_merged_mean)
cov_plot_line <- ggplot(data = cov_merged_mean_dt) + theme_minimal() +
  geom_line(aes(y = value, x = min_count, group = subsampling, color = subsampling)) + labs("Min. Count") +
  ylab("Coverage (%)") + xlab(paste("# Minimum counts for coverage filter, total genes:", length(genes_sel),")")) +
  facet_wrap(~ total_count); cov_plot_line

ggsave(file.path(plot_dir, "coverage_subsampling_repVSsize.png"), cov_plot_line, width = 6, height = 6)



