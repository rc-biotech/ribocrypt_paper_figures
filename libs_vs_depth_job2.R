source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ggplot2)
library(gridExtra)
library(plotly)
library(dplyr)
# devtools::install_github("m-swirski/RiboCrypt")
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
data_dir <- "~/livemount/shared_results/sampling_coverage_multilib/"
collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]


ct_all_h <- countTable(df_all, "mrna", type = "summarized")[filterTranscripts(df_all, longestPerGene = T)]
ct_all_h <- as.data.table(assay(ct_all_h))
libsizes <- colSums(ct_all_h)
libnames <- bamVarName(df_all)
names(libsizes) <- libnames

libsizes <- sort(libsizes)
summary(libsizes) / 1e6

x <- data.table(id = seq_along(libsizes), value = libsizes, libnames = libnames)

subset_sizes <- c(1, 3, 5, 10, 20, 50, 100)
samplings <- 10
#solution <- sum_game(x, wanted = 2e6, draws = subset_sizes[1], attempts = 20000)
total_counts <- c(5, 20, 50, 100, 500, 1000) * 1e6
min_covs <- c(0, 2, 10, 100, 1000)
# genes_3 <- c("RPL7", "TEX101", "ATF4")
# genes_sel <- unique(c(genes_3, sample(collection_genes$gene_symbol,
#                                       size = 4 - length(genes_3), replace = FALSE)))
genes_sel <- collection_genes$gene_symbol %>% unique
genes_sel <- genes_sel[sample(seq_along(genes_sel), size = length(genes_sel), replace = FALSE)]
cov_merged_mean_pergene <- lapply(genes_sel,  function(gene) {
  all_samplings <- lapply(total_counts, function(total_count) {
    gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
    table_long <- RiboCrypt:::load_collection(gene_path_fst)
    table_cds <- subset_fst_by_region(df_all, table_long, id = names(gene_path_fst))
    table_cds[, position := position - min(position) + 1]
    table_cds[, codon := floor((position-1) / 3)]
    table_cds_codons <- table_cds[, .(count = sum(count)), by = .(library, codon)]
    colnames(table_cds_codons)[2] <- "position"
    cov_test_codons <- coverage_repeat_samplings_fst(table_cds_codons, subset_sizes,
                                                     samplings, min_covs, x, total_count)
    cov_test_codons <- as.data.table(cov_test_codons)
    cov_test_codons <- melt(cov_test_codons)
    cov_test_codons$depth <- as.character(total_count)
    cov_test_codons[,nlibs := rep(subset_sizes, .N/length(subset_sizes))]
    return(cov_test_codons)

  })
  all_samplings <- rbindlist(all_samplings)
  fwrite(all_samplings, paste0(data_dir,gene, ".csv"))
  return(all_samplings)
})
