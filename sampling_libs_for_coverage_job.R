source("~/livemount/shared_scripts/fst_metabrowser_utils.R")
library(ggplot2)
library(gridExtra)
library(plotly)
library(dplyr)
# devtools::install_github("m-swirski/RiboCrypt")
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
data_dir <- "~/livemount/shared_results/sampling_coverage_randomlibs/"

collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]



## coverage % number of libraries merged

# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- 20
n_genes <- 20
min_covs <- c(0, 2, 10, 100, 1000) # Minimum count to filter

# # genes <- c("RPL7", "TEX101", "ATF4")
genes_sel <-  collection_genes$gene_symbol %>% unique
genes_sel <- genes_sel[sample(seq_along(genes_sel), size = length(genes_sel), replace = FALSE)]
libs_sampling_files <- system(paste0("ls ", data_dir, "*csv"),intern = TRUE)

genes_computed <- gsub(".*/","",libs_sampling_files) %>% sub(".csv","",.)
genes_sel <- genes_sel[-which(genes_sel %in% genes_computed)]
# genes_sel <- "ENSG00000278817"
cov_all_list <- lapply(genes_sel, function(gene) {try({
  gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
  table_long <- RiboCrypt:::load_collection(gene_path_fst)
  table_cds <- subset_fst_by_region(df_all, table_long, id = names(gene_path_fst))
  table_cds[, position := position - min(position) + 1]
  table_cds[, codon := floor((position-1) / 3)]
  table_cds_codons <- table_cds[, .(count = sum(count)), by = .(library, codon)]
  colnames(table_cds_codons)[2] <- "position"
  cov_test_codons <- coverage_repeat_samplings_fst(table_cds_codons, subset_sizes, samplings, min_covs)
  cov_test_codons <- as.data.table(cov_test_codons)
  cov_test_codons <- melt(cov_test_codons)
  cov_test_codons[,nlibs := rep(subset_sizes, .N/length(subset_sizes))]
  fwrite(cov_test_codons, paste0(data_dir,gene, ".csv"))
  return(cov_test_codons)
}
)
})
