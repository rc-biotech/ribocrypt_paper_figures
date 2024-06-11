# Single gene metabrowser viewer
library(ORFik)
library(data.table)
library(ggplot2)
library(gridExtra)
library(plotly)

df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")

collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]

load_allsamples_fst <- function(gene, df_all, collection_genes, use.names = TRUE) {
  collection_gene <- collection_genes[gene_symbol %in% gene,]
  collection_gene <- collection_gene[!duplicated(gene_symbol),]
  tx <- collection_gene$value
  id <- tx
  gene <- collection_gene$gene_symbol
  names(id) <- gene
  message("ID: ", id, if (!is.null(names(id))) " (", names(id), ")")
  collection_folder <- file.path(resFolder(df_all), "collection_tables/")
  table_path <- paste0(collection_folder, id, ".fst")

  if (!file.exists(dirname(table_path)))
    stop("There is no collection fst tables directory for this organism,",
         " see vignette for more information on how to make these.")
  if (!file.exists(table_path)) stop("Gene has no precomputed table, try another one!")
  if (use.names) names(table_path) <- id
  return(table_path)
}

subset_fst_by_region <- function(df_all, table, id,
                                 gene_mrna = loadRegion(df_all, names.keep = id),
                                 subset = loadRegion(df_all,part = "cds", names.keep = id),
                                 flank_cutoff = 0) {
  stopifnot(flank_cutoff >= 0)
  stopifnot(length(subset) == 1)
  stopifnot(length(gene_mrna) == 1)
  stopifnot(is(table, "data.table"))

  subset_txcoord <- pmapToTranscriptF(subset, gene_mrna)
  subset_pos <- seq(as.integer(start(subset_txcoord)) + flank_cutoff, as.integer(end(subset_txcoord)) - flank_cutoff)
  is_long_format <- all(c("library", "position") %in% colnames(table))
  if (is_long_format) {
    return(table[position %in% subset_pos,])
  }
  return(table[subset_pos,])
}

sample_libs_from_fst <- function(table,
                                 subsets = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                             1500, 2000, 2500)) {
  is_long_format <- all(c("library", "position") %in% colnames(table))
  if (is_long_format) {
    nlibs <- length(unique(table$library))
    subsets <- lapply(subsets, function(s) {
      dt <- table[as.integer(library) %in% sample(seq(nlibs), s, replace = FALSE),
                  .(counts = sum(count)), by = position]
      dt[, size := s]
      dt[, tx_norm_counts := counts / sum(counts)]
    })
  } else {
    nlibs <- ncol(table)
    subsets <- lapply(subsets, function(s) {
      dt <- data.table(counts = rowSums(table[, sample(seq(nlibs), s, replace = FALSE), with = FALSE]),
                       position = seq.int(nrow(table)), size = s)
      dt[, tx_norm_counts := counts / sum(counts)]
    })
  }
  return(subsets)
}

correlation_repeat_samplings_fst <- function(table, subset_sizes = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                                                     1500, 2000, 2500),
                                             samplings = 100, method = "spearman",
                                             verbose = TRUE) {

  samplings_res <- lapply(seq(samplings), function(x) {
    if (verbose) message("Sampling: ", x)
    subsets <- sample_libs_from_fst(table, subset_sizes)
    cor_test <- matrix(unlist(lapply(subsets, function(s) s$tx_norm_counts)),
                       ncol = length(subset_sizes))
    colnames(cor_test) <- subset_sizes
    cor(cor_test, method = method)
  })
  cor_test <- do.call(cbind, lapply(samplings_res, function(x) x[,length(subset_sizes)]))
  message("Done, summary correlation mean statistc for subsets:")
  print(rowMeans(cor_test, na.rm = T))
  return(cor_test)
}

subsampling_fst_animate <- function(table, subset_sizes = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                                            1500, 2000, 2500), normalize = TRUE) {
  all_subsets <- rbindlist(sample_libs_from_fst(table, subset_sizes))
  y_counts <- if(normalize) {
    all_subsets$tx_norm_counts
  } else all_subsets$counts
  plot <- ggplot(data = all_subsets, mapping = aes(x = position, y = y_counts, frame = size)) +
    geom_line() + ggtitle(gene)
  anime <- ggplotly(plot) %>% animation_slider(currentvalue = list(prefix = "Samples ", font = list(color="red"))) %>%
    RiboCrypt:::lineDeSimplify()
  return(anime)
}

merged_cor_subsamplings <- function(list_cor) {
  if (length(names(list_cor)) != length(list_cor)) stop("length(names(list_cor)) != length(list_cor)")
  list_cor_merged <- lapply(seq(length(list_cor)), function(x) {
    cor_test_final <- as.data.table(list_cor[[x]])
    colnames(cor_test_final) <- rep(names(list_cor[x]), ncol(cor_test_final))
    cor_test_final <- cbind(subsampling = rownames(list_cor[[x]]), cor_test_final)
    cor_test_final <- suppressWarnings(melt(cor_test_final))
  })
  cor_merged <- rbindlist(list_cor_merged)
  cor_merged[, subsampling := factor(subsampling, levels = unique(cor_merged$subsampling))]
  return(cor_merged)
}

# Parameters
subset_sizes <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
samplings <- 100

for (gene in collection_genes$gene_symbol) {

gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
table_all_coord <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                        normalization = "tpm", value.var = "count",
                                                        kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
table_cds_nt <- subset_fst_by_region(df_all, table_all_coord, id = names(gene_path_fst))

# Correlation tests (nt level)
cor_test_nt <- correlation_repeat_samplings_fst(table_cds_nt, subset_sizes, samplings)

# Correlation tests (codon level)
table_long <- RiboCrypt:::load_collection(gene_path_fst)
table_cds <- subset_fst_by_region(df_all, table_long, id = names(gene_path_fst))
table_cds[, position := position - min(position) + 1]
table_cds[, codon := floor((position-1) / 3)]
table_cds_codons <- table_cds[, .(count = sum(count)), by = .(library, codon)]
colnames(table_cds_codons)[2] <- "position"
cor_test_codons <- correlation_repeat_samplings_fst(table_cds_codons, subset_sizes, samplings)

result <- merged_cor_subsamplings(list(nt = cor_test_nt, codon = cor_test_codons))

fwrite(result, paste0("~/livemount/shared_results/correlation_all_genes/",gene,".csv"))
}
