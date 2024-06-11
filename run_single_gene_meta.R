library(ORFik)
library(data.table)
df_all <- read.experiment("all_samples-Homo_sapiens", validate = FALSE)
m <- fread("~/livemount/Bio_data/NGS_pipeline/metadata_rc.csv")

collection_genes <- RiboCrypt:::get_gene_name_categories_collection(df_all)
collection_genes[, gene_symbol := gsub("-.*", "", label)]

load_allsamples_fst <- function(gene, df_all, collection_genes) {
  collection_gene <- collection_genes[gene_symbol %in% gene,]
  collection_gene <- collection_gene[!duplicated(gene_symbol),]
  id <- collection_gene$value
  gene <- collection_gene$gene_symbol
  names(id) <- gene
  message("ID: ", id, if (!is.null(names(id))) " (", names(id), ")")
  collection_folder <- file.path(resFolder(df_all), "collection_tables/")
  table_path <- paste0(collection_folder, id, ".fst")
  
  if (!file.exists(dirname(table_path)))
    stop("There is no collection fst tables directory for this organism,",
         " see vignette for more information on how to make these.")
  if (!file.exists(table_path)) stop("Gene has no precomputed table, try another one!")
  
  return(table_path)
}

gene <- c("TEX101")
gene <- c("OAZ1")
gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)

gene <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                             normalization = "tpm",
                                             kmer = 9, metadata = m, min_count = 300, as_list = TRUE)
p <- RiboCrypt:::get_meta_browser_plot(gene$table, clusters = 5,
                                       color_theme = "Matrix (black,green,red)")
p


library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(shiny)
ht = draw(p)

ui = fluidPage(
  InteractiveComplexHeatmapOutput()
)

server = function(input, output, session) {
  makeInteractiveComplexHeatmap(input, output, session, ht)
}

shiny::shinyApp(ui, server)
