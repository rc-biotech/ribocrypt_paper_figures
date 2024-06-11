library(RiboCrypt); library(data.table)
res_folder <- "~/livemount/shared_results/ribocrypt_paper/rds_temp_objects/"
df_merged_human <- read.experiment("all_merged-Homo_sapiens")
symbols <- symbols(df_merged_human)
RFP <- fimport(filepath(df_merged_human, "cov"))
cds <- loadRegion(df_merged_human, part = "cds")
uorfs <- loadRegion(df_merged_human, part = "uorfs")
TIS_uorfs <- split(startSites(uorfs, TRUE, TRUE, TRUE), seq_along(uorfs))
names(TIS_uorfs) <- names(uorfs)
TIS <- split(startSites(cds, TRUE, TRUE, TRUE), names(cds))

flank <- extendLeaders(TIS, 100)
flank <- extendTrailers(flank, 100)

seqs <- txSeqsFromFa(flank, df_merged_human@fafile)

all_starts <- Biostrings::vmatchPattern("ATG", seqs)

all_starts_map <- IRangesList(all_starts)
names(all_starts_map) <- seq_along(all_starts_map)
all_starts_genomic <- pmapFromTranscriptF(all_starts_map, flank, removeEmpty = TRUE)



cov <- coveragePerTiling(all_starts_genomic, RFP, as.data.table = TRUE)[position == 1,]
cov[, id := names(all_starts_genomic)]
cov[, rel_pos := (start(unlist(all_starts)) - 101)]
cov[, is_a_TIS := all_starts_genomic %over% TIS]
cov[, is_a_uORF := all_starts_genomic %over% TIS_uorfs]
cov[, max := max(count), by = id]
cov[, is_max := count == max]



cov[is_a_TIS & !is_a_uORF & is_max,]
cov[is_a_TIS & !is_a_uORF,]
cov[is_a_TIS & !is_a_uORF & !is_max,]

subset_cand <- cov[max > 100  & is_a_TIS & !is_a_uORF,]
subset <- subset_cand[count == 0 & is_max == FALSE,]
print(paste0(round((nrow(subset) / nrow(subset_cand)*100), 2), "% (", nrow(subset), "/", nrow(subset_cand), ")"))
data.table::fwrite(cov, file.path(res_folder, "predicted_TIS_all.csv"))
data.table::fwrite(cov, file.path(res_folder, "predicted_wrong_TIS.csv"))

# index <- 170
indices <- seq(250, 350, by = 10)
for (index in indices) {
  selected_id <- symbols[ensembl_tx_name %in% subset$id,][index,]; selected_id
  browseRC(selected_id$external_gene_name, selected_id$ensembl_gene_id, selected_id$ensembl_tx_name,
           libraries = c("RNA", "RFP"))
  Sys.sleep(2)
}



#' Browse a gene on Ribocrypt webpage
#'
#' Can also disply local RiboCrypt app
#' @param symbol gene symbol, default NULL
#' @param gene_id gene symbol, default NULL
#' @param tx_id gene symbol, default NULL
#' @param exp experiment name, default "all_merged-Homo_sapiens_modalities"
#' @param libraries NULL, default to first in experiment, c("RFP","RNA") would add RNA to default.
#' @param host url, default "https://ribocrypt.org". Set to localhost for local version.
#' @param plot_on_start logical, default TRUE. Plot gene when opening browser.
#' @param frames_type "columns"
#' @param kmer default 1 (kmer window to smear out reads)
#' @param browser getOption("browser")
#' @return browseURL, opens browse with page
#' @export
browseRC <- function(symbol = NULL, gene_id = NULL, tx_id = NULL, exp = "all_merged-Homo_sapiens_modalities",
                     libraries=NULL,
                     host = "https://ribocrypt.org", plot_on_start = TRUE,
                     frames_type = "columns", kmer=1,
                     browser = getOption("browser")) {
  if (is.null(symbol) & is.null(gene_id)) stop("At least on of symbol and gene_id must be defined!")
  exp <- paste0("&dff=", exp)
  frames_type <- paste0("&frames_type=", frames_type)
  kmer <- paste0("&kmer=", kmer)
  prefix_url <- paste0(host, "/?", exp, "&frames_type=columns&kmer=1")
  gene <- paste0(if (!is.na(symbol) & !is.null(symbol)) {paste0(symbol, "-")} else NULL, gene_id)
  if (!is.null(libraries)) libraries <- paste0("&library=", paste(libraries, collapse = ","))
  if (!is.null(tx_id)) tx_id <- paste0("&tx=", tx_id)
  if (!is.null(libraries)) libraries
  select <- paste0("&gene=", gene, tx_id, libraries)
  plot_on_start <- paste0("&go=", as.logical(plot_on_start))
  full_url <- paste0(prefix_url, select, plot_on_start)
  # full_url
  browseURL(full_url, browser = browser)
}
