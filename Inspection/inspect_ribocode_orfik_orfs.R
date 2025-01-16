## Inspect data
library(RiboCrypt)
library(data.table)
# df <- read.experiment("all_merged-Homo_sapiens_2024_8") # On local computer
df <- read.experiment("all_merged-Homo_sapiens") # On server
mrna <- loadRegion(df, "mrna")
cds <- loadRegion(df, "cds")
RFP_covrle <- fimport(filepath(df, "cov"))
#result_dir <- "~/Desktop/benchmark_ORFik_vs_RiboCode/" # On local computer
result_dir <- "~/livemount/shared_data" # On server

# uORFomePipe
uorfs <- readRDS(file.path(result_dir, "uORF_prediction_union_of_ORFik_RiboCode.rds"))
fp_uorfs <- readRDS(file.path(result_dir, "uORF_prediction_RiboCode_not_ORFik.rds"))
fn_uorfs <- readRDS(file.path(result_dir, "uORF_prediction_ORFik_not_RiboCode.rds"))

orf_txs <- unique(txNames(uorfs))
uorfs_true <- uorfs
names(uorfs_true) <- paste0("tuORF", seq_along(uorfs_true))
names(fp_uorfs) <- paste0("fpORF", seq_along(fp_uorfs))
names(fn_uorfs) <- paste0("fnORF", seq_along(fn_uorfs))
uorfs <- GRangesList() # This can also be the full set, so keep for now


for (i in orf_txs) {
  print(multiOmicsPlot_ORFikExp(mrna[i], cds[i],
                                custom_regions = c(uorfs_true, fp_uorfs, fn_uorfs),
                                reads = list(RFP_WT = RFP_covrle),
                                df = df[1,], viewMode = "tx",
                                frames_type = "columns", display_sequence = TRUE))
  readline(prompt="Press 'Enter' for next:")
}
