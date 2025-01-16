source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/ORF_statistics_functions.R")

ref_dir <- ORFik::config()["ref"]
all_exp <- list.experiments(validate = FALSE, pattern = "all_merged-", libtypeExclusive = "RFP")
all_exp <- all_exp[libtypes == "RFP"]
organisms <- all_exp$name
organisms <- organisms[grep("_modalities$|HEK293|_04|Homo_sapiens", organisms, invert = TRUE)]
org <- "all_merged-Homo_sapiens_04_oct_2024_all"
sizes <- c(50, 100)
for (org in organisms) {
  for (size in sizes) {
    df <- read.experiment(org, validate = FALSE)
    org_short <- gsub(" ", "_", tolower(organism(df)))
    message("- ", org_short, " ", size)
    org_dir <- file.path(ref_dir, org_short, "predicted_translons", "ORFik")

    org_dir_out <- file.path(org_dir, "intergenic")
    if (file.exists(file.path(org_dir_out, paste0("intergenic_", size, ".bed12")))) next
    if (!dir.exists(org_dir_out)) dir.create(org_dir_out, recursive = TRUE)


    txdb <- loadTxdb(df)
    tx <- unlistGrl(extendTrailers(extendLeaders(flankPerGroup(loadRegion(txdb, "tx")), 3000), 3000))
    strand(tx) <- "*"
    genome <- GRanges(seqinfo(tx))
    intergenic <- setdiff(genome, tx)

    tiles <- unlistGrl(tile(intergenic, width = size))
    tiles_sample <- tiles[sample(length(tiles), 5000)]
    tiles_grl <- split(tiles_sample, seq_along(tiles_sample))

    export.bed12(tiles_grl, file = file.path(org_dir_out, paste0("intergenic_", size, ".bed12")))
  }
}

