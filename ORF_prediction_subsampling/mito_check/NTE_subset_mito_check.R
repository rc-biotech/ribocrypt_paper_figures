library(ORFik)

exp <- "all_merged-Homo_sapiens"
df <- read.experiment(exp)
table <- readRDS("~/livemount/Bio_data/references/homo_sapiens/predicted_translons/ORFik/predicted_translons_all_merged-Homo_sapiens_RFP_prediction_table.rds")
predicted <- table$predicted == TRUE
NTEs_hits <- table$type == "NTE"

grl <- readRDS("~/livemount/Bio_data/references/homo_sapiens/predicted_translons/ORFik/predicted_translons_all_merged-Homo_sapiens_RFP_candidates.rds")
length(grl)
grl <- grl[predicted & NTEs_hits]
length(grl)

txdb <-loadTxdb(df)
cds <- loadRegion(txdb, "cds")
mrna <- loadRegion(txdb, "mrna")

start <- countOverlaps(grl, cds, type = "start")
grl <- grl[start == 0]
cds_match <- cds[names(grl)]
stopifnot(length(grl) == length(cds_match))

ir_grl <- pmapToTranscriptF(grl, mrna[names(grl)])
ir_cds <- pmapToTranscriptF(cds_match, mrna[names(grl)])
diff <- start(unlist(ir_cds)) - start(unlist(ir_grl))
summary(diff)

grl <- grl[diff > 15]
length(grl)
diff <- diff[diff > 15]
summary(diff)
seqs <- txSeqsFromFa(grl, df, is.sorted = TRUE)
aa <- translate(seqs)
stopifnot(length(grl) == length(seqs))
symbols <- symbols(df)
symbols <- symbols[ensembl_tx_name %in% names(grl),]

res <- data.table::data.table(tx_id = names(grl), extension = diff, sequence_DNA = as.character(seqs), sequence_AA = as.character(aa))
dt <- data.table::merge.data.table(res, symbols, by.x = "tx_id", by.y = "ensembl_tx_name", all.x = TRUE, all.y = FALSE, sort = FALSE)
data.table::setcolorder(dt, c(1,2,5,6,7,3,4))

output_dir <- file.path("~/livemount/shared_results/ribocrypt_paper/mito_check/", exp)
fst::write_fst(dt, file.path(output_dir, "NTE_for_mito_check.fst"))
saveRDS(grl, file.path(output_dir, "NTE_for_mito_check_ranges.rds"))
