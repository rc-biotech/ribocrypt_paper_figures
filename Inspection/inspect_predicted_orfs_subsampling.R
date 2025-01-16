## Inspect data
library(RiboCrypt)
library(data.table)

organisms <- c("all_merged-Homo_sapiens", "all_merged-Saccharomyces_cerevisiae")
orf_type <- c("uORFs", "cds")[1]
df <- read.experiment(organisms[1], validate = FALSE) # On server
mrna <- loadRegion(df, "mrna")
leaders <- loadRegion(df, "leaders")
cds <- loadRegion(df, "cds")
leader_cds <- ORFik:::addCdsOnLeaderEnds(leaders, cds)
# RFP_covrle <- fimport(filepath(df, "cov"))
org <- gsub(" ", "_", organism(df))

RFP_covrle <- qs::qread(file.path("~/livemount/shared_results/ribocrypt_paper/rds_temp_objects/covRLE_whole_library_subsamplings", org, "covRLE_samp_2783_num_1.qs"))
seqinfo(RFP_covrle@forward) <- seqinfo(df)
seqinfo(RFP_covrle@reverse) <- seqinfo(df)
#result_dir <- "~/Desktop/benchmark_ORFik_vs_RiboCode/" # On local computer
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results" # On server

# uORFomePipe
orfs <- qs::qread(file.path(result_dir, org, paste(org, orf_type, "candidate_grl.qs", sep = "_")))
prediction_table <- setDT(fst::read_fst(file.path(result_dir, org, paste(org, orf_type, "pred_by_classifiers.fst", sep = "_"))))
stopifnot(nrow(prediction_table) %% length(orfs) == 0)
names(orfs) <- prediction_table[sampling_size == 1 & sampling_repeat == 1 & classifier == "RiboCode",]$tx_id
orf_txs <- unique(txNames(orfs))

truth_table_RiboC <- prediction_table[classifier == "RiboCode",]
truth_table_ORFquant <- prediction_table[classifier == "ORFquant",]
truth_table_ORFik <- prediction_table[classifier == "ORFik",]
stopifnot(all((names(orfs) == as.character(truth_table_ORFik$tx_id))))

prediction_per <- data.table(RiboCode = truth_table_RiboC$TP, ORFquant = truth_table_ORFquant$TP, ORFik = truth_table_ORFik$TP,
                             sampling_size = truth_table_ORFik$sampling_size, sampling_repeat = truth_table_ORFik$sampling_repeat,
                             tx_id = truth_table_ORFik$tx_id)
prediction_subset <- prediction_per[as.integer(as.character(sampling_size)) == 2783,]
stopifnot(nrow(prediction_subset) == length(orfs))
stopifnot(all((names(orfs) == as.character(prediction_subset$tx_id))))

prediction_subset <- prediction_subset[tx_id %in% names(mrna),]
ribocode_orfs <- orfs[prediction_subset$RiboCode]
ORFquant_orfs <- orfs[prediction_subset$ORFquant]
ORFik_orfs <- orfs[prediction_subset$ORFik]
names(ribocode_orfs) <- paste0("RC_ORF", seq_along(ribocode_orfs))
names(ORFquant_orfs) <- paste0("OC_ORF", seq_along(ORFquant_orfs))
names(ORFik_orfs) <- paste0("O_ORF", seq_along(ORFik_orfs))

# names(uorfs_true) <- paste0("tuORF", seq_along(uorfs_true))
# names(fp_uorfs) <- paste0("fpORF", seq_along(fp_uorfs))
# names(fn_uorfs) <- paste0("fnORF", seq_along(fn_uorfs))
# uorfs <- GRangesList() # This can also be the full set, so keep for now
# TODO: checklist (ENST00000371084)
# to_check <- unique(prediction_subset[(rowSums(prediction_subset[, 1:3]) == 1) & prediction_subset$ORFquant == TRUE,]$tx_id)
to_check <- unique(prediction_subset[rowSums(prediction_subset[, 1:3]) == 2,]$tx_id)
length(to_check)
for (i in to_check[-(1:200)]) {
  # if (countOverlaps(cds[i], c(ribocode_orfs, ORFquant_orfs, ORFik_orfs), type = "equal") != 2) next
  print(multiOmicsPlot_ORFikExp(leader_cds[i], cds[i],
                                custom_regions = c(ribocode_orfs, ORFquant_orfs, ORFik_orfs),
                                reads = list(RFP_WT = RFP_covrle), trailer_extension = 30, leader_extension = 30,
                                df = df[1,], viewMode = "tx", custom_motif = c("HTG"),
                                frames_type = "columns", display_sequence = TRUE))
  readline(prompt="Press 'Enter' for next:")
}
