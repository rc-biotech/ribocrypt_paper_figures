library(ORFik)
library(data.table)
library(ggplot2)
library(massiveNGSpipe)
library(cowplot)

df <- read.experiment("all_merged-Homo_sapiens")
symbols <- symbols(df)
symbols[, uniprot_id := NULL]
RFP <- fimport(filepath(df, "cov"))
bw_paths <- filepath(df, "bigwig")
mane <- data.table::fread("~/livemount/shared_data/canonical_transcripts.txt")[[1]]
cds <- loadRegion(df, "cds")
mrna <- loadRegion(df, "mrna")
leaders <- loadRegion(df, "leaders")
leaders <- leaders[lengths(leaders) == 1]
introns <- loadRegion(df, "introns")
introns <- introns[lengths(introns) > 0]

coding_tx <- names(cds)
name_subset <- coding_tx[coding_tx %in% mane & coding_tx %in% names(introns) & coding_tx %in% names(leaders)]
cds <- cds[names(cds) %in% name_subset]
mrna <- mrna[names(mrna) %in% name_subset]
introns <- introns[names(introns) %in% name_subset]

stopifnot(length(cds) == length(mrna)); length(cds)
stopifnot(length(cds) == length(introns)); length(introns)

introns_copy <- introns
strand_bool <- strandBool(introns)
introns_copy[strand_bool] <- heads(introns[strand_bool], 1)
introns_copy[!strand_bool] <- tails(introns[!strand_bool], 1)
stopifnot(all(lengths(introns_copy) == 1))

intron_starts <- startSites(introns_copy, TRUE, TRUE, TRUE)
cds_starts <- startSites(cds, is.sorted = TRUE)

diff <- rep(as.integer(NA), length(intron_starts))
diff[strand_bool] <-  start(intron_starts[strand_bool]) - cds_starts[strand_bool]
diff[!strand_bool] <- cds_starts[!strand_bool] - end(intron_starts[!strand_bool])
stopifnot(!anyNA(diff))

summary(diff)
summary(diff[strand_bool])
# Update
intron_starts <- intron_starts[diff > 0]
introns_copy <- introns_copy[diff > 0]
cds_starts <- cds_starts[diff > 0]
cds <- cds[diff > 0]
mrna <- mrna[diff > 0]
diff <- diff[diff > 0]
stopifnot(all(diff > 0))

phase <- abs(diff %% -3)

windows_all <- windowPerGroup(intron_starts, tx = introns_copy, upstream = -(phase + 3), downstream = 300)
summary(widthPerGroup(windows_all))


windows <- windows_all[widthPerGroup(windows_all) > 30]
summary(widthPerGroup(windows))

seqs <- txSeqsFromFa(windows, df, TRUE, TRUE)
subseq(seqs, 1, 3) <- "ATG"
orfs <- findORFs(seqs, startCodon = "ATG")
orfs <- orfs[start(orfs) == 1]
length(orfs)

orfs_gr <- ORFik:::mapToGRanges(windows, orfs, groupByTx = FALSE, grl_is_sorted = TRUE)
length(orfs_gr)
orfScores <- orfScore(orfs_gr, RFP, is.sorted = TRUE)
filter <- orfScores$ORFScores > 5 & orfScores$frame_zero_RP > 300
candidates <- orfScores[filter,]
candidates

candidates_gr <- orfs_gr[filter]
candidates_gr
cov <- coveragePerTiling(candidates_gr, RFP, as.data.table = TRUE)
cov[, frame := position %% 3]
cov <- cov[, sum(count > 0), by = .(genes, frame)]
tx_ids <- data.table(CTE_id = names(candidates_gr), tx_ids = txNames(candidates_gr), CTE_length = widthPerGroup(candidates_gr))

merge <- data.table::merge.data.table(tx_ids, symbols, by.x = "tx_ids", by.y = "ensembl_tx_name", all.x = TRUE, all.y = FALSE, sort = FALSE)

final <- cbind(merge, candidates)
final[, ORF_sum := rowSums(final[, .(frame_zero_RP, frame_one_RP, frame_two_RP)])]
final[, ORF_ratio := ORF_sum / CTE_length]
final[, ORF_F1_ratio :=  frame_zero_RP/ CTE_length]
final[, cds_sum := countOverlaps(cds[tx_ids], RFP)]
final[, cds_length := widthPerGroup(cds[tx_ids])]
final[, cds_ratio := cds_sum / cds_length]
final[, ORF_cds_ratio := (ORF_ratio / cds_ratio)*100]
final[, frame_zero_usage := 100*(cov[frame == 0]$V1 / (CTE_length/3))]
final[, frame_one_usage := 100*(cov[frame == 1]$V1 / (CTE_length/3))]
final[, frame_two_usage := 100*(cov[frame == 2]$V1 / (CTE_length/3))]
nrow(final)
final_candidates <- final
nrow(final_candidates)
final_candidates_filtered <- final_candidates[ORF_sum > 300 & frame_zero_usage > 50,]
final_gr <- candidates_gr[final$ORF_sum > 300 & final$frame_zero_usage > 50]
stopifnot(length(final_gr) == nrow(final_candidates_filtered))
saveRDS(final_gr, "~/livemount/shared_results/ribocrypt_paper/intron_candidates.rds")
export.bed12(final_gr, "~/livemount/shared_results/ribocrypt_paper/intron_candidates.bed12")
fwrite(final_candidates_filtered, file = "~/livemount/shared_results/ribocrypt_paper/intron_candidates.csv")
discordr::send_webhook_file("~/livemount/shared_results/ribocrypt_paper/intron_candidates.csv")


final_candidates_filtered <- fread("~/livemount/shared_results/ribocrypt_paper/intron_candidates.csv")
candidates <- as.integer(nrow(final_candidates_filtered))
melt <- melt(final_candidates_filtered[, c("CTE_length", "ORF_cds_ratio")])
melt[variable == "ORF_cds_ratio", variable := "Intron expression level"]
melt[variable == "CTE_length", variable := "Introic extension length"]

empty_x_theme <- theme(axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())
plot_length <- ggplot(melt[variable == "Introic extension length",], aes(y = value)) +
  geom_boxplot(outliers = FALSE) + theme_classic() + ylab("Introic extension length") + empty_x_theme
plot_length
plot_expression_level <- ggplot(melt[variable == "Intron expression level",], aes(y = value)) +
  geom_boxplot(outliers = FALSE) + theme_classic() + ylab("Intron/CDS expression level (%)") +
  ylim(-1, 100) + empty_x_theme
plot_expression_level
plot_counts <- ggplot(data.frame(candidates = candidates), aes(y = candidates)) +
  geom_boxplot() + theme_classic() + ylab("Predicted expressed introns") + xlab(NULL)
+ scale_y_discrete(limits = candidates) + empty_x_theme
plot_counts
title <- cowplot::ggdraw() + draw_label(paste0(type, " analysis (", org, ")"), fontface='bold')

plot <- cowplot::plot_grid(NULL, title, NULL, plot_counts, plot_length, plot_expression_level,
                           nrow = 2, rel_heights = c(0.1, 0.9)); plot

massiveNGSpipe::plot_all_versions(plot, "~/livemount/shared_results/ribocrypt_paper/intron_plot",
                                  send_to_google_drive = TRUE, send_to_discord = TRUE)
discordr::send_webhook_message(paste("Files are also on google drive",
                                     massiveNGSpipe:::google_drive_dir_links()))
