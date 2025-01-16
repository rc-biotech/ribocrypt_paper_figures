library(RiboCrypt); library(data.table)
res_folder <- "~/livemount/shared_results/ribocrypt_paper/rds_temp_objects/"
all_TIS_file <- file.path(res_folder, "predicted_TIS_all.csv")
df_merged_human <- read.experiment("all_merged-Homo_sapiens", validate = FALSE)
symbols <- symbols(df_merged_human)

if (!file.exists(all_TIS_file)) {
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
  TIS_cov <- cov[is_a_TIS == T, .(count, id)]
  colnames(TIS_cov)[1] <- "TIS_count"
  cov <- data.table::merge.data.table(cov, TIS_cov, by = "id", all.x = TRUE, sort = FALSE)
  cov[, in_frame := (rel_pos %% 3) == 0]
  data.table::fwrite(cov, all_TIS_file)
} else cov <- data.table::fread(all_TIS_file)






subset_cand <- cov[max > 100  & is_a_TIS & !is_a_uORF,]
subset <- subset_cand[count == 0 & is_max == FALSE & in_frame == TRUE,]
print(paste0(round((nrow(subset) / nrow(subset_cand)*100), 2), "% (", nrow(subset), "/", nrow(subset_cand), ")"))

predicted_TIS <- cov[!is_a_uORF & is_max & max > 100  & !(rel_pos != 0 & is_a_TIS),]
predicted_TIS[,type := ifelse(rel_pos > 0, "NTE", ifelse(rel_pos == 0, "annotated", "NTT"))]
predicted_TIS <- predicted_TIS[type == "annotated" | (rel_pos != 0 & in_frame & (TIS_count == 0)),]
plot_TIS_type <- ggplot(data = predicted_TIS, aes(x = type)) +
  geom_bar() + scale_y_continuous(trans = "log10") +  theme_minimal() +
  ylab("Count (Log 10)") + xlab("Translon Type") +theme(legend.position = "bottom"); plot_TIS_type
plot_TIS_distance <- ggplot(data = predicted_TIS[type != "annotated",], aes(x = rel_pos)) +
  geom_bar(width=1.5) + theme_minimal() +
  ylab("Count") + xlab("Distance from annotated TIS") + theme(legend.position = "bottom"); plot_TIS_distance

plot_TIS <- cowplot::plot_grid(plot_TIS_type, plot_TIS_distance, nrow = 1); plot_TIS
ggsave(file.path(plot_folder, "TIS_distance.jpg"), plot_TIS, width = 8, height = 8)
ggsave(file.path(plot_folder, "TIS_distance.svg"), plot_TIS, width = 8, height = 8)



## Check
# cov[is_a_TIS & !is_a_uORF & is_max,]
# cov[is_a_TIS & !is_a_uORF,]
# cov[is_a_TIS & !is_a_uORF & !is_max,]
# index <- 170
# indices <- seq(250, 350, by = 10)
indices <- seq(4, 100, by = 10)
for (index in indices) {
  selected_id <- symbols[ensembl_tx_name %in% subset$id,][index,]; selected_id
  browseRC(selected_id$external_gene_name, selected_id$ensembl_gene_id, selected_id$ensembl_tx_name,
           libraries = c("RNA", "RFP"))
  Sys.sleep(2)
}
