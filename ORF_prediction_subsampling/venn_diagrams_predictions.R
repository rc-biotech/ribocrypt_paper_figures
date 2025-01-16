library(VennDiagram)
library(RColorBrewer)
library(data.table)
library(massiveNGSpipe)

plot_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results"
organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
ORF_subsets <- c("uORFs", "cds")
org <- organisms[1]
ORF_subset <- ORF_subsets[1]
orf_length_filter <- 18

all_files <- c()
dt_all_sets <- data.table()
dt_all_sets_summary <- data.table()
dt_all_sets_mean <- data.table()
for (org in organisms) {
  org_short <- gsub("all_samples-", "", org)
  message(org_short)
  main_result_dir <- file.path(result_dir, org_short)

  for (ORF_subset in ORF_subsets) {
    message(ORF_subset)

    grl <- qs::qread(paste(file.path(main_result_dir, org_short), ORF_subset, "candidate_grl.qs", sep = "_"))
    ids <- which(widthPerGroup(grl) >= orf_length_filter)

    fst <- setDT(fst::read_fst(file.path(main_result_dir, paste(org_short, ORF_subset, "pred_by_classifiers.fst", sep = "_"))))
    fst[, id := seq.int(.N), by = .(classifier, sampling_size, sampling_repeat)]
    stopifnot(max(fst$id) == length(grl))
    fst <- fst[id %in% ids]
    fst_template <- fst[classifier == "RiboCode"]

    classifiers <- unique(fst$classifier)
    tps <- setDT(lapply(classifiers, function(x) fst[classifier == x,]$TP))
    colnames(tps) <- classifiers
    total_classifier <- colSums(tps)

    dt_tps <- copy(tps)
    colnames(dt_tps) <- LETTERS[1:3]
    dt_tps[, group := fifelse(A & B & C, "intersection_all",
                              fifelse(A & B & !C, "A & B",
                                      fifelse(A & C & !B, "A & C",
                                              fifelse(B & C & !A, "B & C",
                                                      fifelse(A & !B & !C, "unique A",
                                                              fifelse(B & !A & !C, "unique B",
                                                                      fifelse(C & !A & !B, "unique C", "None")))))))]
    dt_tps[, group := sub("C", classifiers[3], group)]
    dt_tps[, group := sub("A", classifiers[1], group)]
    dt_tps[, group := sub("B", classifiers[2], group)]
    dt_tps_with_meta <- cbind(dt_tps, org = org_short, ORF_subset = ORF_subset, fst_template[, .(sampling_size, sampling_repeat, id)])

    dt_tps_with_meta[, sampling_size := factor(sampling_size, levels = as.character(sort(as.numeric(unique(sampling_size)))))]
    dt_tps_with_meta[, sampling_repeat := factor(sampling_repeat, levels = as.character(sort(as.numeric(unique(sampling_repeat)))))]

    dt_tps_summary <- as.data.table(table(dt_tps_with_meta[, .(group, sampling_size, sampling_repeat)]))
    max_sampling_1_repeat <- dt_tps_summary[sampling_size == "2783" & sampling_repeat == "1",]$N
    dt_tps_summary[sampling_size == "2783" & sampling_repeat != "1",
                   N := rep(max_sampling_1_repeat, length.out = nrow(dt_tps_summary[sampling_size == "2783" & sampling_repeat != "1"]))]

    dt_tps_summary_meta <- cbind(dt_tps_summary, org = org_short, ORF_subset = ORF_subset)
    dt_tps_mean <- dt_tps_summary_meta[, .(mean = mean(N)), by = .(sampling_size, group)]
    dt_tps_mean_meta <- cbind(dt_tps_mean, org = org_short, ORF_subset = ORF_subset)

    dt_all_sets <- rbindlist(list(dt_all_sets, dt_tps_with_meta))
    dt_all_sets_summary <- rbindlist(list(dt_all_sets_summary, dt_tps_summary_meta))
    dt_all_sets_mean <- rbindlist(list(dt_all_sets_mean, dt_tps_mean_meta))
  }
}


# stopifnot(all(file.exists(all_files)))
# sapply(all_files, function(x) browseURL(x))

# All Truth tables (large)
fst::write_fst(dt_all_sets, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction.fst"))
# Summarized all 20 repeat samplings
fwrite(dt_all_sets_summary, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction_summary.csv"))
fst::write_fst(dt_all_sets_summary, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction_summary.fst"))
# Average model
fwrite(dt_all_sets_mean, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction_summary_mean.csv"))
fst::write_fst(dt_all_sets_mean, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction_summary_mean.fst"))

discordr::send_webhook_file(file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction_summary_mean.fst"))
discordr::send_webhook_message(paste("Here are the overlap mean stats for Eivind above, as fst file, with venn diagram categories and average counts, containing columns:",
                                     paste(colnames(dt_tps_mean_meta), collapse = " ")))
message("Done")
