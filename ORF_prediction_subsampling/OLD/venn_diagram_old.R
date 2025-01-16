library(VennDiagram)


# dt_venn_main <- dt_final_copy[sampling_size == "2783" & sampling_repeat %in% c("1"),]
# overlap_per <- dt_venn_main[, c("classifier", "sampling_size", "sampling_repeat", "TP")]
#
# overlap_all <- dt_venn_main$overlap_all_models[1]
# dt_venn <- data.table::melt(dt_venn_main[,.(classifier, RiboCode_ORFquant_overlap, ORFik_RiboCode_overlap,
#                                        ORFik_ORFquant_overlap, sampling_size, sampling_repeat, species, ORF_type)],
#                                      id.vars = c("classifier", "sampling_size", "sampling_repeat", "species", "ORF_type"))
# models <- strsplit(gsub("_ratio", "", dt_venn$variable), "_")
# models <- data.table(model_one = sapply(models, function(x) x[1]), model_two = lapply(models, function(x) x[2]))
# dt_venn <- cbind(dt_venn, models)
# # dt_venn <- dt_venn[classifier == model_one, ]
# dt_venn <- dt_venn[classifier == model_one | classifier == model_two, ]
# dt_venn[, compare_to := paste0(classifier, "_", variable)]
# dt_venn[, value := value - overlap_all]
#
# merge <- cbind(merge.data.table(dt_venn, overlap_per, by = "classifier"), overlap_all)
# merge[, TP := TP - value - overlap_all]
# classifiers <- unique(merge$classifier)
# list <- lapply(classifiers, function(x) paste(x, seq(dt_venn_main[classifier == x,]$TP)))

org <- organisms[1]
ORF_subset <- ORF_subsets[1]

all_files <- c()
dt_all_sets <- data.table()
dt_all_sets_id <- data.table()
for (org in organisms) {
  org_short <- gsub("all_samples-", "", org)
  message(org_short)
  dfc <- read.experiment(org, validate = F)
  main_result_dir <- file.path("~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results", org_short)

  for (ORF_subset in ORF_subsets) {
    message(ORF_subset)

    fst_all <- setDT(fst::read_fst(file.path(main_result_dir, paste(org_short, ORF_subset, "pred_by_classifiers.fst", sep = "_"))))
    sampling_sizes <- as.character(c(1, 3, 10, 500, nrow(dfc)))
    sampling_size <- sampling_sizes[1]
    for (sampling_size in sampling_sizes) {
      message("- ", sampling_size)
      s <- sampling_size
      fst <- fst_all[sampling_size == s & sampling_repeat %in% c("1"),]
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
      dt_tps[, id := .I]
      dt_tps <- dt_tps[, .(id, group, org = org_short, ORF_subset = ORF_subset, sampling_size = sampling_size)]
      dt_all_sets_id <- rbindlist(list(dt_all_sets_id, dt_tps))

      # library(RColorBrewer)
      # myCol <- brewer.pal(length(classifiers), "Pastel2")
      # overlap <- sapply(list(c(1,2), c(1,3), c(2,3)), function(x) {
      #   id1 <- x[1]
      #   id2 <- x[2]
      #   id3 <- which(!(seq(3) %in% x))
      #
      #   res <- colSums(tps[, ..id1] & tps[, ..id2] & !tps[, ..id3])
      #   names(res) <- paste0(classifiers[id1], "_", classifiers[id2])
      #   res
      # })
      # overlap_per <- sapply(classifiers, function(x) {
      #   hits <- grep(x, names(overlap))
      #   res <- sum(overlap[hits])
      # })
      # overlap_all <- c(all = sum(rowSums(tps) == 3))
      # unique <- total_classifier - overlap_per - overlap_all
      # list <- lapply(seq(length(unique)), function(i) {
      #   x <- unique[i]
      #   paste("unique_", names(x), seq(x))}
      # )
      # lengths(list)
      # list <- lapply(seq(length(unique)), function(i) {
      #   x <- classifiers[i]
      #   hits <- grep(x, names(overlap))
      #   res <- c(list[[i]], paste("overlap_", names(overlap)[hits[1]], seq((overlap)[hits[1]])))
      #   res <- c(res, paste("overlap_", names(overlap)[hits[2]], seq((overlap)[hits[2]])))}
      # )
      # lengths(list)
      # list <- lapply(seq(length(unique)), function(i) {
      #   x <- overlap_all
      #   c(list[[i]], paste("overlap_all", seq(x)))}
      # )
      #
      # total_classifier
      # lengths(list)
      # stopifnot(lengths(list) == total_classifier)
      # dt <- data.table(size_of_set = c(unique, overlap, overlap_per, total_classifier, overlap_all), classifier_id = names(c(unique, overlap, overlap_per, total_classifier, overlap_all)),
      #                  set_type_id = c(rep("unique_set", 3), rep("intersection_2_per", 3), rep("intersection_2_total", 3), rep("total_classifier_set", 3), rep("intersection_3_all")),
      #                  org = org_short, ORF_subset = ORF_subset, sampling_size = sampling_size)
      # dt_all_sets <- rbindlist(list(dt_all_sets, dt))
      # Chart
      # if (interactive()) try(dev.off(dev.list()["RStudioGD"]))
      # invisible(capture.output(venn <- venn.diagram(
      #   x = list,
      #   category.names = classifiers,
      #   filename = NULL,
      #   output=TRUE,
      #   disable.logging = TRUE,
      #   # Circles
      #   lwd = total_classifier/100,
      #   lty = 'blank',
      #   fill = myCol,
      #   # Numbers
      #   cex = .6,
      #   fontface = "bold",
      #   fontfamily = "sans",
      #   # Set names
      #   main = paste("Overlap ORF prediction", org_short, ORF_subset),
      #   sub = paste("Libraries merged:", sampling_size),
      #   sub.cex = 0.5,
      #   cat.cex = 0.6,
      #   cat.fontface = "bold",
      #   cat.default.pos = "outer",
      #   cat.pos = c(-27, 27, 135),
      #   cat.dist = c(0.055, 0.055, 0.085),
      #   cat.fontfamily = "sans",
      #   rotation = 1
      # )))
      # grid.draw(venn)
      #
      # file_prefix <- file.path(plot_dir, paste("Overlap_ORF_prediction", org_short, ORF_subset, "nlibs", sampling_size, sep = "_"))
      #
      # plot_all_versions(venn, file_prefix,
      #                   width = 4, height = 4, send_to_google_drive = F, send_to_discord = F)
      #
      # # discordr::send_webhook_file(file)
      # file <- paste0(file_prefix, ".jpg")
      # all_files <- c(all_files, file)
    }
  }
}

stopifnot(all(file.exists(all_files)))
sapply(all_files, function(x) browseURL(x))
fwrite(dt_all_sets, file.path(plot_dir, "venn_sets_Overlap_ORF_prediction.csv"))
discordr::send_webhook_file(file.path(plot_dir, "venn_sets_Overlap_ORF_prediction.csv"))
fwrite(dt_all_sets_id, file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction.csv"))
discordr::send_webhook_file(file.path(plot_dir, "venn_sets_by_id_Overlap_ORF_prediction.csv"))
