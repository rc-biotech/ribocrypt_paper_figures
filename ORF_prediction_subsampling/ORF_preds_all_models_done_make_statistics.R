library(ORFik)
library(data.table)
library(ggplot2)
library(cowplot)
library(fst)
library(qs)
library(massiveNGSpipe)

map_classifier_to_grl <- function(starts, ends, tx_ids, mrna) {
  if (length(starts) != length(ends) || length(starts) != length(tx_ids)) {
    stop("start, ends and tx_ids must be equal length")
  }
  if (length(starts) == 0) return(GRangesList())
  aa <- IRanges(starts, ends)
  names(aa) <- tx_ids
  names(aa) <- chmatch(names(aa), unique(names(aa)))
  orfs <- pmapFromTranscriptF(aa, mrna[unique(tx_ids)], removeEmpty = TRUE)
  stopifnot(all(names(orfs) == tx_ids))
  return(orfs)
}

truth_table <- function (predicted, true, verbose = TRUE) {
  stopifnot(length(predicted) == length(true))
  predicted <- as.logical(predicted)
  true <- as.logical(true)
  tb <- data.table(TP = predicted & true, FP = predicted &
                     !true, FN = !predicted & true, TN = !predicted & !true)
  if (verbose) {
    message("Table (counts):")
    print(colSums(tb))
    message("Table (%):")
    print(round((colSums(tb)/sum(colSums(tb))) * 100, 1))
  }
  return(tb)
}

message(Sys.time())
samplings_all <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"
covRLE_subsampling_dir_all <- "~/livemount/shared_results/ribocrypt_paper//rds_temp_objects/covRLE_whole_library_subsamplings"
organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
remake <- F
dt_final_all <- data.table()

org <- organisms[1]
ORF_subsets <- c("uORFs", "cds")
ORF_subset <- ORF_subsets[1]
for (org in organisms) {
  message(org)
  message("-------------------------------")
  df <- read.experiment(gsub("all_samples", "all_merged", org), validate = F)
  dfc <- read.experiment(org, validate = F)
  org_short <- gsub("all_samples-", "", org)
  main_result_dir <- file.path("~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results", org_short)
  samplings <- unique(c(samplings_all[samplings_all < nrow(dfc)], nrow(dfc)))
  if (org_short == "Homo_sapiens") {
    gtf <- "~/livemount/shared_data/Homo_sapiens_mane/Homo_sapiens_GRCh38_101_subset300_mane.gtf"
    txdb <- paste0(gtf, ".db")
  } else {
    txdb <- loadTxdb(df)
    gtf <- ORFik:::getGtfPathFromTxdb(txdb)
  }
  mrna <- loadRegion(txdb, "mrna")

  for (ORF_subset in ORF_subsets) {
    message("- ", ORF_subset)
    if (ORF_subset == "uORFs") {
      orf_subset <- c("uoORF", "uORF")
      orf_subset_ribocode <- c("Overlap_uORF", "uORF")
    } else if (ORF_subset == "cds") {
      orf_subset <- orf_subset_ribocode <- c("NTE", "annotated")
    }


    dt_final_path <- file.path(main_result_dir, paste(org_short, ORF_subset, "TruthTable_by_classifiers.fst", sep = "_"))
    if (!file.exists(dt_final_path) | remake) {
      # ORF result dirs
      candidate_path <- file.path(main_result_dir, paste(org_short, ORF_subset, "candidate_grl.qs", sep = "_"))
      ORFik_dir <- file.path(main_result_dir, "ORFik")
      ORFquant_dir <- file.path(main_result_dir, "ORFquant")
      RiboCode_dir <- file.path(main_result_dir, "RiboCode")

      ## Set True uORFs and Candidate ORFs from (all_merge prediction merged from all classifiers)
      # Load ORFik results
      file_suffix_last <- paste("ORFik_results", "samp", nrow(dfc),"num", 1,
                                                    name(df), ORFik:::name_decider(df, naming = "full"), sep = "_")
      all_merged_ORFik <- file.path(ORFik_dir, file_suffix_last)
      orfs_table_all <- readRDS(paste0(all_merged_ORFik, "_prediction_table.rds"))
      orfs_cand_all <- readRDS(paste0(all_merged_ORFik, "_candidates.rds"))

      res_ORFik <- orfs_cand_all[orfs_table_all$predicted & orfs_table_all$type %in% orf_subset]
      uorfs_ORFik_tx_coord <- pmapToTranscriptF(res_ORFik, mrna[names(res_ORFik)])
      # Load ORFquant results
      file_suffix_last <- paste0("ORFquant_results_", "samp_", nrow(df),"_num_", 1, ".fst")
      all_merged_ORFquant <- file.path(ORFquant_dir, file_suffix_last)
      res_ORFquant <- setDT(fst::read_fst(all_merged_ORFquant))
      ORFq_uorf_table <- res_ORFquant[ORF_type %in% orf_subset & significant == TRUE,]
      ## Load RiboCode results
      file_suffix_last <- paste("RiboCode_results", "samp", nrow(dfc),"num", 1, sep = "_")
      file_suffix_last <- file.path(file_suffix_last, "results_collapsed.txt")
      all_merged_RiboCode <- file.path(RiboCode_dir, file_suffix_last)
      res_RiboC <- fread(all_merged_RiboCode, header = TRUE)
      res_RiboC <- res_RiboC[ORF_type %in% orf_subset_ribocode,]

      # Make union set
      candidate_uorfs_generic <- orfs_cand_all[orfs_table_all$type %in% orf_subset] # Start with ORFik matches, append..
      starts <- c(as.numeric(start(uorfs_ORFik_tx_coord)), ORFq_uorf_table$start, res_RiboC$ORF_tstart)
      ends <- c(as.numeric(end(uorfs_ORFik_tx_coord)), ORFq_uorf_table$end, res_RiboC$ORF_tstop)
      tx_ids <- c(names(uorfs_ORFik_tx_coord), as.character(ORFq_uorf_table$tx_id), res_RiboC$transcript_id)
      all_pred_all_classifiers <- paste(starts, ends, tx_ids, sep = "_")
      uniques <- !duplicated(all_pred_all_classifiers)
      orfs_merged <- longestORFs(uniqueGroups(map_classifier_to_grl(starts[uniques], ends[uniques], tx_ids[uniques], mrna)))
      candidate_tx_ids <- c(names(candidate_uorfs_generic), tx_ids[uniques][as.numeric(names(orfs_merged))])
      candidate_uorfs_generic <- longestORFs(uniqueGroups(c(candidate_uorfs_generic, orfs_merged)))
      candidate_names <- as.factor(candidate_tx_ids[as.numeric(names(candidate_uorfs_generic))])
      is_true_uorf <-  countOverlaps(candidate_uorfs_generic, orfs_merged, type = "equal") == 1
      stopifnot(max(countOverlaps(candidate_uorfs_generic, candidate_uorfs_generic, type = "equal")) == 1)
      qs::qsave(candidate_uorfs_generic, file = candidate_path)
      # For each covrle
      covRLE_subsampling_dir <- file.path(covRLE_subsampling_dir_all, org_short)
      format_grep <- "\\.qs$"
      rds_paths <- list.files(covRLE_subsampling_dir, format_grep, full.names = TRUE)
      length(rds_paths)
      rds_path <- rds_paths[1]

      dt_final <- data.table()
      all_truth_tables <- data.table()
      all_truth_tables_path <- file.path(main_result_dir, paste(org_short, ORF_subset, "pred_by_classifiers.fst", sep = "_"))
      # Loop over result for each subsampling and aggregate statistics for all models
      if (!file.exists(all_truth_tables_path) | remake) {
        for (rds_path in rds_paths) {
          # Prediction output paths
          out_file_base <- gsub("covRLE_", "", gsub(format_grep, "", basename(rds_path)))
          message("- ", out_file_base)
          out_files_ORFik <- file.path(ORFik_dir, paste("ORFik_results", out_file_base, name(df),
                                                        ORFik:::name_decider(df, naming = "full"), sep = "_"))
          out_files_RiboCode <-  file.path(RiboCode_dir, paste("RiboCode_results", out_file_base, sep = "_"), "results_collapsed.txt")
          out_files_ORFquant <- file.path(ORFquant_dir, paste("ORFquant_results", paste0(out_file_base, ".fst"), sep = "_"))
          ## Load ORFik results
          uorfs_table <- readRDS(paste0(out_files_ORFik, "_prediction_table.rds"))
          orfs_cand <- readRDS(paste0(out_files_ORFik, "_candidates.rds"))
          orfs_ORFik <- uniqueGroups(orfs_cand[uorfs_table$predicted & uorfs_table$type %in% orf_subset])

          # Load ORFquant results
          res_ORFquant <- setDT(fst::read_fst(out_files_ORFquant))
          ORFq_uorf_table <- res_ORFquant[ORF_type %in% orf_subset & significant == TRUE,]
          orfs_ORFq <- uniqueGroups(map_classifier_to_grl(ORFq_uorf_table$start, ORFq_uorf_table$end, ORFq_uorf_table$tx_id, mrna))

          ## Load RiboCode results
          res_RiboC <- fread(out_files_RiboCode)
          res_RiboC <- res_RiboC[ORF_type %in% orf_subset_ribocode,]
          orfs_RiboC <- uniqueGroups(map_classifier_to_grl(res_RiboC$ORF_tstart, res_RiboC$ORF_tstop, res_RiboC$transcript_id, mrna))

          # Truth tables ORFquant
          orfq_orfs_are_active <- countOverlaps(candidate_uorfs_generic, orfs_ORFq, type = "equal") == 1
          truth_table_ORFquant <- truth_table(orfq_orfs_are_active, is_true_uorf, verbose = FALSE)
          # Truth tables RiboCode
          riboc_orfs_are_active <- countOverlaps(candidate_uorfs_generic, orfs_RiboC, type = "equal") == 1
          truth_table_RiboC <- truth_table(riboc_orfs_are_active, is_true_uorf, verbose = FALSE)
          # Truth tables ORFik
          ORFik_orfs_are_active <- countOverlaps(candidate_uorfs_generic, orfs_ORFik, type = "equal") == 1
          truth_table_ORFik <- truth_table(ORFik_orfs_are_active, is_true_uorf, verbose = FALSE)
          # Combine
          sampling <- strsplit(gsub("_num", "", gsub("samp_","", out_file_base)), "_")[[1]]
          truth_tables <- data.table(tx_id = candidate_names,
                                     codons = widthPerGroup(candidate_uorfs_generic, FALSE)/3,
                                     classifier = rep(c("RiboCode","ORFquant", "ORFik"), each = (length(candidate_uorfs_generic))),
                                     rbindlist(list(truth_table_RiboC, truth_table_ORFquant, truth_table_ORFik)),
                                     sampling_size = sampling[1],
                                     sampling_repeat = sampling[2])

          all_truth_tables <- rbindlist(list(all_truth_tables, truth_tables))
        }
        fst::write_fst(all_truth_tables, all_truth_tables_path)
      } else all_truth_tables <- setDT(fst::read_fst(all_truth_tables_path))

      # all_truth_tables <- all_truth_tables[sampling_size == "1" & sampling_repeat %in% c("1", "10"),]

      truth_table_RiboC <- all_truth_tables[classifier == "RiboCode",]
      truth_table_ORFquant <- all_truth_tables[classifier == "ORFquant",]
      truth_table_ORFik <- all_truth_tables[classifier == "ORFik",]
      truth_table_combined <- cbind(truth_table_RiboC, ORFquant = truth_table_ORFquant$TP, ORFik = truth_table_ORFik$TP)
      truth_table_combined[, RiboCode := TP]; truth_table_combined[, TP := NULL];truth_table_combined[, classifier := NULL]
      truth_table_combined[, `:=`(RiboCode_ORFquant_overlap = RiboCode & ORFquant,
                                  ORFik_RiboCode_overlap = ORFik & RiboCode,
                                  ORFik_ORFquant_overlap = ORFik & ORFquant,
                                  overlap_all_models = ORFik & RiboCode & ORFquant)]
      dt <- melt.data.table(truth_table_combined, measure.vars  = c("ORFik", "ORFquant", "RiboCode"),
                            value.name = "TP", variable.name = "classifier")

      dt <- dt[, .(TP = sum(TP),
                  RiboCode_ORFquant_overlap = sum(RiboCode_ORFquant_overlap),
                  ORFik_RiboCode_overlap = sum(ORFik_RiboCode_overlap),
                  ORFik_ORFquant_overlap = sum(ORFik_ORFquant_overlap),
                  overlap_all_models = sum(overlap_all_models)), by = .(classifier, sampling_size, sampling_repeat)]
      # dt_metrics <- data.table(all_truth_tables,
      #                          RiboCode_ORFquant_overlap = truth_table_RiboC$TP & truth_table_ORFquant$TP,
      #                          ORFik_RiboCode_overlap = truth_table_RiboC$TP & truth_table_ORFik$TP,
      #                          ORFik_ORFquant_overlap = truth_table_ORFquant$TP & truth_table_ORFik$TP,
      #                          overlap_all_models = truth_table_RiboC$TP & truth_table_ORFquant$TP & truth_table_ORFik$TP)


      # Summary statistics
      # By classifier table

      # dt <- rbindlist(lapply(unique(dt_metrics$classifier), function(x)
      #   dt_metrics[classifier == x, .(TP = sum(TP), FP = sum(FP), FN = sum(FN), TN = sum(TN),
      #                                RiboCode_ORFquant_overlap = sum(RiboCode_ORFquant_overlap),
      #                                ORFik_RiboCode_overlap = sum(ORFik_RiboCode_overlap),
      #                                ORFik_ORFquant_overlap = sum(ORFik_ORFquant_overlap),
      #                                overlap_all_models = sum(overlap_all_models)), by = .(classifier, sampling_size, sampling_repeat)]))
      # dt[, `:=` (precision = TP / (TP + FP), recall = TP / (TP + FN), accuracy = (TP + TN) / rowSums(dt[, c("TP", "FP", "FN", "TN"), with = FALSE]),
      #            mcc = coverageSim::mcc(TP, FP, TN, FN))]
      dt[is.na(dt)] <- 0
      # dt[, c("precision", "recall", "accuracy", "mcc")] <- round(dt[, colnames(dt) %in% c("precision", "recall", "accuracy", "mcc"), with = FALSE], 3)
      dt[, `:=`(species = as.factor(org_short), ORF_type = as.factor(ORF_subset))]
      stopifnot(!any(dt$overlap_all_models > dt$TP)) # Sanity test
      fst::write_fst(dt, dt_final_path)

    } else dt <- fst::read_fst(dt_final_path)
    dt_final_all <- rbindlist(list(dt_final_all, dt))
  }
  message("Done (", org_short, ")")
}

dt_final_copy <- dt_final_all

dt_final_copy[, sampling_size := factor(sampling_size, levels = sort(unique(as.numeric(as.character(sampling_size)))), ordered = TRUE)]

translon_plot_orfs <- ggplot(data = dt_final_copy, aes(x = sampling_size, y = TP, fill = classifier)) +
  geom_boxplot(aes(color = classifier)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", size = 0.15,
               position = position_dodge(width = 0.75)) +
  theme_minimal() + ylab(paste("Predicted", unique(dt_final_copy$ORF_type))) +  xlab(paste("Number of Libraries (20 repeat samplings)")) +
  theme(legend.position="top") + guides(color="none") +
  guides(fill=guide_legend(title="Classifier")) +
  facet_wrap(ORF_type ~ species, scales = "free"); plot(translon_plot_orfs)


dt_final_copy[, `:=`(RiboCode_other = RiboCode_ORFquant_overlap + ORFik_RiboCode_overlap - overlap_all_models,
                     ORFik_other = ORFik_ORFquant_overlap + ORFik_RiboCode_overlap - overlap_all_models,
                     ORFquant_other = RiboCode_ORFquant_overlap + ORFik_ORFquant_overlap - overlap_all_models)]
dt_final_copy[, `:=`(RiboCode_ratio = round(100 * RiboCode_other / (TP), 2),
                     ORFik_ratio = round(100 * ORFik_other / (TP), 2),
                     ORFquant_ratio = round(100 * ORFquant_other / (TP), 2))]

# dt_final_overlap <- data.table::melt(dt_final_copy[,.(classifier, ORFik_RiboCode_ratio, ORFik_ORFquant_ratio,
#                                                       RiboCode_ORFquant_ratio, sampling_size, sampling_repeat, species)],
#                                      id.vars = c("classifier", "sampling_size", "sampling_repeat", "species"))
dt_final_overlap <- data.table::melt(dt_final_copy[,.(classifier, RiboCode_ratio, ORFik_ratio,
                                                      ORFquant_ratio, sampling_size, sampling_repeat, species, ORF_type)],
                                     id.vars = c("classifier", "sampling_size", "sampling_repeat", "species", "ORF_type"))
models <- strsplit(gsub("_ratio", "", dt_final_overlap$variable), "_")
models <- data.table(model_one = sapply(models, function(x) x[1]))
dt_final_overlap[, variable := factor(variable, levels = c("ORFik_ratio", "ORFquant_ratio", "RiboCode_ratio"))]
# models <- data.table(model_one = sapply(models, function(x) x[1]), model_two = lapply(models, function(x) x[2]))
dt_final_overlap <- cbind(dt_final_overlap, models)
dt_final_overlap <- dt_final_overlap[classifier == model_one, ]
# dt_final_overlap <- dt_final_overlap[classifier == model_one | classifier == model_two, ]
# dt_final_overlap[, compare_to := paste0(classifier, "_", variable)]
# dt_final_overlap[value > 1, value := 1]


translon_plot_overlap <- ggplot(data = dt_final_overlap, aes(x = sampling_size, y = value, fill = variable)) +
  geom_boxplot(aes(color = variable)) +
  stat_summary(fun = median, geom = "crossbar", width = 0.75, color = "black", size = 0.15,
               position = position_dodge(width = 0.75)) +
  theme_minimal() + ylab("Predicted ORFs found in another model (%)") +  xlab(paste("Number of Libraries")) +
  theme(legend.position="top") + guides(color="none") +
  guides(fill=guide_legend(title="Classifier")) +
  facet_wrap(ORF_type ~ species, scales = "free"); plot(translon_plot_overlap)



plot_all_versions(translon_plot_overlap, file.path(plot_dir, "Predicted_ORF_overlap_all_subsamplings"),
                  send_to_google_drive = TRUE, send_to_discord = TRUE)
plot_all_versions(translon_plot_orfs, file.path(plot_dir, "Predicted_ORF_count_all_subsamplings"),
                  send_to_google_drive = TRUE, send_to_discord = TRUE)

