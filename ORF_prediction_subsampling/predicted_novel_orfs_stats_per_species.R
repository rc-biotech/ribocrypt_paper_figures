library(ORFik); library(data.table); library(fst); library(ggplot2)
discordr::set_default_discord_connection(discordr::import_discord_connections("~/livemount/.cache/discordr/r_to_discord_config")[[1]])
plot_folder <- "~/livemount/shared_results/ribocrypt_paper/"
res_folder <- "~/livemount/shared_results/ribocrypt_paper/rds_temp_objects/"
models <- c("ORFik", "RiboCode", "ORFquant")
res_rds_file <- file.path(res_folder, paste0(models, "_metrics_per_species.fst"))
names(res_rds_file) <- models
translon_types <- c("uORF", "uoORF","doORF","dORF")
start_codons <- "ATG|TTG|CTG"
file_name <- "predicted_translons_with_sequence"
file_names <- paste0(file_name, c(".fst", ".csv", "_ranges.rds"))
names(file_names) <- c("fst", "csv", "ranges")
all_exp <- list.experiments(pattern = "all_merged", validate = FALSE)
all_merged_exps <- all_exp[grep("all_merged", name),][libtypes == "RFP",]
stopifnot(!any(duplicated(all_merged_exps$name)))
stopifnot(length(all_merged_exps) > 0)
message("Making translon metrics per species")

get_predictions_all_species <- function(all_merged_exps) {

}
### ORFik
# model <- "ORFik"
# if (!file.exists(res_rds_file[model])) {
#   total_table <- rbindlist(lapply(all_merged_exps$name, function(species) {
#     message(species)
#     df <- read.experiment(species, validate = FALSE)
#     ribo_dir <- riboORFsFolder(df, QCfolder(df))
#     output_dir <- file.path(dirname(df@fafile), "predicted_translons")
#     fst_output_file <- file.path(output_dir, file_names["fst"])
#     if (!file.exists(fst_output_file)) {
#       message("- Species have no translon file!")
#       stop()
#     }
#     try({
#       table <- as.data.table(fst::read_fst(file.path(output_dir, file_names["fst"]), columns = c("type", "length", "start_codons", "F1", "F2", "F3","library")))
#     })
#   }))
#   fst::write_fst(total_table, res_rds_file[model])
# } else total_table <- setDT(fst::read_fst(res_rds_file[model]))
model <- "ORFik"
if (!file.exists(res_rds_file[model])) {
  total_table <- rbindlist(lapply(all_merged_exps$name, function(species) {
    message(species)
    df <- read.experiment(species, validate = FALSE)
    species_short <- gsub("all_merged-", "",species)
    model_folder <- file.path(dirname(df@fafile), "predicted_translons", model)

    out_file_ORFik <- list.files(model_folder, pattern = "_table.rds$", full.names = TRUE)
    stopifnot(length(out_file_ORFik) == 1)
    if (!file.exists(out_file_ORFquant)) {
      message("- Species have no translon file!")
      stop()
    }
    dt_done <- data.table()
    try({
      dt <- readRDS(out_file_ORFik)[, c("type", "length", "start_codons", "F1", "F2", "F3")]
      dt[, library := species_short]
      dt_done <- dt
    })
    return(dt_done)
  }))
  message("Summarizing metrics")

  total_table[, Periodicity_quality := round(F1 / (F1 + F2 + F3), 2)]
  total_table[, `:=`(F1 = NULL, F2 = NULL, F3 = NULL)]
  total_table[, `:=`(type = as.factor(type), start_codons = as.factor(start_codons))]
  total_table[, length := log10(length)]
  total_table[, library := as.factor(library)]
  total_table[, library_id := as.factor(as.character(as.numeric(library)))]
  total_table[, number_of_translons := log10(.N), by = library]
  total_table[type == "Overlap_uORF", type := "uoORF"]
  total_table[type == "Overlap_dORF", type := "doORF"]
  total_table[, Classifier := as.factor(model)]
  fst::write_fst(total_table, res_rds_file[model])
} else total_table <- setDT(fst::read_fst(res_rds_file[model]))
### RiboCode
model <- "RiboCode"
if (!file.exists(res_rds_file[model])) {
  total_table <- rbindlist(lapply(all_merged_exps$name, function(species) {
    message(species)
    df <- read.experiment(species, validate = FALSE)
    species_short <- gsub("all_merged-", "",species)
    model_folder <- file.path(dirname(df@fafile), "predicted_translons", model)
    out_prefix_model <- file.path(model_folder, paste0(model, "_predicted_translons_", species_short))
    out_file_RiboCode <- paste0(out_prefix_model, "_collapsed.txt")
    if (!file.exists(out_file_RiboCode)) {
      message("- Species have no translon file!")
      stop()
    }
    try({
      dt <- fread(out_file_RiboCode)[, c("ORF_type", "ORF_length", "start_codon","Psites_sum_frame0", "Psites_sum_frame1", "Psites_sum_frame2")]
      colnames(dt) <- c("type", "length", "start_codons", "F1", "F2", "F3")
      dt[, library := species_short]
    })
  }))
  message("Summarizing metrics")

  total_table[, Periodicity_quality := round(F1 / (F1 + F2 + F3), 2)]
  total_table[, `:=`(F1 = NULL, F2 = NULL, F3 = NULL)]
  total_table[, `:=`(type = as.factor(type), start_codons = as.factor(start_codons))]
  total_table[, length := log10(length)]
  total_table[, library := as.factor(library)]
  total_table[, library_id := as.factor(as.character(as.numeric(library)))]
  total_table[, number_of_translons := log10(.N), by = library]
  total_table[type == "Overlap_uORF", type := "uoORF"]
  total_table[type == "Overlap_dORF", type := "doORF"]
  total_table[, Classifier := as.factor(model)]
  fst::write_fst(total_table, res_rds_file[model])
} else total_table <- setDT(fst::read_fst(res_rds_file[model]))

model <- "ORFquant"
if (!file.exists(res_rds_file[model])) {
  total_table <- rbindlist(lapply(all_merged_exps$name, function(species) {
    message(species)
    df <- read.experiment(species, validate = FALSE)
    species_short <- gsub("all_merged-", "",species)
    model_folder <- file.path(dirname(df@fafile), "predicted_translons", model)
    out_prefix_model <- file.path(model_folder, paste0(model, "_predicted_translons_", species_short))
    out_file_ORFquant <- paste0(out_prefix_model, ".fst")
    if (!file.exists(out_file_ORFquant)) {
      message("- Species have no translon file!")
      stop()
    }
    dt_done <- data.table()
    try({
      dt <- setDT(fst::read_fst(out_file_ORFquant))
      ir <- IRanges(dt$start, width = 3)
      irl <- split(ir, dt$tx_id)
      tx_names <- names(irl)
      names(irl) <- seq_along(irl)
      start_codons <- ORFik::txSeqsFromFa(pmapFromTranscriptF(irl, loadRegion(df)[tx_names], removeEmpty = TRUE), df@fafile)
      dt_mix <- dt[, c("ORF_type")]
      dt <- cbind(dt_mix, length = dt$end - dt$start, start_codons = as.factor(as.character(start_codons)), dt[, c("frame0", "frame1", "frame2")])
      colnames(dt) <- c("type", "length", "start_codons", "F1", "F2", "F3")
      dt[, library := species_short]
      dt_done <- dt
    })
    return(dt_done)
  }))
  message("Summarizing metrics")

  total_table[, Periodicity_quality := round(F1 / (F1 + F2 + F3), 2)]
  total_table[, `:=`(F1 = NULL, F2 = NULL, F3 = NULL)]
  total_table[, `:=`(type = as.factor(type), start_codons = as.factor(start_codons))]
  total_table[, length := log10(length)]
  total_table[, library := as.factor(library)]
  total_table[, library_id := as.factor(as.character(as.numeric(library)))]
  total_table[, number_of_translons := log10(.N), by = library]
  total_table[type == "Overlap_uORF", type := "uoORF"]
  total_table[type == "Overlap_dORF", type := "doORF"]
  total_table[, Classifier := as.factor(model)]
  fst::write_fst(total_table, res_rds_file[model])
} else total_table <- setDT(fst::read_fst(res_rds_file[model]))


total_table <- rbindlist(lapply(res_rds_file, function(x) fst::read_fst(x)))


species <- gsub("_", " ", unique(total_table$library))
kingdoms_temp <- sapply(species, function(x) {
  res <- biomartr:::ensembl_assembly_hits(x)
  if (is.logical(res)) return(as.character(NA))
  return(res$division[1])
  })
kingdoms <- gsub("Ensembl", "", kingdoms_temp)
kingdoms[species %in% c("Aedes aegypti", "Saccharomyces cerevisiae")] <- c("Vertebrates", "Fungi")
kingdoms[kingdoms %in% "Vertebrates"] <- "v/Animalia"
total_table[, kingdom := as.factor(kingdoms[match(library, unique(library))])]
kingdom_dt <- data.table(species = names(kingdoms), kingdoms)
kingdom_dt_path <- file.path(res_folder, "kingdom_information.csv")
fwrite(kingdom_dt, file = kingdom_dt_path)
googledrive::drive_upload(kingdom_dt_path,
                          "https://drive.google.com/drive/folders/1CAN2_v5TbnfTt6plG-ZkEkhpOvdfMz4o",
                          name = basename(kingdom_dt_path))
# BY species
# table_melt <- melt(total_table[, .(length, Periodicity_quality, library, library_id)], id.vars = c("library", "library_id"))
# levels(table_melt$variable) <- c("Length", "Periodicity quality")
# plot_numeric <- ggplot(data = table_melt, aes(y = value, x = library_id, fill = library)) +
#   geom_boxplot() + facet_wrap(~ variable, scales = "free") + theme_minimal() +
#   theme(legend.position = "none"); plot_numeric
#
# table_melt <- melt(total_table[, .(type, start_codons, library, library_id)], id.vars = c("library", "library_id"))
# table_melt <- table_melt[, .(count = .N ), by = .(library, variable, value)]
# table_melt <- table_melt[, value_total := round((count / sum(count))*100, 2), by = .(library, variable)]
# levels(table_melt$variable) <- c("Translon Type", "Start Codon")
# plot_fctr <- ggplot(data = table_melt, aes(y = value_total, x = value, fill = library)) +
#   geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~ variable, scales = "free") + theme_minimal() +
#   ylab("Usage (%)") + xlab("") +theme(legend.position = "bottom"); plot_fctr
#
# plot <- cowplot::plot_grid(plot_numeric, plot_fctr, ncol = 1); plot
# ggsave(file.path(plot_folder, "novel_orf_stats_by_species.jpg"), plot, width = 8, height = 8)
# ggsave(file.path(plot_folder, "novel_orf_stats_by_species.svg"), plot, width = 8, height = 8)
# By kingdom
message("Plotting")
table_melt <- data.table::melt.data.table(total_table[, .(length, Periodicity_quality, kingdom, library, library_id, Classifier)], id.vars = c("kingdom", "Classifier", "library", "library_id"))
levels(table_melt$variable) <- c("Length", "Periodicity quality")
plot_numeric <- ggplot(data = table_melt, aes(y = value, x = kingdom, fill = Classifier)) +
  geom_boxplot(outliers = FALSE) + facet_wrap(~ variable, scales = "free") + theme_minimal() + ylab("") +
  theme(legend.position = "top"); plot_numeric

table_melt <- data.table::melt.data.table(total_table[, .(type, start_codons, kingdom, library, library_id, Classifier)], id.vars = c("kingdom", "Classifier","library", "library_id"))
table_melt <- table_melt[!(variable == "start_codons" & !(value %in% c("ATG", "GTG", "TTG", "CTG"))), ]
table_melt <- table_melt[!(variable == "type" & !(value %in% c("annotated", "uORF", "uoORF", "NTE"))), ]
table_melt <- table_melt[, .(count = .N ), by = .(Classifier, kingdom, variable, value)]
table_melt <- table_melt[, value_total := round((count / sum(count))*100, 2), by = .(Classifier, kingdom, variable)]
levels(table_melt$variable) <- c("Translon Type", "Start Codon")
plot_fctr <- ggplot(data = table_melt, aes(y = value_total, x = value, fill = kingdom)) +
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap(variable ~ Classifier, scales = "free") + theme_minimal() +
  ylab("Usage (%)") + xlab("") +theme(legend.position = "bottom"); plot_fctr

plot <- cowplot::plot_grid(plot_numeric, plot_fctr, ncol = 1); plot
ggsave(file.path(plot_folder, "novel_orf_stats_by_kingdom.jpg"), plot, width = 8, height = 8)
ggsave(file.path(plot_folder, "novel_orf_stats_by_kingdom.svg"), plot, width = 8, height = 8)
googledrive::drive_upload(file.path(plot_folder, "novel_orf_stats_by_kingdom.svg"),
                       "https://drive.google.com/drive/folders/1CAN2_v5TbnfTt6plG-ZkEkhpOvdfMz4o",
                       name = "novel_orf_stats_by_kingdom.svg")
discordr::send_webhook_ggplot(filename = tempfile(pattern = "discordr",
                                                  fileext = ".jpg"))

