library(ORFik); library(data.table); library(fst); library(ggplot2)

plot_folder <- "~/livemount/shared_results/ribocrypt_paper/"
res_folder <- "~/livemount/shared_results/ribocrypt_paper/rds_temp_objects/"
res_rds_file <- file.path(res_folder, "metrics_per_species.fst")
translon_types <- c("uORF", "uoORF","doORF","dORF")
start_codons <- "ATG|TTG|CTG"
file_name <- "predicted_translons_with_sequence"
file_names <- paste0(file_name, c(".fst", ".csv", "_ranges.rds"))
names(file_names) <- c("fst", "csv", "ranges")
all_exp <- list.experiments(validate = FALSE)
all_merged_exps <- all_exp[grep("all_merged", name),][libtypes == "RFP",]
stopifnot(!any(duplicated(all_merged_exps$name)))
message("Making translon metrics per species")
if (!file.exists(res_rds_file)) {
  total_table <- rbindlist(lapply(all_merged_exps$name, function(species) {
    message(species)
    df <- read.experiment(species)
    ribo_dir <- riboORFsFolder(df, QCfolder(df))
    output_dir <- file.path(dirname(df@fafile), "predicted_translons")
    fst_output_file <- file.path(output_dir, file_names["fst"])
    if (!file.exists(fst_output_file)) {
      message("- Species have no translon file!")
      stop()
    }
    try({
      table <- as.data.table(fst::read_fst(file.path(output_dir, file_names["fst"]), columns = c("type", "length", "start_codons", "F1", "F2", "F3","library")))
    })
  }))
  fst::write_fst(total_table, res_rds_file)
} else total_table <- as.data.table(fst::read_fst(res_rds_file))

message("Summarizing metrics")
total_table[, Periodicity_quality := round(F1 / (F1 + F2 + F3), 2)]
total_table[, `:=`(F1 = NULL, F2 = NULL, F3 = NULL)]
total_table[, `:=`(type = as.factor(type), start_codons = as.factor(start_codons))]
total_table[, length := log10(length)]
total_table[, library := as.factor(gsub("all_merged-", "", library))]
total_table[, library_id := as.factor(as.character(as.numeric(library)))]
total_table[, number_of_translons := log10(.N), by = library]

kingdoms <- sapply(unique(total_table$library), function(x) biomartr:::ensembl_assembly_hits(gsub("_", " ", x))$division[1])
kingdoms <- gsub("Ensembl", "", kingdoms)
kingdoms[1] <- "Vertebrates"
total_table[, kingdom := as.factor(kingdoms[match(library, unique(library))])]

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
table_melt <- melt(total_table[, .(length, Periodicity_quality, kingdom, library, library_id)], id.vars = c("kingdom", "library", "library_id"))
levels(table_melt$variable) <- c("Length", "Periodicity quality")
plot_numeric <- ggplot(data = table_melt, aes(y = value, x = kingdom, fill = kingdom)) +
  geom_boxplot() + facet_wrap(~ variable, scales = "free") + theme_minimal() + ylab("") +
  theme(legend.position = "none"); plot_numeric

table_melt <- melt(total_table[, .(type, start_codons, kingdom, library, library_id)], id.vars = c("kingdom", "library", "library_id"))
table_melt <- table_melt[, .(count = .N ), by = .(kingdom, variable, value)]
table_melt <- table_melt[, value_total := round((count / sum(count))*100, 2), by = .(kingdom, variable)]
levels(table_melt$variable) <- c("Translon Type", "Start Codon")
plot_fctr <- ggplot(data = table_melt, aes(y = value_total, x = value, fill = kingdom)) +
  geom_bar(stat="identity", position=position_dodge()) + facet_wrap(~ variable, scales = "free") + theme_minimal() +
  ylab("Usage (%)") + xlab("") +theme(legend.position = "bottom"); plot_fctr

plot <- cowplot::plot_grid(plot_numeric, plot_fctr, ncol = 1); plot
ggsave(file.path(plot_folder, "novel_orf_stats_by_kingdom.jpg"), plot, width = 8, height = 8)
ggsave(file.path(plot_folder, "novel_orf_stats_by_kingdom.svg"), plot, width = 8, height = 8)

