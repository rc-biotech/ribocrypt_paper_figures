library(massiveNGSpipe)
library(ggplot2)
library(rtracklayer)

# Load and clone jacks library
fetch_github <- FALSE
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir_exists <- dir.exists(file.path(result_dir, "translon-conservation/"))
if (!dir_exists | fetch_github) {
  message("Fetching github repo of phylo data")
  if (dir_exists) system(paste("rm -r", file.path(result_dir, "translon-conservation/")))
  system(paste("cd", result_dir, "&& git clone git@github.com:JackCurragh/translon-conservation.git"))
}
dt_all_all <- data.table()
classifications <- c("", "untranslated_")
classification <- classifications[2]
for (classification in classifications) {
  pred_folder <- file.path(result_dir, "translon-conservation", "results",
                           paste0(classification, "intronic_extended/bed/"))
  files <- list.files(pred_folder, "20241202")[1:3]
  names(files) <- c("phylo_csf", "phylo_p", "phylo_best_30nt")
  dt_all <- data.table()
  for (type in names(files)) {
    message(type)
    phylo <- import(file.path(pred_folder, files[type]))
    names(phylo) <- sub("_.*", "", phylo$name)
    phylo$status <- sub(".*_1_", "", phylo$name)

    to_keep <- names(phylo[phylo$status == "untranslated"])
    phylo <- phylo[names(phylo) %in% to_keep]
    to_keep <- names(phylo[phylo$status == "translated"])
    phylo <- phylo[names(phylo) %in% to_keep]
    length(phylo[phylo$status == "translated"])

    tr <- phylo[phylo$status == "translated"]
    tr <- split(tr, names(tr))
    tru <- phylo[phylo$status == "untranslated"]
    tru <- split(tru, names(tru))

    unmatching_introns <- names(which(lengths(tr) != lengths(tru)))
    tr <- tr[!(names(tr) %in% unmatching_introns)]
    tru <- tru[!(names(tru) %in% unmatching_introns)]
    stopifnot(names(tr) == names(tru))
    length(tr)

    tr_score <- tr@unlistData$score
    tru_score <- tru@unlistData$score
    score <- (tr_score+ 1) / (tru_score + 1)
    predicted <- ifelse(classification == "", "Predicted", "Non_Predicted")
    dt_all <- rbindlist(list(dt_all, data.table(score = score, tr_score, tru_score, type = type, predicted = predicted)))
  }
  dt_all_all <- rbindlist(list(dt_all_all, dt_all))
}


ggplot(dt_all, aes(type, score)) + geom_violin(alpha = 0.8) + geom_boxplot(outliers = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") + theme_classic() +
  ylim(0, 2.5)

dt_all_noscore <- copy(dt_all_all)
dt_all_noscore[, score := NULL]
colnames(dt_all_noscore) <- c("Candidate_part", "Rest_of_intron", "type", "predicted")
dt_all_noscore[, predicted := factor(predicted, unique(predicted))]
dt_all_noscore <- melt.data.table(dt_all_noscore)
dt_all_noscore[, Group := paste0(variable, "_", predicted)]
dt_all_noscore[, Group := factor(Group, unique(Group))]
dt_all_noscore <- dt_all_noscore[type != "phylo_best_30nt",]
dt_all_noscore[, type := factor(type, unique(type))]
levels(dt_all_noscore$type) <- c("Codon conservation", "Amino acid conservation")
phylo_plot <- ggplot(dt_all_noscore, aes(value, colour = Group)) + geom_density() +
  facet_wrap(~ type) + theme_classic() + xlab("score") + theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2)); phylo_plot

# massiveNGSpipe:::plot_all_versions(phylo_plot, file.path(result_dir, "intron_phylo_density"), send_to_discord = TRUE)
