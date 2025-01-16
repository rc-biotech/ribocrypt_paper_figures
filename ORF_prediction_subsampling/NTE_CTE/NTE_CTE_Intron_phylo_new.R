library(massiveNGSpipe)
library(ggplot2)
library(rtracklayer)

df <- read.experiment("all_merged-Homo_sapiens_04_oct_2024_all")
mrna <- loadRegion(df, "mrna")
symbols <- symbols(df)
symbols[, uniprot_id := NULL]

# Load and clone jacks library
fetch_github <- F
localization_dir <- "~/livemount/shared_results/ribocrypt_paper/localization"
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/CTE_NTE_analysis"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir_exists <- dir.exists(file.path(result_dir, "translon-conservation/"))
if (!dir_exists | fetch_github) {
  message("Fetching github repo of phylo data")
  if (dir_exists) system(paste("rm -r", file.path(result_dir, "translon-conservation/")))
  system(paste("cd", result_dir, "&& git clone git@github.com:JackCurragh/translon-conservation.git"))
}
scores <- c("PhyloCSF", "PhyloP")
scores_dir <- paste0(scores, "/")
results_paths <- file.path(result_dir, "translon-conservation", "updated_results", scores_dir)
names(results_paths) <- scores
path <- results_paths[2]
type_list <- lapply (scores, function(score) {
  path <- results_paths[score]
  types_all <- list.files(path, recursive = FALSE, full.names = FALSE)
  types_all <- types_all[grep("_old_|CDS_PhyloP", types_all, invert = TRUE)]
  score <- ifelse(score == "PhyloCSF", "phyloCSF", score)
  types <- gsub(paste0("_", score, ".*"), "", types_all)
  types <- gsub("_5k", "", types)
  stopifnot(length(unique(types)) == length(types))
  pred_files <- file.path(path, types_all)
  stopifnot(all(file.exists(pred_files)))
  names(pred_files) <- types
  pred_files <- pred_files[order(types)]
  return(pred_files)
})
names(type_list) <- scores

groups <- names(type_list[[1]])
groups_dt <- data.table(groups, background = FALSE)
groups_dt[groups %in% c("leaders", "trailers", "nontranslated_introns"), background := TRUE]
group_compare <- groups
group_compare[groups == "leaders"] <- "NTEs"
group_compare[groups == "trailers"] <- "CTEs"
group_compare[groups %in% c("nontranslated_introns", "translated_introns")] <- "Introns"
groups_dt[, group_compare := group_compare]
groups_dt <- cbind(groups_dt, setDT(type_list))
groups <- split(groups, group_compare)

dt_all <- data.table()
score <- scores[1]
unique_groups <- unique(names(groups))
g <- unique_groups[1]
for (score in scores) {
  for (g in unique_groups) {
    t <- g
    message(t, " - ", score)
    file <- groups_dt[group_compare == g & background == FALSE,..score][[1]]
    phylo <- import(file, format = "bed")
    length(phylo)
    stopifnot(length(phylo) > 0)
    if (nrow(groups_dt[group_compare == g,]) > 1) {
      file_un <- groups_dt[group_compare == g & background == TRUE,..score][[1]]
      phylo_un <- import(file_un, format = "bed")
      phylo$status <- "translated"
      phylo_un$status <- "untranslated"
      phylo <- suppressWarnings(c(phylo, phylo_un))
      names(phylo) <- gsub("_.*", "", phylo$name)

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
      ratio <- (tr_score+ 1) / (tru_score + 1)
    } else {
      names(phylo) <- gsub("_.*|\\..*", "", phylo$name)
      ratio <- NA
      tr <- phylo
      tr_score <- tr$score
      tru_score <- NA
    }

    dt_all <- rbindlist(list(dt_all, data.table(score = ratio, tr_score, tru_score,
                                                type = score, region = t, tx = names(tr))))
  }
}



# ggplot(dt_all, aes(type, score, fill = region)) + geom_violin(alpha = 0.8) + geom_boxplot(outliers = FALSE) +
#   geom_hline(yintercept = 1, linetype = "dashed", col = "red") + theme_classic() +
#   ylim(0, 2.5)

dt_all_noscore <- copy(dt_all)
dt_all_noscore[, score := NULL]
colnames(dt_all_noscore) <- c("predicted_part", "whole_UTR", "type", "region", "tx_id")
dt_all_noscore <- melt.data.table(dt_all_noscore)
dt_all_noscore <- dt_all_noscore[type != "phylo_best_30nt",]
dt_all_noscore[, type := factor(type, unique(type))]
levels(dt_all_noscore$type) <- c("Codon conservation", "Nucleotide conservation")
dt_all_noscore[, merged := paste(region, variable, sep = "_")]
phylo_plot <- ggplot(dt_all_noscore, aes(value, colour = variable)) + geom_density() +
  facet_wrap(region ~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot

phylo_plot_merged <- ggplot(dt_all_noscore, aes(y = value, x = region, fill = variable)) + geom_boxplot() +
  facet_wrap( ~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom") + ylim(c(-1000, 1000)); phylo_plot_merged


phylo_plot_NTE <- ggplot(dt_all_noscore[grepl("NTE", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_NTE

phylo_plot_CTE <- ggplot(dt_all_noscore[grepl("CTE", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_CTE

phylo_plot_intronic <- ggplot(dt_all_noscore[grepl("Introns", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_intronic
