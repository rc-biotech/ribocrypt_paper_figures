library(massiveNGSpipe)
library(ggplot2)
library(rtracklayer)

df <- read.experiment("all_merged-Homo_sapiens_04_oct_2024_all")
mrna <- loadRegion(df, "mrna")
symbols <- symbols(df)
symbols[, uniprot_id := NULL]

# Load and clone jacks library

localization_dir <- "~/livemount/shared_results/ribocrypt_paper/localization"
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
plot_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/CTE_NTE_analysis"

fetch_github <- T
github_forks_repo <- "~/livemount/forks"
conservation_dir <- file.path(github_forks_repo, "translon-conservation/")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir_exists <- dir.exists(conservation_dir)
if (!dir_exists | fetch_github) {
  message("Fetching github repo of phylo data")
  if (dir_exists) system(paste("rm -r", paste0(conservation_dir, "/")))
  system(paste("cd", github_forks_repo, "&& git clone git@github.com:JackCurragh/translon-conservation.git"))
}

scores <- c("PhyloCSF", "PhyloP")
scores_dir <- paste0(scores, "/")
results_paths <- file.path(conservation_dir, "updated_results", scores_dir)
names(results_paths) <- scores
path <- results_paths[2]
type_list <- lapply (scores, function(score) {
  path <- results_paths[score]
  types_all <- list.files(path, recursive = FALSE, full.names = FALSE)
  types_all <- types_all[grep("_old_|CDS_PhyloP|CTEs|trailers_p|trailers_P|introns_", types_all, invert = TRUE)]
  score <- ifelse(score == "PhyloCSF", "phyloCSF", score)
  types <- gsub(paste0("_", score, ".*"), "", types_all)
  types <- gsub("_5k", "", types)
  types <- gsub("_new", "", types)
  types <- gsub("_equal_length.bed6", "", types)
  types <- gsub("intron_candidates", "translated_introns", types)
  types <- gsub("intron_not_predicted", "nontranslated_introns", types)
  if (length(unique(types)) != length(types)) browser()
  stopifnot(length(unique(types)) == length(types))
  pred_files <- file.path(path, types_all)
  stopifnot(all(file.exists(pred_files)))
  names(pred_files) <- types
  pred_files <- pred_files[order(types)]
  return(pred_files)
})
names(type_list) <- scores

groups <- names(type_list[[1]])
groups_dt <- data.table(groups, background = FALSE, non_predicted = FALSE)
groups_dt[groups %in% c("leaders", "trailers", "nontranslated_introns"), background := TRUE]
groups_dt[groups %in% c("nontranslated_introns"), non_predicted := TRUE]
group_compare <- groups
group_compare[groups == "leaders"] <- "NTEs"
group_compare[groups %in% c("CTE", "trailers")] <- "CTEs"
group_compare[groups %in% c("nontranslated_introns", "translated_introns")] <- "Introns"
groups_dt[, group_compare := group_compare]
groups_dt <- cbind(groups_dt, setDT(type_list))
groups <- split(groups, group_compare)

dt_all <- dt_regions <- data.table()
score <- scores[1]
unique_groups <- unique(names(groups))
region <- unique_groups[1]
for (score in scores) {
  for (region in unique_groups) {
    message(region, " - ", score)
    file <- groups_dt[group_compare == region & background == FALSE, ..score][[1]]
    phylo <- import(file, format = "bed")
    # if (region == "Introns") browser()
    stopifnot(!is.null(phylo$name) & !is.null(phylo$score))
    tr_N <- length(phylo); tr_N
    tru_N <- overlap_N <- overlap_N_all_len <- 0
    stopifnot(length(phylo) > 0)

    has_background <- sum(groups_dt[group_compare == region,]$background) > 0 & region != "Introns"
    has_nonpredicted <- sum(groups_dt[group_compare == region,]$non_predicted) > 0
    if (has_background | has_nonpredicted) {
      file_un <- groups_dt[group_compare == region & background == TRUE,..score][[1]]
      phylo_un <- import(file_un, format = "bed")
      stopifnot(!is.null(phylo_un$name) & !is.null(phylo_un$score))
      tru_N <- length(phylo_un)
      phylo$status <- "translated"
      phylo_un$status <- "untranslated"
      names(phylo) <- phylo$name
      names(phylo_un) <- phylo_un$name
      if (has_background) {
        trailing_numbers_tr <- length(grep("_\\d+$", phylo$name[1])) > 0
        trailing_numbers_tru <- length(grep("_\\d+$", phylo_un$name[1])) > 0
        trailing_both <- trailing_numbers_tr & trailing_numbers_tru
        phylo <- suppressWarnings(c(phylo, phylo_un))

        if (!trailing_both) {
          names(phylo) <- gsub("_.*", "", names(phylo))
        }

        to_keep <- names(phylo[phylo$status == "untranslated"])
        phylo <- phylo[names(phylo) %in% to_keep]
        to_keep <- names(phylo[phylo$status == "translated"])
        phylo <- phylo[names(phylo) %in% to_keep]
        length(phylo[phylo$status == "translated"])
        if (length(phylo[phylo$status == "translated"]) == 0) {
          warning("No matching tx_ids found for translated/untranslated for region(score):
                ", region, "(", score, ")")
          next
        }
        tr <- phylo[phylo$status == "translated"]
        overlap_N_all_len <- length(tr); length(tr)
        tr <- firstExonPerGroup(split(tr, names(tr)))
        tru <- phylo[phylo$status == "untranslated"]
        tru <- firstExonPerGroup(split(tru, names(tru)))
        unmatching_introns <- names(which(lengths(tr) != lengths(tru)))
        tr <- tr[!(names(tr) %in% unmatching_introns)]
        tru <- tru[!(names(tru) %in% unmatching_introns)]
        stopifnot(names(tr) == names(tru))
      } else {
        tr <- phylo[phylo$status == "translated"]
        overlap_N_all_len <- NA
        tr <- firstExonPerGroup(split(tr, names(tr)))
        tru <- phylo_un[phylo_un$status == "untranslated"]
        tru <- firstExonPerGroup(split(tru, names(tru)))
        min_size <- min(length(tr), length(tru))
        stopifnot(min_size > 0)
        tr <- tr[seq(min_size)]
        tru <- tru[seq(min_size)]
      }
      overlap_N <- length(tr); length(tr)
      tr_width <- widthPerGroup(tr, FALSE)
      tru_width <- widthPerGroup(tru, FALSE)
      tr_score <- tr@unlistData$score
      tru_score <- tru@unlistData$score
      if (region == "Introns") { # Make sure background set never overlaps cds, set those to NA
        pred <- qs::qread("~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/intron_candidates_intron_1_to_7_ranges.qs")
        back <- qs::qread("~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/intron_all_not_predicted_equal_size_1_to_7_ranges.qs")
        pred <- pred[txNames(pred) %in% names(phylo)]
        back <- back[names(back) %in% names(phylo)]
        back <- back[countOverlaps(back, pred) > 0]
        back <- back[countOverlaps(back, pred) > 0]
        cds <- loadRegion(df, "cds")
        whole_intron_overlaps_any_cds <- names(back[countOverlaps(back, cds, minoverlap = 2) > 0])
        tru_score[names(tru) %in% whole_intron_overlaps_any_cds] <- NA
      }

      ratio <- (tr_score+ 1) / (tru_score + 1)
      tx_ids <- names(tr@unlistData)
    } else {
      names(phylo) <- gsub("_.*|\\..*", "", phylo$name)
      ratio <- NA
      tr <- phylo
      tr_score <- tr$score
      tru_score <- NA
      tx_ids <- names(phylo)
      tr_width <- as.integer(width(tr))
      tru_width <- NA
      tr_N <- length(phylo)
    }
    dt_regions <- rbindlist(list(dt_regions, data.table(region, score, tr_N, tru_N, overlap_N_all_len, overlap_N)))
    dt_all <- rbindlist(list(dt_all, data.table(score = ratio, tr_score, tru_score,
                                                type = score, region, tx = tx_ids,
                                                tr_width, tru_width)))
  }
}

dt_regions


# ggplot(dt_all, aes(type, score, fill = region)) + geom_violin(alpha = 0.8) + geom_boxplot(outliers = FALSE) +
#   geom_hline(yintercept = 1, linetype = "dashed", col = "red") + theme_classic() +
#   ylim(0, 2.5)

dt_all_noscore <- copy(dt_all)
dt_all_noscore[, score := NULL]; dt_all_noscore[, tr_width := NULL]; dt_all_noscore[, tru_width := NULL]
colnames(dt_all_noscore) <- c("predicted_part", "whole_UTR", "type", "region", "tx_id")
dt_all_noscore <- melt.data.table(dt_all_noscore)
dt_all_noscore <- dt_all_noscore[type != "phylo_best_30nt",]
dt_all_noscore[, type := factor(type, unique(type))]
levels(dt_all_noscore$type) <- c("Codon conservation", "Nucleotide conservation")
dt_all_noscore[, merged := paste(region, variable, sep = "_")]
# Show total numbers per region:
dt_all_noscore[type == "Codon conservation" & variable == "predicted_part"][, .N, by = region]



phylo_plot <- ggplot(dt_all_noscore, aes(value, colour = variable)) + geom_density() +
  facet_wrap(region ~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot

library(ggplot2)
library(scales)

# Custom pseudo-log transformation
pseudo_log <- trans_new(
  name = "pseudo_log",
  transform = function(x) sign(x) * log1p(abs(x)),  # Transformation formula
  inverse = function(x) sign(x) * (exp(abs(x)) - 1)  # Inverse transformation
)

phylo_plot_merged_pseudo_scale <- ggplot(dt_all_noscore, aes(y = value, x = region, fill = variable)) + geom_boxplot(outliers = FALSE) +
  facet_wrap( ~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom") +
  scale_y_continuous(trans = pseudo_log); phylo_plot_merged_pseudo_scale

phylo_plot_merged <- ggplot(dt_all_noscore, aes(y = value, x = region, fill = variable)) + geom_boxplot(outliers = FALSE) +
  facet_wrap( ~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_merged


phylo_plot_NTE <- ggplot(dt_all_noscore[grepl("NTE", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_NTE

phylo_plot_CTE <- ggplot(dt_all_noscore[grepl("CTE", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_CTE

phylo_plot_intronic <- ggplot(dt_all_noscore[grepl("Introns", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type, scales = "free") + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_intronic

PHYLO_IS_DONE <- TRUE # Check flag
# Intron sanity checks:
# Length of Top 10 worst Intron score regions
width(pred[paste0(dt_all_noscore[variable == "predicted_part" & type == "Codon conservation" &
                                   tx_id %in% whole_intron_overlaps_any_cds & grepl("Introns", region),]
                  [order(value, decreasing = F)][1:10,]$tx_id, "_1")])
# Length of Top 10 best Intron score regions
width(pred[paste0(dt_all_noscore[variable == "predicted_part" & type == "Codon conservation" &
                                   tx_id %in% whole_intron_overlaps_any_cds & grepl("Introns", region),]
                  [order(value, decreasing = T)][1:10,]$tx_id, "_1")])

# Correlation checks:
dt_all_sub <- dt_all[!(region %in% c("intergenic_100", "intergenic_50"))]

dt_all_sub[!is.na(tr_width), .(cor_xy = round(cor(tr_score, tr_width, use = "pairwise.complete.obs", method = "spearman"), 2)), by = .(type, region)]
dt_all_sub[!is.na(tru_width), .(cor_xy = round(cor(tru_score, tru_width, use = "pairwise.complete.obs", method = "spearman"), 2)), by = .(type, region)]
