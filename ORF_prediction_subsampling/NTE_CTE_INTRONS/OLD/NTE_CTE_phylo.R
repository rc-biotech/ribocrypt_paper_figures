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

types_all <- list.dirs(file.path(result_dir, "translon-conservation", "results/"), recursive = FALSE, full.names = FALSE)
regions <- c("nte_leader","cte_trailer")
stopifnot(all(regions %in% types_all))
pred_folders <- file.path(result_dir, "translon-conservation", "results", regions, "bed/")
names(pred_folders) <- regions
region <- regions[1]

dt_all <- data.table()
for (region in regions) {
  pred_folder <- pred_folders[region]
  files <- list.files(pred_folder, "_translated")
  files_untrans <- list.files(pred_folder, "_untranslated")
  scores <- c("phylo_csf", "phylo_p", "phylo_best_30nt")
  names(files) <- names(files_untrans) <- scores
  score <- scores[1]
  for (score in scores) {
    message(region, " - ", score)
    phylo <- import(file.path(pred_folder, files_untrans[score]))
    length(phylo)
    stopifnot(length(phylo) > 0)
    phylo_un <- import(file.path(pred_folder, files[score]))
    phylo$status <- "translated"
    phylo_un$status <- "untranslated"
    phylo <- c(phylo, phylo_un)
    names(phylo) <- phylo$name

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
    dt_all <- rbindlist(list(dt_all, data.table(score = score, tr_score, tru_score,
                                                type = type, region, tx = names(tr))))
  }
}


ggplot(dt_all, aes(type, score, fill = region)) + geom_violin(alpha = 0.8) + geom_boxplot(outliers = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", col = "red") + theme_classic() +
  ylim(0, 2.5)

dt_all_noscore <- copy(dt_all)
dt_all_noscore[, score := NULL]
colnames(dt_all_noscore) <- c("predicted_part", "whole_UTR", "type", "region", "tx_id")
dt_all_noscore <- melt.data.table(dt_all_noscore)
dt_all_noscore <- dt_all_noscore[type != "phylo_best_30nt",]
dt_all_noscore[, type := factor(type, unique(type))]
levels(dt_all_noscore$type) <- c("Codon conservation", "Nucleotide conservation")
phylo_plot <- ggplot(dt_all_noscore, aes(value, colour = variable)) + geom_density() +
  facet_wrap(region ~ type) + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot

phylo_plot_NTE <- ggplot(dt_all_noscore[grepl("nte", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type) + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_NTE

phylo_plot_CTE <- ggplot(dt_all_noscore[grepl("cte", region)], aes(value, colour = variable)) + geom_density() +
  facet_wrap(~ type) + theme_classic() + xlab("score") + theme(legend.position = "bottom"); phylo_plot_CTE
#
# massiveNGSpipe:::plot_all_versions(phylo_plot, file.path(plot_dir, "NTE_CTE_phylo_density"), send_to_discord = TRUE)
# counts <- dt_all[, .N, by = .(type, region)][c(1,4)]
#
# count_plot <- ggplot(counts, aes(N, x = region)) + geom_boxplot(outliers = FALSE) + theme_classic()
# massiveNGSpipe:::plot_all_versions(count_plot, file.path(plot_dir, "NTE_CTE_counts"), send_to_discord = TRUE)
#
#
# RFP <- fimport(filepath(df, "cov"))
# gr_NTE <- import(list.files(pred_folders, "_translated", full.names = T)[3:4][2])
# gr_NTE <- split(gr_NTE, gr_NTE$name)
# stopifnot(max(lengths(gr_NTE)) == 1)
#
# extend <- 18
# gr_to_plot <- windowPerGroup(gr_NTE, mrna[names(gr_NTE)], 18, 45)
# gr_to_plot <- extendLeaders(gr_to_plot, extend)
# seqlevelsStyle(gr_to_plot) <- seqlevelsStyle(df)[1]
# dt <- coveragePerTiling(gr_to_plot, RFP, as.data.table = TRUE, withFrames = TRUE)
# cov_final_filtered <- dt[position <= 45,]
# plot_with_cds <- pSitePlot(cov_final_filtered, facet = TRUE, title = "Coverage per Intron index (with CDS flank)")
# plot_with_cds <- plot_with_cds + theme_classic() + xlab("Position of Intron (from first in frame nt relative to CDS)") + ylab("Transcript Normalized counts")
# plot_with_cds
#
# dt_all_phylo <- dt_all[type == "phylo_csf"]
# dt_all_phylo[, `:=`(score_csf = tr_score, score_p = dt_all[type == "phylo_p"]$tr_score,
#                     score_window =dt_all[type == "phylo_best_30nt"]$tr_score)]
# dt_all_phylo[, `:=`(tr_score = NULL, tru_score = NULL, score = NULL, type = NULL)]
# dt_local <- rbindlist(lapply(list.files(localization_dir, full.names = T), function(x) fread(x)), fill = TRUE, idcol="file")
# ids <- sub("_all.*", "", list.files(localization_dir))
# dt_local[, file := factor(file, labels = ids)]
# dt_local[file == "cte", tx := Protein_ID]
# dt_local[, Protein_ID := NULL]
# dt_local[, region := file]; dt_local$file <- NULL
# levels(dt_local$region) <- c("cte_trailer", "nte_leader")
# cols_numeric <- which(sapply(dt_local, is.numeric))  # Identify numeric columns
# dt_local[, (names(dt_local)[cols_numeric]) := lapply(.SD, abs), .SDcols = cols_numeric]
# dt_merged_pl <- merge.data.table(dt_all_phylo, dt_local, by = c("region", "tx"), all = FALSE)
#
# dt_merged_pl_nte <- dt_merged_pl[region == "nte_leader"]
# dt_merged_pl_cte <- dt_merged_pl[region == "cte_trailer"]
# cor_table <- cor(dt_merged_pl_nte[, !c("region", "tx")], use = "pairwise.complete.obs", method = "spearman")
# abs(cor_table) > 0.1
#
# cor_table <- cor(dt_merged_pl_nte[tr_score > quantile(score, 0.99), !c("region", "tx")], use = "pairwise.complete.obs")
# abs(cor_table) > 0.2
#
# cor_table <- cor(dt_merged_pl_cte[, !c("region", "tx")], use = "pairwise.complete.obs")
# abs(cor_table) > 0.2
#
#
#
# lm_nte <- lm(score_csf + score_p + score_window ~ ., data = dt_merged_pl_nte[, !c("region", "tx")])
# lm_cte <- lm(score_csf + score_p + score_window ~ ., data = dt_merged_pl_cte[, !c("region", "tx", "mitoprot", "signalP")])
# summary(lm_nte)
# summary(lm_cte)
