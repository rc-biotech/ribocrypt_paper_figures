library(data.table)
library(ggplot2)
all_exp <- list.experiments(validate = FALSE, "~/Bio_data/ORFik_experiments")
res_folder <- "~/livemount/shared_results/ribocrypt_paper/"
metadata <- fread("~/livemount/Bio_data/NGS_pipeline/FINAL_LIST.csv")
metadata$INHIBITOR[metadata$INHIBITOR == "chx"] <- "Cyclohexamide"
metadata$INHIBITOR[metadata$INHIBITOR == "frozen"] <- "Frozen"
metadata$INHIBITOR[metadata$INHIBITOR == "puro"] <- "Puromycin"
metadata$INHIBITOR[metadata$INHIBITOR == "harr"] <- "Harringtonin"
metadata$INHIBITOR[metadata$INHIBITOR == "ltm"] <- "Lactimidomycin"
metadata$CELL_LINE[metadata$CELL_LINE == "C57BL/6"] <- ""
metadata_used <- fread("~/livemount/Bio_data/NGS_pipeline/FINAL_LIST.csv")
all_libs <- 13023 # The new stats, update when we have the csv

# top5 <- function(x, other_id) {
#   top5_table <- head(sort(table(x[x != "" & x != "NONE"]), decreasing = TRUE), 5)
#   other <- length(x) - sum(top5_table)
#   names(other) <- paste("Other", other_id)
#   return(c(other, top5_table))
# }
top5 <- function(x, other_id, top = 5) {
  full_table <- sort(table(x[x != "" & x != "NONE"]), decreasing = TRUE)
  top5 <- full_table[seq(top)]
  other <- sum(full_table) - sum(top5)
  names(other) <- paste("Other", other_id)
  return(c(other, top5))
}
all_ids <- c("Organism (# Libraries)", "Inhibitor", "Tissues", "Cell Lines", "Filter pass rate", "File breakdown", "CDS' expressed(%) (> 100 RPFs)", "CDS RPFs (log10)")
top5_orgs <- top5(metadata$ScientificName, "Organisms")
top5_inhib <- top5(metadata$INHIBITOR, "Inhibitors")

top5_tissue <- top5(metadata$TISSUE, "Tissues")
top5_cell_line <- top5(metadata$CELL_LINE, "Cell Lines")
passes_quality <- c("All libraries" = all_libs, "High quality libraries" = nrow(metadata))
file_breakdown <- c(FASTA = all_libs, BAM = all_libs, OFST = all_libs, Pshifted_OFST =nrow(metadata),
                    Pshifted_BIGWIG = nrow(metadata))

top <- gsub(" ", "_", names(top5_orgs[-1]))
all_allmerged <- all_exp[grep("all_merged", name)][libtypes == "RFP",][]
all_merged <- grep(paste(top, collapse = "|"), all_allmerged$name, value = T)
rpf_stats <- function(x) {
  df <- read.experiment(x, validate = FALSE, "~/Bio_data/ORFik_experiments")
  txkeep <- filterTranscripts(df, 0,1,0)
  counts <- suppressWarnings(countTable(df, region = "cds"))
  counts <- counts[rownames(counts) %in% txkeep]
  data.table(sum(counts),  round(sum(counts > 100) / nrow(counts), 2)*100)
}
# expression <- rbindlist(lapply(all_allmerged$name, function(x) rpf_stats(x)))

expression <- rbindlist(lapply(all_merged, function(x) rpf_stats(x)))

ordering <- order(expression$V1, decreasing = TRUE)
expression <- expression[ordering,]
psites <- log10(expression[,1][[1]])
cds_ratio_100 <- expression[,2][[1]]
names(psites) <-  names(cds_ratio_100) <- gsub("_", " ", gsub("all_merged-", "", all_merged[ordering]))


all_tables_list <- list(top5_orgs, top5_inhib, top5_tissue, top5_cell_line, passes_quality, file_breakdown, cds_ratio_100, psites)


final <- lapply(seq_along(all_tables_list),
                function(x) {
                  y <- all_tables_list[[x]]
                  data.table(value = y, name = names(y), id = all_ids[x])
                })
final <- rbindlist(final)
final[, id := factor(id, levels = unique(id))]
final[, name := factor(name, levels = rev(unique(name)))]

final_alt <- final[!(id %in% c("Filter pass rate", "File breakdown")), ]
final_alt[name == "Other Inhibitors", name := "Other/missing"]
final_alt <- final_alt[!(grepl("other", name, ignore.case = TRUE) & name != "Other Organisms"),]
ggplot(final_alt, aes(x=name, y=value)) +
  geom_bar(aes(fill = id), stat = "identity") + facet_wrap(~ id, scales = "free", ncol = 2) + theme_minimal() +
  coord_flip() + xlab("") + ylab("") +  theme(legend.position = "none")
ggsave(filename = file.path(res_folder, "database_stats_alt.png"), width = 10, heigh = 5, dpi = 600, units = "in")
ggsave(filename = file.path(res_folder, "database_stats_alt.jpg"), width = 10, heigh = 5, dpi = 600, units = "in")
ggsave(filename = file.path(res_folder, "database_stats_alt.svg"), width = 10, heigh = 5, dpi = 600, units = "in")

# OLD
ggplot(final, aes(x=name, y=value)) +
  geom_bar(aes(fill = id), stat = "identity") + facet_wrap(~ id, scales = "free") + theme_minimal() +
  coord_flip() + xlab("") + ylab("") +  theme(legend.position = "none")
ggsave(filename = file.path(res_folder, "database_stats.png"), width = 10, heigh = 5, dpi = 600, units = "in")
ggsave(filename = file.path(res_folder, "database_stats.pdf"), width = 10, heigh = 5, dpi = 600, units = "in")
ggsave(filename = file.path(res_folder, "database_stats.svg"), width = 10, heigh = 5, dpi = 600, units = "in")
