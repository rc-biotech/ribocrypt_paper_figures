library(ORFik)

df <- read.experiment("all_merged-Homo_sapiens")
table <- setDT(fst::read_fst("~/livemount/Bio_data/references/homo_sapiens/predicted_translons/mane_isoforms_only/predicted_translons_with_sequence.fst"))
predicted <- table$predicted == TRUE
type_hits <- table$type %in% c("uORF", "uoORF")
length_valid <- table$length > 15*3
subset <- predicted & type_hits & length_valid
grl <- readRDS("~/livemount/Bio_data/references/homo_sapiens/predicted_translons/mane_isoforms_only/predicted_translons_with_sequence_ranges.rds")
stopifnot(nrow(table) == length(grl))
length(grl)
grl <- grl[subset]
length(grl)

dt <- table[subset,]
fst::write_fst(dt, "~/livemount/shared_results/ribocrypt_paper/mito_check/uORFs_for_mito_check.fst")
saveRDS(grl, "~/livemount/shared_results/ribocrypt_paper/mito_check/uORFs_for_mito_check_ranges.rds")
