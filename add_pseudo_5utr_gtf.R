library(ORFik)
df <- read.experiment("/home/rstudio/Bio_data/ORFik_experiments/all_merged-Caenorhabditis_elegans.csv", validate = FALSE)
org_short <- gsub(" ", "_", organism(df))
gtf_path <- ORFik:::getGtfPathFromTxdb(loadTxdb(df))
gtf_out_path <- file.path("~/livemount/shared_data", paste0(org_short, "_mane"),
                          paste0(ORFik:::remove.file_ext(gtf_path, TRUE), "_mane",".",tools::file_ext(gtf_path)))
file.exists(gtf_out_path)



