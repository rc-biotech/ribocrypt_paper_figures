## Run all classifiers (libraries are preload upstream)
library(ORFik)
library(ORFik) # Needed for ORFik:::take_Fvals_spect
library(data.table)


message("Starting classifiers:")
start_time <- Sys.time()
message(start_time)


# Main paths
result_dir_all <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"
covRLE_subsampling_dir_all <- "~/livemount/shared_results/ribocrypt_paper//rds_temp_objects/covRLE_whole_library_subsamplings"
## Create dirs
dir.create(result_dir_all, recursive = TRUE, showWarnings = FALSE)

# Parameters
export_metrics_table <- FALSE
# Include NTE to make it fair comparison
ORF_categories_to_keep <- c("uORF", "uoORF", "annotated", "NTE") # "NTE", "NTT", "internal", "doORF", "dORF"
name_of_result <- "uORF"
longestORF = TRUE
startCodon = paste(startDefinition(1), "GTG",sep = "|")
stopCodon = stopDefinition(1)
minimumLength = 0


remake <- TRUE
organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
org <- organisms[1]
for (org in organisms) {
  message(org)
  message("-------------------------------")
  df <- read.experiment(gsub("all_samples", "all_merged", org), validate = F)
  dfc <- read.experiment(org, validate = F)
  org_short <- gsub("all_samples-", "", org)
  result_dir <- file.path(result_dir_all, org_short)
  ORFik_dir <- file.path(result_dir, "ORFik/")
  file_suffix_last <- paste("ORFik_results", "samp", nrow(dfc),"num", 1,
                             name(df), ORFik:::name_decider(df, naming = "full"),
                            "prediction.rds", sep = "_")

  file_last <- file.path(ORFik_dir, file_suffix_last)
  this_org_is_done <- file.exists(file_last)
  if (this_org_is_done & !remake) {
    message("done (", org_short, ")")
    message("-------------------------------")
    next
  }
  dir.create(ORFik_dir, recursive = TRUE, showWarnings = FALSE)
  # Load ranges
  if (org_short == "Homo_sapiens") {
    gtf <- "~/livemount/shared_data/Homo_sapiens_mane/Homo_sapiens_GRCh38_101_subset300_mane.gtf"
    txdb <- paste0(gtf, ".db")
  } else {
    txdb <- loadTxdb(df)
    gtf <- ORFik:::getGtfPathFromTxdb(txdb)
  }

  faFile <- df@fafile
  seqinfos <- seqinfo(df)
  mrna <- loadRegion(txdb, "mrna")
  cds <- loadRegion(txdb, "cds")[names(mrna)]
  orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df, TRUE), longestORF = longestORF, startCodon = startCodon,
                                  stopCodon = stopCodon, minimumLength = minimumLength)
  orfs_gr_all <- ORFik:::categorize_and_filter_ORFs(orf_candidate_ranges,
                                                    ORF_categories_to_keep, cds, mrna)
  # Add lost CDS back
  cds_not_part_of_ranges <- ORFik:::removeMetaCols(cds[countOverlaps(cds, orfs_gr_all, type = "equal") == 0])
  mcols(cds_not_part_of_ranges) <- DataFrame(category = rep("annotated", length(cds_not_part_of_ranges)))
  orfs_gr_all <- c(orfs_gr_all, cds_not_part_of_ranges)
  print(table(mcols(orfs_gr_all)$category))


  # For each covrle
  format_grep <- "\\.qs$"
  covRLE_subsampling_dir <- file.path(covRLE_subsampling_dir_all, org_short)
  covrle_paths <- list.files(covRLE_subsampling_dir, format_grep, full.names = TRUE)
  length(covrle_paths)
  covrle_path <- covrle_paths[1]

  for (covrle_path in covrle_paths) {
    out_file_base <- paste0("ORFik_results", gsub("covRLE", "", gsub(format_grep, "", basename(covrle_path))))
    message(out_file_base)
    out_file <- file.path(ORFik_dir, paste(out_file_base, name(df), ORFik:::name_decider(df, naming = "full"),
                                           "prediction.rds", sep = "_"))

    if (file.exists(out_file) & !remake) next
    start_time <- Sys.time()
    P_sites_all_cov <- qs::qread(covrle_path, nthreads = 5)
    seqinfo(P_sites_all_cov@forward) <- seqinfos
    seqinfo(P_sites_all_cov@reverse) <- seqinfos

    ORFik::detect_ribo_orfs(df, out_folder = ORFik_dir, prefix_result = out_file_base,
                            mrna = mrna, cds = cds,
                            libraries = list(P_sites_all_cov),
                            orf_candidate_ranges = orf_candidate_ranges,
                            orfs_gr = orfs_gr_all,
                            ORF_categories_to_keep = ORF_categories_to_keep)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }
}

