## Run all classifiers (libraries are preload upstream)
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/predictor_functions.R")

message("Starting classifiers:")
start_time <- Sys.time()
message(start_time)


# Main paths
result_dir_all <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"
covRLE_subsampling_dir_all <- "~/livemount/shared_results/ribocrypt_paper//rds_temp_objects/covRLE_whole_library_subsamplings"

# Parameters
export_metrics_table <- FALSE
ORF_categories_to_keep <- c("uORF", "uoORF", "annotated") # "NTE", "NTT", "internal", "doORF", "dORF"
name_of_result <- "uORF"
longestORF = TRUE
startCodon = paste(startDefinition(1), "GTG",sep = "|")
stopCodon = stopDefinition(1)
minimumLength = 0



## Create dirs
dir.create(result_dir_all, recursive = TRUE, showWarnings = FALSE)

remake <- TRUE
organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
org <- organisms[1]
for (org in organisms) {
  message(org)
  message("-------------------------------")
  df <- read.experiment(org, validate = F)
  org_short <- gsub("all_samples-", "", org)
  result_dir <- file.path(result_dir_all, org_short)
  ORFquant_dir <- file.path(result_dir, "ORFquant/")
  file_suffix_last <- paste0("ORFquant_results_", "samp_", nrow(df),"_num_", 1, ".fst")
  file_last <- file.path(ORFquant_dir, file_suffix_last)
  this_org_is_done <- file.exists(file_last)
  if (this_org_is_done & !remake) {
    message("done (", org_short, ")")
    message("-------------------------------")
    next
  }
  dir.create(ORFquant_dir, recursive = TRUE, showWarnings = FALSE)
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
  mrna <- loadRegion(txdb, "mrna")# ["ENST00000372324"]
  cds <- loadRegion(txdb, "cds")[names(mrna)]
  orf_candidate_ranges = findORFs(seqs = txSeqsFromFa(mrna, df, TRUE), longestORF = longestORF, startCodon = startCodon,
                                  stopCodon = stopCodon, minimumLength = minimumLength)
  orfs_gr_all <- ORFik:::categorize_and_filter_ORFs(orf_candidate_ranges,
                                                    ORF_categories_to_keep, cds, mrna)
  # Add lost CDS back
  cds_not_part_of_ranges <- ORFik:::removeMetaCols(cds[countOverlaps(cds, orfs_gr_all, type = "equal") == 0])
  mcols(cds_not_part_of_ranges) <- DataFrame(category = rep("annotated", length(cds_not_part_of_ranges)))
  orfs_gr_all <- c(orfs_gr_all, cds_not_part_of_ranges)


  # For each covrle
  format_grep <- "\\.qs$"
  covRLE_subsampling_dir <- file.path(covRLE_subsampling_dir_all, org_short)
  covrle_paths <- list.files(covRLE_subsampling_dir, format_grep, full.names = TRUE)
  length(covrle_paths)
  covrle_path <- covrle_paths[1]

  for (covrle_path in covrle_paths) {
    out_file_base <- paste0("ORFquant_results", gsub("covRLE", "", gsub(format_grep, ".fst", basename(covrle_path))))
    message(out_file_base)
    out_file <- file.path(ORFquant_dir, out_file_base)

    if (file.exists(out_file) & !remake) next
    start_time <- Sys.time()
    P_sites_all_cov <- qs::qread(covrle_path, nthreads = 5)
    seqinfo(P_sites_all_cov@forward) <- seqinfos
    seqinfo(P_sites_all_cov@reverse) <- seqinfos

    fst::write_fst(run_ORFquant_fast(orfs_gr_all, P_sites_all_cov, mrna), out_file)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }
}

