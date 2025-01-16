
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/predictor_functions.R")

## First time install -> NOTE: Your need to replace code there, with my custom RiboCode code if you use different docker
#system(paste(conda, "create -n 'ribocode_env' python=3.7 ipython"))
#system(paste(conda, "install -n 'ribocode_env' -c bioconda ribocode"))
conda <- "~/miniconda3/bin/conda"
conda_works <- system(paste(conda, "info"), intern = TRUE)
conda_works <- is.null(attr(conda_works, "status"))
if (!conda_works) stop("Conda is was not found in path specified!")
call_prefix <- paste(conda, 'run -n ribocode_env')




message("Starting classifiers:")
message(Sys.time())

# Main paths
result_dir_all <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"
covRLE_subsampling_dir_all <- "~/livemount/shared_results/ribocrypt_paper//rds_temp_objects/covRLE_whole_library_subsamplings"

remake <- TRUE
remake_hd5 <- FALSE
organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
org <- organisms[1]
for (org in organisms) {
  message(org)
  message("-------------------------------")
  df <- read.experiment(gsub("all_samples", "all_merged", org), validate = F)
  dfc <- read.experiment(org, validate = F)
  org_short <- gsub("all_samples-", "", org)

  ## Create dirs
  result_dir <- file.path(result_dir_all, org_short)
  riboc_folder <- file.path(result_dir, "RiboCode/")

  file_suffix_last <- paste("RiboCode_results", "samp", nrow(dfc),"num", 1, sep = "_")
  file_suffix_last <- file.path(file_suffix_last, "results_collapsed.txt")

  file_last <- file.path(riboc_folder, file_suffix_last)
  this_org_is_done <- file.exists(file_last)
  if (this_org_is_done & !remake) {
    message("done (", org_short, ")")
    message("-------------------------------")
    next
  }
  dir.create(riboc_folder, recursive = TRUE, showWarnings = FALSE)
  # Load ranges
  if (org_short == "Homo_sapiens") {
    gtf <- "~/livemount/shared_data/Homo_sapiens_mane/Homo_sapiens_GRCh38_101_subset300_mane.gtf"
    txdb <- paste0(gtf, ".db")
  } else {
    gtf <- "~/livemount/shared_data/Saccharomyces_cerevisiae_mane/Saccharomyces_cerevisiae.R64-1-1.108_ensembl_mane.gtf"
    txdb <- paste0(gtf, ".db")
  }

  faFile <- df@fafile
  seqinfos <- seqinfo(df)
  tx <- loadRegion(txdb, "tx")
  config_path <- RiboCode_make_annotation_and_toy_config(riboc_folder, gtf, faFile)



  # For each covrle
  format_grep <- "\\.qs$"
  covRLE_subsampling_dir <- file.path(covRLE_subsampling_dir_all, org_short)
  covrle_paths <- list.files(covRLE_subsampling_dir, format_grep, full.names = TRUE)
  length(covrle_paths)
  covrle_path <- covrle_paths[1]

  for (covrle_path in covrle_paths) {
    out_file_base <- paste0("RiboCode_results", gsub("covRLE", "", gsub(format_grep, "", basename(covrle_path))))
    message(out_file_base)
    out_dir <- file.path(riboc_folder, out_file_base)
    out_prefix <- file.path(out_dir, "results")
    result_table_path <- paste0(out_prefix, "_collapsed.txt")
    if (file.exists(result_table_path) & !remake) next
    start_time <- Sys.time()
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


    message("- Loading dataset..")
    P_sites_all_cov <- qs::qread(covrle_path, nthreads = 5)
    seqinfo(P_sites_all_cov@forward) <- seqinfos
    seqinfo(P_sites_all_cov@reverse) <- seqinfos
    run_RiboCode_fast(tx, P_sites_all_cov, out_dir,
                      out_prefix, config_path,
                      conda, call_prefix,
                      verbose = FALSE, remake_hd5 = FALSE)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }
}

