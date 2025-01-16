## Run all classifiers (libraries are preload upstream)
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
start_time <- Sys.time()
message(start_time)


# Main paths
result_dir_all <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/"

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
result_dir_all <- ORFik::config()["ref"]

file_name <- "predicted_translons_with_sequence"
file_names <- paste0(file_name, c(".fst", ".csv", "_ranges.rds"))
names(file_names) <- c("fst", "csv", "ranges")
all_exp <- list.experiments(validate = FALSE)
all_merged_exps <- all_exp[grep("all_merged", name),][libtypes == "RFP",]
stopifnot(!any(duplicated(all_merged_exps$name)))
organisms <- all_merged_exps$name
organisms <- organisms[!(organisms %in% c("all_merged-Mus_musculus"))]
#organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
org <- organisms[1]
remake <- FALSE
for (org in organisms) {
  message("-------------------------------")
  message(org)
  message("-------------------------------")
  df <- read.experiment(gsub("all_samples", "all_merged", org), validate = F)
  # dfc <- read.experiment(org, validate = F)
  org_short <- gsub("all_merged-", "", org)
  result_dir <- file.path(result_dir_all, tolower(org_short), "predicted_translons")
  ORFik_dir <- file.path(result_dir, "ORFik/")
  riboc_folder <- file.path(result_dir, "RiboCode/")
  ORFquant_dir <- file.path(result_dir,  "ORFquant/")

  libnames <- ORFik:::name_decider(df, naming = "full")
  out_prefix_ORFik <-"predicted_translons"
  out_prefix_ORFik_full <- file.path(ORFik_dir, paste0(out_prefix_ORFik, "_", name(df), "_", libnames))
  out_file_ORFik <- paste0(out_prefix_ORFik_full, "_prediction_table.rds")
  out_file_orfquant <- file.path(ORFquant_dir, paste0("ORFquant_predicted_translons_", org_short, ".fst"))
  out_prefix_RiboCode <- file.path(riboc_folder, paste0("RiboCode_predicted_translons_", org_short))
  out_file_RiboCode <- paste0(out_prefix_RiboCode, "_collapsed.txt")

  file_last <- file.path(out_file_RiboCode)
  this_org_is_done <- file.exists(file_last)
  if (this_org_is_done & !remake) {
    message("done (", org_short, ")")
    message("-------------------------------")
    next
  }
  dir.create(ORFik_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(riboc_folder, recursive = TRUE, showWarnings = FALSE)
  dir.create(ORFquant_dir, recursive = TRUE, showWarnings = FALSE)
  # Load ranges
  gtf_path <- ORFik:::getGtfPathFromTxdb(loadTxdb(df))
  file_ext <- tools::file_ext(gtf_path)
  gtf_out_path <- file.path("~/livemount/shared_data/pseudo_5utrs_gtfs", paste0(org_short, "_mane"),
                            paste0(ORFik:::remove.file_ext(gtf_path, TRUE), "_mane",".", "gtf"))
  if (file_ext %in% c("gff")) {
    gtf_out_path <- gtf_path
  } else {
    if (!file.exists(gtf_out_path)) {
      txdb <- loadTxdb(df)
      leaders <- loadRegion(txdb, "leaders")
      cds <- loadRegion(txdb, "cds")
      if (length(cds) == 0)
        stop("Can not add pseudo 5' UTRs to a genome without Coding sequences!")
      percentage <- round((length(leaders) / length(cds))*100, 1)
      if (percentage < 30) {
        message("-- Adding pseudo utrs")
        add_pseudo_leaders_gtf(gtf_path, gtf_out_path, extension = 51, df = df, maximum_tx_length = 18000)
      } else {
        add_pseudo_leaders_gtf(gtf_path, gtf_out_path, extension = 0, df = df, maximum_tx_length = 18000)
      }
    }
  }


  txdb_path <- paste0(gtf_out_path, ".db")
  txdb <- loadTxdb(txdb_path)


  faFile <- df@fafile
  seqinfos <- seqinfo(df)
  ORFik_and_ORFquant_not_done <- !file.exists(out_file_ORFik) | !file.exists(out_file_orfquant) | remake
  if (ORFik_and_ORFquant_not_done) {
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
  }


  P_sites_all_cov <- readRDS(filepath(df, "cov"))

  if (!file.exists(out_file_ORFik) | remake) { # ORFik
    message("-- ORFik")
    start_time <- Sys.time()
    ORFik::detect_ribo_orfs(df, out_folder = ORFik_dir, prefix_result = out_prefix_ORFik,
                            mrna = mrna, cds = cds,
                            libraries = list(P_sites_all_cov),
                            orf_candidate_ranges = orf_candidate_ranges,
                            orfs_gr = orfs_gr_all,
                            ORF_categories_to_keep = ORF_categories_to_keep)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }

  if (!file.exists(out_file_orfquant) | remake) { # ORFquant
    message("-- ORFquant")
    start_time <- Sys.time()

    fst::write_fst(run_ORFquant_fast(orfs_gr_all, P_sites_all_cov, mrna), out_file_orfquant)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }

  if (!file.exists(out_file_RiboCode) | remake) { # ORFquant
    message("-- RiboCode")
    start_time <- Sys.time()
    config_path <- RiboCode_make_annotation_and_toy_config(riboc_folder, gtf_out_path, faFile)
    tx <- loadRegion(txdb, "tx")
    run_RiboCode_fast(tx, P_sites_all_cov, riboc_folder,
                      out_prefix_RiboCode, config_path,
                      conda, call_prefix,
                      verbose = TRUE, remake_hd5 = FALSE)
    time_dif <- Sys.time() - start_time
    message(paste("--", round(time_dif, 2), attr(time_dif, "units")))
  }
}
