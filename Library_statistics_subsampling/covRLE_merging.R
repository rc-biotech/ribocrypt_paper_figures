library(ORFik)
library(BiocParallel)
library(Rcpp)
library(qs)
Rcpp::sourceCpp("~/livemount/shared_scripts/ribocrypt_paper_figures/sum_rle_lists_single_thread.cpp")

Reduce_Rle_fast_no_comp <- function(cov) {
  stopifnot(is(cov, "list"))
  if (length(cov) == 0) return(list())
  stopifnot(is(cov[[1]], "RleList"))
  if (length(unique(lengths(cov))) != 1) stop("All RleLists must have equal lengths")

  values <- lapply(cov, function(x) runValue(x))
  lengths <- lapply(cov, function(x) runLength(x))
  chr_names <- names(cov[[1]])

  res <- sumMultipleRleListsOptimized2(values, lengths)

  res_RleList <- lapply(res, function(x) Rle(x$values, x$lengths))
  res_RleList <- as(res_RleList, "List")
  names(res_RleList) <- chr_names
  return(res_RleList)
}

Reduce_Rle_fast <- function(cov) {
  Rcpp::sourceCpp("~/livemount/shared_scripts/ribocrypt_paper_figures/sum_rle_lists_single_thread.cpp")
  stopifnot(is(cov, "list"))
  if (length(cov) == 0) return(list())
  stopifnot(is(cov[[1]], "RleList"))
  if (length(unique(lengths(cov))) != 1) stop("All RleLists must have equal lengths")

  values <- lapply(cov, function(x) runValue(x))
  lengths <- lapply(cov, function(x) runLength(x))
  chr_names <- names(cov[[1]])

  res <- sumMultipleRleListsOptimized2(values, lengths)

  res_RleList <- lapply(res, function(x) Rle(x$values, x$lengths))
  res_RleList <- as(res_RleList, "List")
  names(res_RleList) <- chr_names
  return(res_RleList)
}

send_to_discord <- !interactive()
if (send_to_discord) {
  message("Sending results to Discord")
  discordr::set_default_discord_connection(discordr::import_discord_connections("~/livemount/.cache/discordr/r_to_discord_config")[[1]])
}

plot_dir <- "~/livemount/shared_results/ribocrypt_paper/"
rds_dir <- file.path(plot_dir, "rds_temp_objects")
res_dir_all <- file.path(rds_dir, "covRLE_whole_library_subsamplings")
res_sampling_info_dir_all <- file.path(res_dir_all, "subsamplings_info")
dir.create(res_sampling_info_dir_all, showWarnings = FALSE, recursive = TRUE)
if (send_to_discord) {
  discordr::send_webhook_message("Starting covRLE pipeline..")
}
samplings_all <- c(1, 3, 6, 10, 30, 100, 500, 1000, 1500, 2000, 2500, 2783)
num_samplings <- seq(20)

organisms <- c("all_samples-Homo_sapiens", "all_samples-Saccharomyces_cerevisiae")
org <- organisms[1]
for (org in organisms) {
  message(org)
  message("-------------------------------")
  dfc <- read.experiment(org, validate = F)
  org_short <- gsub("all_samples-", "", org)
  res_dir <- file.path(res_dir_all, org_short)
  res_sampling_info_dir <- file.path(res_sampling_info_dir_all, org_short)
  file_suffix_last <- paste0("samp_", nrow(dfc),"_num_", 1, ".qs")
  file_last <- file.path(res_dir, paste0("covRLE_", file_suffix_last))
  this_org_is_done <- file.exists(file_last)
  if (this_org_is_done) {
    message("done (", organism(dfc), ")")
    next
  }

  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res_sampling_info_dir, recursive = TRUE, showWarnings = FALSE)

  all_paths <- filepath(dfc, "cov", base_folders = libFolder(dfc, mode = "all"), suffix_stem = c("_pshifted", ""))
  first_start_time <- Sys.time()


  all_paths_qs_object <- file.path(res_dir_all, paste0("all_covrles_", org_short, ".qs"))
  if (file.exists(all_paths_qs_object)) {
    message("Loading all libraries as covRLE.qs")
    all_covs <- qs::qread(all_paths_qs_object, nthreads = 5)
  } else {
    message("Creating all libraries as covRLE.qs")
    all_covs <- lapply(all_paths, function(x) readRDS(x))
    qs::qsave(all_covs, all_paths_qs_object, nthreads = 5)
  }
  samplings <- unique(c(samplings_all[samplings_all < nrow(dfc)], nrow(dfc)))
  large_genome <- sum(seqlengths(seqinfo(dfc))) > 309975071 # (1/10 of human)
  bp <- BiocParallel::MulticoreParam(exportglobals = FALSE)
  sp <- BiocParallel::SnowParam(workers = 2, exportvariables = TRUE, exportglobals = FALSE)
  res <- try({
    for (i in samplings) {
      message(i)
      for (n in num_samplings) {
        message("- ", n)
        file_suffix <- paste0("samp_", i,"_num_", n, c(".rds", ".qs"))
        file <- file.path(res_dir, paste0("covRLE_", file_suffix))
        names(file_suffix) <- names(file) <- c("rds", "qs")

        if (any(file.exists(file))) {
          if (file.exists(file["rds"]) & !file.exists(file["qs"])) {
            message("--- converting rds to qs")
            qs::qsave(readRDS(file["rds"]), file["qs"], nthreads = 5)
            file.remove(file["rds"])
          }
          if (i == max(samplings)) {break} else next
        }


        samp <- sample(x = seq(nrow(dfc)), size = i, replace = FALSE)
        saveRDS(samp, file.path(res_sampling_info_dir, file_suffix["rds"]))
        # dfc_sub <- dfc[samp,]
        start_time <- Sys.time()
        message("-- Loading libs...")
        if (exists("all_covs") & length(all_covs) == length(all_paths)) {
          covs <- all_covs[samp]
        } else {
          covs <- lapply(all_paths[samp], function(x) fimport(x))
        }
        covs_split <- list(lapply(covs, function(x) f(x)),
                           lapply(covs, function(x) r(x)))
        rm(covs)

        message("--- Computing sum of covRles..")

        more_than_6_libs_do_parallel <- i > 6 & large_genome
        if (more_than_6_libs_do_parallel) {
          len <- i
          s <- ifelse(len > 30, 5, 3); splits <- ceiling(seq(len) / (len/s))
          sps <- BiocParallel::SnowParam(workers = s, exportvariables = TRUE, exportglobals = FALSE)
          if (s == 1) {
            cov_sum <- bplapply(covs_split, function(cov_x, Reduce_Rle_fast) Reduce_Rle_fast(cov_x),
                                BPPARAM = sp, Reduce_Rle_fast = Reduce_Rle_fast)
          } else {
            cov_sum <- bplapply(covs_split, function(cov_x, bp, splits) Reduce_Rle_fast(bplapply(split(cov_x, splits), function(x) Reduce_Rle_fast(x), BPPARAM = bp)), BPPARAM = sp, bp = bp, splits = splits)
          }
        } else {
          cov_sum <- lapply(covs_split, function(cov_x) Reduce_Rle_fast_no_comp(cov_x))
        }
        message("---- Saving object..")
        qs::qsave(covRle(cov_sum[[1]], cov_sum[[2]]),
                file["qs"], nthreads = 5)

        time_dif <- Sys.time() - start_time
        time_dif_nice <- paste("--", round(time_dif, 2), attr(time_dif, "units"))
        if (send_to_discord) {
          discordr::send_webhook_message(paste0("covRLE pipeline (", organism(dfc),"): size (", i, ")", " sampling (", n, "), done.. \n-time used: ", time_dif_nice))
        }
        print(time_dif_nice)
        if (i == max(samplings)) break
      }
    }
  })

  if (is(res, "try-error")) {
    time_dif <- Sys.time() - first_start_time
    time_dif_nice <- paste("--", round(time_dif, 2), attr(time_dif, "units"))
    if (!exists("i")) i <- 1
    if (!exists("n")) n <- 1
    res <- system("free -h", intern = TRUE)
    mem_info_clean <- strsplit(gsub(" +", ",", res[2]), ",")[[1]]
    free_memory <- paste("\n-free memory:", mem_info_clean[4], "/", mem_info_clean[2])
    error <- paste0("Error, covRLE pipeline (", organism(dfc),") failed at: size (", i, ")", " sampling (", n,
                    ") \n-time used total: ", time_dif_nice,
                    free_memory)
    discordr::send_webhook_message(error)
    stop(cat(error))
  }
  message("done (", organism(dfc), ")")
  discordr::send_webhook_message("done (", organism(dfc), ")")
}
message("done")


