library(ORFik)
library(ORFquant) # Needed for ORFquant:::take_Fvals_spect
library(data.table)
library(rhdf5)
library(hdf5r)

#' Create internal hdf5 file to bypass bam load
#'
#' Internally runs coveragePerTiling and creates a matrix, where each row is a gene
#' Since genes are not of equal size, the longest gene defines number of columns,
#' and all short genes are filled with 0's on the right side. Therefor if you cut
#' out extreme long genes, you will save a lot of space.
make_ribocode_hd5 <- function(tx, P_sites_all_cov, out_dir, remake_hd5 = TRUE) {
  message("- HDF5 file..")
  stopifnot(!anyNA(names(tx)))
  stopifnot(length(x) > 0)
  h5file <- file.path(out_dir, "RFP_1_psites.hd5")
  if (file.exists(h5file) & !remake_hd5) {
    message("Using existing HDF5 file")
    return(invisible(NULL))
  } else suppressWarnings(file.remove(h5file))


  message("-- making psite coverage for HDF5")
  psites <- IntegerList(coveragePerTiling(tx, P_sites_all_cov))
  names(psites) <- NULL
  psites_number <- sum(sum(psites))
  psites <- as(psites, "list")

  message("-- Writing HDF5 file")
  h5file_link <- H5File$new(h5file, mode = "w")
  vlen_int_type <- H5T_VLEN$new(h5types$H5T_NATIVE_INT) # Define the VLEN datatype for integers

  h5file_link[["transcript_ids"]] <- names(tx)
  h5file_link$create_dataset("p_sites", robj =  psites, dtype = vlen_int_type, gzip_level = 6, dims = length(psites))
  h5file_link$create_dataset("psites_number", robj =  psites_number, dtype = h5types$int64_t, gzip_level = 6, dims = 1)
  h5file_link$close_all()

  suppressWarnings(rm(psites))
  return(invisible(NULL))
}

#'  Run RiboCode prediction with hd5 speedup
#'
#' First time install -> NOTE: Your need to replace code there, with my custom RiboCode code if you use different docker
#' system(paste(conda, "create -n 'ribocode_env' python=3.7 ipython"))
#' system(paste(conda, "install -n 'ribocode_env' -c bioconda ribocode"))
#' @param tx GRangesList, All transcripts from loadRegion(txdb, "tx"),
#' used to calculate transcript coverage saved to the hdf5 file, ignored if
#' hdf5 file exists and remake_hd5 is FALSE.
#' @param P_sites_all_cov coverage (preferably as a covRle, since it is much faster)
#' @param out_dir character, output directory for all files.
#' @param
run_RiboCode_fast <- function(tx, P_sites_all_cov, out_dir,
                              out_prefix = file.path(out_dir, "results"),
                              config_path,
                              conda = "~/miniconda3/bin/conda",
                              call_prefix = paste(conda, 'run -n ribocode_env'),
                              alt_start_codons = c("CTG", "TTG", "GTG"),
                              verbose = TRUE,
                              remake_hd5 = TRUE,
                              longest_orf = "no",
                              threads = 2) {
  stopifnot(longest_orf %in% c("yes", "no"))
  make_ribocode_hd5(tx, P_sites_all_cov, out_dir, remake_hd5)
  # Predict
  message("-- Prediction..")
  detect_step <- paste("RiboCode",
                       "-a", out_dir,
                       "-c",  config_path,
                       paste("-l", longest_orf),
                       paste("-A", paste(alt_start_codons, collapse = ",")),
                       paste("-m", threads),
                       "-o", out_prefix)
  Ribo_code_predict <- paste("cd", out_dir, "&&" , call_prefix, detect_step)
  res <- system(Ribo_code_predict, ignore.stdout = !verbose, ignore.stderr = !verbose)
  stopifnot(res == 0)
  return(invisible(res))
}

RiboCode_make_annotation_and_toy_config <- function(riboc_folder, gtf, faFile,
                                                    conda = "~/miniconda3/bin/conda",
                                                    call_prefix = paste(conda, 'run -n ribocode_env')) {
  # First time init
  config_path <- file.path(riboc_folder, "config.txt")
  init_not_done_before <- !file.exists(config_path)
  if (init_not_done_before) {
    message("RiboCode --")

    prepare_step <- paste("prepare_transcripts",
                          "-g", gtf,
                          "-f", faFile,
                          "-o", riboc_folder)
    res <- system(paste(call_prefix, prepare_step))
    stopifnot(res == 0)
    toy_bam_file <- paste0(tempdir(), "test.bam")
    saveRDS(toy_bam_file, toy_bam_file)
    config.RC <- data.table(SampleName = "RFP_1",
                            AlignmentFile = toy_bam_file, # This path is not used, we use hdf5 file
                            Stranded = "yes",
                            plen = "27",
                            psite =  "0")
    colnames(config.RC) <- NULL
    fwrite(config.RC, config_path, sep = "\t")
    message("-- Annotation init done")
  }
  return(config_path)
}

#' Multitaper Fourier transform, with 24 tapers, with resulting F test for significance.
run_ORFquant_fast <- function(orfs_gr_all, P_sites_all_cov, mrna,
                              minimum_coverage = 10, active_psites_min = 4,
                              frame0_ratio_min = 0.5, tapers=24, bw=12) {
  orfs_gr <- orfs_gr_all[countOverlaps(orfs_gr_all, P_sites_all_cov) > minimum_coverage]
  length(orfs_gr_all); length(orfs_gr)
  # sum(psite) > 0 (4 from previous)
  # Positions with psites > 2 (4 from previous)
  # In frame reads > 50%
  # pval < 0.05
  cov_dt <- coveragePerTiling(orfs_gr, P_sites_all_cov, as.data.table = TRUE, withFrames = TRUE)
  rm(P_sites_all_cov)
  cov_dt_frame <- cov_dt[, (frame_sum_gene = sum(count)), by = .(genes, frame)]
  cov_dt_agg <- cov_dt[, .(sum_gene = sum(count)), by = genes]
  cov_dt_agg[, `:=`(frame0 = cov_dt_frame[frame == 0]$V1, frame1 = cov_dt_frame[frame == 1]$V1, frame2 = cov_dt_frame[frame == 2]$V1)]
  cov_dt_agg[, `:=`(frame0_ratio = frame0 / sum_gene, frame1_ratio = frame1 / sum_gene, frame2_ratio = frame2 / sum_gene)]
  cov_dt_agg[, `:=`(frame0_ratio = round(frame0_ratio, 2), frame1_ratio = round(frame1_ratio, 2), frame2_ratio = round(frame2_ratio, 2))]
  cov_dt_agg <- setDT(cbind(cov_dt_agg))
  cov_dt_agg[, active_psites := cov_dt[, sum(count > 0), by = genes]$V1]
  cov_dt_agg[is.na(cov_dt_agg)] <- 0
  cov_dt_agg[, count_filter := FALSE]
  cov_dt_agg[sum_gene > 4 & active_psites > active_psites_min & frame0_ratio > frame0_ratio_min, count_filter := TRUE]
  table(cov_dt_agg$count_filter)
  cov_dt_agg_passed <- cov_dt_agg[count_filter == TRUE,]

  cov_dt_res <- cov_dt[genes %in% cov_dt_agg_passed$genes,]
  cov_dt_res[, psit := count]
  length(unique(cov_dt_res$genes))

  orfs_gr_passed <- orfs_gr[unique(cov_dt_agg_passed$genes)]
  ir_passed <- ranges(pmapToTranscriptF(orfs_gr_passed, mrna[names(orfs_gr_passed)]))
  cov_dt_agg_passed <- cbind(start = as.integer(start(ir_passed)), end = as.integer(end(ir_passed)), cov_dt_agg_passed, tx_id = as.factor(names(orfs_gr_passed)))
  cov_dt_agg_passed[, stop_group := paste0(cov_dt_agg_passed$end, "_", cov_dt_agg_passed$tx_id)]
  cov_dt_agg_passed[, stop_group_duplicate := duplicated(stop_group)]
  length(unique(cov_dt_agg_passed$genes))
  sum(!cov_dt_agg_passed$stop_group_duplicate)

  cov_dt_agg_passed_stop_unique <- cov_dt_agg_passed[stop_group_duplicate == FALSE,]
  message("- Candidate ORFs: ", nrow(cov_dt_agg_passed_stop_unique))
  cov_dt_agg_passed_stop_unique[, pval := cov_dt_res[genes %in% cov_dt_agg_passed_stop_unique$genes,][
    , pf(q=ORFquant:::take_Fvals_spect(x = psit,n_tapers = tapers,time_bw = bw,
                                       slepians_values = multitaper::dpss(n= if(length(psit)<25){length(psit)+(50-length(psit))} else {length(psit)},k=tapers,nw=bw))[1],
         df1=2,df2=(2*24)-2,lower.tail=F), by = genes]$V1]
  cov_dt_agg_passed_stop_unique[, significant := pval < 0.05]
  message("- Significant ORFs: ", sum(cov_dt_agg_passed_stop_unique$significant))
  cov_dt_agg_passed_stop_unique[, ORF_type := mcols(orfs_gr[cov_dt_agg_passed_stop_unique$genes])$category]
}

#' Add pseudo leaders
#' RiboCode needs new gtf if no leaders are present
#' @param gtf_path the path to original gtf
#' @param gtf_out_path the path to new gtf, will auto make new directory.
#' @param df ORFik experiment, default NULL. If defined will make optimized txdb
#' @param maximum_tx_length numeric, default NULL.
#'  The maximum length of gene to keep (default is keep all)
#' @return path to new gtf
add_pseudo_leaders_gtf <- function(gtf_path, gtf_out_path, extension = 51,
                                   df = NULL, maximum_tx_length = NULL) {
  file_ext <- tools::file_ext(gtf_path)
  stopifnot(file_ext %in% c("gff3", "gtf"))
  gtf <- rtracklayer::import(gtf_path)
  if (length(gtf) == 0) {
    stop("GTF is of size 0, check that this is correct!")
  }
  stopifnot(ncol(mcols(gtf)) > 0)
  if (file_ext == "gff3") {
    message("Converting gff3 to gtf format")
    gtf <- convert_gff3_to_gtf(gtf)
  }


  stopifnot(all(c("gene","transcript", "exon", "CDS") %in% gtf$type))
  gtf_original <- gtf
  gtf_original_length <- length(gtf_original)

  removed <- removed_gtf <- 0
  filter_on_maxlength <- !is.null(maximum_tx_length) && is.numeric(maximum_tx_length)
  if (filter_on_maxlength) {
    tx <- gtf[gtf$type == "transcript"]
    exons <- gtf[gtf$type == "exon"]
    exons_dt <- data.table(width = width(exons), tx = exons$transcript_id)
    tx_dt <- exons_dt[, .(tx_width = sum(width)), by = tx]
    stopifnot(length(tx_dt$tx) == length(tx))
    stopifnot(all(tx_dt$tx == tx$transcript_id))

    original_num_tx <- length(tx)
    tx <- tx[tx_dt$tx_width <= maximum_tx_length]
    removed <- original_num_tx - length(tx)
    message("Transcripts removed over max length: ", removed, " (", round(100*removed/original_num_tx, 2), "%)")
    gtf <- gtf[gtf$gene_id %in% tx$gene_id]
    removed_gtf <- gtf_original_length - length(gtf)
  }
  if (extension != 0) {
    gtf$index <- seq(length(gtf))
    cds_transcripts <- unique(gtf[gtf$type == "CDS"]$transcript_id)
    stopifnot("Can not add pseudo 5' UTRs to a genome without Coding sequences!" =
                length(cds_transcripts) != 0)
    tx <- gtf[gtf$transcript_id %in% cds_transcripts & gtf$type == "transcript"]
    gene <- gtf[gtf$gene_id %in% tx$gene_id & gtf$type == "gene"]
    exons <- gtf[gtf$transcript_id %in% cds_transcripts & gtf$type == "exon"]
    exons_grl <- split(exons, exons$transcript_id)
    exons_grl <- exons_grl[tx$transcript_id]
    exons_grl_extend <- extendLeaders(exons_grl, extension = extension)
    exons_extend <- unlistGrl(exons_grl_extend)
    names(exons_extend) <- NULL
    new_starts <- startSites(exons_grl_extend, is.sorted = TRUE, keep.names = TRUE)


    max(table(tx$gene_id))
    dt <- data.table(tx_id = tx$transcript_id, gene_id = tx$gene_id,
                     start = as.integer(end(tx)), end = as.integer(end(tx)),
                     strand = strandBool(tx), selected = FALSE)
    dt[strand == TRUE, selected := start == max(start), by = gene_id]
    dt[strand == FALSE, selected := end == min(end), by = gene_id]

    dt <- dt[selected == TRUE,] # Longest upstream isoform
    dt <- dt[!duplicated(gene_id), ] # Only first if multiple match
    stopifnot(nrow(dt) == length(gene$gene_id))
    stopifnot(all(dt$gene_id == gene$gene_id))

    strand_bool_tx <- strandBool(tx)
    strand_bool_gene <- strandBool(gene)
    start(tx[strand_bool_tx]) <- new_starts[tx[strand_bool_tx]$transcript_id]
    end(tx[!strand_bool_tx]) <- new_starts[tx[!strand_bool_tx]$transcript_id]
    start(gene[strand_bool_gene]) <- new_starts[dt[strand == TRUE,]$tx_id]
    end(gene[!strand_bool_gene]) <- new_starts[dt[strand == FALSE,]$tx_id]

    # Append back

    gtf[exons_extend$index] <- exons_extend
    gtf[tx$index] <- tx
    gtf[gene$index] <- gene
    gtf$index <- NULL
  }

  stopifnot(length(gtf) == (gtf_original_length - removed_gtf))
  stopifnot(ncol(mcols(gtf)) == ncol(mcols(gtf_original)))
  stopifnot(nrow(mcols(gtf)) == (nrow(mcols(gtf_original)) - removed_gtf))

  dir.create(dirname(gtf_out_path), showWarnings = FALSE, recursive = TRUE)
  rtracklayer::export(gtf, gtf_out_path, format = "gtf")
  make_txdb <- !is.null(df) & is(df, "experiment")
  if (make_txdb) {
    makeTxdbFromGenome(gtf_out_path, df@fafile, organism(df), optimize = TRUE)
    txdb <- loadTxdb(paste0(gtf_out_path, ".db"))
    message("Widths of new leaders:")
    print(head(table(widthPerGroup(loadRegion(txdb, "leaders"))), 10))
    message("Widths of original leaders:")
    print(head(table(widthPerGroup(loadRegion(df, "leaders"))), 10))
  }
  return(gtf_out_path)
}

convert_gff3_to_gtf <- function(gtf) {
  if (!("chromosome" %in% unique(gtf$type))) stop("gff3 format must have type chromosome!")
  stopifnot(c("type", "Parent", "Name") %in% colnames(mcols(gtf)))

  gtf[gtf$type %in% c("chromosome", "biological_region")] <- NULL
  gtf_original_length <- length(gtf)
  stopifnot(!any(strand(gtf) == "*"))
  gtf$description <- NULL
  dt <- as.data.table(mcols(gtf))
  dt[, type := as.character(type)]
  dt[type %in% grep("gene$", unique(type), value = TRUE), type := "gene"]

  to_become_tx <- c("mRNA", "pseudogenic_transcript", "rRNA", "ncRNA", "tRNA", "snRNA", "snoRNA","RNase_P_RNA", "SRP_RNA", "RNase_MRP_RNA")
  dt[type %in% to_become_tx, type := "transcript"]

  dt[,`:=`(type = as.factor(type), tag = as.character(NA),
           gene_biotype = biotype, transcript_biotype = biotype,
           exon_number = rank)]

  dt[, id := cumsum(dt$type == "gene")]
  dt[, id_tx := cumsum(dt$type == "transcript")]
  dt[1, id_tx := 1]
  dt[, id_gene_index := which(type == "gene")[id]]
  dt[type == "gene", transcript_biotype := as.character(NA)]
  dt[, `:=`(gene_source = source, transcript_source = source)]
  dt[type == "gene", transcript_source := as.character(NA)]
  dt[grep("gene:", Parent), gene_id := gsub("gene:", "", grep("gene:", Parent, value = TRUE))]
  dt[grep("transcript:", Parent), transcript_id := gsub("transcript:", "", grep("transcript:", Parent, value = TRUE))]
  dt[, gene_name := Name[id_gene_index]]
  dt[grep(":exon:[0-9]+$", gene_name), gene_name := as.character(NA)]
  gene_ids <- unique(dt$gene_id)
  gene_ids <- gene_ids[!is.na(gene_ids)]
  tx_ids <- unique(dt$transcript_id)
  tx_ids <- tx_ids[!is.na(tx_ids)]
  dt[, gene_id := gene_ids[id]]
  dt[, transcript_id := tx_ids[id_tx]]
  dt[, gene_biotype := gene_biotype[id_gene_index]]
  dt[, transcript_biotype := gene_biotype]
  dt[type == "gene", transcript_biotype := as.character(NA)]
  dt[type == "gene", transcript_id := as.character(NA)]
  dt[, transcript_name := gene_name]

  # delete gff3 unique columns
  cols_to_remove <- c("ID", "Alias", "biotype", "description", "logic_name", "Parent",
                      "Name", "constitutive", "ensembl_end_phase", "ensembl_phase", "rank", "external_name", "Is_circular")
  cols_to_remove <- cols_to_remove[cols_to_remove %in% colnames(dt)]
  dt[, cols_to_remove] <- NULL
  # Reorder to gtf format
  gtf_col_order <-  c("source", "type", "score", "phase", "gene_id", "gene_name", "gene_source",
                      "gene_biotype", "transcript_id", "transcript_name", "transcript_source",
                      "transcript_biotype", "tag", "exon_number", "exon_id", "protein_id")
  gtf_col_order <- gtf_col_order[gtf_col_order %in% colnames(dt)]
  dt <- dt[, ..gtf_col_order]
  stopifnot(nrow(dt) == length(gtf))
  mcols(gtf) <- DataFrame(dt)
  stopifnot(length(gtf) == gtf_original_length)
  no_duplicated_tx_ids_for_transcripts <- !any(duplicated(gtf[gtf$type == "transcript"]$transcript_id))
  stopifnot(no_duplicated_tx_ids_for_transcripts)
  no_duplicated_gene_ids_for_genes <- !any(duplicated(gtf[gtf$type == "gene"]$gene_id))
  stopifnot(no_duplicated_gene_ids_for_genes)
  return(gtf)
}
