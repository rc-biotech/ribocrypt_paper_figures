library(ORFik)
library(data.table)
library(ggplot2)
library(plotly)

df <- read.experiment("all_samples-Homo_sapiens_04oct_2024", validate = FALSE)
nrow(df)
paths <- filepath(df, type = "cov", base_folders = libFolder(df, mode = "all"), suffix_stem = c("", "_pshifted"))
symbols <- symbols(df)
symbols <- symbols[!duplicated(ensembl_gene_id)]
cds <- loadRegion(df, "cds")
mrna <- loadRegion(df, "mrna")
countTable <- countTable(df, type = "summarized")
faFile <- FaFile(df@fafile)
ids <- names(cds)[names(cds) %in% symbols$ensembl_tx_name]
cds <- cds[ids]
mrna <- mrna[ids]
countTable <- countTable[ids]
dim(countTable)
srr <- runIDs(df)

dt_all <- lapply(seq_along(paths), function(i) {
  message(i)
  dt <- try(codon_usage(fimport(paths[i]), cds = cds, mrna = mrna, filter_table = assay(countTable[,i]), faFile = faFile))
  if (is(dt, "try-error")) return(data.table())
  dt[, variable := srr[i]]
  return(dt)
})
dt_all_dt <- rbindlist(dt_all, fill = TRUE)
fst::write_fst(dt_all_dt, "~/livemount/shared_results/ribocrypt_paper/codon_usage_homo_sapiens_all_libs.fst")
fwrite(dt_all_dt, "~/livemount/shared_results/ribocrypt_paper/codon_usage_homo_sapiens_all_libs.csv")
dt_all_dt <- fst::read_fst("~/livemount/shared_results/ribocrypt_paper/codon_usage_homo_sapiens_all_libs.fst")
dt <- copy(dt_all_dt)
length(unique(dt$seqs))
dt <- dt[seqs %in% dt_all_dt[, sum(N.total), by = seqs][V1 > 15e3,]$seqs, ]
length(unique(dt$seqs))
levels(dt$seqs) <- sort(levels(dt$seqs))
plot <- ggplot(dt, aes(x = seqs, y = relative_to_max_score)) + geom_boxplot(outliers = FALSE) + facet_wrap(~ type, scales = "free") + coord_flip() + theme_bw()
ggsave(plot, filename = "~/livemount/shared_results/ribocrypt_paper/codon_usage_homo_sapiens_all_libs.png",
       widt = 6, height = 8, dpi = 300)
p <- plotly_build(plot)

for(i in 1:length(p$x$data)) {
  p$x$data[[i]]$marker$opacity = 0
}
p

