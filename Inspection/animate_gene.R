animate_gene <- function(gene, df_all, collection_genes, subregion = NA, subset_sizes = c(1, 3, 6, 10, 30, 100, 500, 1000,
                                                                                          1500, 2000, 2500), normalize = TRUE) {
  gene_path_fst <- load_allsamples_fst(gene, df_all, collection_genes)
  table <- RiboCrypt:::compute_collection_table(gene_path_fst, NULL, df_all, "TISSUE",
                                                          normalization = "tpm", value.var = "count",
                                                          kmer = 1, metadata = m, min_count = 0, as_list = TRUE)$table
  # if (!is.na(subregion)[1]) {
  #   table <- table[subregion,]
  # }
  all_subsets <- rbindlist(sample_libs_from_fst(table, subset_sizes))
  all_subsets[,frame := as.factor((position-1) %% 3)]
  if(normalize) {
     all_subsets[, counts := tx_norm_counts]
   }
  if (!is.na(subregion)[1])  all_subsets <- all_subsets[position %in% subregion]
  all_subsets[counts > 5000, counts := 5000]
  plot <- ggplot(data = all_subsets, mapping = aes(x = position, y = counts, fill = frame)) + theme_bw() +
    geom_col() + transition_states(size, transition_length = 2, state_length = 1) +
    labs(title = "no. of libraries: {closest_state}") +
    theme(axis.title = element_text(size = 36),
          axis.text = element_text(size = 22),
          title = element_text(size = 45),
          legend.text = element_text(size = 36)) +
    enter_fade() +
    exit_shrink() +
    ease_aes('sine-in-out')
    # plot2 <- ggplot() + geom_rect(aes(xmin = min(all_subsets$position),ymin=0, xmax = max(all_subsets$position), ymax = 1))
  animation <- animate(plot, renderer = gifski_renderer())
  return(animation)
}

slc35a4_animation <- animate_gene("SLC35A4", df_all, collection_genes,subregion = 900:1500, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500), normalize = FALSE)
atf4_animation <- animate_gene("ATF4", df_all, collection_genes,subregion = 800:1000, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)
bag1_animation <- animate_gene("BAG1", df_all, collection_genes,subregion = 1:300, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)
NIPAL3_animation <- animate_gene("NIPAL3", df_all, collection_genes,subregion = 1:300, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)
nras_animation <- animate_gene("NRAS", df_all, collection_genes,subregion = 1:500, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)
nxt1_animation <- animate_gene("NXT1", df_all, collection_genes,subregion = 1:500, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)

nras_animation_full <- animate_gene("NRAS", df_all, collection_genes,subregion = NA, subset_sizes = c(1,3,5,10,20,30,50,100,300,500,750,1000,1500,2000,2500, 2700), normalize = FALSE)

all_subsets

p <- ggplot(all_subsets[size == 2700][600:1500],aes(x = position, fill = frame, y= counts)) +
  geom_col()
ggplotly(p)
