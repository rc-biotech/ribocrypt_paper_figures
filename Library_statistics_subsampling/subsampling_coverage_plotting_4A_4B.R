 sampling_libs_dir <- "livemount/shared_results/ribocrypt_paper/sampling_coverage_randomlibs/"
 sampling_depth_dir <- "livemount/shared_results/ribocrypt_paper/sampling_coverage_multilib/"
 libs_sampling_files <- system(paste0("ls ", sampling_libs_dir, "*csv"),intern = TRUE)
 all_libs_samplings <- lapply(libs_sampling_files, fread)
 names(all_libs_samplings) <- gsub(".*/","",libs_sampling_files) %>% sub(".csv","",.)
 all_libs_samplings[[1]]
 samplings_merged <- rbindlist(all_libs_samplings, idcol = "gene")
 samplings_merged[, min_count := as.numeric(sub("min_count_","",variable))]
 samplings_merged[min_count %in% c(0,2), min_count := min_count + 1]
 samplings_merged[, min_count := as.factor(min_count)][,nlibs := as.factor(nlibs)]
 samplings_merged_mean <- samplings_merged[, .(value = mean(value)), by = .(min_count, nlibs)]
 cov_plot_line <- ggplot(data = samplings_merged_mean) + theme_minimal() +
     geom_line(aes(y = value, x = nlibs, group = min_count, 
                   # linetype = min_count, 
                   color = min_count), linewidth = 2) + labs("Min. Count") +
     ylab("Coverage (%)") + xlab(paste("# Libraries (of", samplings, "subsamplings)")) +
     theme(axis.text = element_text(size = 14),
                     axis.title = element_text(size = 18),
                     legend.text = element_text(size = 12),
                     legend.title = element_text(size = 16)); cov_plot_line
 cov_plot_line %>% ggplotly %>% plotly::config(
   toImageButtonOptions = list(
     format = "svg"))
 
 depth_sampling_files <- system(paste0("ls ", sampling_depth_dir, "*csv"),intern = TRUE)
 
   all_depth_samplings <- lapply(depth_sampling_files, fread)
 names(all_depth_samplings) <- gsub(".*/","",depth_sampling_files) %>% sub(".csv","",.)
 all_depth_samplings[[1]]

 
 depth_merged <- rbindlist(all_depth_samplings, idcol = "gene")
 depth_merged[, min_count := as.numeric(sub("min_count_","",variable))]
 depth_merged[min_count %in% c(0,2), min_count := min_count + 1]
 depth_merged[, min_count := as.factor(min_count)][,nlibs := as.factor(nlibs)]
 depth_merged_mean <- depth_merged[, .(value = mean(value, na.rm = TRUE), med_value = median(value, na.rm = TRUE)), by = .(min_count, nlibs,depth)]
 depth_merged_mean[,million_reads := as.factor(depth / 10^6)]
 
   cov_plot_line <- ggplot(data = depth_merged_mean) + theme_minimal() +
     geom_line(aes(y = value, x = million_reads, group = nlibs, color = nlibs)) + labs("Min. Count") +
     ylab("Coverage (%)") + xlab(paste("# Million reads, total genes:", length(unique(depth_merged$gene)))) + ggtitle("minimum count") +
     facet_wrap(~ min_count); cov_plot_line

   cov_plot_line <- ggplot(data = depth_merged_mean) + theme_minimal() +
     geom_line(aes(y = value, x = min_count, group = nlibs, color = nlibs)) + labs("Min. Count") +
     ylab("Coverage (%)") + xlab(paste("# minimum count, total genes:", length(unique(depth_merged$gene)))) + ggtitle("Million reads") +
     facet_wrap(~ million_reads); cov_plot_line
 
   cov_plot_line <- ggplot(data = depth_merged_mean) + theme_minimal() +
     geom_line(aes(y = value, x = nlibs, group = million_reads, color = million_reads)) + labs("Min. Count") +
     ylab("Coverage (%)") + xlab(paste("# Number of libraries, total genes:", length(unique(depth_merged$gene)))) + ggtitle("minimum count") +
     facet_wrap(~ min_count); cov_plot_line
   
   cov_plot_line <- ggplot(data = depth_merged_mean[min_count == 10][million_reads %in% c(50,500,1000)][!(million_reads == 1000 & nlibs %in% c(1,3))][!(nlibs %in% c(1,3))]) + theme_minimal() +
     geom_line(aes(y = value, x = nlibs, group = million_reads, color = million_reads, linetype = million_reads), linewidth = 2) + labs("Min. Count") +
     ylab("Coverage (%)") + xlab(paste("# Number of libraries, total genes:", length(unique(depth_merged$gene)))) +  theme(axis.text = element_text(size = 14),
                                                                                                                                                      axis.title = element_text(size = 18),
                                                                                                                                                      legend.text = element_text(size = 12),
                                                                                                                                                      legend.title = element_text(size = 16)); cov_plot_line
   cov_plot_line %>% ggplotly %>% plotly::config(
     toImageButtonOptions = list(
       format = "svg"))
   