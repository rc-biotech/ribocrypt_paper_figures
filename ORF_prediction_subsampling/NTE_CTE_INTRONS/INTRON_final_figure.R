# Source sub scripts
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/intron_results/"
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/introns/introns_exons_pair_analysis.R")
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/introns/introns_exons_pair_analysis_tables.R")
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/introns/intron_phylo_analysis.R")


# Combined plot with phylo
all_plots <- gridExtra::arrangeGrob(metacoverage_plot, phylo_plot, ncol = 1)
all_plots <- gridExtra::grid.arrange(intronic_stats_plot, all_plots, nrow = 1)

massiveNGSpipe:::plot_all_versions(all_plots, file.path(result_dir, "intron_figure_final"),
                                   send_to_discord = TRUE, width = 11)
# library(plotly)
# ggplot(cov_final_filtered_intron) +
#   geom_line(aes(x = position, y = score, color = fraction))
#
# cov_final[,intron := ifelse(fraction > 2,"3+",fraction)]
# cov_intron_plot <- cov_final[,.(score = mean(score)),by = .(intron, position)]
#
# intp <- ggplot(cov_intron_plot) +
#   geom_line(aes(x = position, y = score, color = intron ))
# ggplotly(intp)
# intp <- ggplot(cov_final) +
#   geom_line(aes(x = position, y = score, color = as.factor(fraction )))
# ggplotly(intp)
