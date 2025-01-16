
result_dir <- "~/livemount/shared_results/ribocrypt_paper/ORF_prediction_results/CTE_NTE_analysis/"
#source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/NTE_CTE/CTE_subset.R")
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/NTE_CTE/CTE_subset_tables.R")
source("~/livemount/shared_scripts/ribocrypt_paper_figures/ORF_prediction_subsampling/NTE_CTE/NTE_CTE_phylo.R")

# Combined plot with phylo
all_plots <- gridExtra::arrangeGrob(metacoverage_plot_CTE, phylo_plot_CTE, ncol = 1)
all_plots <- gridExtra::grid.arrange(stats_plot_CTE, all_plots, nrow = 1)

# Send to discord and save to disc
massiveNGSpipe:::plot_all_versions(all_plots, file.path(result_dir, "CTE_figure_final"),
                                   send_to_discord = TRUE, width = 11)
