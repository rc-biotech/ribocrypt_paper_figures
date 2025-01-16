
source("~/livemount/shared_scripts/ribocrypt_paper_figures/predicted_novel_orfs_stats_per_species.R")
source("~/livemount/shared_scripts/ribocrypt_paper_figures/TIS_is_correct_check.R")

# Merge with novel orf plot
plot_TIS_novel <- cowplot::plot_grid(plot_TIS, plot, ncol = 1); plot_TIS_novel
svg_name <- "Figure5_TIS_distance_novel.svg"
jpg_name <- "Figure5_TIS_distance_novel.jpg"
ggsave(file.path(plot_folder, "Figure5_TIS_distance_novel.jpg"), plot_TIS_novel, width = 8, height = 8)
ggsave(file.path(plot_folder, svg_name), plot_TIS_novel, width = 8, height = 8)
googledrive::drive_put(file.path(plot_folder, svg_name),
                       "https://drive.google.com/drive/folders/13o2M5KtqBKICo_gFh66anp3_38vSsbrg",
                       name = svg_name)
googledrive::drive_put(file.path(plot_folder, jpg_name),
                       "https://drive.google.com/drive/folders/13o2M5KtqBKICo_gFh66anp3_38vSsbrg",
                       name = jpg_name)
