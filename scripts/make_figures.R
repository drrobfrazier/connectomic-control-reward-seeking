# Connectomic Control of Maladaptive Reward Seeking
# Recreate summary figures from archived result tables.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(forcats)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 0) stop("Run this script with Rscript scripts/make_figures.R")
script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
repo_root <- dirname(dirname(script_path))
results_dir <- file.path(repo_root, "results")
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

matched <- read_csv(file.path(results_dir, "matched_dan_null_10000.csv"), show_col_types = FALSE)
structural <- read_csv(file.path(results_dir, "structural_candidate_ranking.csv"), show_col_types = FALSE)
lesion <- read_csv(file.path(results_dir, "virtual_lesion_summary.csv"), show_col_types = FALSE)

obs_shared <- 2636
obs_all4 <- 1181
obs_balance <- 0.18564486164746236
obs_syn <- 120654

save_plot <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 320)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height)
}

base_theme <- theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p1 <- ggplot(matched, aes(shared_targets)) +
  geom_histogram(bins = 50) + geom_vline(xintercept = obs_shared, linetype = "dashed") +
  labs(title = "Matched DAN null: shared downstream targets", x = "Shared downstream targets", y = "Permutation count") + base_theme
save_plot(p1, "matched_null_shared_targets")

p2 <- ggplot(matched, aes(all4_shared)) +
  geom_histogram(bins = 50) + geom_vline(xintercept = obs_all4, linetype = "dashed") +
  labs(title = "Common coverage by all four punishment DANs", x = "Shared targets contacted by all four punishment DANs", y = "Permutation count") + base_theme
save_plot(p2, "all4_punishment_convergence")

p3 <- ggplot(matched, aes(mean_balance)) +
  geom_histogram(bins = 50) + geom_vline(xintercept = obs_balance, linetype = "dashed") +
  labs(title = "Reward-punishment polarization", x = "Mean balance score", y = "Permutation count") + base_theme
save_plot(p3, "reward_punishment_polarization")

p4 <- ggplot(matched, aes(shared_total_synapses)) +
  geom_histogram(bins = 50) + geom_vline(xintercept = obs_syn, linetype = "dashed") +
  labs(title = "Synaptic mass of the shared network", x = "Total synapses within shared targets", y = "Permutation count") + base_theme
save_plot(p4, "shared_synaptic_mass")

top20 <- structural %>%
  arrange(desc(interface_score)) %>%
  slice_head(n = 20) %>%
  mutate(label = paste0(type, " (", body_post, ")"), label = fct_reorder(label, interface_score))
p5 <- ggplot(top20, aes(interface_score, label)) +
  geom_segment(aes(x = 0, xend = interface_score, y = label, yend = label)) +
  geom_point() +
  labs(title = "Top structural interface candidates", x = "Composite structural score", y = NULL) + base_theme
save_plot(p5, "top_structural_candidates", 10, 7)

p6 <- ggplot(lesion, aes(median_reward_delta_rdi, median_conflict_delta_rdi, size = robust_lesion_score, label = type)) +
  geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(alpha = 0.8) + geom_text_repel(max.overlaps = Inf) +
  labs(title = "Virtual lesion effects", x = "Median reward-state ΔRDI", y = "Median conflict-state ΔRDI", size = "Robustness score") + base_theme
save_plot(p6, "virtual_lesion_effects", 10, 7)

p7 <- lesion %>% arrange(desc(robust_lesion_score)) %>% mutate(type = fct_reorder(type, robust_lesion_score)) %>%
  ggplot(aes(robust_lesion_score, type)) + geom_col() +
  labs(title = "Virtual lesion robustness ranking", x = "Robustness score", y = NULL) + base_theme
save_plot(p7, "virtual_lesion_ranking", 9, 6)

cat("Figures written to:", fig_dir, "\n")
