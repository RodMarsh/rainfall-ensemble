#!/usr/bin/env Rscript
# Diagnostic figures for ensemble perturbation structure.
# Simulates directly from parameters — no forcing data needed.
#
# Outputs: docs/figures/{ensemble_parameters,perturbation_distribution,spatial_correlation}.png

# Set working directory to repo root (parent of scripts/)
repo_root <- tryCatch(
  normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")),
  error = function(e) getwd()
)
setwd(repo_root)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(yaml)
})

source("R/forcing_uncertainty.R")
source("R/uncertainty_config.R")

# Config
cfg_raw   <- yaml::read_yaml("conf/ensemble_parameters.yml")
unc_cfg   <- normalize_uncertainty(cfg_raw$uncertainty)
out_dir   <- "docs/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

month_labels <- c("Jan","Feb","Mar","Apr","May","Jun",
                   "Jul","Aug","Sep","Oct","Nov","Dec")

# Colour palette
col_rain <- "#2166AC"
col_pet  <- "#B2182B"
col_tail <- "#E08214"
col_fill <- "#377EB8"

# Distance matrix
ref_path <- "data/metadata/spatial_reference.rds"
if (file.exists(ref_path)) {
  ref     <- readRDS(ref_path)
  dist_km <- ref$dist_km
  n_cats  <- nrow(ref$cats_xy)
  cat_ids <- ref$cats_xy$subcatchment_id_std
  message("Loaded spatial reference: ", n_cats, " subcatchments")
} else {
  # Synthetic fallback: 37 subcatchments spread over ~140 km
  message("No spatial_reference.rds found — using synthetic 37-subcatchment grid")
  n_cats  <- 37L
  set.seed(42)
  xy <- data.frame(x = runif(n_cats, 0, 100), y = runif(n_cats, 0, 140))
  dist_km <- as.matrix(dist(xy))
  cat_ids <- sprintf("cat_%02d", seq_len(n_cats))
  dimnames(dist_km) <- list(cat_ids, cat_ids)
}

# ----
cholP <- precompute_chol_by_month(dist_km, unc_cfg$sigma_P, unc_cfg$L_P_km)
cholE <- precompute_chol_by_month(dist_km, unc_cfg$sigma_E, unc_cfg$L_E_km)

# Fig 1: monthly parameter profiles
message("Rendering ensemble_parameters.png ...")

param_df <- tibble(
  month = rep(1:12, 6),
  value = c(unc_cfg$sigma_P, unc_cfg$phi_P, unc_cfg$L_P_km,
            unc_cfg$sigma_E, unc_cfg$phi_E, unc_cfg$p_tail),
  param = rep(c("sigma[P]~(mm)", "phi[P]", "L[P]~(km)",
                "sigma[E]~(mm)", "phi[E]", "p[tail]"), each = 12),
  variable = rep(c("Rainfall","Rainfall","Rainfall","PET","PET","Rainfall"), each = 12)
)
param_df$param <- factor(param_df$param,
  levels = c("sigma[P]~(mm)", "phi[P]", "L[P]~(km)",
             "sigma[E]~(mm)", "phi[E]", "p[tail]"))

p1 <- ggplot(param_df, aes(month, value, colour = variable)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  facet_wrap(~ param, scales = "free_y", labeller = label_parsed, ncol = 3) +
  scale_x_continuous(breaks = 1:12, labels = substr(month_labels, 1, 1)) +
  scale_colour_manual(values = c("Rainfall" = col_rain, "PET" = col_pet)) +
  labs(
    title = "Monthly ensemble perturbation parameters",
    subtitle = paste0("Exponential spatial covariance, AR(1) temporal persistence, ",
                      "heavy-tail mixture (k = ", unc_cfg$k_tail, ")"),
    x = "Month", y = NULL, colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 13)
  )

ggsave(file.path(out_dir, "ensemble_parameters.png"), p1,
       width = 10, height = 6, dpi = 300, bg = "white")


# Fig 2: perturbation distributions by month
message("Rendering perturbation_distribution.png ...")

# Simulate 20 members x 3 years of daily perturbation fields
sim_dates <- seq(as.Date("2020-01-01"), as.Date("2022-12-31"), by = "day")
n_members <- 20L
base_seed <- 12345L

all_eps <- lapply(seq_len(n_members), function(mem) {
  seed_i <- base_seed + 1000L * mem
  eps_P <- simulate_field_fast(sim_dates, cholP, unc_cfg$phi_P,
                               p_tail_vec = unc_cfg$p_tail,
                               k_tail = unc_cfg$k_tail,
                               seed = seed_i, recenter_month = TRUE)
  eps_E <- simulate_field_fast(sim_dates, cholE, unc_cfg$phi_E,
                               p_tail_vec = NULL, k_tail = 1.0,
                               seed = seed_i + 1L, recenter_month = TRUE)
  tibble(
    date   = rep(sim_dates, n_cats),
    cat_id = rep(seq_len(n_cats), each = length(sim_dates)),
    month  = as.integer(format(rep(sim_dates, n_cats), "%m")),
    eps_P  = as.vector(eps_P),
    eps_E  = as.vector(eps_E),
    member = mem
  )
}) |> bind_rows()

# Summary stats per month
eps_summary <- all_eps |>
  group_by(month) |>
  summarise(
    P_mean = mean(eps_P),
    P_sd   = sd(eps_P),
    P_p05  = quantile(eps_P, 0.05),
    P_p25  = quantile(eps_P, 0.25),
    P_p50  = median(eps_P),
    P_p75  = quantile(eps_P, 0.75),
    P_p95  = quantile(eps_P, 0.95),
    E_mean = mean(eps_E),
    E_sd   = sd(eps_E),
    E_p05  = quantile(eps_E, 0.05),
    E_p25  = quantile(eps_E, 0.25),
    E_p50  = median(eps_E),
    E_p75  = quantile(eps_E, 0.75),
    E_p95  = quantile(eps_E, 0.95),
    .groups = "drop"
  )

# Panel A: Rainfall perturbation violin
pA <- ggplot(all_eps |> mutate(month_f = factor(month, labels = month_labels)),
             aes(x = month_f, y = eps_P)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_violin(fill = col_rain, colour = NA, alpha = 0.35, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.2,
               fill = "white", colour = col_rain, linewidth = 0.35) +
  labs(title = "Rainfall perturbation (mm/day)",
       x = NULL, y = "Perturbation (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# Panel B: PET perturbation violin
pB <- ggplot(all_eps |> mutate(month_f = factor(month, labels = month_labels)),
             aes(x = month_f, y = eps_E)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_violin(fill = col_pet, colour = NA, alpha = 0.35, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.3, outlier.alpha = 0.2,
               fill = "white", colour = col_pet, linewidth = 0.35) +
  labs(title = "PET perturbation (mm/day)",
       x = NULL, y = "Perturbation (mm)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

# Panel C: Inter-subcatchment SD (spatial spread) by month
spatial_spread <- all_eps |>
  group_by(member, date, month) |>
  summarise(sd_P = sd(eps_P), sd_E = sd(eps_E), .groups = "drop")

pC <- spatial_spread |>
  pivot_longer(cols = c(sd_P, sd_E),
               names_to = "variable", values_to = "sd") |>
  mutate(
    variable = ifelse(variable == "sd_P", "Rainfall", "PET"),
    month_f  = factor(month, labels = month_labels)
  ) |>
  ggplot(aes(x = month_f, y = sd, fill = variable)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.15, linewidth = 0.3,
               position = position_dodge(width = 0.75), width = 0.6) +
  scale_fill_manual(values = c("Rainfall" = col_rain, "PET" = col_pet)) +
  labs(title = "Spatial spread: inter-subcatchment SD per day",
       x = NULL, y = "SD across subcatchments (mm)", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11))

p2 <- (pA | pB) / pC +
  plot_annotation(
    title = "Distribution of ensemble perturbations across subcatchments",
    subtitle = paste0("Simulated from ", n_members, " members x 3 years x ",
                      n_cats, " subcatchments. Monthly recentring ensures zero mean."),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey30")
    )
  )

ggsave(file.path(out_dir, "perturbation_distribution.png"), p2,
       width = 12, height = 9, dpi = 300, bg = "white")


# Fig 3: spatial correlation decay
message("Rendering spatial_correlation.png ...")

# Empirical pairwise correlations, pooled across 3 members
months_show <- c(1, 4, 7, 10)  # Jan, Apr, Jul, Oct
month_names <- c("1" = "Jan (wet)", "4" = "Apr (transition)",
                 "7" = "Jul (dry)", "10" = "Oct (onset)")

cor_data <- lapply(months_show, function(mo) {
  # Pool 3 members for more stable correlations
  eps_pool <- lapply(1:3, function(mem) {
    seed_i <- base_seed + 1000L * mem
    eps <- simulate_field_fast(sim_dates, cholP, unc_cfg$phi_P,
                               p_tail_vec = unc_cfg$p_tail,
                               k_tail = unc_cfg$k_tail,
                               seed = seed_i, recenter_month = TRUE)
    idx <- which(as.integer(format(sim_dates, "%m")) == mo)
    eps[idx, , drop = FALSE]
  }) |> do.call(what = rbind)

  C <- cor(eps_pool)
  # Extract all unique pairs
  pairs <- expand.grid(i = seq_len(n_cats), j = seq_len(n_cats)) |>
    filter(i < j)
  tibble(
    month  = mo,
    dist   = dist_km[cbind(pairs$i, pairs$j)],
    corr   = C[cbind(pairs$i, pairs$j)]
  )
}) |> bind_rows()

# Theoretical curves
dist_seq <- seq(0, max(dist_km), length.out = 200)
theory <- lapply(months_show, function(mo) {
  tibble(
    month = mo,
    dist  = dist_seq,
    corr  = exp(-dist_seq / unc_cfg$L_P_km[mo])
  )
}) |> bind_rows()

p3 <- ggplot() +
  geom_point(data = cor_data |> mutate(month_lbl = month_names[as.character(month)]),
             aes(x = dist, y = corr), alpha = 0.25, size = 1, colour = col_rain) +
  geom_line(data = theory |> mutate(month_lbl = month_names[as.character(month)]),
            aes(x = dist, y = corr), linewidth = 0.8, colour = "black", linetype = "dashed") +
  facet_wrap(~ month_lbl, nrow = 1) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  scale_y_continuous(limits = c(-0.2, 1.05), breaks = seq(-0.2, 1, 0.2)) +
  labs(
    title = "Spatial correlation of rainfall perturbations",
    subtitle = paste0("Points: empirical pairwise correlations (",
                      n_cats, " subcatchments, 3 members pooled). ",
                      "Dashed: theoretical exp(-D/L)"),
    x = "Inter-subcatchment distance (km)",
    y = "Correlation"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30")
  )

ggsave(file.path(out_dir, "spatial_correlation.png"), p3,
       width = 12, height = 4.5, dpi = 300, bg = "white")

message("\nAll figures saved to ", out_dir, "/")
