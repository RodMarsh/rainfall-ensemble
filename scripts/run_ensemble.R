#!/usr/bin/env Rscript
# Generate perturbed rainfall/PET ensemble from SILO forcing.
#
# Usage:
#   Rscript scripts/run_ensemble.R
#   Rscript scripts/run_ensemble.R --config conf/ensemble_parameters.yml
#   Rscript scripts/run_ensemble.R --n-members 10 --format csv

suppressPackageStartupMessages({
  library(dplyr)
  library(yaml)
  library(qs)
})

source("R/forcing_uncertainty.R")
source("R/uncertainty_config.R")
source("R/spatial.R")
source("R/forcing.R")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

config_path <- "conf/ensemble_parameters.yml"
n_members_override <- NULL
output_format <- "qs"  # "qs" or "csv"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--n-members" && i < length(args)) {
    n_members_override <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--format" && i < length(args)) {
    output_format <- args[i + 1]; i <- i + 2
  } else {
    message("Unknown argument: ", args[i]); i <- i + 1
  }
}

# ----
message("Reading config: ", config_path)
cfg <- yaml::read_yaml(config_path)
paths <- cfg$paths

n_members <- n_members_override %||% cfg$ensemble$n_members %||% 50L
base_seed <- cfg$ensemble$seed %||% 20251107L
use_parallel <- cfg$ensemble$parallel %||% FALSE
max_workers <- cfg$ensemble$max_workers %||% 4L

message("Ensemble: ", n_members, " members, seed=", base_seed,
        ", parallel=", use_parallel, ", format=", output_format)

message("\nReading forcing data ...")
forcing_raw <- read_silo_forcing(
  dir_rain = paths$forcing_rain,
  dir_pet  = paths$forcing_pet
)
message("  ", nrow(forcing_raw), " forcing rows, ",
        length(unique(forcing_raw$subcatchment_id)), " subcatchments, ",
        "date range: ", min(forcing_raw$date), " to ", max(forcing_raw$date))

# Attach standardised IDs
alias_map <- readr::read_csv(paths$alias_map, show_col_types = FALSE)
areas_tbl <- readr::read_csv(paths$areas_tbl, show_col_types = FALSE)
forcing_std_cat <- attach_ids_and_areas(forcing_raw, alias_map, areas_tbl)
forcing_std <- forcing_std_cat$forcing_daily

# Spatial structure
message("\nBuilding spatial covariance ...")
cats_sf <- read_cats_sf(paths$shapefile)
cats_xy <- build_cats_xy_std(cats_sf, alias_map, id_col = "MUID")
dist_km <- compute_dist_km_metric(cats_xy)
message("  ", nrow(cats_xy), " subcatchments, distance matrix: ",
        nrow(dist_km), "x", ncol(dist_km))

message("\nPrecomputing Cholesky factors ...")
unc_cfg <- normalize_uncertainty(cfg$uncertainty)

cholP_by_mo <- precompute_chol_by_month(dist_km, unc_cfg$sigma_P, unc_cfg$L_P_km)
cholE_by_mo <- precompute_chol_by_month(dist_km, unc_cfg$sigma_E, unc_cfg$L_E_km)
message("  Done (12 months x 2 variables)")

message("\nGenerating ", n_members, " ensemble members ...")

out_dir <- paths$output_dir %||% "outputs/ensemble_members"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

generate_member <- function(mem) {
  t0 <- Sys.time()
  result <- perturb_forcing_member(
    forcing_daily = forcing_std,
    cats_xy       = cats_xy,
    cfg           = unc_cfg,
    member_id     = mem,
    seed          = base_seed,
    cholP_by_mo   = cholP_by_mo,
    cholE_by_mo   = cholE_by_mo,
    use_rho_PE    = !all(is.na(unc_cfg$rho_PE))
  ) |>
    dplyr::rename(
      rainfall_mm = rainfall_mm_p,
      pet_mm      = pet_mm_p,
      forcingMember = member
    ) |>
    dplyr::select(date, subcatchment_id_std, rainfall_mm, pet_mm, forcingMember)

  # Write output
  fname <- sprintf("forcingMember_%03d", as.integer(mem))
  if (output_format == "csv") {
    path <- file.path(out_dir, paste0(fname, ".csv"))
    readr::write_csv(result, path)
  } else {
    path <- file.path(out_dir, paste0(fname, ".qs"))
    qs::qsave(result, path, preset = "balanced")
  }

  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message("  Member ", mem, " -> ", basename(path), " (", elapsed, "s)")
  path
}

if (use_parallel && n_members > 1 && requireNamespace("future.apply", quietly = TRUE)) {
  future::plan(future::multisession, workers = min(max_workers, n_members))
  paths_out <- future.apply::future_lapply(seq_len(n_members), generate_member)
  future::plan(future::sequential)
} else {
  paths_out <- lapply(seq_len(n_members), generate_member)
}

message("Wrote ", length(paths_out), " ensemble members to ", out_dir, "/")
