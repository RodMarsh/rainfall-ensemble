# Subcatchment shapefile → centroids → inter-subcatchment distance matrix.

# Read subcatchment polygons → metric CRS, valid, 2D
read_cats_sf <- function(cats_path, crs_metric = 3577, id_col = "MUID") {
  stopifnot(file.exists(cats_path))
  sf::read_sf(cats_path, quiet = TRUE) |>
    sf::st_make_valid() |>
    sf::st_transform(crs_metric) |>
    sf::st_zm(drop = TRUE) |>
    dplyr::transmute(!!id_col := !!rlang::sym(id_col), geometry = .data$geometry) |>
    dplyr::rename(MUID = !!rlang::sym(id_col))
}

# Subcatchment centroids (metric CRS) with standardised IDs
build_cats_xy_std <- function(cats_sf_metric, alias_map, id_col = "MUID") {
  pts <- sf::st_point_on_surface(cats_sf_metric$geometry)
  coords <- sf::st_coordinates(pts)
  ids <- sf::st_drop_geometry(cats_sf_metric) |>
    dplyr::select(!!rlang::sym(id_col)) |>
    dplyr::rename(subcatchment_id = !!rlang::sym(id_col)) |>
    dplyr::left_join(alias_map, by = "subcatchment_id") |>
    dplyr::mutate(subcatchment_id_std = dplyr::coalesce(
      .data$subcatchment_id_std, .data$subcatchment_id)) |>
    dplyr::select(.data$subcatchment_id_std)

  tibble::tibble(
    subcatchment_id_std = ids$subcatchment_id_std,
    x = coords[, 1],
    y = coords[, 2]
  )
}

# Inter-subcatchment distance matrix (km) from metric centroids
compute_dist_km_metric <- function(cats_xy_metric, crs = 3577L) {
  stopifnot(all(c("subcatchment_id_std", "x", "y") %in% names(cats_xy_metric)))

  pts <- sf::st_as_sf(cats_xy_metric, coords = c("x", "y"), crs = crs, remove = FALSE)
  Dm  <- sf::st_distance(pts)
  Dkm <- as.numeric(units::set_units(Dm, "km"))
  dim(Dkm) <- dim(Dm)
  dimnames(Dkm) <- list(cats_xy_metric$subcatchment_id_std,
                        cats_xy_metric$subcatchment_id_std)
  Dkm
}
