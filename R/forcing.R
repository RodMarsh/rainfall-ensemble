# Read SILO forcing (tab-delimited text, one file per subcatchment) and
# attach standardised IDs and areas.

read_silo_forcing <- function(dir_rain, dir_pet, tz = "Australia/Brisbane") {
  rain_files <- fs::dir_ls(dir_rain, glob = "*Rainfall_*_*.txt")
  pet_files  <- fs::dir_ls(dir_pet,  glob = "*Evap_*_*.txt")

  if (length(rain_files) == 0) stop("No rainfall files found in: ", dir_rain)
  if (length(pet_files)  == 0) stop("No PET files found in: ", dir_pet)

  id_from_rain <- function(path) {
    id <- stringr::str_match(fs::path_file(path), "Rainfall_\\d+-\\d+_(.+)\\.txt$")[,2]
    if (is.na(id)) stop("Could not parse subcatchment_id from rainfall filename: ", path)
    id
  }
  id_from_pet  <- function(path) {
    id <- stringr::str_match(fs::path_file(path), "Evap_\\d+-\\d+_(.+)\\.txt$")[,2]
    if (is.na(id)) stop("Could not parse subcatchment_id from PET filename: ", path)
    id
  }

  parse_date_local <- function(x) {
    as.Date(lubridate::ymd_hms(x, tz = tz, quiet = TRUE))
  }

  read_dhi_rain <- function(path) {
    readr::read_delim(path, delim = "\t", skip = 1,
      col_types = readr::cols(Time = readr::col_character(), Current = readr::col_double()),
      show_col_types = FALSE) |>
      dplyr::mutate(date = parse_date_local(.data$Time)) |>
      dplyr::transmute(
        date = .data$date,
        subcatchment_id = id_from_rain(path),
        rainfall_mm = .data$Current
      ) |>
      dplyr::group_by(.data$subcatchment_id, .data$date) |>
      dplyr::summarise(rainfall_mm = mean(.data$rainfall_mm, na.rm = TRUE), .groups = "drop")
  }

  read_dhi_pet <- function(path) {
    readr::read_delim(path, delim = "\t", skip = 1,
      col_types = readr::cols(Time = readr::col_character(), Current = readr::col_double()),
      show_col_types = FALSE) |>
      dplyr::mutate(date = parse_date_local(.data$Time)) |>
      dplyr::transmute(
        date = .data$date,
        subcatchment_id = id_from_pet(path),
        pet_mm = .data$Current
      ) |>
      dplyr::group_by(.data$subcatchment_id, .data$date) |>
      dplyr::summarise(pet_mm = mean(.data$pet_mm, na.rm = TRUE), .groups = "drop")
  }

  rain_all <- purrr::map_dfr(rain_files, read_dhi_rain)
  pet_all  <- purrr::map_dfr(pet_files,  read_dhi_pet)

  forcing_daily <- list(rain_all, pet_all) |>
    purrr::reduce(dplyr::full_join, by = c("date", "subcatchment_id")) |>
    dplyr::group_by(.data$subcatchment_id) |>
    dplyr::arrange(.data$date, .by_group = TRUE) |>
    tidyr::complete(date = seq(min(.data$date), max(.data$date), by = "day")) |>
    dplyr::ungroup() |>
    dplyr::select(.data$date, .data$subcatchment_id, .data$rainfall_mm, .data$pet_mm)

  # Fail fast on missing
  miss <- forcing_daily |>
    dplyr::filter(is.na(.data$rainfall_mm) | is.na(.data$pet_mm))
  if (nrow(miss) > 0) {
    print(miss |> dplyr::slice_head(n = 10))
    stop("Missing forcing detected (see first 10 rows).")
  }

  forcing_daily
}

# Attach standardised IDs and areas to forcing data
attach_ids_and_areas <- function(forcing_daily, alias_map, areas_tbl) {
  stopifnot(all(c("subcatchment_id","subcatchment_id_std") %in% names(alias_map)))
  stopifnot(all(c("subcatchment_id_std","area_km2") %in% names(areas_tbl)))

  forcing_ids <- forcing_daily |> dplyr::distinct(subcatchment_id)
  id_map <- alias_map |>
    dplyr::right_join(forcing_ids, by = "subcatchment_id", relationship = "one-to-one") |>
    dplyr::mutate(subcatchment_id_std = dplyr::coalesce(subcatchment_id_std, subcatchment_id))

  cat_catalogue <- id_map |>
    dplyr::distinct(subcatchment_id_std, .keep_all = TRUE) |>
    dplyr::left_join(areas_tbl, by = "subcatchment_id_std", relationship = "many-to-one") |>
    dplyr::arrange(subcatchment_id_std)

  n0 <- nrow(forcing_daily)
  forcing_std <- forcing_daily |>
    dplyr::left_join(id_map |> dplyr::select(subcatchment_id, subcatchment_id_std),
                     by = "subcatchment_id", relationship = "many-to-one")
  stopifnot(nrow(forcing_std) == n0)

  list(forcing_daily = forcing_std, cat_catalogue = cat_catalogue)
}
