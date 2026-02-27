# AR(1) spatio-temporal perturbation of daily P and PET fields.
#
# eps_t = phi_m * eps_{t-1} + sqrt(1 - phi_m^2) * (C_m %*% Z_t)
#
# C_m = chol of exponential spatial covariance for month m.
# Z_t ~ N(0,I) with optional mixture inflation (p_tail, k_tail).
# See Jeffrey et al. (2001) for parameter basis.


make_cov <- function(D_km, sigma, L_km, nugget = 1e-6) {
  stopifnot(is.matrix(D_km), nrow(D_km) == ncol(D_km))
  stopifnot(is.numeric(sigma), length(sigma) == 1L, is.finite(sigma), sigma >= 0)
  stopifnot(is.numeric(L_km), length(L_km) == 1L, is.finite(L_km), L_km > 0)

  K <- exp(-D_km / L_km)
  Sigma <- (sigma^2) * K
  diag(Sigma) <- diag(Sigma) + nugget^2
  Sigma
}


# Monthly Cholesky factors — returns list(12) of upper-tri factors
precompute_chol_by_month <- function(dist_km, sigma_vec, L_vec) {
  stopifnot(is.matrix(dist_km), nrow(dist_km) == ncol(dist_km))

  s <- as.numeric(unlist(sigma_vec))
  L <- as.numeric(unlist(L_vec))
  if (length(s) != 12L || length(L) != 12L || any(!is.finite(s)) || any(!is.finite(L))) {
    stop("precompute_chol_by_month(): sigma_vec and L_vec must be 12 finite numeric values.")
  }

  lapply(seq_len(12), function(mo) {
    Sigma <- make_cov(dist_km, s[mo], L[mo])
    chol(Sigma)
  })
}


# ----
simulate_field_fast <- function(dates, chol_by_mo, phi_vec,
                                p_tail_vec = NULL, k_tail = 1.0,
                                seed = 1, recenter_month = TRUE) {
  set.seed(seed)
  n  <- ncol(chol_by_mo[[1]])
  Tn <- length(dates)
  eps <- matrix(0.0, nrow = Tn, ncol = n)
  z_prev <- rep(0, n)

  months <- as.integer(format(dates, "%m"))
  for (t in seq_len(Tn)) {
    m <- months[t]
    z  <- drop(chol_by_mo[[m]] %*% stats::rnorm(n))
    if (!is.null(p_tail_vec) && stats::runif(1) < p_tail_vec[m]) z <- k_tail * z
    eps[t, ] <- phi_vec[m] * z_prev + sqrt(1 - phi_vec[m]^2) * z
    z_prev   <- eps[t, ]
  }
  if (recenter_month) {
    for (m in 1:12) {
      idx <- which(months == m)
      if (length(idx)) {
        mu <- colMeans(eps[idx, , drop = FALSE])
        eps[idx, ] <- sweep(eps[idx, , drop = FALSE], 2L, mu, FUN = "-")
      }
    }
  }
  eps
}


# Perturb forcing for a single ensemble member
perturb_forcing_member <- function(forcing_daily, cats_xy, cfg,
                                   member_id, seed,
                                   cholP_by_mo, cholE_by_mo,
                                   rain_threshold = 0.5,
                                   use_rho_PE = TRUE) {

  ids  <- cats_xy$subcatchment_id_std
  days <- sort(unique(as.Date(forcing_daily$date)))

  fx <- forcing_daily |>
    dplyr::filter(.data$subcatchment_id_std %in% ids) |>
    dplyr::mutate(date = as.Date(.data$date)) |>
    dplyr::arrange(.data$date, .data$subcatchment_id_std)

  # AR(1) innovations (P & E), month-specific
  seed_i <- as.integer(seed) + 1000L * as.integer(member_id)
  P_eps <- simulate_field_fast(days, chol_by_mo = cholP_by_mo,
                               phi_vec = cfg$phi_P,
                               p_tail_vec = cfg$p_tail,
                               k_tail = cfg$k_tail,
                               seed = seed_i, recenter_month = TRUE)
  E_eps <- simulate_field_fast(days, chol_by_mo = cholE_by_mo,
                               phi_vec = cfg$phi_E,
                               p_tail_vec = NULL, k_tail = 1.0,
                               seed = seed_i + 1L, recenter_month = TRUE)

  # Optional P-E coupling on spatial means
  if (use_rho_PE && !all(is.na(cfg$rho_PE))) {
    mos <- as.integer(format(days, "%m"))
    for (t in seq_along(days)) {
      rho <- cfg$rho_PE[mos[t]]
      if (!is.na(rho) && abs(rho) > 0) {
        muP <- mean(P_eps[t, ]); sP <- sd(P_eps[t, ])
        muE <- mean(E_eps[t, ]); sE <- sd(E_eps[t, ])
        if (sP > 0 && sE > 0) {
          z1 <- rnorm(1); z2 <- rnorm(1)
          muP2 <- muP + sP * z1
          muE2 <- muE + sE * (rho * z1 + sqrt(1 - rho^2) * z2)
          P_eps[t, ] <- (P_eps[t, ] - muP) + muP2
          E_eps[t, ] <- (E_eps[t, ] - muE) + muE2
        }
      }
    }
  }

  # Map innovations to rows via (date, id) index
  key <- tibble::tibble(
    date = rep(days, each = length(ids)),
    subcatchment_id_std = rep(ids, times = length(days)),
    irow = match(date, days),
    icol = match(subcatchment_id_std, ids),
    eP   = P_eps[cbind(irow, icol)],
    eE   = E_eps[cbind(irow, icol)]
  )

  out <- fx |>
    dplyr::left_join(key[, c("date","subcatchment_id_std","eP","eE")],
                     by = c("date","subcatchment_id_std")) |>
    dplyr::mutate(
      month = as.integer(format(.data$date, "%m")),
      addP  = dplyr::if_else(.data$rainfall_mm > (cfg$pmin_mm[.data$month]), .data$eP, 0),
      rainfall_mm_p = pmax(0, .data$rainfall_mm + .data$addP),
      pet_mm_p      = pmax(0, .data$pet_mm + .data$eE),
      member        = as.integer(member_id)
    ) |>
    dplyr::select(.data$date, .data$subcatchment_id_std,
                  rainfall_mm_p = .data$rainfall_mm_p,
                  pet_mm_p      = .data$pet_mm_p,
                  .data$member)

  out
}
