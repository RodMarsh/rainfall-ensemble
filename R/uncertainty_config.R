# Normalise uncertainty parameters from YAML config to monthly vectors (12).
# Defaults from Jeffrey et al. for Australian tropical/semi-arid conditions.

normalize_uncertainty <- function(u) {
  defaults <- list(
    sigma_P = c(5,5,5,4,3,2,2,2,2,3,4,5),
    phi_P   = c(0.6,0.6,0.6,0.5,0.4,0.3,0.3,0.3,0.3,0.4,0.5,0.6),
    L_P_km  = c(90,90,90,110,130,150,150,150,140,120,100,90),
    p_tail  = c(0.10,0.10,0.10,0.08,0.05,0.03,0.03,0.03,0.04,0.06,0.08,0.10),
    k_tail  = 2.0,
    sigma_E = rep(1.2, 12),
    phi_E   = rep(0.4, 12),
    L_E_km  = rep(200, 12),
    rho_PE  = rep(NA_real_, 12),
    pmin_mm = rep(0.5, 12)
  )

  u <- modifyList(defaults, u %||% list())

  to_num12 <- function(x, allowNA = FALSE) {
    x <- suppressWarnings(as.numeric(unlist(x)))
    if (length(x) == 1L) x <- rep(x, 12)
    if (length(x) != 12L) stop("Uncertainty vector must be 12 numeric values.")
    if (!allowNA && any(!is.finite(x))) stop("Uncertainty vector must be 12 numeric values.")
    x
  }

  list(
    sigma_P = to_num12(u$sigma_P),
    phi_P   = to_num12(u$phi_P),
    L_P_km  = to_num12(u$L_P_km),
    p_tail  = to_num12(u$p_tail),
    k_tail  = as.numeric(u$k_tail %||% 1.0),
    sigma_E = to_num12(u$sigma_E),
    phi_E   = to_num12(u$phi_E),
    L_E_km  = to_num12(u$L_E_km),
    rho_PE  = to_num12(u$rho_PE, allowNA = TRUE),
    pmin_mm = to_num12(u$pmin_mm)
  )
}
