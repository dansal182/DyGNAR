eval_mu_term <- function(t, T0 = 100, type="linear") { # nolint
    x = t / T0 # nolint
    if (type == "linear") { # nolint
        out = 3 * x # nolint
    } else if (type == "doppler") { # nolint # nolint: indentation_linter.
        out = 5 * sqrt(x * (1 - x)) * sin((2.1 * x * pi ) / (x + 0.05)) # nolint
    } else if (type == "sinc") { #  nolint # nolint
        out = sin(10 * (x + 0.05 )) / (x + 0.05) # nolint
    } else if (type == "sin_prod") { # nolint # nolint: indentation_linter.
        out = sin(10 * (x + 0.05 )) * (10 * x + 0.05) # nolint
  } else {
    out = sin(6 * (x - 1/3 + 0.05)) / (x - 1/3 + 0.05) + sin(4 * (x - 2/3 + 0.05)) / (x - 2/3 + 0.05) # nolint
  }
  return(out)
}


eval_seasonal_term <- function(t, period = 1, type = "one") {
    k = 2 * pi / period # nolint
    if (type == "one") {
    out <- 2 * cos(t * k)
    } else if (type == "two") { # nolint: indentation_linter.
    out <- sin(t * k)
  } else if (type == "three") {
    out <- cos(t * k) * sin(t * k) + cos(t * k)
  } else if (type == "four") {
    out <- 3 * (cos(t * k) * sin(t * k) + sin(t * k))
  } else {
    out <- 2 * cos(t * k) * sin(t * k)
  }
  return(out)
}


structural_GNARsim <- function(GNARsim, covariate_groups, predictor_matrix, eta_coeffs, mu_types, pi_types, period = 1) { # nolint
    T0 = nrow(GNARsim) # nolint
    d = ncol (GNARsim) # nolint
    Zdmat <- build_predictor_design(predictor_mat = predictor_matrix, covariate_groups = covariate_groups, vts_dimension = d) # nolint
    Omega <- residuals_matrix(Zdmat %*% eta_coeffs, d) # nolint
    mu_mat <- vapply(1:d, function(x) {return(eval_mu_term(c(1:T0), T0 = T0, type = mu_types[x]))}, rep(0.0, T0)) # nolint
    pi_mat <- vapply(1:d, function(x) {return(eval_seasonal_term(c(1:T0), period = period, type = pi_types[x]))}, rep(0.0, T0)) # nolint
    Ymat <- mu_mat + pi_mat + Omega + GNARsim # nolint
    colnames(Ymat) <- vapply(1:d, function(x) {return(paste0("Y", as.character(x)))}, "") # nolint
    return(Ymat)
}