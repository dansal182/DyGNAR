structural_GNARfit <- function(longitudinal_network_data, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups) { # nolint
    vts <- longitudinal_network_data
     # nolint
}


temporal_trend_spline_smoother <- function(spline_residuals, vts_dimension, plag) { # nolint
    spline_res_mat = residuals_matrix(spline_residuals, vts_dimension = vts_dimension) # nolint
    t1 = plag + 1
    tn = nrow(spline_res_mat) + plag
    time_range = c(t1:tn) # nolint
    mu_hat_list <- lapply(1:vts_dimension, function(x) {smooth.spline(time_range, spline_res_mat[, x], cv = FALSE)}) # nolint
    return(mu_hat_list)
}


seasonal_trend_smoother <- function(seasonal_residuals, vts_dimension, plag, basis_dim=3, period = 2 * pi) { # nolint
    seasonal_res_mat = residuals_matrix(seasonal_residuals, vts_dimension = vts_dimension) # nolint
    t1 = plag + 1
    tn = nrow(seasonal_res_mat) + plag
    time_range = c(t1:tn) # nolint
    fourier_basis <- fda::create.fourier.basis(rangeval=range(time_range), nbasis=basis_dim, period=period) # nolint
    pi_hat_list <- lapply(1:vts_dimension, function(x) {fda::smooth.basis(time_range, seasonal_res_mat[, x], fdParaobj=fourier_basis)}) # nolint
    return(pi_hat_list)
}