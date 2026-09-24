structural_GNARfit <- function(longitudinal_network_data, predictor_design, network, periodic_basis_dim=3, period=2*pi, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups, iter=100, convergence_tolreance=1e-7) { # nolint
    # set up the main variables and design matrices  # nolint: indentation_linter, line_length_linter.
    ymat <- longitudinal_network_data # nolint
    vts_dimension = ncol(ymat) # nolint
    d = vts_dimension # nolint
    ycol_means <- colMeans(ymat)
    plag = max(alpha_order) # nolint
    C = length(covariate_groups) # nolint
    dz = ncol(predictor_design) # nolint
    q = sum(vapply(1:C, function(x) {get_community_parameter_number(alpha_order[x], beta_order[[x]], interaction_groups[[x]])}, 0.0)) # nolint
    time_steps = nrow(ymat) # nolint
    t1 = plag + 1 # nolint
    tn = time_steps # nolint
    ymat <- vapply(1:d, function(x) {return(ymat[, x] - ycol_means[x])}, rep(0.0, tn)) # nolint
    yvec <-  recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(ymat[x, ]))})) # nolint
    # check if the model includes a fixed effects predictor term # nolint
    Zdmat <- build_predictor_design(predictor_mat = predictor_design, covariate_groups = covariate_groups, vts_dimension = d) # nolint
    eta_hat = rep(0, C * dz) # nolint
    # nolint: indentation_linter.
    # Initialise values for the backnet-fitting algorithm # nolint
    theta_hat =  rep(0, q) # nolint
    fourier_basis <- fda::create.fourier.basis(rangeval=range(1, tn), nbasis=periodic_basis_dim, period=period)  # nolint
    Phi_mat <- fda::eval.basis(c(t1:tn), fourier_basis) # nolint
    phi_vec <- fda::eval.basis(c(1:plag), fourier_basis) # nolint
    mu_hat <- rep(0.0, tn * d) # nolint
    pi_hat <- rep(0.0, tn * d) # nolint
    beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0) # nolint
    r_stage_adjacency_matrices <- get_r_stages_adjacency_list(as.matrix(network), max(beta_depths)) # nolint
    t0 = d * plag + 1 # nolint
    td = tn * d # nolint
    xhat <- rep(0.0, td - d * plag) # nolint
    # run the backnet-fitting algorithm until convergence # nolint
    for (j in 1:iter) { # nolint
        # t0:tdprint(length(yvec[(d * plag + 1):(d * tn)])) # nolint
        mu_res <- yvec[t0:td] - (pi_hat[t0:td] + Zdmat[t0:td, ] %*% eta_hat + xhat) # nolint # nolint
        mu_hat_mat <- temporal_trend_spline_smoother(spline_residuals = mu_res, spar_vec = c(0.75, 1, 1.5, 0.75, 1), vts_dimension = d, plag = plag) # nolint
        mu_hat <- recursive_dmat_cols_cpp(lapply(1:tn, function(x) {return(as.matrix(mu_hat_mat[x, ]))})) # nolint
        pi_res <- yvec[t0:td] - (mu_hat[t0:td] + Zdmat[t0:td, ] %*% eta_hat + xhat) # nolint
        pi_hat_mat <- seasonal_trend_smoother(seasonal_residuals = pi_res, vts_dimension = d, eval_fourier_basis = Phi_mat, eval_basis_pre_lag = phi_vec) # nolint
        pi_hat <- recursive_dmat_cols_cpp(lapply(1:tn, function(x) {return(as.matrix(pi_hat_mat[x, ]))})) # nolint
        eta_res <- yvec[t0:td] - (mu_hat[t0:td] + pi_hat[t0:td] + xhat) # nolint
        eta_hat <- least_squares_solver(Zdmat[t0:td, ], eta_res) # nolint
        xhat <- yvec - (mu_hat + pi_hat + Zdmat %*% eta_hat) # nolint
        xhat_vts <- residuals_matrix(xhat, d) # nolint
        xhat_means <- colMeans(xhat_vts) # nolint
        xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
        GNAR_design <- build_GNAR_design(vts = xhat_vts, stages_tensor = r_stage_adjacency_matrices, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = weight_matrix, covariate_groups = covariate_groups, interaction_groups = interaction_groups) # nolint # nolint
        Rmat <- GNAR_design[, 2:(q + 1)] # nolint
        theta_hat <- least_squares_solver(Rmat, GNAR_design[, 1]) # nolint
        xhat <- Rmat %*% theta_hat
        # xhat_vts <- residuals_matrix(xhat, d) # nolint
        # GNAR_design <- build_GNAR_design(vts = xhat_vts, stages_tensor = r_stage_adjacency_matrices, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = weight_matrix, covariate_groups = covariate_groups, interaction_groups = interaction_groups) # nolint # nolint
        # Rmat_aux <- GNAR_design[, 2:(q + 1)] # nolint
        # theta_hat <- least_squares_solver(Rmat_aux, GNAR_design[, 1]) # nolint
        # xhat <- Rmat %*% theta_hat
    } # nolint
    xhat <- yvec - (mu_hat + pi_hat + Zdmat %*% eta_hat) # nolint
    xhat_vts <- residuals_matrix(xhat, d) # nolint
    xhat_means <- colMeans(xhat_vts) # nolint
    xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
    rownames(theta_hat) <- colnames(GNAR_design)[2:(q + 1)]
    out <- list(mu_hat_mat, pi_hat_mat, eta_hat, xhat_vts, theta_hat, eta_res, Zdmat[t0:td, ], GNAR_design[, 1], Rmat)
    names(out) <- c("temporal_trend_fitted_values", "seasonal_trend_fitted_values", "fixed_network_effects", "GNAR_residuals", "GNAR_parameters", "adjusted_response", "predictor_design", "GNAR_response", "GNAR_design") # nolint
    return(out)
}


temporal_trend_spline_smoother <- function(spline_residuals, spar_vec=NULL, vts_dimension, plag) { # nolint
    spline_res_mat = residuals_matrix(spline_residuals, vts_dimension = vts_dimension) # nolint
    t1 = plag + 1 # nolint
    tn = nrow(spline_res_mat) + plag # nolint
    time_range = c(t1:tn) # nolint
    # spar_vec <- c(0.75, 1, 1.5, 0.75, 1) for simulation in the introduction
    if (is.null(spar_vec)) {
        spar_vec <- runif(vts_dimension, 0.6, 1.2) # nolint: indentation_linter, line_length_linter.
    } # nolint: indentation_linter.
    mu_hat_mat <- vapply(1:vts_dimension, function(x) {smooth.spline(time_range, spline_res_mat[, x], spar = spar_vec[x], cv = NA, keep.stuff = FALSE)$y}, rep(0.0, nrow(spline_res_mat))) # nolint
    mu_means <- colMeans(mu_hat_mat) # nolint
    mu_hat_mat <- vapply(1:vts_dimension, function(x) {return(mu_hat_mat[, x] - mu_means[x])}, rep(0.0, nrow(spline_res_mat))) # nolint
    slopes_hat <- vapply(1:vts_dimension, function(x) {return(mu_hat_mat[2, x] - mu_hat_mat[1, x])}, 0.0) # nolint
    mu_boundary_mat <- vapply(1:vts_dimension, function(x) {return(slopes_hat[x] * c(1:plag))}, rep(0.0, plag)) # nolint
    mu_hat_mat <- rbind(mu_boundary_mat, mu_hat_mat) # nolint
    return(mu_hat_mat)
}

temporal_trend_spline_smoother_cv <- function(spline_residuals, spar_vec=NULL, vts_dimension, plag) { # nolint
    spline_res_mat = residuals_matrix(spline_residuals, vts_dimension = vts_dimension) # nolint
    t1 = plag + 1 # nolint
    tn = nrow(spline_res_mat) + plag # nolint
    time_range = c(t1:tn) # nolint
    # spar_vec <- c(0.75, 1, 1.5, 0.75, 1) for simulation in the introduction
    mu_hat_mat <- vapply(1:vts_dimension, function(x) {smooth.spline(time_range, spline_res_mat[, x], cv = TRUE, keep.stuff = FALSE)$y}, rep(0.0, nrow(spline_res_mat))) # nolint
    mu_means <- colMeans(mu_hat_mat) # nolint
    mu_hat_mat <- vapply(1:vts_dimension, function(x) {return(mu_hat_mat[, x] - mu_means[x])}, rep(0.0, nrow(spline_res_mat))) # nolint
    slopes_hat <- vapply(1:vts_dimension, function(x) {return(mu_hat_mat[2, x] - mu_hat_mat[1, x])}, 0.0) # nolint
    mu_boundary_mat <- vapply(1:vts_dimension, function(x) {return(slopes_hat[x] * c(1:plag))}, rep(0.0, plag)) # nolint
    mu_hat_mat <- rbind(mu_boundary_mat, mu_hat_mat) # nolint
    return(mu_hat_mat)
}


seasonal_trend_smoother <- function(seasonal_residuals, vts_dimension, eval_fourier_basis, eval_basis_pre_lag) { # nolint
    seasonal_res_mat = residuals_matrix(seasonal_residuals, vts_dimension = vts_dimension) # nolint
    # t1 = plag + 1 # nolint
    # tn = nrow(seasonal_res_mat) + plag # nolint
    # time_range = c(t1:tn) # nolint
    theta_pi_hat_list <- lapply(1:vts_dimension, function(x) {return(least_squares_solver(eval_fourier_basis, seasonal_res_mat[, x]))}) # nolint
    pi_hat_mat <- vapply(1:vts_dimension, function(x) {return(rbind(eval_basis_pre_lag, eval_fourier_basis) %*% theta_pi_hat_list[[x]])}, rep(0.0, nrow(seasonal_res_mat) + nrow(eval_basis_pre_lag))) # nolint
    pi_means <- colMeans(pi_hat_mat) # nolint
    pi_hat_mat <- vapply(1:vts_dimension, function(x) {return(pi_hat_mat[, x] - pi_means[x])}, rep(0.0, nrow(seasonal_res_mat) + nrow(eval_basis_pre_lag))) # nolint
    return(pi_hat_mat)
}


build_GNAR_design <- function(vts, stages_tensor, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups) { # nolint
    # beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0) # nolint
    # stages_tensor = get_r_stages_adjacency_list(as.matrix(as.matrix(GNARtoigraph(network))), max(beta_depths)) # nolint
    cov_data_blocks = get_covariate_data_blocks(vts, covariate_groups) # nolint
    W_list = get_covariate_weight_matrices(weight_matrix, covariate_groups) # nolint
    max_lag = max(alpha_order) # nolint
    d = ncol(vts) # nolint
    covariate_models <- lapply(seq(1:length(covariate_groups)), function(x) {build_covariate_design(stages_tensor, alpha_order[x], # nolint
                                                                                                  beta_order[[x]], # nolint
                                                                                                  W_list, cov_data_blocks, x,  # nolint
                                                                                                  interaction_groups[[x]])}) # nolint: line_length_linter.
    if (length(unique(alpha_order)) != 1) { # nolint
        covariate_balanced_designs <- lapply(seq(1:length(covariate_models)), function(x) {grab_valid_observations(covariate_models[[x]], # nolint
                                                                                                             alpha_order[x],  # nolint # nolint
                                                                                                             max_lag, d)}) # nolint
    full_model <- merge_covariate_designs(covariate_balanced_designs) # nolint
    } else { # nolint
        full_model <- merge_covariate_designs(covariate_models) # nolint
    } # nolint
    return(as.matrix(full_model))
}


build_predictor_design <- function(predictor_mat, covariate_groups, vts_dimension) { # nolint
    C = length(covariate_groups) # nolint
    d = vts_dimension # nolint
    nsize = nrow(predictor_mat) # nolint
    predictor_group_blocks <- lapply(1:C, function(x) {return(rep(covariate_groups[[x]], d) * diag(nsize) %*% predictor_mat)}) # nolint
    out_mat <- recursive_dmat_rows_cpp(predictor_group_blocks) # nolint
    return(out_mat)# nolint
}


temp_trend_GNARfit <- function(longitudinal_network_data, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups, rmax, iter=100, convergence_tolreance=1e-7) { # nolint
    # set up the main variables and design matrices  # nolint: indentation_linter, line_length_linter.
    ymat <- longitudinal_network_data # nolint
    vts_dimension = ncol(ymat) # nolint
    d = vts_dimension # nolint
    ycol_means <- colMeans(ymat)
    plag = max(alpha_order) # nolint
    C = length(covariate_groups) # nolint
    q = sum(vapply(1:C, function(x) {get_community_parameter_number(alpha_order[x], beta_order[[x]], interaction_groups[[x]])}, 0.0)) # nolint
    time_steps = nrow(ymat) # nolint
    t1 = plag + 1 # nolint
    tn = time_steps # nolint
    ymat <- vapply(1:d, function(x) {return(ymat[, x] - ycol_means[x])}, rep(0.0, tn)) # nolint
    yvec <-  recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(ymat[x, ]))})) # nolint
    theta_hat =  rep(0, q) # nolint
    beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0) # nolint
    r_stage_adjacency_matrices <- get_r_stages_adjacency_list(as.matrix(network), rmax) # nolint
    t0 = d * plag + 1 # nolint
    td = tn * d # nolint
    xhat <- rep(0.0, td - t0 + 1) # nolint
    mu_hat <- rep(0.0, td)
    # run the backnet-fitting algorithm until convergence # nolint
    for (j in 1:iter) { # nolint
        # t0:tdprint(length(yvec[(d * plag + 1):(d * tn)])) # nolint
        #Rmat <- rbind(GNAR_design[, 2:(q + 1)], 0.001 * diag(q))  # nolint
        #yvec_ridge <- c(GNAR_design[, 1], rep(0, q)) # nolint
        mu_res <- yvec[t0:td] - xhat # nolint # nolint
        mu_hat_mat <- temporal_trend_spline_smoother(spline_residuals = mu_res, spar_vec = rep(0.65, d), vts_dimension = d, plag = plag) # nolint
        mu_hat <- recursive_dmat_cols_cpp(lapply(1:tn, function(x) {return(as.matrix(mu_hat_mat[x, ]))})) # nolint
        xhat <- yvec - mu_hat # nolint
        xhat_vts <- residuals_matrix(xhat, d) # nolint
        xhat_means <- colMeans(xhat_vts) # nolint
        xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
        GNAR_design <- build_GNAR_design(vts = xhat_vts, stages_tensor = r_stage_adjacency_matrices, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = weight_matrix, covariate_groups = covariate_groups, interaction_groups = interaction_groups) # nolint # nolint
        Rmat <- GNAR_design[, 2:(q + 1)] # nolint
        yvec_ridge <- GNAR_design[, 1]
        #Rmat <- rbind(GNAR_design[, 2:(q + 1)], 0.01 * diag(q))  # nolint
        #yvec_ridge <- c(GNAR_design[, 1], rep(0, q)) # nolint
        theta_hat <- least_squares_solver(Rmat, yvec_ridge) # nolint
        xhat <- GNAR_design[, 2:(q + 1)] %*% theta_hat
    } # nolint
    # theta_hat <- c(theta_hat)  nolint
    xhat <- yvec - mu_hat # nolint
    xhat_vts <- residuals_matrix(xhat, d) # nolint
    xhat_means <- colMeans(xhat_vts) # nolint
    xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
    uhat <- yvec[t0:td] - (mu_hat[t0:td] + GNAR_design[, 2:(q + 1)] %*% theta_hat)
    uhat_mat <- residuals_matrix(uhat, d)
    rownames(theta_hat) <- get_com_interaction_gnar_names(alpha_order, beta_order, length(covariate_groups), interaction_groups) # nolint
    out <- list(mu_hat_mat, xhat_vts, uhat, uhat_mat, theta_hat)
    names(out) <- c("temporal trend fitted values", "GNAR residuals", "Residuals Vector", "Residuals Matrix", "GNAR parameters") # nolint
    return(out)
}

temp_trend_cv_GNARfit <- function(longitudinal_network_data, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups, rmax, iter=100, convergence_tolreance=1e-7) { # nolint
    # set up the main variables and design matrices  # nolint: indentation_linter, line_length_linter.
    ymat <- longitudinal_network_data # nolint
    vts_dimension = ncol(ymat) # nolint
    d = vts_dimension # nolint
    ycol_means <- colMeans(ymat)
    plag = max(alpha_order) # nolint
    C = length(covariate_groups) # nolint
    q = sum(vapply(1:C, function(x) {get_community_parameter_number(alpha_order[x], beta_order[[x]], interaction_groups[[x]])}, 0.0)) # nolint
    time_steps = nrow(ymat) # nolint
    t1 = plag + 1 # nolint
    tn = time_steps # nolint
    ymat <- vapply(1:d, function(x) {return(ymat[, x] - ycol_means[x])}, rep(0.0, tn)) # nolint
    yvec <-  recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(ymat[x, ]))})) # nolint
    theta_hat =  rep(0, q) # nolint
    beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0) # nolint
    r_stage_adjacency_matrices <- get_r_stages_adjacency_list(as.matrix(network), rmax) # nolint
    t0 = d * plag + 1 # nolint
    td = tn * d # nolint
    xhat <- rep(0.0, td - t0 + 1) # nolint
    mu_hat <- rep(0.0, td)
    # run the backnet-fitting algorithm until convergence # nolint
    for (j in 1:iter) { # nolint
        # t0:tdprint(length(yvec[(d * plag + 1):(d * tn)])) # nolint
        #Rmat <- rbind(GNAR_design[, 2:(q + 1)], 0.001 * diag(q))  # nolint
        #yvec_ridge <- c(GNAR_design[, 1], rep(0, q)) # nolint
        mu_res <- yvec[t0:td] - xhat # nolint # nolint
        mu_hat_mat <- temporal_trend_spline_smoother_cv(spline_residuals = mu_res, spar_vec = NULL, vts_dimension = d, plag = plag) # nolint
        mu_hat <- recursive_dmat_cols_cpp(lapply(1:tn, function(x) {return(as.matrix(mu_hat_mat[x, ]))})) # nolint
        xhat <- yvec - mu_hat # nolint
        xhat_vts <- residuals_matrix(xhat, d) # nolint
        xhat_means <- colMeans(xhat_vts) # nolint
        xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
        GNAR_design <- build_GNAR_design(vts = xhat_vts, stages_tensor = r_stage_adjacency_matrices, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = weight_matrix, covariate_groups = covariate_groups, interaction_groups = interaction_groups) # nolint # nolint
        Rmat <- GNAR_design[, 2:(q + 1)] # nolint
        yvec_ridge <- GNAR_design[, 1]
        #Rmat <- rbind(GNAR_design[, 2:(q + 1)], 0.01 * diag(q))  # nolint
        #yvec_ridge <- c(GNAR_design[, 1], rep(0, q)) # nolint
        theta_hat <- least_squares_solver(Rmat, yvec_ridge) # nolint
        xhat <- GNAR_design[, 2:(q + 1)] %*% theta_hat
    } # nolint
    # theta_hat <- c(theta_hat)  nolint
    xhat <- yvec - mu_hat # nolint
    xhat_vts <- residuals_matrix(xhat, d) # nolint
    xhat_means <- colMeans(xhat_vts) # nolint
    xhat_vts <- vapply(1:d, function(x) {return(xhat_vts[, x] - xhat_means[x])}, rep(0.0, tn)) # nolint
    uhat <- yvec[t0:td] - (mu_hat[t0:td] + GNAR_design[, 2:(q + 1)] %*% theta_hat)
    uhat_mat <- residuals_matrix(uhat, d)
    rownames(theta_hat) <- get_com_interaction_gnar_names(alpha_order, beta_order, length(covariate_groups), interaction_groups) # nolint
    out <- list(mu_hat_mat, xhat_vts, uhat, uhat_mat, theta_hat)
    names(out) <- c("temporal trend fitted values", "GNAR residuals", "Residuals Vector", "Residuals Matrix", "GNAR parameters") # nolint
    return(out)
}

trend_forecast_fit <- function(fitted_trend_vals, trend_lags, realisation_length) { # nolint
  time_range <- seq(realisation_length - (trend_lags - 1), 
                    realisation_length, 1)
    vts_dim <- ncol(fitted_trend_vals)
    yvals_mat <- fitted_trend_vals[time_range, ]
    xvals <- cbind(rep(1, length(time_range)), time_range)
    beta_coefffs_mat <- vapply(1:vts_dim, function(x) { return(as.numeric(least_squares_solver(xvals, yvals_mat[, x])))}, rep(0.0, 2)) # nolint
    return(beta_coefffs_mat)
}

one_step_ahead_forecast <- function(beta_coeffs_mat, fitted_GNAR_term, network, alpha_order, beta_order, weight_matrix, covariate_groups, interaction_groups, theta_hat, time_val) { # nolint
    realisation_length <- nrow(fitted_GNAR_term) # nolint
    pmax <- max(alpha_order)
    time_lags <- seq(realisation_length - (pmax - 1), realisation_length, 1) # nolint
    C = length(covariate_groups) # nolint
    q = sum(vapply(1:C, function(x) {get_community_parameter_number(alpha_order[x], beta_order[[x]], interaction_groups[[x]])}, 0.0)) # nolint
    beta_depths = vapply(seq(1:length(covariate_groups)), function(x) {max(beta_order[[x]])}, 0) # nolint
    r_stage_adjacency_matrices <- get_r_stages_adjacency_list(as.matrix(network), max(beta_depths)) # nolint
    one_step_GNAR_design <- build_GNAR_design(vts = fitted_GNAR_term, stages_tensor = r_stage_adjacency_matrices, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = weight_matrix, covariate_groups = covariate_groups, interaction_groups = interaction_groups) # nolint
    daux <- nrow(one_step_GNAR_design)
    vts_dim <- ncol(fitted_GNAR_term)
    drange <- seq(daux - (vts_dim - 1), daux, 1)
    trend_pred <- t(beta_coeffs_mat) %*% c(1, time_val)
    GNAR_pred <- one_step_GNAR_design[drange, 2:(q + 1)] %*% theta_hat # nolint
    return(trend_pred + GNAR_pred)
}

struct_fit_residuals <- function(longitudinal_network_data, fitted_trend, fitted_GNAR) { # nolint
    ymat <- longitudinal_network_data # nolint
    vts_dimension = ncol(ymat) # nolint
    d = vts_dimension # nolint
    ycol_means <- colMeans(ymat)
    time_steps = nrow(ymat) # nolint
    ymat <- vapply(1:d, function(x) {return(ymat[, x] - ycol_means[x])}, rep(0.0, time_steps)) # nolint
    yvec <- recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(ymat[x, ]))})) # nolint
    muvec <- recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(fitted_trend[x, ]))})) # nolint
    xvec <- recursive_dmat_cols_cpp(lapply(1:time_steps, function(x) {return(as.matrix(fitted_GNAR[x, ]))})) # nolint
    struct_residuals <- yvec - (muvec + xvec)
    return(struct_residuals)
}