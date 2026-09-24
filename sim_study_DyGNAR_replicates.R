#' Replicate Structural GNAR Simulation Experiment Across Time Steps
#'
#' @param tsteps_vec       Integer vector of time step values to iterate over (default: c(50, 100, 500, 1000))
#' @param net              Network object (default: fiveNet)
#' @param alpha_params     Alpha parameters for comm_GNARsim
#' @param beta_params      Beta parameters for comm_GNARsim
#' @param gamma_params     Gamma parameters for comm_GNARsim
#' @param cov_groups       Covariate groups for comm_GNARsim
#' @param interactions     Interactions argument for comm_GNARsim
#' @param real_phi         Named numeric vector of true fixed-effect (phi) coefficients
#' @param real_theta       Named numeric vector of true GNAR (theta) coefficients
#' @param eta_coeffs       Eta coefficients passed to structural_GNARsim
#' @param mu_types         Mu type strings passed to structural_GNARsim
#' @param pi_types         Pi type strings passed to structural_GNARsim
#' @param period           Period argument passed to structural_GNARsim and structural_GNARfit
#' @param covariate_groups Covariate groups passed to structural_GNARsim and structural_GNARfit
#' @param alphaOrder       Alpha order passed to structural_GNARfit (default: c(1, 2))
#' @param betaLags         Beta lags list passed to structural_GNARfit (default: list(c(1), c(1, 1)))
#' @param eta_groups       Eta groups list passed to structural_GNARfit (default: list(c(2), c(1)))
#' @param iter             Number of iterations for structural_GNARfit (default: 10)
#' @param n_nodes          Number of nodes in the network (default: 5)
#' @param max_stage        Maximum stage for weights and corbit (default: 3)
#' @param nacf_lags        Number of lags for corbit_plot_cpp (default: 10)
#' @param level            Confidence level for check_coverage (default: 0.95)
#' @param seed             Random seed for reproducibility (default: NULL)
#'
#' @return A data frame with one row per tsteps value and columns:
#'         tsteps, rmse_response, rmse_theta, rmse_phi,
#'         fraction_theta_covered, fraction_phi_covered, nacf_max_val
#'         
#'
################
################
run_experiment <- function(
    tsteps_vec      = c(50, 100, 500, 1000),
    net             = fiveNet,
    alpha_params,
    beta_params,
    gamma_params,
    cov_groups,
    interactions,
    real_phi,
    real_theta,
    eta_coeffs,
    mu_types,
    pi_types,
    period,
    covariate_groups,
    alphaOrder      = c(1, 2),
    betaLags        = list(c(1), c(1, 1)),
    eta_groups      = list(c(2), c(1)),
    iter            = 30,
    n_nodes         = 5,
    max_stage       = 3,
    nacf_lags       = 10,
    level_theta     = 1 - 0.05 / 9,
    level_phi       = 1 - 0.05 / 4,
    seed            = NULL
) {
  
  # --- Input Validation ---
  if (!is.numeric(tsteps_vec) || any(tsteps_vec <= 0)) {
    stop("`tsteps_vec` must be a vector of positive integers.")
  }
  if (!is.numeric(real_phi) || !is.numeric(real_theta)) {
    stop("`real_phi` and `real_theta` must be numeric vectors.")
  }
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Precompute weight matrix (shared across all tsteps)
  W <- weights_matrix_cpp(net, max_stage)
  
  # --- Iterate Over tsteps Values ---
  results <- lapply(tsteps_vec, function(tsteps) {
    
    message(sprintf("Running experiment for tsteps = %d ...", tsteps))
    
    # 1. Build predictor matrix
    ZZ <- matrix(
      rnorm(tsteps * n_nodes * 2),
      nrow = tsteps * n_nodes,
      ncol = 2
    )
    
    # 2. Simulate community GNAR process
    sim1 <- comm_GNARsim(
      tsteps, net, W, 2,
      alpha_params, beta_params, gamma_params,
      cov_groups, interactions, max_stage
    )
    
    # 3. Simulate structural GNAR process
    struct_sim <- structural_GNARsim(
      GNARsim          = sim1,
      predictor_matrix = ZZ,
      covariate_groups = covariate_groups,
      eta_coeffs       = eta_coeffs,
      mu_types         = mu_types,
      pi_types         = pi_types,
      period           = period
    )
    
    # 4. Fit structural GNAR model
    test_fit <- structural_GNARfit(
      struct_sim, ZZ, net, n_nodes, period,
      alphaOrder, betaLags, W,
      covariate_groups, eta_groups,
      iter = iter
    )
    
    # 5. Fixed effects (phi) regression
    fixed_effects_fit <- lm(
      test_fit$adjusted_response ~ test_fit$predictor_design + 0
    )
    
    # 6. GNAR (theta) regression
    GNAR_fit <- lm(
      test_fit$GNAR_response ~ test_fit$GNAR_design + 0
    )
    
    # 7. Compute metrics
    rmse_response <- sqrt(
      sum(fixed_effects_fit$residuals^2) / length(fixed_effects_fit$residuals)
    )
    
    fraction_phi_covered   <- check_coverage(fixed_effects_fit, real_phi,  level = level_phi)
    fraction_theta_covered <- check_coverage(GNAR_fit,          real_theta, level = level_theta)
    
    rmse_theta <- sqrt(
      #(tsteps - 2) * 2 *
        sum((GNAR_fit$coefficients - real_theta)^2)
    )
    
    rmse_phi <- sqrt(
      #(tsteps - 2) * 2 *
        sum((fixed_effects_fit$coefficients - real_phi)^2)
    )
    
    # 8. Compute NACF and extract max absolute value
    nacf_vals    <- corbit_plot_cpp(
      residuals_matrix(fixed_effects_fit$residuals, n_nodes),
      net, nacf_lags, max_stage, W
    )
    nacf_max_val <- max(abs(nacf_vals))
    
    # 9. Return named result vector
    c(
      tsteps                 = tsteps,
      rmse_response          = rmse_response,
      rmse_theta             = rmse_theta,
      rmse_phi               = rmse_phi,
      fraction_theta_covered = fraction_theta_covered,
      fraction_phi_covered   = fraction_phi_covered,
      nacf_max_val           = nacf_max_val
    )
  })
  
  # --- Combine Into Data Frame ---
  result_df <- as.data.frame(do.call(rbind, results))
  
  message("Done. Returning results data frame.")
  return(result_df)
}

#' Replicate Structural GNAR Experiment and Summarise Results
#'
#' Runs run_experiment() n_reps times for each value of tsteps and returns
#' a summary table of means and standard deviations for all metrics.
#'
#' @param n_reps        Number of replications per tsteps value (default: 1000)
#' @param tsteps_vec    Integer vector of time step values (default: c(50, 100, 500, 1000))
#' @param seed          Base random seed for reproducibility (default: NULL)
#' @param ...           Additional arguments forwarded to run_experiment()
#'
#' @return A data frame with one row per tsteps value containing the mean and
#'         standard deviation of each metric across all replications
replicate_experiment <- function(
    n_reps     = 1000,
    tsteps_vec = c(50, 100, 500, 1000),
    seed       = NULL,
    ...
) {
  
  # --- Input Validation ---
  if (!is.numeric(n_reps) || n_reps < 1) {
    stop("`n_reps` must be a positive integer.")
  }
  if (!is.numeric(tsteps_vec) || any(tsteps_vec <= 0)) {
    stop("`tsteps_vec` must be a vector of positive integers.")
  }
  
  # Metric names (excluding tsteps, which is fixed per group)
  metric_names <- c(
    "rmse_response",
    "rmse_theta",
    "rmse_phi",
    "fraction_theta_covered",
    "fraction_phi_covered",
    "nacf_max_val"
  )
  
  # --- Iterate Over Each tsteps Value ---
  summary_rows <- lapply(tsteps_vec, function(tsteps) {
    
    message(sprintf(
      "\n===== Replicating experiment for tsteps = %d (%d reps) =====",
      tsteps, n_reps
    ))
    
    # Run n_reps replications for this single tsteps value
    reps <- lapply(seq_len(n_reps), function(rep_i) {
      
      # Derive a unique seed per replication if a base seed is provided
      if (!is.null(seed)) set.seed(seed + rep_i)
      
      if (rep_i %% 100 == 0) {
        message(sprintf("  tsteps = %4d | rep %d / %d", tsteps, rep_i, n_reps))
      }
      
      # Run experiment for a single tsteps value, suppress per-step messages
      suppressMessages(
        run_experiment(tsteps_vec = tsteps, ...)
      )
    })
    
    # Stack all replications into one matrix (n_reps x n_metrics)
    rep_matrix <- do.call(rbind, reps)
    
    # Drop the tsteps column before computing statistics
    metrics <- rep_matrix[, metric_names, drop = FALSE]
    
    # Compute mean and SD for each metric
    means <- colMeans(metrics)
    sds   <- apply(metrics, 2, sd)
    
    # Build named result row
    row <- c(tsteps = tsteps)
    for (m in metric_names) {
      row[paste0("mean_", m)] <- means[m]
      row[paste0("sd_",   m)] <- sds[m]
    }
    
    row
  })
  
  # --- Combine Into Summary Data Frame ---
  summary_df <- as.data.frame(do.call(rbind, summary_rows))
  
  # Reorder columns: tsteps | mean_x, sd_x | mean_y, sd_y | ...
  ordered_cols <- c(
    "tsteps",
    as.vector(outer(c("mean_", "sd_"), metric_names, paste0))
  )
  summary_df <- summary_df[, ordered_cols]
  
  message("\nAll replications complete. Returning summary table.")
  return(summary_df)
}
