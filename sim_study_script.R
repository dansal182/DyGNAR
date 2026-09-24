tsteps = 50
ZZ = matrix(rnorm(tsteps * 5 * 2), nrow = tsteps * 5, ncol = 2)
sim1 <-comm_GNARsim(tsteps, fiveNet,  weights_matrix_cpp(fiveNet, 3), 2, alpha_params, beta_params, gamma_params, cov_groups, interactions, 3)
struct_sim_test <- structural_GNARsim(GNARsim = sim1, predictor_matrix = ZZ, covariate_groups = list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), eta_coeffs = c(5, 1, 0, 0), 
                                      mu_types = c("sinc", "doppler", "linear", "sin_prod", "alt"), pi_types = c("one", "two", "three", "four", "five"),
                                        period = 7)
test_fit <- structural_GNARfit(struct_sim_test, ZZ, fiveNet, 5, 7, c(1, 2), list(c(1), c(1, 1)), weights_matrix_cpp(fiveNet, 3), 
                                 list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(2), c(1)), iter = 10)
fixed_effects_fit <- lm(test_fit$adjusted_response ~ test_fit$predictor_design + 0)
rmse_response = sqrt(sum(fixed_effects_fit$residuals^2) / length(fixed_effects_fit$residuals))
fraction_phi_covered = check_coverage(fixed_effects_fit, real_phi, level = 0.95)
GNAR_fit <- lm(test_fit$GNAR_response ~ test_fit$GNAR_design + 0)
fraction_theta_covered <- check_coverage(GNAR_fit, real_theta, level = 0.95)
rmse_theta <- sqrt(length(GNAR_fit$reisduals) * sum(GNAR_fit$coefficients - real_theta)^2)
rmse_phi <- sqrt( length(fixed_effects_fit$residuals) * sum(fixed_effects_fit$coefficients - real_phi)^2)
nacf_vals <- corbit_plot_cpp(residuals_matrix(fixed_effects_fit$residuals, 5), fiveNet, 10, 3, weights_matrix_cpp(fiveNet, 3))
nacf_max_val <- max(abs(nacf_vals))
c(tsteps, rmse_response, rmse_theta, rmse_phi, fraction_theta_covered, fraction_phi_covered, nacf_max_val)
###################
# experiment ends here
#####################
alpha_params <- list(c(0.23), c(0.29, 0.35))
beta_params <- list(list(c(0.47)), list(c(0.15), c(0.17)))
gamma_params <- list(list(list(c(0.0, 0.12))), 
                     list(list(c(0.13, 0.0)), list(c(0.14, 0.0)))
)
real_theta <- c(0.23, 0.47, 0.12, 0.29, 0.15, 0.13, 0.35, 0.17, 0.14)
real_phi <- c(5, 1, 0, 0)
cov_groups <- list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1))
interactions <- list(c(2), c(1))
var_mats <- build_phi_mats(fiveNet,  weights_matrix_cpp(fiveNet, 3), 2, alpha_params, beta_params, gamma_params, cov_groups, interactions, 3)
var_mats
test_sim <- comm_GNARsim(100, fiveNet,  weights_matrix_cpp(fiveNet, 3), 2, alpha_params, beta_params, gamma_params, cov_groups, interactions, 3)

alpha_test <- c(1, 2)
beta_test <- list(c(1), c(1, 1))
test_int_fit <- community_interaction_gnar_fit(test_fit$GNAR_residuals, fiveNet, alpha_test, beta_test, weights_matrix_cpp(fiveNet, 3), cov_groups, 
                                                      interactions)
summary(test_int_fit)
Zmat <- build_predictor_design(ZZ, cov_groups, 5)
test_res <- c(recursive_dmat_cols_cpp(lapply(1:tsteps, function(x) {return(as.matrix(struct_sim_test[x, ]))}))) - (c(recursive_dmat_cols_cpp(lapply(1:tsteps, function(x) {return(as.matrix(test_fit$`temporal trend fitted values`[x, ]))}))) 
                                                                                                                   + c(recursive_dmat_cols_cpp(lapply(1:tsteps, function(x) {return(as.matrix(test_fit$`seasonal trend fitted values`[x, ]))})))
                                                                                                                   + Zmat %*% test_fit$`fixed network effects` 
                                                                                                                   #+ c(recursive_dmat_cols_cpp(lapply(1:tsteps, function(x) {return(as.matrix(test_fit$`GNAR residuals`[x, ]))})))
                                                                                                                   )
hist(test_res, breaks = 10)
nacf_vals <- corbit_plot_cpp(residuals_matrix(test_res, 5), fiveNet, 10, 3, weights_matrix_cpp(fiveNet, 3))
max(abs(nacf_vals))
hist(test_int_fit$residuals, breaks = 30)
nacf_vals <- corbit_plot_cpp(residuals_matrix(test_int_fit$residuals, 5), fiveNet, 10, 3, weights_matrix_cpp(fiveNet, 3))
max(abs(nacf_vals))
test_int_fit <- community_interaction_gnar_fit(test_sim, fiveNet, alpha_test, beta_test, weights_matrix_cpp(fiveNet, 3), cov_groups, 
                                               interactions)

nacf_vals <- corbit_plot_cpp(test_sim, fiveNet, 10, 3, weights_matrix_cpp(fiveNet, 3), partial = 'yes')
max(abs(nacf_vals))


check_coverage <- function(model, real, level = 0.95) {
  ci <- as.data.frame(confint(model, level = level))
  colnames(ci) <- c("lower", "upper")
  names(real) <- rownames(ci)
  # --- Coverage Check ---
  covered <- real >= ci$lower & real <= ci$upper
  # --- Return fraction covered ---
  return(mean(covered))
}
nacf_vals <- corbit_plot_cpp(residuals_matrix(sim1, 5), fiveNet, 10, 3, weights_matrix_cpp(fiveNet, 3))
max(abs(nacf_vals))

#####################
#####################
# Replicate experiments start here
results <- run_experiment(
  tsteps_vec      = c(50, 100, 500),
  net             = fiveNet,
  alpha_params    = alpha_params,
  beta_params     = beta_params,
  gamma_params    = gamma_params,
  cov_groups      = cov_groups,
  interactions    = interactions,
  real_phi        = real_phi,
  real_theta      = real_theta,
  eta_coeffs      = c(5, 1, 0, 0),
  mu_types        = c("sinc", "doppler", "linear", "sin_prod", "alt"),
  pi_types        = c("one", "two", "three", "four", "five"),
  period          = 7,
  covariate_groups = list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1))
)

print(results)

summary_table <- replicate_experiment(
  n_reps           = 1000,
  tsteps_vec       = c(50, 100, 150, 500, 1000),
  seed             = 2026,
  # --- run_experiment() arguments below ---
  net              = fiveNet,
  alpha_params     = alpha_params,
  beta_params      = beta_params,
  gamma_params     = gamma_params,
  cov_groups       = cov_groups,
  interactions     = interactions,
  real_phi         = real_phi,
  real_theta       = real_theta,
  eta_coeffs       = c(5, 1, 0, 0),
  mu_types         = c("sinc", "doppler", "linear", "sin_prod", "alt"),
  pi_types         = c("one", "two", "three", "four", "five"),
  period           = 7,
  covariate_groups = list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1))
)

print(summary_table)

