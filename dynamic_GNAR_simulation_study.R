source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/structural_GNAR_fitting_methods.r")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/struc_GNAR_sim_setup.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
# PNACF level plots
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/nacf_level_plots.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/com_GNARfit_tv_weights.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/weight_adjustment_methods.R")
require(fda)
require(fastCorbit)
Z = matrix(rnorm(50), nrow=10, ncol=5)
cov_groups = list(c(1, 0, 0, 0, 1), c(0, 1, 1, 1, 0))

Zdmat = build_predictor_design(Z, cov_groups, 5)
Zdmat[3:10, ] %*% rnorm(10)
View(t(Zdmat) %*% Zdmat)
###########################
cos_fits = seasonal_trend_smoother(rnorm(995), 5, B[1:199, ], eval.basis(runif(1), create.fourier.basis(c(-pi, pi), nbasis= 3, period = 2 * pi)))
spline_fits = temporal_trend_spline_smoother(rnorm(995), 5, 1)


y = recursive_dmat_cols_cpp(lapply(1:10, function(x) {return(as.matrix(Z[x, ]))}))

mu_vec = recursive_dmat_cols_cpp(lapply(1:200, function(x) {return(as.matrix(spline_fits[x, ]))}))
all.equal(mu_vec[1:5], spline_fits[1, ])

pi_mat <- vapply(1:5, function(x) {return(cos_fits[[x]])}, rep(0, 200))
pi_vec = recursive_dmat_cols_cpp(lapply(1:200, function(x) {return(as.matrix(cos_fits[x, ]))}))

cos_fits <- seasonal_trend_smoother(rnorm(1000), 5, 2, B[1:200, ])

gnar_design_test <- build_GNAR_design(sim0, get_r_stages_adjacency_list(as.matrix(fiveNet), 3), rep(2, 2), list(c(1, 1), c(1, 1)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(0), c(0)))
community_interaction_gnar_fit(sim0, fiveNet, rep(2, 2), list(c(1, 1), c(1, 1)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(0), c(0)))

gnar_design_test[, 2:9]

#######################
tsteps = 200
ZZ = matrix(rnorm(tsteps * 5 * 2), nrow = tsteps * 5, ncol = 2)

sim1 <- sim1 <-comm_GNARsim(tsteps, fiveNet,  weights_matrix_cpp(fiveNet, 3), 2, alpha_params, beta_params, gamma_params, cov_groups, interactions, 3)
struct_sim_test <- structural_GNARsim(GNARsim = sim1, predictor_matrix = ZZ, covariate_groups = list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), eta_coeffs = c(5, 1, 0, 1), 
                                      mu_types = c("sinc", "doppler", "linear", "sin_prod", "alt"), pi_types = c("one", "two", "three", "four", "five"),
                                      period = 7)
###########
test_fit <- structural_GNARfit(struct_sim_test, ZZ, fiveNet, 5, 7, c(1, 2), list(c(1), c(1, 1)), weights_matrix_cpp(fiveNet, 3), list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(0), c(0)), iter = 10)
print(test_fit$`fixed network effects`)
print(test_fit $`GNAR parameters`)
struc_means <- colMeans(struct_sim_test)
struct_sim_test <- vapply(1:5, function(x) {return(struct_sim_test[, x] - struc_means[x])}, rep(0, tsteps))
#####################
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/sim_ts_plots.pdf", width = 30, height = 5)
par(mfrow = c(1, 5))
plot(c(1:tsteps), struct_sim_test[, 1], 'b', lty = 2, main = "Node 1", ylab = "Longitudinal Series", pch = 16, xlab = 't', cex.lab = 1.25)
lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 1], ylab = "", xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:tsteps), test_fit$`seasonal trend fitted values`[3:tsteps, 1] + test_fit$`temporal trend fitted values`[3:tsteps, 1], ylab = '', xlab = "", type = "l", col = 'cyan3', lwd = 2)
plot(c(1:tsteps), struct_sim_test[, 2], 'b', lty = 2, main = "Node 2", ylab = '', pch = 16, xlab = '', cex.lab = 1)
lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 2], ylab = "", xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:tsteps), test_fit$`seasonal trend fitted values`[3:tsteps, 2] + test_fit$`temporal trend fitted values`[3:tsteps, 2], ylab = '', xlab = "", type = "l", col = 'cyan3', lwd = 2)
plot(c(1:tsteps), struct_sim_test[, 3], 'b', lty = 2, main = "Node 3", ylab = '', pch = 16, xlab = '', cex.lab = 2)
lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 3], ylab = "", xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:tsteps), test_fit$`seasonal trend fitted values`[3:tsteps, 3] + test_fit$`temporal trend fitted values`[3:tsteps, 3], ylab = '', xlab = "", type = "l", col = 'cyan3', lwd = 2)
plot(c(1:tsteps), struct_sim_test[, 4], 'b', lty = 2, main = "Node 4", ylab = '', pch = 16, xlab = '', cex.lab = 2)
lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 4], ylab = "", xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:tsteps), test_fit$`seasonal trend fitted values`[3:tsteps, 4] + test_fit$`temporal trend fitted values`[3:tsteps, 4], ylab = '', xlab = "", type = "l", col = 'cyan3', lwd = 2)
plot(c(1:tsteps), struct_sim_test[, 5], 'b', lty = 2, main = "Node 5", ylab = '', pch = 16, xlab = '', cex.lab = 2)
lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 5], ylab = "", xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:tsteps), test_fit$`seasonal trend fitted values`[3:tsteps, 5] + test_fit$`temporal trend fitted values`[3:tsteps, 5], ylab = '', xlab = "", type = "l", col = 'cyan3', lwd = 2)
dev.off()
#par(mfrow = c(1, 1))
#par(mfrow = c(3, 5))
###########################
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/sim_trend_fit_plots.pdf", width = 22, height = 4)
par(mfrow = c(1, 5))
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 1], ylab = "Trend Decomposition", xlab = "", type = "l", col = 'blue', lwd = 2, cex.lab = 1.25)
lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 1], ylab = '', xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 2], ylab = '', xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 2], ylab = '', xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 3], ylab = '', xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 3], ylab = '', xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 4], ylab = '', xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 4], ylab ='', xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 5], ylab = '', xlab = "", type = "l", col = 'blue', lwd = 2)
lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 5], ylab = '', xlab = "", type = "l", col = 'cyan', lwd = 2)
dev.off()
#par(mfrow = c(1, 1))
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/sim_seasonal_fit_plots.pdf", width = 22, height = 4)
par(mfrow = c(1, 5))
end_point = 7 * 10
aux_basis <- create.fourier.basis(c(3, 17), nbasis = 5, period = 7)
eval_seasonal_points <- eval.basis(seq(3, 17, 0.1), aux_basis)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 1], ylab = expression(hat(pi)["1, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 2], ylab = expression(hat(pi)["2, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 3], ylab = expression(hat(pi)["3, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 4], ylab = expression(hat(pi)["4, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 5], ylab = expression(hat(pi)["5, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
dev.off()
######################
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/sim_GNAR_fit_plots.pdf", width = 30, height = 5)
par(mfrow = c(1, 5))
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 1], ylab = "GNAR Term", xlab = "t", type = "b", 
     col = 'forestgreen', lwd = 1, pch = 16, cex.lab = 1.25)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 2], type = 'b', ylab = '', xlab = "t",
     col = 'forestgreen', lwd = 1, pch = 16)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 3], type = 'b', ylab = '', xlab = "t",
     col = 'forestgreen', lwd = 1, pch = 16)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 4], type = 'b', ylab = '', xlab = "t",
     col = 'forestgreen', lwd = 1, pch = 16)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 5], type = 'b', ylab = '', xlab = "t",
     col = 'forestgreen', lwd = 1, pch = 16)
#####################
Zmat1 <- residuals_matrix(ZZ[, 1], 5)
Zmat2 <- residuals_matrix(ZZ[, 2], 5)
#par(mfrow = c(2, 3))
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/sim_adj_response_plots.pdf", width = 30, height = 5)
par(mfrow = c(1, 5))
plot(Zmat1[, 1],  struct_sim_test[, 1] - (test_fit$`temporal trend fitted values`[, 1] + test_fit$`seasonal trend fitted values`[, 1] + test_fit$`GNAR residuals`[, 1]), pch=16,
     xlab=expression(v["t"]), ylab="Adjusted Response", cex.lab = 1.25)+
  lines(Zmat1[, 1], test_fit$`fixed network effects`[3, 1] * Zmat1[, 1], lwd = 3, col = "orange")
plot(Zmat1[, 2],  struct_sim_test[, 2] - (test_fit$`temporal trend fitted values`[, 2] + test_fit$`seasonal trend fitted values`[, 2] + test_fit$`GNAR residuals`[, 2]), pch=16,
     xlab=expression(v["t"]), ylab='')+
  lines(Zmat1[, 2], test_fit$`fixed network effects`[1, 1] * Zmat1[, 2], lwd = 3, col = "orange")
plot(Zmat1[, 3],  struct_sim_test[, 3] - (test_fit$`temporal trend fitted values`[, 3] + test_fit$`seasonal trend fitted values`[, 3] + test_fit$`GNAR residuals`[, 3]), pch=16,
     xlab=expression(v["t"]), ylab='')+
  lines(Zmat1[, 3], test_fit$`fixed network effects`[1, 1] * Zmat1[, 3], lwd = 3, col = "orange")
plot(Zmat1[, 4],  struct_sim_test[, 4] - (test_fit$`temporal trend fitted values`[, 4] + test_fit$`seasonal trend fitted values`[, 4] + test_fit$`GNAR residuals`[, 4]), pch=16,
     xlab=expression(v["t"]), ylab='')+
  lines(Zmat1[, 4], test_fit$`fixed network effects`[1, 1] * Zmat1[, 4], lwd = 3, col = "orange")
plot(Zmat1[, 5],  struct_sim_test[, 5] - (test_fit$`temporal trend fitted values`[, 5] + test_fit$`seasonal trend fitted values`[, 5] + test_fit$`GNAR residuals`[, 5]), pch=16,
     xlab=expression(v["t"]), ylab='')+
  lines(Zmat1[, 5], test_fit$`fixed network effects`[3, 1] * Zmat1[, 5], lwd = 3, col = "orange")
dev.off()
par(mfrow = c(1, 1))
#################
#################
par(mfrow = c(2, 3))
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 1], 'l')
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 2], 'l')
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 3], 'l')
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 4], 'l')
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 5], 'l')
par(mfrow = c(1, 1))
gnar_com_diff_order <- community_gnar_fit(test_fit$`GNAR residuals`, fiveNet, c(1, 2), list(c(1), c(1, 1)),  weights_matrix_cpp(fiveNet, 3), list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
summary(gnar_com_diff_order)
XX = build_GNAR_design(sim1, get_r_stages_adjacency_list(as.matrix(fiveNet), 3), c(1, 2), list(c(1), c(1, 1)),  weights_matrix_cpp(fiveNet, 3), list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(0), c(0)))
least_squares_solver(XX[, 2:7], XX[, 1])
Xhat = build_GNAR_design(test_fit$`GNAR residuals`, get_r_stages_adjacency_list(as.matrix(fiveNet), 3), c(1, 2), list(c(1), c(1, 1)),  weights_matrix_cpp(fiveNet, 3), list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)), list(c(0), c(0)))
least_squares_solver(Xhat[, 2:7], Xhat[, 1])

fixed_effects <- vapply(1:5, function(x) {ifelse(x==1 || x==5, return(cbind(Zmat1[, x], Zmat2[, x]) %*% test_fit$`fixed network effects`[3:4, ]), 
                                                 return(cbind(Zmat1[, x], Zmat2[, x]) %*% test_fit$`fixed network effects`[1:2, ]))}, rep(0.0, tsteps))
end_point = tsteps
struc_means <- colMeans(struct_sim_test)
struct_sim_test <- vapply(1:5, function(x) {return(struct_sim_test[, x] - struc_means[x])}, rep(0, tsteps))
par(mfrow=c(2, 3))
plot(fiveNet)
for (l in 1:5){
    plot(struct_sim_test[, l], ylab = expression(Y["i, t"]), xlab = "t", main = paste0("Time Series: i = ", as.character(l)), pch=16)+
    lines(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, l], type = "l", col = 'blue', lwd = 3)+
    lines(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, l] + test_fit$`temporal trend fitted values`[3:tsteps, l], type = "l", col = 'forestgreen', lwd = 2)+
    points(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, l] + test_fit$`temporal trend fitted values`[3:tsteps, l] + fixed_effects[3:tsteps, l], pch=17, col = 'orange', lwd = 2)
  #lines(c(3:length(test_fit$`GNAR residuals`[, 1])), test_fit$`seasonal trend fitted values`[3:end_point, l] + test_fit$`temporal trend fitted values`[3:tsteps, l] + test_fit$`GNAR residuals`[3:tsteps, l], type = "l", col = 'forestgreen', lwd = 2)
}
par(mfrow=c(1, 1))
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 1], ylab = expression(mu["i, t"]), xlab = "", type = "l", col = 'blue', lwd = 3)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 2], ylab = expression(mu["i, t"]), xlab = "", type = "l", col = 'blue', lwd = 3)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 3], ylab = expression(mu["i, t"]), xlab = "", type = "l", col = 'blue', lwd = 3)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 4], ylab = expression(mu["i, t"]), xlab = "", type = "l", col = 'blue', lwd = 3)
plot(c(3:tsteps), test_fit$`temporal trend fitted values`[3:tsteps, 5], ylab = expression(mu["i, t"]), xlab = "", type = "l", col = 'blue', lwd = 3)
###
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 1], ylab = expression(pi["i, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 2], ylab = expression(pi["i, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 3], ylab = expression(pi["i, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 4], ylab = expression(pi["i, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
plot(c(3:end_point), test_fit$`seasonal trend fitted values`[3:end_point, 5], ylab = expression(pi["i, t"]), xlab = "", type = "l", col = 'cyan', lwd = 2)
###
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 1], ylab = expression(hat(x)["1, t"]), xlab = "t", type = "l", col = 'forestgreen', lwd = 1)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 2], ylab = expression(hat(x)["2, t"]), xlab = "t", type = "l", col = 'forestgreen', lwd = 1)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 3], ylab = expression(hat(x)["3, t"]), xlab = "t", type = "l", col = 'forestgreen', lwd = 1)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 4], ylab = expression(hat(x)["4, t"]), xlab = "t", type = "l", col = 'forestgreen', lwd = 1)
plot(c(1:length(test_fit$`GNAR residuals`[, 1])), test_fit$`GNAR residuals`[, 5], ylab = expression(hat(x)["5, t"]), xlab = "t", type = "l", col = 'forestgreen', lwd = 1)
####
plot(Zmat2[, 1],  struct_sim_test[, 1] - (test_fit$`temporal trend fitted values`[, 1] + test_fit$`seasonal trend fitted values`[, 1] + test_fit$`GNAR residuals`[, 1]), pch=16,
     xlab='', ylab=expression(hat(u)["1, t"]))+
  lines(Zmat2[, 1], test_fit$`fixed network effects`[3, 1] * Zmat2[, 1], lwd = 3, col = "orange")
plot(Zmat1[, 2],  struct_sim_test[, 2] - (mean(struct_sim_test[, 2]) + test_fit$`temporal trend fitted values`[, 2] + test_fit$`seasonal trend fitted values`[, 2] + test_fit$`GNAR residuals`[, 2]), pch=16,
     xlab=expression(z["t"]), ylab=expression(hat(u)["2, t"]))+
  lines(Zmat1[, 2], test_fit$`fixed network effects`[1, 1] * Zmat1[, 2], lwd = 3, col = "orange")
plot(Zmat[, 3],  struct_sim_test[, 3] - (mean(struct_sim_test[, 3]) + test_fit$`temporal trend fitted values`[, 3] + test_fit$`seasonal trend fitted values`[, 3] + test_fit$`GNAR residuals`[, 3]), pch=16,
     xlab=expression(z["t"]), ylab=expression(hat(u)["2, t"]))+
  lines(Zmat[, 3], test_fit$`fixed network effects`[1, 1] * Zmat[, 3], lwd = 3, col = "orange")
plot(Zmat[, 4],  struct_sim_test[, 4] - (mean(struct_sim_test[, 4]) + test_fit$`temporal trend fitted values`[, 4] + test_fit$`seasonal trend fitted values`[, 4] + test_fit$`GNAR residuals`[, 4]), pch=16,
     xlab=expression(z["t"]), ylab=expression(hat(u)["2, t"]))+
  lines(Zmat[, 4], test_fit$`fixed network effects`[1, 1] * Zmat[, 4], lwd = 3, col = "orange")
plot(Zmat[, 5],  struct_sim_test[, 5] - (mean(struct_sim_test[, 5]) + test_fit$`temporal trend fitted values`[, 5] + test_fit$`seasonal trend fitted values`[, 5] + test_fit$`GNAR residuals`[, 5]), pch=16,
     xlab=expression(z["t"]), ylab=expression(hat(u)["2, t"]))+
  lines(Zmat[, 5], test_fit$`fixed network effects`[3, 1] * Zmat[, 5], lwd = 3, col = "orange")
###
par(mfrow=c(1, 1))

test_fit$`fixed network effects`

