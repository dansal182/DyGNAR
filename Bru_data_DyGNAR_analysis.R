Bru_data <- load("~/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/Bru_data.RData")
library(GNAR)
library(igraph)
library(Rcpp)
library(RcppArmadillo)
library(fastCorbit)
library(ggplot2)
library(viridis)
sourceCpp('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/nacf_inner.cpp')
sourceCpp('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/get_k_stage_adjacency_list.cpp')
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/nacf_pnacf_rwrapper.r')
source('/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/weight_adjustment_methods.R')
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
#######################
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/structural_GNAR_fitting_methods.r")
#####################
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/struc_GNAR_sim_setup.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/one-stepcommunity_GNAR_prediction.r", encoding = "UTF-8")
require(fda)
# Compute the network adjacency matrix
S1 <- as.matrix(GNARtoigraph(bur_netGNAR))
# Find the longest, shortest path and rmax
D = distances(graph_from_adjacency_matrix(S1))
D[!is.finite(D)] <- 0
rmax = max(D)
locind_edt <- c(rep(1, 10), rep(2, 10), rep(3, 17), rep(5, 23), rep(6, 11), rep(4, 15), rep(6, 12))
# Set community membership values
loc1_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==1, 1, 0))}, 0)
loc1_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==2, 1, 0))}, 0)
loc2_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==3, 1, 0))}, 0)
loc2_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==4, 1, 0))}, 0)
loc3_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==5, 1, 0))}, 0)
loc3_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==6, 1, 0))}, 0)

hist(recursive_dmat_cols_cpp(lapply(1:43, function(j) {return(as.matrix(bru_vts[2:44, j]))})), breaks = 50, xlab = 'y', main ='Bru VTS')

plot(MXD2[, 78], type = 'l', xlab="t", ylab=expression(mu["i, t"]))
dev.off()
pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/bur_netGNAR_series_plots.pdf")
par(mfrow=c(10, 10), mar=c(1, 1, 1, 1))
for (i in 1:98) {
  plot(MXD2[, i], type='l', xlab='t', ylab=expression(expression(X["i, t"])))
}
par(mfrow=c(1, 1))
dev.off()


bru_data_community_indicators <- list(loc1_sat1, loc1_sat2, loc2_sat1, loc2_sat2, loc3_sat1, loc3_sat2)
W_bru_data <- fastCorbit::weights_matrix_cpp(bur_netGNAR, rmax)
bru_data_interactions <- list(c(2:6), c(1, 3:6), c(1,2, 4:6), c(1, 2, 3, 5, 6), c(1, 2, 3, 4, 6), c(1:5))
alpha_order <- rep(1, 6)
beta_order <- list(c(2), c(2), c(2), c(2), c(3), c(3)) # one for each community

bru_struct_GNARfit <- temp_trend_GNARfit(longitudinal_network_data = 100.0 * MXD2, network = bur_netGNAR, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = W_bru_data, 
                                         covariate_groups = bru_data_community_indicators, interaction_groups = bru_data_interactions, rmax = rmax, iter = 100)
summary(bru_struct_GNARfit)
print(bru_struct_GNARfit$`GNAR parameters`)
sum(bru_struct_GNARfit$`GNAR residuals`[2:44, ]^2)
sum(bru_struct_GNARfit$`Residuals Vector`^2) / (98 * 43 - 90)
sum(abs(bru_struct_GNARfit$`Residuals Vector`)) / (98 * 43 - 90)
sum((bru_vts[2:44, ] - bru_struct_GNARfit$`temporal trend fitted values`[2:44, ])^2)
sum((bru_vts[2:44, ] - bru_struct_GNARfit$`GNAR residuals`[2:44, ])^2)
bru_struct_fit_residuals <- struct_fit_residuals(100 * MXD2, bru_struct_GNARfit$`temporal trend fitted values`, bru_struct_GNARfit$`GNAR residuals`)
bru_fit_res_ts <- residuals_matrix(bru_struct_fit_residuals, 98)
bru_res_alt <- vapply(2:44, function(j) {bru_vts[j, ] - (bru_struct_GNARfit$`temporal trend fitted values`[j, ] + bru_struct_GNARfit$`GNAR residuals`[j, ])}, rep(0, 98))
sum(bru_res_alt^2)
hist(bru_struct_GNARfit$`Residuals Vector`, breaks = 50)
GNAR_AIC(bru_struct_GNARfit$`Residuals Vector`, 98, 34)
GNAR_BIC(bru_struct_GNARfit$`Residuals Vector`, 98, 34)
GNAR_HQC(bru_struct_GNARfit$`Residuals Vector`, 98, 34)
nacf_vals <- corbit_plot_cpp(bru_struct_GNARfit$`Residuals Matrix`, bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))


bru_struct_GNARfit_cv <- temp_trend_cv_GNARfit(longitudinal_network_data = 100.0 * MXD2, network = bur_netGNAR, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = W_bru_data, 
                                         covariate_groups = bru_data_community_indicators, interaction_groups = bru_data_interactions, rmax = rmax, iter = 100)
summary(bru_struct_GNARfit_cv)
print(bru_struct_GNARfit_cv$`GNAR parameters`)
sum(bru_struct_GNARfit_cv$`GNAR residuals`[2:44, ]^2)
sum(bru_struct_GNARfit_cv$`Residuals Vector`^2)
sum(abs(bru_struct_GNARfit_cv$`Residuals Vector`)) / (98 * 43 - 90)
sum((bru_vts[2:44, ] - bru_struct_GNARfit_cv$`temporal trend fitted values`[2:44, ])^2)
sum((bru_vts[2:44, ] - bru_struct_GNARfit_cv$`GNAR residuals`[2:44, ])^2)
bru_struct_fit_residuals_cv <- struct_fit_residuals(100 * MXD2, bru_struct_GNARfit_cv$`temporal trend fitted values`, bru_struct_GNARfit_cv$`GNAR residuals`)
bru_fit_res_ts_cv <- residuals_matrix(bru_struct_fit_residuals_cv, 98)
hist(bru_struct_GNARfit_cv$`Residuals Vector`, breaks = 50)
hist(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(bru_struct_GNARfit_cv$`Residuals Matrix`[1:43, j]))})), breaks = 50)
GNAR_AIC(bru_struct_GNARfit_cv$`Residuals Vector`, 98, 90)
GNAR_BIC(bru_struct_GNARfit_cv$`Residuals Vector`, 98, 90)
GNAR_HQC(bru_struct_GNARfit_cv$`Residuals Vector`, 98, 90)
nacf_vals <- corbit_plot_cpp(bru_struct_GNARfit_cv$`Residuals Matrix`, bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
length(bru_struct_GNARfit_cv$`GNAR parameters`)
sample_kurtosis(bru_struct_GNARfit_cv$`Residuals Vector`)
sample_skewness(bru_struct_GNARfit_cv$`Residuals Vector`)
sample_kurtosis(bru_struct_GNARfit_cv$`Residuals Matrix`)
sample_skewness(bru_struct_GNARfit_cv$`Residuals Matrix`)
hist(sample_kurtosis(bru_struct_GNARfit_cv$`Residuals Matrix`))
hist(sample_skewness(bru_struct_GNARfit_cv$`Residuals Matrix`))


bru_res_global_fit <- global_gnar_fit(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR, 1, 3, W_bru_data)
summary(bru_res_global_fit)
bru_res_global_fit_ts = residuals_matrix(bru_res_global_fit$residuals, 98)
hist(bru_res_global_fit$residuals, breaks = 50)
quantile(bru_res_global_fit$residuals, probs = c(0.0275, 0.975))
GNAR_AIC(bru_res_global_fit$residuals, 98, 2)
GNAR_BIC(bru_res_global_fit$residuals, 98, 2)
GNAR_HQC(bru_res_global_fit$residuals, 98, 2)

bru_res_interaction_fit <- community_interaction_gnar_fit(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR, alpha_order, beta_order, W_bru_data, bru_data_community_indicators, 
                                                          bru_data_interactions)
summary(bru_res_interaction_fit)
sd(bru_res_interaction_fit$residuals)
hist(bru_res_interaction_fit$residuals, breaks = 50)
sum((bru_res_interaction_fit$residuals)^2) / (43 * 98 - 90)
sum(abs(bru_res_interaction_fit$residuals))
sqrt(sum((bru_res_interaction_fit$coefficients)^2))
sum(abs(bru_res_interaction_fit$coefficients))
bru_res_interaction_fit_ts = residuals_matrix(bru_res_interaction_fit$residuals, 98)
GNAR_AIC(bru_res_interaction_fit$residuals, 98, 18)
GNAR_BIC(bru_res_interaction_fit$residuals, 98, 30)
GNAR_HQC(bru_res_interaction_fit$residuals, 98, 30)
nacf_vals <- corbit_plot_cpp(residuals_matrix(bru_res_interaction_fit$residuals, 98), bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
length(bru_res_interaction_fit$coefficients)
sample_kurtosis(bru_struct_GNARfit$`Residuals Vector`)
sample_skewness(bru_res_interaction_fit$residuals)
hist(sample_kurtosis(bru_struct_GNARfit$`Residuals Matrix`))
hist(sample_skewness(bru_res_interaction_fit_ts))



bru_global_fit <- global_gnar_fit(bru_vts, bur_netGNAR, 1, 3, W_bru_data)
summary(bru_global_fit)
bru_global_fit_ts = residuals_matrix(bru_global_fit$residuals, 98)
sum((bru_global_fit$residuals)^2) / (93 * 43 - 4)
sum(abs(bru_global_fit$residuals))
sqrt(sum((bru_global_fit$coefficients)^2))
sum(abs(bru_global_fit$coefficients))
hist(bru_global_fit$residuals, breaks = 50)
GNAR_AIC(bru_global_fit$residuals, 98, 0)
GNAR_BIC(bru_global_fit$residuals, 98, 3)
GNAR_HQC(bru_global_fit$residuals, 98, 3)
nacf_vals <- corbit_plot_cpp(residuals_matrix(bru_global_fit$residuals, 98), bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))


bru_interaction_fit <- community_interaction_gnar_fit(bru_vts, bur_netGNAR, alpha_order, beta_order, W_bru_data, bru_data_community_indicators, 
                                                          bru_data_interactions)
summary(bru_interaction_fit)
bru_interaction_fit_ts = residuals_matrix(bru_interaction_fit$residuals, 98)
hist(bru_interaction_fit$residuals, breaks = 50)
sum((bru_interaction_fit$residuals)^2) / (43 * 98 - 34)
sum(abs(bru_interaction_fit$residuals))
sqrt(sum((bru_interaction_fit$coefficients)^2))
sum(abs(bru_interaction_fit$coefficients))
GNAR_AIC(bru_interaction_fit$residuals, 98, 0)
GNAR_BIC(bru_interaction_fit$residuals, 98, 30)
GNAR_HQC(bru_interaction_fit$residuals, 98, 30)
nacf_vals <- corbit_plot_cpp(residuals_matrix(bru_interaction_fit$residuals, 98), bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))


bru_interaction_fit <- community_interaction_gnar_fit(100.0 * MXD2, bur_netGNAR, alpha_order, beta_order, W_bru_data, bru_data_community_indicators, 
                                                      bru_data_interactions)
summary(bru_interaction_fit)
bru_interaction_fit_ts = residuals_matrix(bru_interaction_fit$residuals, 98)
hist(bru_interaction_fit$residuals, breaks = 50)
sum((bru_interaction_fit$residuals)^2)
sum(abs(bru_interaction_fit$residuals))
GNAR_AIC(bru_interaction_fit$residuals, 98, 0)
GNAR_BIC(bru_interaction_fit$residuals, 98, 30)
GNAR_HQC(bru_interaction_fit$residuals, 98, 30)




pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/bur_netGNAR_fitted_trend_plots.pdf")
par(mfrow=c(10, 10), mar=c(1, 1, 1, 1))
for (i in 1:98) {
  plot(bru_struct_GNARfit$`temporal trend fitted values`[, i], type='l', xlab='t', ylab=expression(expression(mu["i, t"])))
}
par(mfrow=c(1, 1))
dev.off()


pdf("/Users/danielsalnikov/Documents/PhD/my_stuff/papers/bur_netGNAR_GNARresiduals_plots.pdf")
par(mfrow=c(10, 10), mar=c(1, 1, 1, 1))
for (i in 1:98) {
  plot(bru_struct_GNARfit$`GNAR residuals`[, i], type='l', xlab='t', ylab=expression(expression(mu["i, t"])))
}
par(mfrow=c(1, 1))
dev.off()

#tt = data.frame(param=rep(c('alpha', 'beta'), 6), est=round(bru_struct_GNARfit$`GNAR parameters`, 3), community = rep(c('loc1_sat1', 'loc1_sat2', 'loc2_sat1', 'loc2_sat2', 'loc3_sat1', 'loc3_sat2'), 2))

gnar_com_diff_order_fit <- community_gnar_fit(MXD2, bur_netGNAR, alpha_order, beta_order, W_bru_data, bru_data_community_indicators)
GNAR_AIC(gnar_com_diff_order_fit$residuals, 98, 18)
summary(gnar_com_diff_order_fit)


bru_vts <- vapply(1:98, function(x) {return(100.0 * MXD2[, x] - colMeans(100.0 * MXD2)[x])}, rep(0.0, 44))
XX = build_GNAR_design(bru_vts - mu_hat_mat, get_r_stages_adjacency_list(as.matrix(bur_netGNAR), 4),  alpha_order, beta_order,  W_bru_data, bru_data_community_indicators, bru_data_interactions)
least_squares_solver(XX[, 2:19], XX[, 1])

mu_hat_mat <- temporal_trend_spline_smoother(spline_residuals = bru_vts[2:44, ], spar_vec = rep(0.65, 98), vts_dimension = 98, plag = 1)
mu_res_mat <- bru_vts - mu_hat_mat
mu_hat_res <- recursive_dmat_cols_cpp(lapply(1:44, function(x) {return(as.matrix(mu_res_mat[x, ]))}))
sum((mu_hat_res)^2)
GNAR_AIC(mu_hat_res, 98, 0)
#############
mu_hat_mat <- temporal_trend_spline_smoother(spline_residuals = bru_vts[2:44, ], spar_vec = rep(0.65, 98), vts_dimension = 98, plag = 1)
mu_res_mat <- bru_vts - mu_hat_mat
mu_hat_res <- recursive_dmat_cols_cpp(lapply(2:44, function(x) {return(as.matrix(mu_res_mat[x, ]))}))
sum(mu_hat_res^2) / (98 * 43 - 1)
sum(abs(mu_hat_res)) / (98 * 43 - 1)
GNAR_AIC(mu_hat_res, 98, 0)
############
mu_hat_mat <- temporal_trend_spline_smoother_cv(spline_residuals = bru_vts[2:44, ], vts_dimension = 98, plag = 1)
mu_res_mat <- bru_vts[2:44, ] - mu_hat_mat[2:44, ]
mu_hat_res <- recursive_dmat_cols_cpp(lapply(1:43, function(x) {return(as.matrix(mu_res_mat[x, ]))}))

sum(mu_hat_res^2) / (98 * 43 - 1)
sum(abs(mu_hat_res)) / (98 * 43 - 1)
##########
##############
##########
sparse_var_fit <- fitVAR(100 * MXD2, p = 1, nlambda = 100)
theta_mat = sparse_var_fit$A[[1]]
print(sum(sparse_var_fit$residuals[2:44, ]^2) / (98 * 43 - 672))
sum(sparse_var_fit$residuals[2:44, ]^2)
sum(abs(sparse_var_fit$residuals))
sparse_var_fit$lambda
sparse_var_res <- t(vapply(2:44, function(x) {bru_vts[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
hist(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), breaks = 50, xlab = 'u', main ='sparse VAR res')
sqrt(sum(theta_mat^2))
sum((sparse_var_res)^2) / (98 * 43 - 672)
b = theta_mat
b[abs(b) > 0.0] = 1
sum(rowSums(b))
GNAR_AIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_BIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_HQC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
summary(sparse_var_fit)
sum(bru_vts^2)
sum(sparse_var_fit$residuals[2:44, ]^2)
nacf_vals <- corbit_plot_cpp(sparse_var_fit$residuals[2:44, ], bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
VAR_mat_plot(1, 1, theta_mat)
eigen(theta_mat, only.values = TRUE)$values[1]
sample_kurtosis(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
sample_skewness(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
hist(sample_kurtosis(sparse_var_res))
hist(sample_skewness(sparse_var_res))
###############
sparse_var_fit <- fitVAR(bru_struct_GNARfit_cv$`temporal trend fitted values`, p = 1, nlambda = 100)
theta_mat = sparse_var_fit$A[[1]]
print(sum(sparse_var_fit$residuals[2:44, ]^2) / (98 * 43 - 672))
sum(sparse_var_fit$residuals[2:44, ]^2)
sum(abs(sparse_var_fit$residuals))
sparse_var_fit$lambda
sparse_var_res <- t(vapply(2:44, function(x) {bru_vts[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
hist(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), breaks = 50, xlab = 'u', main ='sparse VAR res')
sqrt(sum(theta_mat^2))
sum((sparse_var_res)^2) / (98 * 43 - 672)
b = theta_mat
b[abs(b) > 0.0] = 1
sum(rowSums(b))
GNAR_AIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_BIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_HQC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
summary(sparse_var_fit)
sum(bru_vts^2)
sum(sparse_var_fit$residuals[2:44, ]^2)
nacf_vals <- corbit_plot_cpp(sparse_var_fit$residuals[2:44, ], bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
VAR_mat_plot(1, 1, theta_mat)
eigen(theta_mat, only.values = TRUE)$values[1]
sample_kurtosis(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
sample_skewness(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
hist(sample_kurtosis(sparse_var_res))
hist(sample_skewness(sparse_var_res))
#######################
sparse_var_fit <- fitVAR(bru_struct_GNARfit_cv$`GNAR residuals`, p = 1, nlambda = 100)
theta_mat = sparse_var_fit$A[[1]]
print(sum(sparse_var_fit$residuals[2:44, ]^2) / (98 * 43 - 672))
sum(sparse_var_fit$residuals[2:44, ]^2)
sum(abs(sparse_var_fit$residuals))
sparse_var_fit$lambda
sparse_var_res <- t(vapply(2:44, function(x) {bru_vts[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
hist(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), breaks = 50, xlab = 'u', main ='sparse VAR res')
sqrt(sum(theta_mat^2))
sum((sparse_var_res)^2) / (98 * 43 - 672)
b = theta_mat
b[abs(b) > 0.0] = 1
sum(rowSums(b))
GNAR_AIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_BIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
GNAR_HQC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
summary(sparse_var_fit)
sum(bru_vts^2)
sum(sparse_var_fit$residuals[2:44, ]^2)
nacf_vals <- corbit_plot_cpp(sparse_var_fit$residuals[2:44, ], bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
VAR_mat_plot(1, 1, theta_mat)
eigen(theta_mat, only.values = TRUE)$values[1]
sample_kurtosis(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
sample_skewness(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})))
hist(sample_kurtosis(sparse_var_res))
hist(sample_skewness(sparse_var_res))
###################################
sparse_var_fit <- fitVAR(mu_res_mat, p = 1, nlambda = 100)
theta_mat = sparse_var_fit$A[[1]]
print(sum(sparse_var_fit$residuals^2))
sum(abs(sparse_var_fit$residuals))
sparse_var_fit$lambda
sparse_var_res <- t(vapply(2:44, function(x) {bru_vts[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
sqrt(sum(theta_mat^2))
sum((sparse_var_res)^2) / (98 * 43 - 672)
b = theta_mat
b[abs(b) > 0.0] = 1
q = sum(rowSums(b))
GNAR_AIC(recursive_dmat_cols_cpp(lapply(1:98, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), 98, 696)
hist(recursive_dmat_cols_cpp(lapply(1:43, function(j) {return(as.matrix(sparse_var_fit$residuals[2:44, j]))})), breaks = 50, xlab = 'u', main ='sparse VAR res')
summary(sparse_var_fit)
sum(bru_vts^2)
sum(sparse_var_fit$residuals[2:44, ]^2)
nacf_vals <- corbit_plot_cpp(sparse_var_fit$residuals[2:44, ], bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_vals))
VAR_mat_plot(1, 1, theta_mat)
sparse_var_mu_res <- t(vapply(2:44, function(x) {bru_struct_GNARfit$`temporal trend fitted values`[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
sum((sparse_var_mu_res)^2)
sparse_var_gnar_res <- t(vapply(2:44, function(x) {bru_struct_GNARfit$`GNAR residuals`[x, 1:98] - theta_mat %*% bru_vts[x - 1, 1:98]}, rep(0, 98)))
sum((sparse_var_gnar_res)^2)
###############################
############################
###################
naive_mean <- bru_vts[2:44, ]
naive_mean_vec <- recursive_dmat_cols_cpp(lapply(2:44, function(x) {return(as.matrix(bru_vts[x, ]))}))
sum(naive_mean_vec^2) / (98 * 43 - 1)
sum(abs(naive_mean_vec)) / (98 * 43 - 1)


naive_res <- 100 * t(vapply(2:44, function(x) {MXD2[x, 1:98] - MXD2[x - 1, 1:98]}, rep(0, 98)))
naive_res_vec <- recursive_dmat_cols_cpp(lapply(1:43, function(x) {return(as.matrix(naive_res[x, ]))}))
sum((naive_res)^2) / (98 * 43 - 1)
sum(abs(naive_res)) / (98 * 43 - 1)
##########################################
#########################################
########################################
ar_residuals <- vapply(1:98, function(x){forecast::auto.arima(bru_vts[1:44, x], d = 0, D = 0, max.p = 1, max.q = 0,
                      max.P = 0, max.Q = 0, stationary = TRUE, seasonal = FALSE, ic = 'bic',
                      allowmean = FALSE, allowdrift = FALSE, trace = FALSE)$residuals}, rep(0, 44))
ar_res_vec <- recursive_dmat_cols_cpp(lapply(2:44, function(x) {return(as.matrix(ar_residuals[x, ]))}))
sum(ar_res_vec^2) / (98 * 42)
sum(abs(ar_res_vec)) / (98 * 42)
################################################
solve(t(XX[, 2:19]) %*% XX[, 2:19], t(XX[, 2:19]) %*% XX[, 1])
XXridge = rbind(XX[, 2:13], 0.1 * diag(12))
yvec_ridge = c(XX[, 1], rep(0, 12))
solve(t(XXridge) %*% XXridge, t(XXridge) %*% yvec_ridge)
least_squares_solver(XXridge, yvec_ridge)

beta_coef_mat = trend_forecast_fit(bru_struct_GNARfit$`temporal trend fitted values`, 10, 24)
bru_rstage_adj_list <- get_r_stages_adjacency_list(as.matrix(bur_netGNAR), 4)
lagged_forecasts = cbind(bru_struct_GNARfit$`GNAR residuals`[24, ], matrix(rep(0.0, 98 * 20), nrow = 98, ncol = 20))
lagged_gnar_vals = cbind(bru_struct_GNARfit$`GNAR residuals`[24, ], matrix(rep(0.0, 98 * 20), nrow = 98, ncol = 20))
for (t in 1:20) {
  Rmat = vapply(1:3, function(i) {(bru_rstage_adj_list[[i]] * W_bru_data) %*% lagged_gnar_vals[, t] }, rep(0.0, 98))
  Rmat = cbind(lagged_gnar_vals[, t], Rmat)
  yhat = t(beta_coef_mat) %*% c(1, (t + 24)) + Rmat %*% bru_res_global_fit$coefficients
  lagged_forecasts[, t + 1] = yhat
  lagged_gnar_vals[, t + 1] = Rmat %*% bru_res_global_fit$coefficients
}


lagged_forecast_ar = cbind(100.0 * MXD2[24, ], matrix(rep(0.0, 98 * 20), nrow = 98, ncol = 20))
for (t in 1:20) {
  Rmat = vapply(1:2, function(i) {(bru_rstage_adj_list[[i]] * W_bru_data) %*% lagged_forecast_ar[, t] }, rep(0.0, 98))
  Rmat = cbind(lagged_forecast_ar[, t], Rmat)
  yhat = Rmat %*% bru_global_fit$coefficients
  lagged_forecast_ar[, t + 1] = yhat 
}

