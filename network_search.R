Bru_data <- load("~/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/Bru_data.RData")
library(GNAR)
library(igraph)
library(Rcpp)
library(RcppArmadillo)
library(fastCorbit)
library(glasso)
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
source("/Users/danielsalnikov/Documents/PhD/code_scripts/latex_output_and_utilities/latex_outputs_utilities.R")
require(fda)

# Compute the network adjacency matrix
locind_edt <- c(rep(1, 10), rep(2, 10), rep(3, 17), rep(5, 23), rep(6, 11), rep(4, 15), rep(6, 12))
# Set community membership values
loc1_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==1, 1, 0))}, 0)
loc1_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==2, 1, 0))}, 0)
loc2_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==3, 1, 0))}, 0)
loc2_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==4, 1, 0))}, 0)
loc3_sat1 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==5, 1, 0))}, 0)
loc3_sat2 <- vapply(1:98, function(x) {return(ifelse(locind_edt[x]==6, 1, 0))}, 0)
############################################################################
bru_data_community_indicators <- list(loc1_sat1, loc1_sat2, loc2_sat1, loc2_sat2, loc3_sat1, loc3_sat2)
bru_data_interactions <- list(c(2), c(1), c(4), c(3), c(6), c(5))
alpha_order <- rep(1, 6)
beta_order <- list(c(2), c(2), c(2), c(3), c(3), c(3))

UTM_coords <- matrix(c(X, Y), nrow = 98, ncol =2)
UTM_dist <- function(x1, x2) {
  a = as.numeric(t(x1 - x2) %*% (x1 - x2))
  return(sqrt(a))
}
UTM_dist_matrix <- matrix(rep(0.0, 98*98), nrow = 98, ncol = 98)
#UTM_dist_matrix <- vapply(1:98, function(i){UTM_dist_inner(UTM_coords[i, ])}, rep(0.0, 98))
for (i in 1:98) {
  UTM_dist_matrix[i, ] = vapply(1:98, function(j){UTM_dist(UTM_coords[i, ], UTM_coords[j, ])}, 0.0)
}
mean_dist = median(UTM_dist_matrix)
UTM_kernel = exp(- (UTM_dist_matrix / (2 * 200))^2)
net_mat <- glasso(s = UTM_kernel, rho = 0.12, penalize.diagonal = FALSE)
b =net_mat$wi
b[abs(b) > 0.0] = 1
diag(b) = 0
net_result = graph_from_adjacency_matrix(b, 'undirected')
bur_netGNAR <- igraphtoGNAR(net_result)
S1 <- as.matrix(GNARtoigraph(bur_netGNAR))
# Find the longest, shortest path and rmax
D = distances(graph_from_adjacency_matrix(S1))
D[!is.finite(D)] <- 0
rmax = max(D)
W_bru_data <- fastCorbit::weights_matrix_cpp(bur_netGNAR, rmax)

bru_struct_GNARfit <- temp_trend_GNARfit(longitudinal_network_data = 100.0 * MXD2, network = bur_netGNAR, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = W_bru_data, 
                                         covariate_groups = bru_data_community_indicators, interaction_groups = bru_data_interactions, rmax = rmax, iter = 10)

bru_res_interaction_fit <- community_interaction_gnar_fit(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR, alpha_order, beta_order, W_bru_data, bru_data_community_indicators, 
                                                          bru_data_interactions)
sum((bru_res_interaction_fit$residuals)^2)
summary(bru_res_interaction_fit)
sqrt(sum((bru_res_interaction_fit$coefficients)^2))
sum(abs(bru_res_interaction_fit$coefficients))
bru_res_interaction_fit_ts = residuals_matrix(bru_res_interaction_fit$residuals, 98)

#####################################
bru_data = 100.0 * MXD2
rho_optima = 0.1
RSS_optima = sum((bru_res_interaction_fit$residuals)^2)
for (h in seq(0.19, 0.29, 0.01)) {
  net_mat <- glasso(s = UTM_kernel, rho = h, penalize.diagonal = FALSE)
  b = net_mat$wi
  b[abs(b) > 0.0] = 1
  diag(b) = 0
  net_result = graph_from_adjacency_matrix(b, 'undirected')
  bur_netGNAR0 <- igraphtoGNAR(net_result)
  S1 <- as.matrix(GNARtoigraph(bur_netGNAR0))
  # Find the longest, shortest path and rmax
  D = distances(graph_from_adjacency_matrix(S1))
  D[!is.finite(D)] <- 0
  rmax0 = max(D)
  W_bru_data0 <- fastCorbit::weights_matrix_cpp(bur_netGNAR0, rmax0)
  bru_struct_GNARfit <- temp_trend_GNARfit(longitudinal_network_data = 100.0 * MXD2, network = bur_netGNAR0, alpha_order = alpha_order, beta_order = beta_order, weight_matrix = W_bru_data0, 
                                           covariate_groups = bru_data_community_indicators, interaction_groups = bru_data_interactions, rmax = rmax0, iter = 100)
  
  bru_res_interaction_fit <- community_interaction_gnar_fit(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR0, alpha_order, beta_order, W_bru_data0, bru_data_community_indicators, 
                                                            bru_data_interactions)
  rss_aux = sum((bru_res_interaction_fit$residuals)^2)
  print(rss_aux)
  if (rss_aux < RSS_optima) {
    RSS_optima = rss_aux
    rho_optima = h
  }
}
print(RSS_optima)
print(rho_optima)

