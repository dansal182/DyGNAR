# Build Bru network based on spatial correlation and a Kriging Graphical Lasso approach
Bru_data <- load("~/Documents/PhD/my_stuff/papers/communal_gnar_costa_rica/dynamic_GNAR_revised/Bru_data.RData")
library(GNAR)
library(igraph)
library(Rcpp)
library(RcppArmadillo)
library(fastCorbit)
library(glasso)
library(sparsevar)
sourceCpp('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/nacf_inner.cpp')
sourceCpp('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/get_k_stage_adjacency_list.cpp')
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source('/Users/danielsalnikov/Documents/PhD/code_scripts/fast_nacf_cpp_files/nacf_pnacf_rwrapper.r')
source('/Users/danielsalnikov/Documents/PhD/code_scripts/network_autocorrelation/corbit_scripts/weight_adjustment_methods.R')
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/structural_GNAR_fitting_methods.r")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/struc_GNAR_sim_setup.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_interactions_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/one-stepcommunity_GNAR_prediction.r", encoding = "UTF-8")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/community_gnar_fitting_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/communal_gnar_simulation_methods.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/global_gnar_fitting.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/communal_gnar/com_GNARfit_tv_weights.R")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNAR_simulation")
source("/Users/danielsalnikov/Documents/PhD/code_scripts/vs_code_editing/community_GNARfit.r")
require(fda)
############################################
############################################
############################################
# Start by computing distance using longitude and latitude.
plot(bru2)
plot(X, Y)
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
View(UTM_dist_matrix)
mean_dist = median(UTM_dist_matrix)
UTM_kernel = exp(- (UTM_dist_matrix / (2 * 200))^2)
hist(UTM_dist_matrix[UTM_dist_matrix > 0], breaks = 50)
net_mat <- glasso(s = UTM_kernel, rho = 0.1, penalize.diagonal = FALSE)
View(net_mat$wi)
b =net_mat$wi
b[abs(b) > 0.0] = 1
diag(b) = 0
net_result = graph_from_adjacency_matrix(b, 'undirected')
V(net_result)$Xcoord <- X
V(net_result)$Ycoord <- Y
V(net_result)$color <- rep('black', 98)
node_labels <- vapply(1:98, function(i) {return(paste0('(', as.character(i), ')'))}, '')
group_legend_labels <- c(vapply(1:3, function(i) {c(paste0('Location: ', as.character(i), ' & Satellite: 1'), 
                                                    paste0('Location: ', as.character(i), ' & Satellite: 2'))}, c('', '')))
location <- matrix(c(V(net_result)$Xcoord, V(net_result)$Ycoord), ncol = 2)
plot.igraph(net_result, vertex.label='', vertex.label.cex = .5, vertex.color = 'black', edge.width = 0.5, vertex.label.color = 'black', 
            vertex.shape = 'circle', vertex.size = 500, mark.col = c('#ADD5D7', '#B3DBB6', '#FDB06C', '#7871A8', '#AACAF5', '#F6C5E2'),
            mark.groups = list(c(1:10), c(11:20), c(21:37), c(72:86), c(38:60), c(c(61:71), c(87:98))), mark.border = NA,
            layout = location, rescale=FALSE, asp = 0, xlim = range(V(net_result)$Xcoord), ylim = range(V(net_result)$Ycoord))
title(main = "\n Brumadinho Dam Location Network", cex.main=3)
legend("bottomright",
       legend = group_legend_labels,
       pch = 16,
       text.col = 'black',
       pt.cex = 2,
       cex = 2,
       bty = 'o', 
       horiz = FALSE,
       col = c('#ADD5D7', '#B3DBB6', '#FDB06C', '#7871A8', '#AACAF5', '#F6C5E2')
       )
# layout = location, rescale=FALSE, asp = 0, xlim = range(V(net_result)$Xcoord), ylim = range(V(net_result)$Ycoord),
# Load necessary libraries
# install.packages(c("sf", "ggplot2", "ggspatial"))
library(sf)
library(ggplot2)
library(ggspatial)
library(prettymapr)
# 1. Define your UTM coordinates
# Easting: 591817, Northing: 7775020
# Zone 23S (EPSG: 31983 for SIRGAS 2000 / UTM zone 23S)
coords_df <- data.frame(
  id = "Target Location",
  x = X,
  y = Y
)

# 2. Convert to an 'sf' object 
# We set the CRS to 31983, which is the standard for UTM 23S in Brazil
point_sf <- st_as_sf(coords_df, coords = c("x", "y"), crs = 31983)

# 3. Create the plot
ggplot() +
  # Add a basemap layer (requires internet connection)
  annotation_map_tile(type = "osm", zoom = 18) + 
  # Plot the point
  geom_sf(data = point_sf, color = "black", size = 4, shape = 16) +
  # Formatting
  ggtitle("Ground Displacements at the Brumadinho Dam by Location in Belo Horizonte, Brazil (UTM 23S) \n")+
  theme_minimal()+
  theme(plot.title = element_text(color="black", size=22, face="bold"),
        axis.text.x = element_text(color="black", size=18, face="bold"),
        axis.text.y = element_text(color="black", size=18, face="bold"))


#######################
bur_netGNAR <- igraphtoGNAR(net_result)

nacf_raw_data <- corbit_plot_cpp(MXD2, bur_netGNAR, 10, 4, W_bru_data)
max(abs(nacf_raw_data))
corbit_plot_cpp(MXD2, bur_netGNAR, 10, 4, W_bru_data, partial = "yes")

corbit_plot_cpp(MXD2, bru2, 10, 6, fastCorbit::weights_matrix_cpp(bru2, 6))
corbit_plot_cpp(MXD2, bru2, 10, 6, fastCorbit::weights_matrix_cpp(bru2, 6), partial = "yes")

# Community-wise decomposition for each vector time series
loc1_sat1_nts <- t(vapply(seq(1:44), function(x) {return(loc1_sat1 * MXD2[x, ])}, rep(0.0, 98)))
loc1_sat2_nts <- t(vapply(seq(1:44), function(x) {return(loc1_sat2 * MXD2[x, ])}, rep(0.0, 98)))
loc2_sat1_nts <- t(vapply(seq(1:44), function(x) {return(loc2_sat1 * MXD2[x, ])}, rep(0.0, 98)))
loc2_sat2_nts <- t(vapply(seq(1:44), function(x) {return(loc2_sat2 * MXD2[x, ])}, rep(0.0, 98)))
loc3_sat1_nts <- t(vapply(seq(1:44), function(x) {return(loc3_sat1 * MXD2[x, ])}, rep(0.0, 98)))
loc3_sat2_nts <- t(vapply(seq(1:44), function(x) {return(loc3_sat2 * MXD2[x, ])}, rep(0.0, 98)))


# R-Corbit plots
# NACF
# PNACF
r_corbit_plot_cpp(list(loc1_sat1_nts, loc1_sat2_nts, loc2_sat1_nts, loc2_sat2_nts, loc3_sat1_nts, loc3_sat2_nts), list(bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR),
                  max_lag = 10, max_stage = 4, list(W_bru_data, W_bru_data, W_bru_data, W_bru_data, W_bru_data, W_bru_data), c("loc1_sat1", "loc1_sat2", "loc2_sat1", "loc2_sat2", "loc3_sat1", "loc3_sat2"), 
                  same_net = "no", partial = "no")

r_corbit_plot_cpp(list(loc1_sat1_nts, loc1_sat2_nts, loc2_sat1_nts, loc2_sat2_nts, loc3_sat1_nts, loc3_sat2_nts), list(bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR, bur_netGNAR),
                  max_lag = 5, max_stage = 4, list(W_bru_data, W_bru_data, W_bru_data, W_bru_data, W_bru_data, W_bru_data), c("loc1_sat1", "loc1_sat2", "loc2_sat1", "loc2_sat2", "loc3_sat1", "loc3_sat2"), 
                  same_net = "no", partial = "yes")

# Plots after removing temporal trend
corbit_plot_cpp(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR, 10, 4, W_bru_data)
corbit_plot_cpp(bru_struct_GNARfit$`GNAR residuals`, bur_netGNAR, 10, 4, W_bru_data, partial = "yes")

#votes_vts_matrix_diff <- votes_vts_matrix[2:12, ] - votes_vts_matrix[1:11, ]
#colnames(votes_vts_matrix_diff) <- usa_net_data$Tag
#red_vts_diff <- t(vapply(seq(1:11), function(x) {return(red_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
#blue_vts_diff <- t(vapply(seq(1:11), function(x) {return(blue_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))
#swing_vts_diff <- t(vapply(seq(1:11), function(x) {return(swing_states * votes_vts_matrix_diff[x, ])}, rep(0, 51)))

# Community-wise decomposition for each vector time series
loc1_sat1_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, 1:10]
loc1_sat2_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, 2:20]
loc2_sat1_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, 21:37]
loc2_sat2_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, 72:86]
loc3_sat1_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, 38:60]
loc3_sat2_trend <- bru_struct_GNARfit$`temporal trend fitted values`[, c(61:71, 87:98)]
par(mfrow=c(2, 3))
#########
mean_loc1_sat1_trend <- rowMeans(loc1_sat1_trend)
loc1_sat1_quantile <- vapply(1:44, function(x) {quantile(loc1_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat1_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 1", xlab = '', cex.main = 2, cex.lab = 2, ylab = "Trend Estimate" )+
  lines(c(1:44), c(loc1_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 1:10) {
      #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
      #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    points(c(1:44), 100.0 * (loc1_sat1_nts[1:44, j] - mean(loc1_sat1_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
############################
mean_loc2_sat1_trend <- rowMeans(loc2_sat1_trend)
loc2_sat1_quantile <- vapply(1:44, function(x) {quantile(loc2_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 21:37) {
    points(c(1:44), 100.0 * (loc2_sat1_nts[1:44, j] - mean(loc2_sat1_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
      #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
      #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
########
mean_loc3_sat1_trend <- rowMeans(loc3_sat1_trend)
loc3_sat1_quantile <- vapply(1:44, function(x) {quantile(loc3_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc3_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 38:60) {
    points(c(1:44), 100.0 * (loc3_sat1_nts[1:44, j] - mean(loc3_sat1_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1:44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
      #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
      #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
################
mean_loc1_sat2_trend <- rowMeans(loc1_sat2_trend)
loc1_sat2_quantile <- vapply(1:44, function(x) {quantile(loc1_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat2_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 2", ylab = "Trend Estimate", xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc1_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 11:20) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc1_sat2_nts[1:44, j] - mean(loc1_sat2_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
########
mean_loc2_sat2_trend <- rowMeans(loc2_sat2_trend)
loc2_sat2_quantile <- vapply(1:44, function(x) {quantile(loc2_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat2_trend[1:44], type ='l', ylim = c(-2, 2.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 72:86) {
    points(c(1:44), 100.0 * (loc2_sat2_nts[1:44, j] - mean(loc2_sat2_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
      #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
      #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
#############
mean_loc3_sat2_trend <- rowMeans(loc3_sat2_trend)
loc3_sat2_quantile <- vapply(1:44, function(x) {quantile(loc3_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat2_trend[1:44], type ='l', ylim = c(-1.5, 1.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 2", ylab = expression(hat(mu["t"])), xlab = "Time: t" , cex.main = 2, cex.lab = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in c(61:71, 87:98)) {
      #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
      #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc3_sat2_nts[1:44, j] - mean(loc3_sat2_nts[1:44, j])) -  bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
par(mfrow = c(1, 1))
###################
###################
#################
# plotting the resoiduals and model diagnostics
##############
#############
blue_pal <- paletteer::paletteer_c("ggthemes::Classic Blue", 98)
orange_pal <- paletteer::paletteer_c("ggthemes::Classic Orange", 98)
gray_pal = paletteer::paletteer_c("ggthemes::Classic Gray", 98)
viridis_pal <- paletteer::paletteer_c("viridis::viridis", 98)

################
##############
residual_nts <- residuals_matrix(bru_res_interaction_fit$residuals, 98)
plot(rep(0, 43), type='b', col='black', pch = 16, lwd = 3, ylim = c(-3, 3))+
  points(c(1:44), residual_nts[, 1], col = blue_pal[1], pch = 16)
  for (j in 2:98){
    points(c(1:44), residual_nts[, j], col = blue_pal[j], pch = 16)
  }
  lines(c(1:44), rep(0, 43), col = 'black', lwd = 3)+
    lines(c(1:44), rep(1.5, 43), col = 'orange', lwd = 3)+
    lines(c(1:44), rep(-1.5, 43), col = 'orange', lwd = 3)

  
corbit_plot_cpp(residual_nts, bur_netGNAR, 10, 4, W_bru_data)
corbit_plot_cpp(residual_nts, bur_netGNAR, 20, 4, W_bru_data, partial = "yes")


plot(colMeans(lagged_forecasts), type='b', col = 'cyan', pch = 16, lwd = 5, ylim = c(-3.5, 3.5))
for (j in 1:98){
  points(c(1:21), lagged_forecasts[j, ], col = blue_pal[j], pch = 16)
  points(c(1:21), lagged_forecast_ar[j, ], col = orange_pal[j], pch = 17)
  lines(c(1:21), 100 * MXD2[24:44, j], type = 'h')
}
lines(c(1:21), rep(-1, 21) + colMeans(lagged_forecasts), col = 'cyan', lwd = 3, type = 'b')
lines(c(1:21), rep(1, 21) + colMeans(lagged_forecasts), col = 'cyan', lwd = 3, type = 'b')
############################
############################
############################
par(mfrow =c(2, 3))
plot(mean_loc1_sat1_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc1_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
############################
plot(mean_loc2_sat1_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc2_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
########
plot(mean_loc3_sat1_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc3_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
################
plot(mean_loc1_sat2_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc1_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
########
plot(mean_loc2_sat2_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc2_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
#############
plot(mean_loc3_sat2_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t" )+
  lines(c(1:44), c(loc3_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)
par(mfrow = c(1, 1))
###################
###################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-3, 3))
l1_norm = 0
for (j in 1:98) {
  points(c(1:43), c(bru_res_global_fit_ts[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_res_global_fit_ts[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
  l1_norm = l1_norm + sum(abs(least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_res_global_fit_ts[, j])))
}
print(l1_norm)
##############
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-2, 2))
l1_norm = 0
for (j in 1:98) {
  points(c(1:43), c(bru_res_interaction_fit_ts[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_res_interaction_fit_ts[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
  l1_norm = l1_norm + sum(abs(least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_res_interaction_fit_ts[, j])))
}
print(paste0("Sum of residual linear trend coefficients: ", as.character(l1_norm)))
#####################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-2, 2))
l1_norm = 0
for (j in 1:98) {
  points(c(1:43), c(bru_struct_GNARfit_cv$`Residuals Matrix`[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_struct_GNARfit_cv$`Residuals Matrix`[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
  l1_norm = l1_norm + sum(abs(least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_struct_GNARfit_cv$`Residuals Matrix`[, j])))
}
print(paste0("Sum of residual linear trend coefficients: ", as.character(l1_norm)))
#################
plot(rep(0, 44), type = 'l', col = 'blue', ylim = c(-3, 3))
for (j in 1:98) {
  points(c(1:44), c(bru_fit_res_ts[, j]), pch = 16)
  lines(c(1:44), vapply(1:44, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 44), c(1:44)), bru_fit_res_ts[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
}
############################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-3, 3))
for (j in 1:98) {
  points(c(1:43), c(bru_global_fit_ts[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_global_fit_ts[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
}
################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-3, 3))
l1_norm = 0
for (j in 1:98) {
  points(c(1:43), c(bru_interaction_fit_ts[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_interaction_fit_ts[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
  l1_norm = l1_norm + sum(abs(least_squares_solver(cbind(rep(1, 43), c(1:43)), bru_interaction_fit_ts[, j])))
}
print(l1_norm)
################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-2, 2))
for (j in 1:98) {
  points(c(1:43), c(sparse_var_fit$residuals[2:44, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), sparse_var_fit$residuals[2:44, j])}, 0.0), col = orange_pal[j],
        lty = 2)
}
print(l1_norm)
################
plot(rep(0, 43), type = 'l', col = 'blue', ylim = c(-2, 2))
l1_norm = 0
for (j in 1:98) {
  points(c(1:43), c(sparse_var_res[, j]), pch = 16)
  lines(c(1:43), vapply(1:43, function(i) {c(1, i) %*% least_squares_solver(cbind(rep(1, 43), c(1:43)), sparse_var_res[, j])}, 0.0), col = orange_pal[j],
        lty = 2)
  l1_norm = l1_norm + sum(abs(least_squares_solver(cbind(rep(1, 43), c(1:43)), sparse_var_res[, j])))
}
print(l1_norm)
################
##################
##################
par(mfrow=c(2, 3))
#########
mean_loc1_sat1_trend <- rowMeans(loc1_sat1_trend)
loc1_sat1_quantile <- vapply(1:44, function(x) {quantile(loc1_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat1_trend[1:44], type ='l', ylim = c(-1, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 1", xlab = '', cex.main = 2, cex.lab = 2, ylab = "Trend Estimate" )+
  lines(c(1:44), c(loc1_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 1:10) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    points(c(1:44), 100.0 * (loc1_sat1_nts[1:44, j] - mean(loc1_sat1_nts[1:44, j])) - (bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + bru_struct_GNARfit$`GNAR residuals`[1:44, j]),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
############################
mean_loc2_sat1_trend <- rowMeans(loc2_sat1_trend)
loc2_sat1_quantile <- vapply(1:44, function(x) {quantile(loc2_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 21:37) {
    points(c(1:44), 100.0 * (loc2_sat1_nts[1:44, j] - mean(loc2_sat1_nts[1:44, j])) - (bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + bru_struct_GNARfit$`GNAR residuals`[1:44, j]),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
########
mean_loc3_sat1_trend <- rowMeans(loc3_sat1_trend)
loc3_sat1_quantile <- vapply(1:44, function(x) {quantile(loc3_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc3_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 38:60) {
    points(c(1:44), 100.0 * (loc3_sat1_nts[1:44, j] - mean(loc3_sat1_nts[1:44, j])) - (bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + bru_struct_GNARfit$`GNAR residuals`[1:44, j]),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
################
mean_loc1_sat2_trend <- rowMeans(loc1_sat2_trend)
loc1_sat2_quantile <- vapply(1:44, function(x) {quantile(loc1_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat2_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 2", ylab = "Trend Estimate", xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc1_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 11:20) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc1_sat2_nts[1:44, j] - mean(loc1_sat2_nts[1:44, j])) - (bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + bru_struct_GNARfit$`GNAR residuals`[1:44, j]),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
########
mean_loc2_sat2_trend <- rowMeans(loc2_sat2_trend)
loc2_sat2_quantile <- vapply(1:44, function(x) {quantile(loc2_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat2_trend[1:44], type ='l', ylim = c(-2, 2.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 72:86) {
    points(c(1:44), 100.0 * (loc2_sat2_nts[1:44, j] - mean(loc2_sat2_nts[1:44, j])) - bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] - bru_struct_GNARfit$`GNAR residuals`[1:44, j],
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
#############
mean_loc3_sat2_trend <- rowMeans(loc3_sat2_trend)
loc3_sat2_quantile <- vapply(1:44, function(x) {quantile(loc3_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat2_trend[1:44], type ='l', ylim = c(-1.5, 1.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 2", ylab = expression(hat(mu["t"])), xlab = "Time: t" , cex.main = 2, cex.lab = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in c(61:71, 87:98)) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc3_sat2_nts[1:44, j] - mean(loc3_sat2_nts[1:44, j])) - (bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + bru_struct_GNARfit$`GNAR residuals`[1:44, j]), col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
par(mfrow = c(1, 1))
#########################
########################
###### Plots for DGNAR fit with chosen spline lambdas
par(mfrow=c(2, 3))
#########
mean_loc1_sat1_trend <- rowMeans(loc1_sat1_trend)
loc1_sat1_quantile <- vapply(1:44, function(x) {quantile(loc1_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat1_trend[1:44], type ='l', ylim = c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 1", xlab = '', cex.main = 2, cex.lab = 2, ylab = "Trend Estimate" )+
  lines(c(1:44), c(loc1_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 1:10) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    points(c(1:44), 100.0 * (loc1_sat1_nts[1:44, j] - mean(loc1_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
############################
mean_loc2_sat1_trend <- rowMeans(loc2_sat1_trend)
loc2_sat1_quantile <- vapply(1:44, function(x) {quantile(loc2_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 21:37) {
    points(c(1:44), 100.0 * (loc2_sat1_nts[1:44, j] - mean(loc2_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
########
mean_loc3_sat1_trend <- rowMeans(loc3_sat1_trend)
loc3_sat1_quantile <- vapply(1:44, function(x) {quantile(loc3_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat1_trend[1:44], type ='l', ylim = c(-1.5, 1.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc3_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 38:60) {
    points(c(1:44), 100.0 * (loc3_sat1_nts[1:44, j] - mean(loc3_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1:44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
################
mean_loc1_sat2_trend <- rowMeans(loc1_sat2_trend)
loc1_sat2_quantile <- vapply(1:44, function(x) {quantile(loc1_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat2_trend[1:44], type ='l', ylim = c(-1.5, 1), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 2", ylab = "Trend Estimate", xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc1_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 11:20) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc1_sat2_nts[1:44, j] - mean(loc1_sat2_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
########
mean_loc2_sat2_trend <- rowMeans(loc2_sat2_trend)
loc2_sat2_quantile <- vapply(1:44, function(x) {quantile(loc2_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat2_trend[1:44], type ='l', ylim = c(-2, 2.5), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 72:86) {
    points(c(1:44), 100.0 * (loc2_sat2_nts[1:44, j] - mean(loc2_sat2_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
#############
mean_loc3_sat2_trend <- rowMeans(loc3_sat2_trend)
loc3_sat2_quantile <- vapply(1:44, function(x) {quantile(loc3_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat2_trend[1:44], type ='l', ylim = c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 2", ylab = expression(hat(mu["t"])), xlab = "Time: t" , cex.main = 2, cex.lab = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in c(61:71, 87:98)) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc3_sat2_nts[1:44, j] - mean(loc3_sat2_nts[1:44, j])), col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
par(mfrow = c(1, 1))
##########
#########################
########################
###### Plots for DyGNAR fit with cv spline lambdas
par(mfrow=c(2, 3))
#########
mean_loc1_sat1_trend <- rowMeans(loc1_sat1_trend)
loc1_sat1_quantile <- vapply(1:44, function(x) {quantile(loc1_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat1_trend[1:44], type ='l', ylim = c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 1", xlab = '', cex.main = 2, cex.lab = 2, ylab = "" )+
  lines(c(1:44), c(loc1_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 1:10) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    #lines(c(1:44), .63 + bru_struct_GNARfit$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)+
    points(c(1:44), 100.0 * (loc1_sat1_nts[1:44, j] - mean(loc1_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
############################
mean_loc2_sat1_trend <- rowMeans(loc2_sat1_trend)
loc2_sat1_quantile <- vapply(1:44, function(x) {quantile(loc2_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat1_trend[1:44], type ='l', ylim =  c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 21:37) {
    points(c(1:44), 100.0 * (loc2_sat1_nts[1:44, j] - mean(loc2_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
########
mean_loc3_sat1_trend <- rowMeans(loc3_sat1_trend)
loc3_sat1_quantile <- vapply(1:44, function(x) {quantile(loc3_sat1_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat1_trend[1:44], type ='l', ylim =  c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 1", ylab = expression(hat(mu)["t"]), xlab = "", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc3_sat1_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat1_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 38:60) {
    points(c(1:44), 100.0 * (loc3_sat1_nts[1:44, j] - mean(loc3_sat1_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
################
mean_loc1_sat2_trend <- rowMeans(loc1_sat2_trend)
loc1_sat2_quantile <- vapply(1:44, function(x) {quantile(loc1_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc1_sat2_trend[1:44], type ='l', ylim =  c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 1 Satellite 2", ylab = "", xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc1_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc1_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 11:20) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc1_sat2_nts[1:44, j] - mean(loc1_sat2_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
########
mean_loc2_sat2_trend <- rowMeans(loc2_sat2_trend)
loc2_sat2_quantile <- vapply(1:44, function(x) {quantile(loc2_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc2_sat2_trend[1:44], type ='l', ylim =  c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 2 Satellite 2", ylab = expression(hat(mu)["t"]), xlab = "Time: t", cex.main = 2, cex.lab = 2 )+
  lines(c(1:44), c(loc2_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc2_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in 72:86) {
    points(c(1:44), 100.0 * (loc2_sat2_nts[1:44, j] - mean(loc2_sat2_nts[1:44, j])),
           col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
    #lines(c(1:44), -.63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
  }
#############
mean_loc3_sat2_trend <- rowMeans(loc3_sat2_trend)
loc3_sat2_quantile <- vapply(1:44, function(x) {quantile(loc3_sat2_trend[x, ], probs = c(0.0, 1))}, rep(0.0, 2))
plot(mean_loc3_sat2_trend[1:44], type ='l', ylim =  c(-2, 2), col = 'forestgreen', lwd=3, pch = 16,
     main = "Location 3 Satellite 2", ylab = expression(hat(mu["t"])), xlab = "Time: t" , cex.main = 2, cex.lab = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[1, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  lines(c(1:44), c(loc3_sat2_quantile[2, 1:44]), col = 'blue', lwd = 2, type = 'l', lty = 2)+
  for (j in c(61:71, 87:98)) {
    #lines(c(1:44), -.63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    #lines(c(1:44), .63 + bru_struct_GNARfit_cv$`temporal trend fitted values`[1:44, j] + 
    #        bru_struct_GNARfit$`GNAR residuals`[1:44, j], col = 'orange', type = 'h', lwd = 1)
    points(c(1:44), 100.0 * (loc3_sat2_nts[1:44, j] - mean(loc3_sat2_nts[1:44, j])), col = 'black', pch = 16, cex = 1)
    lines(c(1:44), bru_struct_GNARfit_cv$`temporal trend fitted values`[1: 44, j], col = viridis_pal[j], lwd = 1, type = 'l', lty = 1)
  }
par(mfrow = c(1, 1))
