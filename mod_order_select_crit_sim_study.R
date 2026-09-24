set.seed(2024)
sim1 <- gnar_sim_community(100, fiveNet, W, c(1, 2), c(1, 1), list(c(0.23), c(0.20, 0.18)), list(list(c(0.47)), list(c(0.30), c(0.27))), 
                           list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
gnar_com_diff_order <- community_gnar_fit(sim1, fiveNet, c(1, 2), list(c(1), c(1, 1)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
p = 10
lags = c(1:p)
crit_mat = matrix(rep(0, p * 3), nrow = p, ncol = 3)
for (m in lags) {
    gnar_com_diff_order <- community_gnar_fit(sim1, fiveNet, c(m, m), list(rep(1, m), rep(1, m)), W, list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))
    fit_res = gnar_com_diff_order$residuals
    q = length(gnar_com_diff_order$coefficients)
    crit_mat[m, 1]=GNAR_AIC(fit_res, 5, q)
    crit_mat[m, 2]=GNAR_BIC(fit_res, 5, q)
    crit_mat[m, 3]=GNAR_HQC(fit_res, 5, q)
}
View(crit_mat)
gnar_com_diff_order$coefficients
GNARcovariance_estimator_cpp(residuals_matrix(fit_res, 5))
summary(gnar_com_diff_order)
u = residuals_matrix(fit_res, 5)

n = 1000
X = matrix(rnorm(100*100, 1, 1), nrow = 100, ncol = 100)
Z = 1 / 10 * t(X) %*% X
Smat <- bdiag(lapply(1:(n/100), function(x) {return(Z)}))
u = c(chol(as.matrix(Smat)) %*% rnorm(n))
x = rnorm(n)
xx = matrix(c(rep(1:n), x), nrow = n, ncol = 2)
y = 1 + 2 * x + u
plot(x, y)
least_squares_solver(xx, y)
gls_estimator_cpp(xx, y, Z, n/100)
 
