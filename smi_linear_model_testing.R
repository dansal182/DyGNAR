T0 = 1000
omega = runif(1000)
fixed_effect = 2 * omega
trend_component = vapply(c(1:T0), function(x){return(log(x) + cos(x))}, 0.0)
ar_component =  gnar_sim_community(T0, fiveNet, W, c(1, 2), c(1, 1), list(c(0.23), c(0.20, 0.18)), list(list(c(0.47)), list(c(0.30), c(0.27))), 
                                   list(c(0, 1, 1, 1, 0), c(1, 0, 0, 0, 1)))[, 1]


Ysim = vapply(c(1:T0), function(x){return(fixed_effect[x] + trend_component[x] + ar_component[x])}, 0.0)
plot(c(1:T0), Ysim, 'l')

summary(lm(Ysim - trend_component - ar_component ~ omega))
for (j in 1:1000) {
  
}


View(Omat)
train_points = 1000
T0 = train_points
omega = runif(T0)
fixed_effect = 2 * omega
#Omat = matrix(c(rep(1, T0), omega), nrow = T0, ncol = 2)
a0 = rnorm(1)
phi = 0.4
ar_component = rep(0, T0)
for (t in 1:T0){
  if (t == 1) {
    ar_component[t] = phi * a0 + rnorm(1)
  } else {
    ar_component[t] = phi * ar_component[t - 1] + rnorm(1)
  }
}
acf(ar_component)
summary(lm(ar_component[2:T0] ~ ar_component[1:(T0-1)] + 0))
plot(c(1:T0), ar_component, 'l')
y = f(c(1:train_points)/T0, 'doppler') + 2 * cos(c(1:train_points)) + fixed_effect + ar_component
plot(c(1:T0), y, 'l')
for(j in 1:1000) {
  if (j == 1){
    eta_res = y[2:T0]
  } else {
    eta_res = y[2:T0] - mu_hat$y - phi * ar_component[1:(T0-1)] - 2 * cos(c(2:T0))
  }
  eta_hat = lm(eta_res ~ omega[2:T0] + 0)$coefficients
  #eta_hat[1] = eta_hat[1] -0.30
  mu_res = y[2:T0] - eta_hat * omega[2:T0] - phi * ar_component[1:(T0-1)] - 2 * cos(c(2:T0))
  mu_hat = smooth.spline(c(2:train_points), mu_res, cv = FALSE)
}
eta_hat
plot(mu_hat$x, mu_hat$y, 'l')
plot(omega,  y - mu_hat$y, 'l')
summary(lm(y - f(c(1:train_points)/T0, 'sin_prod') ~ omega + 0))
plot(c(1:train_points), f(c(1:train_points)/ T0, 'doppler'), 'l')
ar_residuals = y[2:T0] - mu_hat$y - eta_hat * omega[2:T0] - 2 * cos(c(2:T0))
plot(c(2:T0), ar_residuals, 'l')
acf(ar_residuals, 10)
summary(lm(ar_residuals[2:T0] ~ ar_residuals[1:(T0-1)] + 0))
omega_sort = sort(omega)
plot(omega_sort[2:T0], eta_res, 'l')


plot(c(1:T0), y, 'l')
  lines(mu_hat$y, col='blue', lwd=3)
  lines(2 * cos(c(2:T0)) + mu_hat$y, col='cyan')
  lines(2 * cos(c(2:T0)) + mu_hat$y + eta_hat * omega, col='orange')

par(mfrow=c(5, 1))
  plot(c(1:T0), y, 'l')
  plot(mu_hat$x, mu_hat$y, col='blue', lwd=3, 'l')
  plot(c(1:T0), 2 * cos(c(1:T0)), col='cyan', 'l')
  plot(omega[2:T0], eta_res, col='grey')+
    lines(omega_sort[2:T0], eta_hat * omega_sort[2:T0], col='orange', lwd=3)
  plot(c(2:T0), ar_residuals, col='forestgreen', 'l')
par(mfrow=c(1, 1))

plot(omega, y)
  lines(omega_sort, eta_hat * omega_sort)
plot(omega[2:T0], eta_res)  
  lines(omega[2:T0], eta_hat * omega[2:T0], col='orange')  
plot(c(1:T0), eta_hat * omega, 'l')  

types = c("sin_prod", "sinc", "doppler", "sinc", "sin_prod")
a <- vapply(1:5, function(x) {return(f(runif(200), types[x]))}, rep(0.0, 200))
