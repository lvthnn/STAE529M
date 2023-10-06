library(invgamma)
library(latex2exp)
library(rjags)

set.seed(42)

# true parameter values
mu_t     <- 10
delta_t  <- 1
sigma_t  <- 2
n <- m <- 50

# priors
pri    <- 100^2
a <- b <- 0.01

Y1 <- rnorm(n, mu_t, sigma_t)
Y2 <- rnorm(m, mu_t + delta_t, sigma_t)
Y <- c(Y1, Y2)

# Gibbs sampler
S <- 1e+4
theta <- matrix(NA, nrow = S, ncol = 3)
colnames(theta) <- c("mu", "delta", "sigma^2")

theta[1,] <- c(mean(Y), 0, var(Y))

for (s in 2:S) {
 
  tau   <- 1 / theta[s - 1, "sigma^2"]
  delta <- theta[s - 1, "delta"]

  # mu | Y, delta, sigma2
  denom <- (n + m) * tau + 1/pri
  num   <- tau * sum(Y[1:n]) + tau * sum(Y[(n+1):(n+m)] - delta)
  mu_s  <- rnorm(1, num/denom, 1/sqrt(denom))

  # delta | Y, mu, sigma2
  denom <- tau * m + 1/pri
  num   <- tau * sum(Y[(n+1):(n+m)] - mu_s)
  delta_s <- rnorm(1, num/denom, 1/sqrt(denom))

  # sigma2 | Y, mu, delta
  A <- a + (n + m) / 2
  B <- sum((Y[1:n] - mu_s)^2)/2 + sum((Y[(n+1):(n+m)] - mu_s - delta_s)^2)/2 + b
  sigma2_s <- rinvgamma(1, A, B)

  theta[s,] <- c(mu_s, delta_s, sigma2_s)

}

# summaries
apply(theta, 2, mean)
apply(theta, 2, quantile, c(0.025, 0.975))

# visualisations
labs <- paste0("$\\", colnames(theta), "$")
mas  <- apply(theta, 2, function(t) cumsum(t) / seq_along(t))
true <- c(mu_t, delta_t, sigma_t^2)

pdf("trace_plots_vars_ex_5.pdf", width = 12)
par(mfrow = c(1, 3))
for (i in 1:3) {
  plot(theta[,i], type = 'l', ylab = TeX(labs[i]))
  lines(mas[,i], type = 'l', col = 'red')
  abline(a = true[i], b = 0, col = 'blue')
}
dev.off()

pdf("density_vars_ex_5.pdf", width = 12)
par(mfrow = c(1, 3))
for (i in 1:3) {
  plot(density(theta[,i]), main = NA, xlab = TeX(labs[i]))
}
dev.off()


# JAGS simulation

model_string <- textConnection("model{
  # Likelihood 

  for (i in 1:n) {
    Y[i] ~ dnorm(mu, tau)
  }

  for (i in (n+1):(n+m)) {
    Y[i] ~ dnorm(mu + delta, tau)
  }

  # Priors
  
  mu     ~ dnorm(0, 0.0001)
  delta  ~ dnorm(0, 0.0001)
  tau    ~ dgamma(0.1, 0.1)
  sigma2 <- 1/tau
}")

inits <- list(mu = mean(Y), delta = 0, tau = 1/var(Y))
data  <- list(n = n, m = m, Y = Y)

model <- jags.model(model_string, inits = inits, data = data, n.chains = 2, quiet = TRUE)

update(model, 2000, progress.bar = "none")

params <- c("mu", "delta", "sigma2")

samples <- coda.samples(model, variable.names = params, n.iter = 10000, progress.bar = "none")
