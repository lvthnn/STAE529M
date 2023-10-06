library(rjags)

n <- 10
a <- 1
Y <- 1:10

data <- list(n = n, Y = Y, a = a)

model_string <- textConnection("model{
  for (i in 1:n) {
    Y[i] ~ dnorm(0, tau[i])
  }

  for (i in 1:n) {
    tau[i] ~ dgamma(a, b)
    sigma[i] <- 1/tau[i]
  }

  b ~ dgamma(1, 1)
}")

inits <- list(tau = rep(1/var(Y), 10), b = 1)

model <- jags.model(model_string, inits = inits, data = data, n.chains = 2, quiet = TRUE)

update(model, 2000, progress.bar = "none")

params <- c("sigma", "b")

samples <- coda.samples(model, variable.names = params, n.iter = 1e+4, progress.bar = "none")
