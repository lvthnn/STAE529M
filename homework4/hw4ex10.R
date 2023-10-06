library(rjags)

set.seed(42)

Y    <- c(64, 13, 33, 18, 30, 20)
dT   <- length(Y)

data  <- list(Y = Y, dT = dT)

model_string <- textConnection("model{
  for (t in 1:dT) {
    Y[t] ~ dpois(exp(alpha + beta * t))
  }

  alpha ~ dnorm(0, 0.01)
  beta  ~ dnorm(0, 0.01)

}")

inits <- list(alpha = 1, beta = 0)

model <- jags.model(model_string, inits = inits, data = data, n.chains = 2, quiet = TRUE)

update(model, 1e+4, progress.bar = "none")

samples <- coda.samples(model, variable.names = c("alpha", "beta"), n.iter = 3e+4, progress.bar = "none")

summary(samples)

plot(samples)

geweke.diag(samples[[1]]) # good 
gelman.diag(samples) # good
